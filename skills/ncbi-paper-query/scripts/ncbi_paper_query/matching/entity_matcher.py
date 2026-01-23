"""
NCBI Paper Query - Entity Matcher

Semantic entity matching using Claude CLI for intelligent interpretation.
Falls back to rule-based matching when Claude CLI is unavailable.
"""

import re
import json
import hashlib
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Union

from ..config import CLAUDE_CLI_AVAILABLE

logger = logging.getLogger(__name__)


class EntityMatcher:
    """
    Semantic entity matching using Claude CLI for intelligent interpretation.
    Falls back to rule-based matching when Claude CLI is unavailable.
    """

    # Cache file for Claude responses to avoid redundant CLI calls
    CACHE_FILE = ".entity_match_cache.json"

    def __init__(self, use_claude: bool = True, cache_dir: str = None):
        """
        Initialize EntityMatcher.

        Args:
            use_claude: Whether to use Claude CLI for semantic matching (default True)
            cache_dir: Directory for cache file (default current directory)
        """
        self.use_claude = use_claude and CLAUDE_CLI_AVAILABLE
        self.cache_dir = Path(cache_dir) if cache_dir else Path(".")
        self.cache = self._load_cache()

        if self.use_claude:
            logger.info("Claude CLI-based entity matching enabled")
        else:
            if not CLAUDE_CLI_AVAILABLE:
                logger.info("Claude CLI not found, using rule-based matching")

        # Fallback: organism patterns for rule-based matching
        self.organism_patterns = {
            "human": r"\b(human|homo sapiens|patient|subjects?|clinical|trial)\b",
            "mouse": r"\b(mouse|mice|mus musculus|murine|c57bl|balb)\b",
            "rat": r"\b(rat|rats|rattus norvegicus|wistar|sprague)\b",
        }

    def _load_cache(self) -> Dict:
        """Load cached Claude responses."""
        cache_path = self.cache_dir / self.CACHE_FILE
        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    return json.load(f)
            except Exception:
                return {}
        return {}

    def _save_cache(self):
        """Save Claude responses to cache."""
        cache_path = self.cache_dir / self.CACHE_FILE
        try:
            with open(cache_path, 'w') as f:
                json.dump(self.cache, f, indent=2)
        except Exception as e:
            logger.debug(f"Failed to save entity match cache: {e}")

    def _get_cache_key(self, text: str, disease: str, tissue: str,
                       organisms: List[str], cell_type: str) -> str:
        """Generate cache key from input parameters."""
        # Use hash of text to keep keys manageable
        text_hash = hashlib.md5(text.encode()).hexdigest()[:16]
        org_str = ",".join(sorted(organisms)) if organisms else ""
        key_str = f"{text_hash}|{disease}|{tissue}|{org_str}|{cell_type or ''}"
        return hashlib.md5(key_str.encode()).hexdigest()

    def _call_claude_cli(self, text: str, disease: str, tissue: str,
                         organisms: List[str], cell_type: str) -> Dict:
        """
        Use Claude CLI to semantically match publication against search criteria.

        Returns dict with match results and reasoning.
        """
        org_list = ", ".join(organisms) if organisms else "any"
        cell_str = cell_type if cell_type else "not specified"

        prompt = f"""Analyze this scientific publication text and determine if it matches the search criteria.

SEARCH CRITERIA:
- Target Disease/Indication: {disease}
- Target Tissue/Organ: {tissue}
- Target Organism(s): {org_list}
- Target Cell Type: {cell_str}

PUBLICATION TEXT:
{text[:3000]}

INSTRUCTIONS:
Evaluate semantic relevance, not just keyword matching. Consider:
1. Disease: Include related conditions, subtypes, phenotypes, complications, and therapeutic areas
   - Example: "fibrotic stricture" is related to Crohn's disease
   - Example: "anti-TNF therapy" implies IBD treatment context
   - Example: "Montreal classification B2" indicates Crohn's stricturing phenotype
2. Tissue: Include anatomical regions, cell layers, and related structures
   - Example: "colonic mucosa" matches "intestine"
   - Example: "ileal Peyer's patches" matches "intestine"
3. Organism: Look for species mentions, model organisms, patient populations
4. Cell Type: Match specific cell types, subtypes, and related populations

Respond ONLY with this exact JSON format (no other text):
{{
    "disease_match": true,
    "disease_found": "specific disease term found or inferred",
    "disease_confidence": 0.95,
    "disease_reasoning": "brief explanation",
    "tissue_match": true,
    "tissue_found": "specific tissue term found or inferred",
    "tissue_confidence": 0.90,
    "tissue_reasoning": "brief explanation",
    "organism_match": true,
    "organism_found": "specific organism found",
    "organism_confidence": 1.0,
    "cell_type_match": true,
    "cell_type_found": "specific cell type found or null",
    "cell_type_confidence": 1.0,
    "overall_relevance": 0.92,
    "relevance_reasoning": "brief overall assessment"
}}

Be generous with matches - if the paper is clearly relevant to the therapeutic area, match it even without exact keyword matches."""

        try:
            # Call Claude CLI with the prompt
            result = subprocess.run(
                ["claude", "-p", prompt, "--output-format", "text"],
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                logger.warning(f"Claude CLI returned error: {result.stderr}")
                return None

            response_text = result.stdout.strip()

            # Extract JSON from response (handle markdown code blocks)
            if "```json" in response_text:
                json_str = response_text.split("```json")[1].split("```")[0]
            elif "```" in response_text:
                json_str = response_text.split("```")[1].split("```")[0]
            else:
                # Try to find JSON object in response
                start_idx = response_text.find('{')
                end_idx = response_text.rfind('}')
                if start_idx != -1 and end_idx != -1:
                    json_str = response_text[start_idx:end_idx + 1]
                else:
                    json_str = response_text

            result = json.loads(json_str.strip())
            return result

        except subprocess.TimeoutExpired:
            logger.warning("Claude CLI call timed out")
            return None
        except json.JSONDecodeError as e:
            logger.warning(f"Failed to parse Claude CLI response as JSON: {e}")
            return None
        except Exception as e:
            logger.warning(f"Claude CLI call failed: {e}")
            return None

    def expand_terms(self, term: str, term_type: str = "disease") -> List[str]:
        """Expand a search term - kept for backward compatibility."""
        # With Claude, we don't need extensive synonym lists
        # Just return the original term
        return [term]

    def interpret_diseases(self, user_input: str) -> List[str]:
        """
        Interpret user disease input and expand to related terms.

        PRIMARY: Claude CLI-based semantic interpretation
        FALLBACK: Comma-split + fuzzy string match (when Claude unavailable)

        Args:
            user_input: User-provided disease string (may be comma-separated)

        Returns:
            List of expanded disease terms for OR logic query
        """
        # Check cache first
        cache_key = f"interpret_diseases:{hashlib.md5(user_input.encode()).hexdigest()}"
        if cache_key in self.cache:
            logger.debug(f"Using cached disease interpretation for: {user_input[:50]}...")
            return self.cache[cache_key]

        # Try Claude CLI first
        if self.use_claude:
            result = self._claude_interpret_diseases(user_input)
            if result:
                self.cache[cache_key] = result
                self._save_cache()
                return result

        # Fallback to fuzzy matching
        logger.info("Using fuzzy string matching for disease interpretation")
        return self._fuzzy_interpret_diseases(user_input)

    def _claude_interpret_diseases(self, user_input: str) -> Optional[List[str]]:
        """Use Claude CLI to semantically interpret disease input."""
        prompt = f"""Interpret this disease/indication input and expand to a list of related terms for a PubMed search.

USER INPUT: {user_input}

INSTRUCTIONS:
1. Parse any comma-separated diseases (these use OR logic)
2. For each disease, expand to include:
   - Official disease name
   - Common synonyms and abbreviations
   - Major subtypes if relevant
   - Related conditions if clearly implied
3. Correct obvious typos (e.g., "Crohns" → "Crohn's disease")
4. Do NOT over-expand (keep list focused on what user likely wants)

Examples:
- "IBD" → ["inflammatory bowel disease", "Crohn's disease", "ulcerative colitis", "colitis"]
- "Crohns, UC" → ["Crohn's disease", "ulcerative colitis"]
- "pulmonary fibrosis" → ["pulmonary fibrosis", "idiopathic pulmonary fibrosis", "IPF", "lung fibrosis"]
- "RA" → ["rheumatoid arthritis", "RA"]

Respond ONLY with a JSON array of disease terms (no other text):
["disease term 1", "disease term 2", ...]"""

        try:
            result = subprocess.run(
                ["claude", "-p", prompt, "--output-format", "text"],
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode != 0:
                logger.warning(f"Claude CLI returned error: {result.stderr}")
                return None

            response_text = result.stdout.strip()

            # Extract JSON array from response
            if "```json" in response_text:
                json_str = response_text.split("```json")[1].split("```")[0]
            elif "```" in response_text:
                json_str = response_text.split("```")[1].split("```")[0]
            else:
                # Try to find JSON array in response
                start_idx = response_text.find('[')
                end_idx = response_text.rfind(']')
                if start_idx != -1 and end_idx != -1:
                    json_str = response_text[start_idx:end_idx + 1]
                else:
                    json_str = response_text

            diseases = json.loads(json_str.strip())
            if isinstance(diseases, list) and len(diseases) > 0:
                logger.info(f"Claude interpreted diseases: {diseases}")
                return diseases
            return None

        except subprocess.TimeoutExpired:
            logger.warning("Claude CLI call timed out for disease interpretation")
            return None
        except json.JSONDecodeError as e:
            logger.warning(f"Failed to parse Claude disease interpretation as JSON: {e}")
            return None
        except Exception as e:
            logger.warning(f"Claude CLI disease interpretation failed: {e}")
            return None

    def _fuzzy_interpret_diseases(self, user_input: str) -> List[str]:
        """
        Fallback fuzzy interpretation when Claude is unavailable.

        Uses comma-split + common disease expansions.
        """
        # Known disease expansions for common abbreviations
        disease_expansions = {
            "ibd": ["inflammatory bowel disease", "Crohn's disease", "ulcerative colitis", "colitis"],
            "crohn": ["Crohn's disease"],
            "crohns": ["Crohn's disease"],
            "crohn's": ["Crohn's disease"],
            "uc": ["ulcerative colitis"],
            "ra": ["rheumatoid arthritis"],
            "ms": ["multiple sclerosis"],
            "sle": ["systemic lupus erythematosus", "lupus"],
            "lupus": ["systemic lupus erythematosus", "lupus"],
            "ipf": ["idiopathic pulmonary fibrosis", "pulmonary fibrosis"],
            "copd": ["chronic obstructive pulmonary disease", "COPD"],
            "nash": ["non-alcoholic steatohepatitis", "NASH", "fatty liver disease"],
            "nafld": ["non-alcoholic fatty liver disease", "NAFLD", "fatty liver disease"],
            "ad": ["Alzheimer's disease", "Alzheimer disease"],
            "pd": ["Parkinson's disease", "Parkinson disease"],
            "als": ["amyotrophic lateral sclerosis", "ALS", "Lou Gehrig's disease"],
            "t1d": ["type 1 diabetes", "T1D", "juvenile diabetes"],
            "t2d": ["type 2 diabetes", "T2D", "diabetes mellitus type 2"],
        }

        # Split by comma
        raw_terms = [t.strip() for t in user_input.split(",") if t.strip()]
        expanded_terms = []

        for term in raw_terms:
            term_lower = term.lower().strip()

            # Check for known expansion
            if term_lower in disease_expansions:
                expanded_terms.extend(disease_expansions[term_lower])
            else:
                # Try fuzzy match against known terms
                best_match = None
                best_ratio = 0.0
                try:
                    from rapidfuzz import fuzz
                    for key in disease_expansions.keys():
                        ratio = fuzz.ratio(term_lower, key) / 100.0
                        if ratio > best_ratio and ratio > 0.8:
                            best_ratio = ratio
                            best_match = key
                except ImportError:
                    # If rapidfuzz not available, use simple substring match
                    for key in disease_expansions.keys():
                        if term_lower in key or key in term_lower:
                            best_match = key
                            break

                if best_match:
                    expanded_terms.extend(disease_expansions[best_match])
                    logger.info(f"Fuzzy matched '{term}' to '{best_match}'")
                else:
                    # No match found, use original term
                    expanded_terms.append(term)

        # Remove duplicates while preserving order
        seen = set()
        unique_terms = []
        for t in expanded_terms:
            t_lower = t.lower()
            if t_lower not in seen:
                seen.add(t_lower)
                unique_terms.append(t)

        logger.info(f"Fuzzy interpreted diseases: {unique_terms}")
        return unique_terms

    def interpret_tissues(self, user_input: str) -> List[str]:
        """
        Interpret user tissue input and expand to related anatomical terms.

        PRIMARY: Claude CLI-based semantic interpretation
        FALLBACK: Predefined anatomical synonym mappings

        Args:
            user_input: User-provided tissue string (may be comma-separated)

        Returns:
            List of expanded tissue terms for OR logic query
        """
        # Check cache first
        cache_key = f"interpret_tissues:{hashlib.md5(user_input.encode()).hexdigest()}"
        if cache_key in self.cache:
            logger.debug(f"Using cached tissue interpretation for: {user_input[:50]}...")
            return self.cache[cache_key]

        # Try Claude CLI first
        if self.use_claude:
            result = self._claude_interpret_tissues(user_input)
            if result:
                self.cache[cache_key] = result
                self._save_cache()
                return result

        # Fallback to predefined mappings
        logger.info("Using predefined mappings for tissue interpretation")
        return self._fuzzy_interpret_tissues(user_input)

    def _claude_interpret_tissues(self, user_input: str) -> Optional[List[str]]:
        """Use Claude CLI to semantically interpret tissue input."""
        prompt = f"""Interpret this tissue/organ input and expand to a list of related anatomical terms for a PubMed search.

USER INPUT: {user_input}

INSTRUCTIONS:
1. Parse any comma-separated tissues (these use OR logic)
2. For each tissue, expand to include:
   - The official anatomical name
   - Common synonyms and abbreviations
   - Anatomical subregions that would be relevant
   - Related terms commonly used in scientific literature
3. Keep the list focused and clinically relevant
4. Do NOT over-expand to unrelated organs

Examples:
- "intestine" → ["intestine", "colon", "ileum", "jejunum", "duodenum", "cecum", "rectum", "gut", "bowel", "gastrointestinal", "enteric", "colonic", "intestinal mucosa"]
- "lung" → ["lung", "pulmonary", "bronchial", "alveolar", "airway", "respiratory"]
- "liver" → ["liver", "hepatic", "hepatocyte", "biliary"]
- "skin" → ["skin", "dermal", "epidermis", "cutaneous", "keratinocyte"]
- "brain" → ["brain", "cerebral", "cortex", "hippocampus", "neuronal", "CNS"]
- "colon" → ["colon", "colonic", "colorectal", "large intestine", "sigmoid", "cecum", "rectum"]

Respond ONLY with a JSON array of tissue terms (no other text):
["tissue term 1", "tissue term 2", ...]"""

        try:
            result = subprocess.run(
                ["claude", "-p", prompt, "--output-format", "text"],
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode != 0:
                logger.warning(f"Claude CLI returned error for tissue interpretation: {result.stderr}")
                return None

            response_text = result.stdout.strip()

            # Extract JSON array from response
            if "```json" in response_text:
                json_str = response_text.split("```json")[1].split("```")[0]
            elif "```" in response_text:
                json_str = response_text.split("```")[1].split("```")[0]
            else:
                # Try to find JSON array in response
                start_idx = response_text.find('[')
                end_idx = response_text.rfind(']')
                if start_idx != -1 and end_idx != -1:
                    json_str = response_text[start_idx:end_idx + 1]
                else:
                    json_str = response_text

            tissues = json.loads(json_str.strip())
            if isinstance(tissues, list) and len(tissues) > 0:
                logger.info(f"Claude interpreted tissues: {tissues}")
                return tissues
            return None

        except subprocess.TimeoutExpired:
            logger.warning("Claude CLI call timed out for tissue interpretation")
            return None
        except json.JSONDecodeError as e:
            logger.warning(f"Failed to parse Claude tissue interpretation as JSON: {e}")
            return None
        except Exception as e:
            logger.warning(f"Claude CLI tissue interpretation failed: {e}")
            return None

    def _fuzzy_interpret_tissues(self, user_input: str) -> List[str]:
        """
        Fallback tissue interpretation using predefined anatomical mappings.
        """
        # Predefined tissue expansions for common anatomical terms
        tissue_expansions = {
            "intestine": ["intestine", "intestinal", "colon", "colonic", "ileum", "ileal", "jejunum",
                         "duodenum", "cecum", "rectum", "rectal", "gut", "bowel", "gastrointestinal",
                         "GI tract", "enteric", "enterocyte", "mucosa"],
            "colon": ["colon", "colonic", "colorectal", "large intestine", "sigmoid", "cecum",
                     "rectum", "rectal", "bowel"],
            "small intestine": ["small intestine", "ileum", "ileal", "jejunum", "duodenum", "enterocyte"],
            "gut": ["gut", "intestine", "intestinal", "colon", "colonic", "gastrointestinal", "enteric"],
            "lung": ["lung", "pulmonary", "bronchial", "alveolar", "airway", "respiratory", "bronchus"],
            "liver": ["liver", "hepatic", "hepatocyte", "hepatocellular", "biliary"],
            "kidney": ["kidney", "renal", "nephron", "glomerular"],
            "skin": ["skin", "dermal", "dermis", "epidermis", "epidermal", "cutaneous", "keratinocyte"],
            "brain": ["brain", "cerebral", "cortex", "cortical", "hippocampus", "hippocampal",
                     "neuronal", "CNS", "central nervous system"],
            "heart": ["heart", "cardiac", "cardiomyocyte", "myocardial", "cardiovascular"],
            "muscle": ["muscle", "muscular", "skeletal muscle", "myocyte", "myofiber"],
            "adipose": ["adipose", "fat", "adipocyte", "white adipose", "brown adipose"],
            "pancreas": ["pancreas", "pancreatic", "islet", "beta cell"],
            "spleen": ["spleen", "splenic"],
            "thymus": ["thymus", "thymic"],
            "bone marrow": ["bone marrow", "marrow", "hematopoietic"],
            "lymph node": ["lymph node", "lymphoid", "lymphatic"],
            "blood": ["blood", "peripheral blood", "PBMC", "serum", "plasma"],
        }

        # Split by comma
        raw_terms = [t.strip().lower() for t in user_input.split(",") if t.strip()]
        expanded_terms = []

        for term in raw_terms:
            # Check for exact match
            if term in tissue_expansions:
                expanded_terms.extend(tissue_expansions[term])
            else:
                # Check for partial matches
                matched = False
                for key, expansions in tissue_expansions.items():
                    if term in key or key in term:
                        expanded_terms.extend(expansions)
                        matched = True
                        break

                if not matched:
                    # No match found, use original term
                    expanded_terms.append(term)

        # Remove duplicates while preserving order
        seen = set()
        unique_terms = []
        for t in expanded_terms:
            t_lower = t.lower()
            if t_lower not in seen:
                seen.add(t_lower)
                unique_terms.append(t)

        logger.info(f"Interpreted tissues: {unique_terms}")
        return unique_terms

    def score_publication(
        self,
        pub: Dict,
        disease: Union[str, List[str]],
        tissues: List[str],
        organisms: List[str],
        cell_type: str = None
    ) -> Dict:
        """
        Score a publication for relevance using Claude semantic understanding.

        Args:
            pub: Publication dict with 'title' and 'abstract'
            disease: Target disease/indication (string or list for OR logic)
            tissues: List of target tissues/organs (OR logic)
            organisms: List of target organisms
            cell_type: Optional target cell type

        Returns:
            Dict with match results and relevance score
        """
        title = pub.get('title', '')
        abstract = pub.get('abstract', '')
        text = f"{title}\n\n{abstract}"

        # Normalize tissues to list
        if isinstance(tissues, str):
            tissues = [tissues]

        # Normalize disease to list for OR logic
        if isinstance(disease, str):
            diseases = [disease]
        else:
            diseases = disease

        # Check cache first (use sorted tissues and diseases for consistent cache key)
        tissue_str = "|".join(sorted(tissues))
        disease_str = "|".join(sorted(diseases))
        cache_key = self._get_cache_key(text, disease_str, tissue_str, organisms or [], cell_type)
        if cache_key in self.cache:
            logger.debug(f"Using cached match result for {pub.get('pmid', 'unknown')}")
            return self.cache[cache_key]

        # Use Claude CLI if available
        if self.use_claude:
            # Join tissues and diseases for Claude prompt (OR logic)
            tissue_prompt = " or ".join(tissues)
            disease_prompt = " or ".join(diseases) if len(diseases) > 1 else diseases[0]
            claude_result = self._call_claude_cli(text, disease_prompt, tissue_prompt, organisms or [], cell_type)

            if claude_result:
                result = {
                    "matched": claude_result.get("disease_match", False) and claude_result.get("tissue_match", False),
                    "matched_disease": claude_result.get("disease_found") if claude_result.get("disease_match") else None,
                    "matched_tissue": claude_result.get("tissue_found") if claude_result.get("tissue_match") else None,
                    "matched_organism": claude_result.get("organism_found") if claude_result.get("organism_match") else None,
                    "matched_cell_type": claude_result.get("cell_type_found") if claude_result.get("cell_type_match") else None,
                    "relevance_score": claude_result.get("overall_relevance", 0.0),
                    "disease_reasoning": claude_result.get("disease_reasoning", ""),
                    "tissue_reasoning": claude_result.get("tissue_reasoning", ""),
                    "relevance_reasoning": claude_result.get("relevance_reasoning", ""),
                    "match_method": "claude"
                }

                # Cache the result
                self.cache[cache_key] = result
                self._save_cache()

                return result

        # Fallback to rule-based matching
        return self._rule_based_match(text, diseases, tissues, organisms, cell_type)

    def _rule_based_match(
        self,
        text: str,
        diseases: Union[str, List[str]],
        tissues: List[str],
        organisms: List[str],
        cell_type: str = None
    ) -> Dict:
        """Fallback rule-based matching when Claude is unavailable."""
        text_lower = text.lower()

        # Normalize tissues to list
        if isinstance(tissues, str):
            tissues = [tissues]

        # Normalize diseases to list for OR logic
        if isinstance(diseases, str):
            diseases = [diseases]

        # Disease match (OR logic - any disease matches)
        disease_match = False
        disease_found = None
        for disease in diseases:
            if disease.lower() in text_lower:
                disease_match = True
                disease_found = disease
                break
        disease_score = 1.0 if disease_match else 0.0

        # Tissue match (OR logic - any tissue matches)
        tissue_match = False
        tissue_found = None
        for tissue in tissues:
            if tissue.lower() in text_lower:
                tissue_match = True
                tissue_found = tissue
                break
        tissue_score = 1.0 if tissue_match else 0.0

        # Organism match
        org_match = False
        org_found = None
        org_score = 0.0
        for org in (organisms or []):
            org_lower = org.lower()
            pattern = self.organism_patterns.get(org_lower, rf"\b{org_lower}\b")
            if re.search(pattern, text_lower, re.IGNORECASE):
                org_match = True
                org_found = org
                org_score = 1.0
                break

        # Cell type match
        cell_match = True
        cell_found = None
        cell_score = 1.0
        if cell_type:
            cell_match = cell_type.lower() in text_lower
            cell_found = cell_type if cell_match else None
            cell_score = 1.0 if cell_match else 0.0

        # Calculate overall score
        weights = {"disease": 0.4, "tissue": 0.3, "organism": 0.2, "cell_type": 0.1}
        total_score = (
            disease_score * weights["disease"] +
            tissue_score * weights["tissue"] +
            org_score * weights["organism"] +
            cell_score * weights["cell_type"]
        )

        return {
            "matched": disease_match and tissue_match,
            "matched_disease": disease_found,
            "matched_tissue": tissue_found,
            "matched_organism": org_found,
            "matched_cell_type": cell_found,
            "relevance_score": total_score,
            "match_method": "rule_based"
        }
