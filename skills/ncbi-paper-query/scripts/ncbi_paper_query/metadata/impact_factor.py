"""
NCBI Paper Query - Impact Factor Lookup Module

Looks up journal impact factors from multiple sources including
reference files, cache, and web search (bioxbio.com, Scimago).
"""

import re
import json
import logging
from pathlib import Path
from typing import Dict, Optional
from functools import lru_cache

import requests
from bs4 import BeautifulSoup

from ..config import SKILL_DIR, logger


class ImpactFactorLookup:
    """Look up journal impact factors from multiple sources."""

    CACHE_FILE = ".if_cache.json"
    SCIMAGO_URL = "https://www.scimagojr.com/journalsearch.php"

    # Preprint servers (no IF)
    PREPRINT_SERVERS = {
        "biorxiv", "medrxiv", "arxiv", "chemrxiv", "preprints",
        "ssrn", "research square", "authorea"
    }

    # Impact factors are loaded from:
    # 1. journal_if.json (user-maintainable reference file)
    # 2. .if_cache.json (web search cache)
    # 3. Web search (bioxbio.com)
    IF_REFERENCE_FILE = "journal_if.json"

    def __init__(self, cache_dir: str = "."):
        self.cache_path = Path(cache_dir) / self.CACHE_FILE
        # Look for journal_if.json in SKILL_DIR/scripts first, then cache_dir
        skill_ref_path = SKILL_DIR / "scripts" / self.IF_REFERENCE_FILE
        local_ref_path = Path(cache_dir) / self.IF_REFERENCE_FILE

        # Also check for data directory within package
        package_ref_path = Path(__file__).parent.parent / "data" / self.IF_REFERENCE_FILE

        if skill_ref_path.exists():
            self.reference_path = skill_ref_path
        elif package_ref_path.exists():
            self.reference_path = package_ref_path
        else:
            self.reference_path = local_ref_path

        self.cache = self._load_cache()
        self.reference_if = self._load_reference_if()
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        })

    def _load_reference_if(self) -> Dict:
        """Load impact factors from reference file (user-maintainable)."""
        if self.reference_path.exists():
            try:
                with open(self.reference_path, 'r') as f:
                    data = json.load(f)
                    # Normalize keys to lowercase
                    return {k.lower().strip(): v for k, v in data.items()}
            except Exception as e:
                logger.warning(f"Failed to load IF reference file: {e}")
        return {}

    def _load_cache(self) -> Dict:
        """Load cached impact factors."""
        if self.cache_path.exists():
            try:
                with open(self.cache_path, 'r') as f:
                    return json.load(f)
            except:
                pass
        return {}

    def _save_cache(self):
        """Save impact factor cache."""
        try:
            with open(self.cache_path, 'w') as f:
                json.dump(self.cache, f)
        except Exception as e:
            logger.warning(f"Failed to save IF cache: {e}")

    def is_preprint(self, journal: str) -> bool:
        """Check if journal is a preprint server."""
        journal_lower = journal.lower()
        return any(pp in journal_lower for pp in self.PREPRINT_SERVERS)

    @lru_cache(maxsize=1000)
    def get_impact_factor(self, journal: str) -> Optional[float]:
        """Get impact factor for a journal.

        Lookup order:
        1. Reference file (journal_if.json) - user-maintainable
        2. Cache file (.if_cache.json) - web search results
        3. Web search (bioxbio.com)
        """
        if not journal:
            return None

        # Check if preprint
        if self.is_preprint(journal):
            return None

        cache_key = journal.lower().strip()

        # 1. Check reference file first (user-maintainable, most reliable)
        if cache_key in self.reference_if:
            ref_value = self.reference_if[cache_key]
            if ref_value is not None and self._validate_if_value(ref_value):
                return ref_value

        # 2. Check web search cache
        if cache_key in self.cache:
            cached_value = self.cache[cache_key]
            if cached_value is not None and self._validate_if_value(cached_value):
                return cached_value
            # Invalid cached value, remove it
            del self.cache[cache_key]
            self._save_cache()

        # 3. Try web search (bioxbio.com)
        try:
            if_value = self._fetch_via_web_search(journal)
            if if_value is not None and self._validate_if_value(if_value):
                self.cache[cache_key] = if_value
                self._save_cache()
                return if_value
        except Exception as e:
            logger.debug(f"Web search failed for {journal}: {e}")

        # 4. Try Scimago as fallback (often blocked)
        try:
            if_value = self._fetch_from_scimago(journal)
            if if_value is not None and self._validate_if_value(if_value):
                self.cache[cache_key] = if_value
                self._save_cache()
                return if_value
        except Exception as e:
            logger.debug(f"Scimago fetch failed for {journal}: {e}")

        return None

    def _validate_if_value(self, value: float) -> bool:
        """
        Validate that a value looks like a real impact factor.

        Rejects:
        - Values < 0.1 (too small to be a real IF)
        - Values > 200 (unrealistically high)
        - Values in range 1900-2100 (likely years, not IF values)
        """
        if value < 0.1 or value > 200:
            return False
        # Reject if it looks like a year (1900-2100)
        if 1900 <= value <= 2100:
            return False
        return True

    def _fetch_via_web_search(self, journal: str) -> Optional[float]:
        """Fetch impact factor via web search (bioxbio.com as reliable source)."""
        try:
            # Try multiple URL variations for bioxbio.com
            journal_variations = [
                journal.replace(' ', '-').replace('&', 'AND').upper(),
                journal.replace(' ', '-').replace('&', '-AND-').upper(),
                journal.replace(' ', '-').upper(),
                # Remove common prefixes/suffixes
                re.sub(r'\s*\(.*?\)\s*', '', journal).replace(' ', '-').upper(),
                # Abbreviation patterns
                ''.join(word[0] for word in journal.split() if word[0].isupper()).upper() if len(journal.split()) > 1 else journal.upper(),
            ]

            for journal_url in journal_variations:
                if not journal_url:
                    continue
                search_url = f"https://www.bioxbio.com/if/html/{journal_url}.html"
                try:
                    response = self.session.get(search_url, timeout=10)
                    if response.status_code != 200:
                        continue

                    soup = BeautifulSoup(response.text, 'html.parser')

                    # Verify this is the right journal (check title matches)
                    title = soup.title.string if soup.title else ""
                    if title and journal.lower()[:10] not in title.lower():
                        continue

                    # bioxbio.com uses tables with IF data
                    # Structure: Year | Impact Factor (IF) | Total Articles | Total Cites
                    tables = soup.find_all('table')
                    for table in tables:
                        rows = table.find_all('tr')
                        for row in rows[1:]:  # Skip header row
                            cells = row.find_all(['td', 'th'])
                            if len(cells) >= 2:
                                year_cell = cells[0].get_text(strip=True)
                                if_cell = cells[1].get_text(strip=True)
                                # Look for recent year IF values (2022-2025)
                                if any(str(y) in year_cell for y in [2022, 2023, 2024, 2025]):
                                    try:
                                        if_value = float(if_cell)
                                        if self._validate_if_value(if_value):
                                            return if_value
                                    except ValueError:
                                        continue

                    # Fallback: try regex on page text
                    page_text = soup.get_text()
                    if_patterns = [
                        r'(\d+\.?\d*)\s*(?:Impact\s*Factor|IF)',
                        r'Impact\s*Factor[:\s]+(\d+\.?\d*)',
                    ]
                    for pattern in if_patterns:
                        match = re.search(pattern, page_text, re.IGNORECASE)
                        if match:
                            value = float(match.group(1))
                            if self._validate_if_value(value):
                                return value

                except Exception:
                    continue

            return None

        except Exception as e:
            logger.debug(f"Web search IF fetch failed: {e}")
            return None

    def _fetch_from_scimago(self, journal: str) -> Optional[float]:
        """Fetch SJR score from Scimago (used as proxy for IF)."""
        try:
            params = {"q": journal}
            response = self.session.get(self.SCIMAGO_URL, params=params, timeout=10)

            if response.status_code != 200:
                return None

            soup = BeautifulSoup(response.text, 'html.parser')

            # Scimago search results page structure: look for journal entries
            # The search results contain journal cards with SJR values

            # Try to find the search result entries
            search_results = soup.find_all('a', class_='jrnlname')
            if not search_results:
                # Try alternative selectors
                search_results = soup.find_all('div', class_='search_results')

            # Look for SJR value in the page
            # Scimago displays SJR in various formats - try multiple patterns
            sjr_patterns = [
                r'SJR\s*[:\s]*(\d+\.?\d*)',
                r'(\d+\.?\d*)\s*SJR',
                r'sjr[:\s]*(\d+\.?\d*)',
            ]

            page_text = soup.get_text()
            for pattern in sjr_patterns:
                match = re.search(pattern, page_text, re.IGNORECASE)
                if match:
                    try:
                        sjr_value = float(match.group(1))
                        # SJR is typically 0-10, but IF can be higher
                        # Scale SJR to approximate IF (rough conversion)
                        if sjr_value > 0:
                            # SJR of ~3 often corresponds to IF of ~10
                            estimated_if = round(sjr_value * 3.5, 1)
                            # Validate the estimated IF value
                            if self._validate_if_value(estimated_if):
                                return estimated_if
                    except (ValueError, IndexError):
                        continue

            return None

        except Exception as e:
            logger.debug(f"Scimago fetch failed: {e}")
            return None
