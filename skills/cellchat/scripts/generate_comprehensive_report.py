#!/usr/bin/env python3
"""
CellChat Comprehensive Analysis Report Generator with ULTRATHINK Deep Biological Interpretation
================================================================================================

This script generates comprehensive markdown reports from CellChat analysis outputs with
embedded ULTRATHINK prompts for Claude Code to provide deep biological interpretations.

ULTRATHINK FRAMEWORK:
Claude Code must act as a professional biologist with expertise in molecular biology,
disease biology (context-specific), and computational biology when interpreting results.

For EACH section, apply:
1. MOLECULAR MECHANISM ANALYSIS - What interactions drive patterns?
2. BIOLOGICAL PROCESS IMPLICATIONS - How do changes relate to pathophysiology?
3. CELL TYPE FUNCTIONAL INTERPRETATION - What are cell roles in this context?
4. THERAPEUTIC/TRANSLATIONAL ASSESSMENT - What are intervention opportunities?
5. EXPERIMENTAL VALIDATION PRIORITIES - What requires further investigation?

WEB SEARCH INTEGRATION:
Use WebFetch to retrieve context from appropriate databases:
- Literature: PubMed, Google Scholar
- Gene/Protein: UniProt, GeneCards, NCBI Gene
- Pathways: KEGG, Reactome, STRING
- Disease: OMIM, DisGeNET, GWAS Catalog
- Drug/Target: DrugBank, ChEMBL, OpenTargets
- Cell/Tissue: Human Protein Atlas, GTEx, CellMarker

Usage:
    python generate_comprehensive_report.py -i <cellchat_results_dir> [-o <output_report.md>]

Input: Directory containing CellChat analysis outputs
Output: cellchat_comprehensive_report.md with ULTRATHINK interpretation prompts
"""

import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime
import json
import re


# =============================================================================
# ULTRATHINK BIOLOGICAL KNOWLEDGE DATABASE
# =============================================================================
# This knowledge base provides context-aware biological interpretations.
# Claude Code should EXTEND this with WebFetch for study-specific context.

PATHWAY_KNOWLEDGE = {
    'COLLAGEN': {
        'full_name': 'Collagen Signaling',
        'category': 'ECM',
        'description': 'Collagen signaling mediates ECM structure, mechanical signaling, and tissue stiffness through integrin and DDR receptors.',
        'general_role': 'Central to tissue architecture, wound healing, and fibrotic processes. Collagens provide structural scaffold and activate mechanotransduction pathways.',
        'key_genes': ['COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL6A1', 'ITGA1', 'ITGA2', 'ITGB1', 'DDR1', 'DDR2'],
        'therapeutic_targets': ['DDR1/2 inhibitors', 'Integrin antagonists', 'LOX/LOXL2 inhibitors'],
        'interpretation_up': 'Increased COLLAGEN signaling indicates active ECM remodeling. In disease contexts, may suggest tissue repair, fibrosis, or desmoplastic response.',
        'interpretation_down': 'Decreased COLLAGEN signaling may indicate reduced ECM production, tissue resolution, or matrix degradation phases.'
    },
    'LAMININ': {
        'full_name': 'Laminin Signaling',
        'category': 'ECM',
        'description': 'Laminin-integrin signaling maintains basement membrane integrity and epithelial-stromal communication.',
        'general_role': 'Essential for basement membrane structure, cell polarity, and epithelial organization. Disruption indicates tissue remodeling.',
        'key_genes': ['LAMA1', 'LAMA2', 'LAMA3', 'LAMB1', 'LAMB2', 'LAMC1', 'ITGA3', 'ITGA6', 'ITGB4'],
        'therapeutic_targets': ['Integrin α6β4 modulators'],
        'interpretation_up': 'Enhanced LAMININ signaling suggests active basement membrane remodeling and epithelial-stromal communication.',
        'interpretation_down': 'Reduced LAMININ signaling indicates basement membrane disruption or loss of epithelial organization.'
    },
    'TENASCIN': {
        'full_name': 'Tenascin Signaling',
        'category': 'ECM',
        'description': 'Tenascins are matricellular proteins upregulated in tissue injury, inflammation, and remodeling.',
        'general_role': 'Damage-associated matricellular protein that amplifies inflammation and modulates cell behavior. Normally low in adult tissue but induced during injury/disease.',
        'key_genes': ['TNC', 'TNXB', 'TNR', 'ITGAV', 'ITGB3', 'ITGB6', 'TLR4'],
        'therapeutic_targets': ['TNC-targeting antibodies', 'Integrin αvβ3/αvβ6 inhibitors'],
        'interpretation_up': 'Elevated TENASCIN signaling strongly indicates active tissue injury and remodeling. Context-dependent pro- or anti-inflammatory effects.',
        'interpretation_down': 'Reduced TENASCIN signaling suggests decreased injury response or transition to quiescent state.'
    },
    'FN1': {
        'full_name': 'Fibronectin Signaling',
        'category': 'ECM',
        'description': 'Fibronectin is a key ECM glycoprotein mediating cell adhesion, migration, and wound healing.',
        'general_role': 'Forms provisional matrix during tissue repair, serves as scaffold for cell migration and ECM assembly. EDA splice variant particularly important in disease.',
        'key_genes': ['FN1', 'ITGA5', 'ITGB1', 'ITGAV', 'ITGB3', 'ITGB5'],
        'therapeutic_targets': ['Integrin α5β1 inhibitors', 'FN-EDA splice variant targeting'],
        'interpretation_up': 'Increased FN1 signaling indicates active wound healing and ECM assembly. May suggest ongoing tissue remodeling.',
        'interpretation_down': 'Decreased FN1 signaling suggests reduced provisional matrix formation or resolved repair process.'
    },
    'TGFb': {
        'full_name': 'TGF-beta Signaling',
        'category': 'Secreted Signaling',
        'description': 'TGF-beta is a master regulator of cell differentiation, ECM production, and tissue homeostasis.',
        'general_role': 'Pleiotropic cytokine controlling cell proliferation, differentiation, apoptosis, and ECM synthesis. Context-dependent tumor suppressor or promoter.',
        'key_genes': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD7'],
        'therapeutic_targets': ['Anti-TGF-β antibodies', 'TGFBR1 inhibitors', 'Pirfenidone'],
        'interpretation_up': 'Increased TGF-β signaling indicates active cell differentiation programs. May drive ECM production, EMT, or immune modulation depending on context.',
        'interpretation_down': 'Reduced TGF-β signaling may indicate decreased differentiation drive or loss of homeostatic regulation.'
    },
    'PDGF': {
        'full_name': 'Platelet-Derived Growth Factor Signaling',
        'category': 'Secreted Signaling',
        'description': 'PDGF drives mesenchymal cell proliferation, migration, and activation.',
        'general_role': 'Potent mitogen for fibroblasts, pericytes, and smooth muscle cells. Critical for wound healing and vascular development.',
        'key_genes': ['PDGFA', 'PDGFB', 'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB'],
        'therapeutic_targets': ['Imatinib', 'Nintedanib'],
        'interpretation_up': 'Enhanced PDGF signaling indicates active mesenchymal cell proliferation and migration.',
        'interpretation_down': 'Reduced PDGF signaling suggests decreased mesenchymal activation or quiescent state.'
    },
    'WNT': {
        'full_name': 'WNT Signaling',
        'category': 'Secreted Signaling',
        'description': 'WNT signaling regulates tissue regeneration, stem cell maintenance, and developmental patterning.',
        'general_role': 'Normally quiescent in adult tissue but reactivated during regeneration/disease. Controls cell fate decisions and proliferation.',
        'key_genes': ['WNT1', 'WNT2', 'WNT3A', 'WNT5A', 'WNT7B', 'FZD1', 'FZD2', 'LRP5', 'LRP6', 'CTNNB1', 'SFRP1', 'SFRP2'],
        'therapeutic_targets': ['Porcupine inhibitors', 'β-catenin antagonists'],
        'interpretation_up': 'Increased WNT signaling indicates developmental pathway reactivation. May suggest regenerative response or aberrant proliferation.',
        'interpretation_down': 'Reduced WNT signaling suggests decreased regenerative/stem cell signaling. SFRP expression may actively suppress pathway.'
    },
    'BMP': {
        'full_name': 'Bone Morphogenetic Protein Signaling',
        'category': 'Secreted Signaling',
        'description': 'BMP signaling regulates cell differentiation and can counteract TGF-β effects.',
        'general_role': 'Important for tissue patterning and homeostasis. BMP7 has tissue-protective effects in many contexts. GREM1/GREM2 are BMP antagonists.',
        'key_genes': ['BMP2', 'BMP4', 'BMP5', 'BMP7', 'BMPR1A', 'BMPR1B', 'BMPR2', 'GREM1', 'GREM2', 'NOG'],
        'therapeutic_targets': ['Recombinant BMP7', 'BMP antagonist neutralization'],
        'interpretation_up': 'Increased BMP signaling may represent protective counter-regulation or morphogenetic activity.',
        'interpretation_down': 'Reduced BMP signaling, especially with BMP antagonist expression, removes protective signaling.'
    },
    'NOTCH': {
        'full_name': 'NOTCH Signaling',
        'category': 'Cell-Cell Contact',
        'description': 'Notch signaling controls cell fate decisions through direct cell-cell contact.',
        'general_role': 'Juxtacrine signaling pathway regulating differentiation, proliferation, and cell fate. Context-dependent effects on tissue homeostasis.',
        'key_genes': ['NOTCH1', 'NOTCH2', 'NOTCH3', 'JAG1', 'JAG2', 'DLL1', 'DLL4', 'HES1', 'HEY1'],
        'therapeutic_targets': ['Gamma-secretase inhibitors', 'NOTCH antibodies'],
        'interpretation_up': 'Enhanced NOTCH signaling indicates active contact-dependent cell fate decisions.',
        'interpretation_down': 'Reduced NOTCH signaling suggests decreased contact-mediated differentiation signals.'
    },
    'VEGF': {
        'full_name': 'Vascular Endothelial Growth Factor Signaling',
        'category': 'Secreted Signaling',
        'description': 'VEGF signaling drives angiogenesis and vascular permeability.',
        'general_role': 'Master regulator of blood vessel formation. Critical for tissue vascularization during development, wound healing, and disease.',
        'key_genes': ['VEGFA', 'VEGFB', 'VEGFC', 'KDR', 'FLT1', 'NRP1', 'NRP2'],
        'therapeutic_targets': ['Bevacizumab', 'VEGFR inhibitors'],
        'interpretation_up': 'Increased VEGF signaling suggests active angiogenesis and vascular remodeling.',
        'interpretation_down': 'Reduced VEGF signaling indicates decreased angiogenic drive.'
    },
    'CXCL': {
        'full_name': 'CXC Chemokine Signaling',
        'category': 'Secreted Signaling',
        'description': 'CXC chemokines recruit and position immune and stromal cells.',
        'general_role': 'Control immune cell trafficking and tissue positioning. Context-dependent pro- or anti-inflammatory effects.',
        'key_genes': ['CXCL1', 'CXCL2', 'CXCL8', 'CXCL12', 'CXCL14', 'CXCR1', 'CXCR2', 'CXCR4'],
        'therapeutic_targets': ['CXCR4 inhibitors', 'CXCL8 neutralization'],
        'interpretation_up': 'Enhanced CXCL signaling indicates active cell recruitment and chemotactic signaling.',
        'interpretation_down': 'Reduced CXCL signaling suggests decreased chemotactic activity.'
    },
    'CCL': {
        'full_name': 'CC Chemokine Signaling',
        'category': 'Secreted Signaling',
        'description': 'CC chemokines mediate monocyte/macrophage and lymphocyte recruitment.',
        'general_role': 'Major regulators of myeloid cell trafficking. CCL2/MCP-1 is the primary monocyte chemoattractant.',
        'key_genes': ['CCL2', 'CCL3', 'CCL4', 'CCL5', 'CCL7', 'CCR2', 'CCR5'],
        'therapeutic_targets': ['CCR2/5 inhibitors', 'CCL2 neutralization'],
        'interpretation_up': 'Increased CCL signaling indicates active monocyte/macrophage recruitment.',
        'interpretation_down': 'Decreased CCL signaling suggests reduced myeloid cell recruitment.'
    },
    'MIF': {
        'full_name': 'Macrophage Migration Inhibitory Factor',
        'category': 'Secreted Signaling',
        'description': 'MIF is a pro-inflammatory cytokine that inhibits macrophage migration and promotes inflammation.',
        'general_role': 'Pleiotropic cytokine with roles in inflammation, cell proliferation, and metabolism. Counter-regulates glucocorticoid effects.',
        'key_genes': ['MIF', 'CD74', 'CXCR2', 'CXCR4'],
        'therapeutic_targets': ['MIF inhibitors', 'CD74 targeting'],
        'interpretation_up': 'Increased MIF signaling indicates pro-inflammatory activation and may promote cell survival/proliferation.',
        'interpretation_down': 'Decreased MIF signaling suggests reduced inflammatory drive.'
    },
    'MHC-I': {
        'full_name': 'MHC Class I Signaling',
        'category': 'Cell-Cell Contact',
        'description': 'MHC-I presents intracellular peptides for immune surveillance.',
        'general_role': 'Essential for cytotoxic T cell recognition. Altered expression affects immune surveillance.',
        'key_genes': ['HLA-A', 'HLA-B', 'HLA-C', 'B2M'],
        'therapeutic_targets': ['Checkpoint modulators'],
        'interpretation_up': 'Enhanced MHC-I signaling suggests increased immune surveillance or interferon exposure.',
        'interpretation_down': 'Reduced MHC-I may indicate immune evasion mechanisms.'
    },
    'MHC-II': {
        'full_name': 'MHC Class II Signaling',
        'category': 'Cell-Cell Contact',
        'description': 'MHC-II presents extracellular antigens to CD4+ T cells.',
        'general_role': 'Professional antigen presentation. Aberrant expression on non-immune cells indicates inflammatory activation.',
        'key_genes': ['HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'CD74'],
        'therapeutic_targets': ['Immunomodulators'],
        'interpretation_up': 'Increased MHC-II indicates interferon-gamma exposure or antigen-presenting capacity.',
        'interpretation_down': 'Reduced MHC-II suggests decreased T cell-stromal crosstalk.'
    },
    'THBS': {
        'full_name': 'Thrombospondin Signaling',
        'category': 'ECM',
        'description': 'Thrombospondins are matricellular proteins that activate TGF-β and regulate angiogenesis.',
        'general_role': 'THBS1 is a major activator of latent TGF-β. Regulates cell-matrix interactions and angiogenesis.',
        'key_genes': ['THBS1', 'THBS2', 'THBS4', 'CD36', 'CD47'],
        'therapeutic_targets': ['TGF-β activation blockers'],
        'interpretation_up': 'Enhanced THBS signaling indicates active TGF-β activation and matricellular signaling.',
        'interpretation_down': 'Reduced THBS signaling suggests decreased TGF-β activation through this pathway.'
    },
    'SEMA': {
        'full_name': 'Semaphorin Signaling',
        'category': 'Secreted Signaling',
        'description': 'Semaphorins regulate cell guidance, angiogenesis, and immune function.',
        'general_role': 'Axon guidance molecules with broader roles in cell positioning, angiogenesis, and immune regulation.',
        'key_genes': ['SEMA3A', 'SEMA3C', 'SEMA4D', 'SEMA7A', 'NRP1', 'NRP2', 'PLXNA1'],
        'therapeutic_targets': ['Semaphorin-targeting antibodies'],
        'interpretation_up': 'Increased SEMA signaling indicates active cell guidance and positioning during tissue remodeling.',
        'interpretation_down': 'Reduced SEMA signaling suggests decreased guidance cue activity.'
    },
    'NRG': {
        'full_name': 'Neuregulin Signaling',
        'category': 'Secreted Signaling',
        'description': 'Neuregulins signal through ErbB receptors to regulate cell differentiation and survival.',
        'general_role': 'Growth factors with tissue-protective effects. Important in cardiac, neural, and epithelial biology.',
        'key_genes': ['NRG1', 'NRG2', 'NRG3', 'NRG4', 'ERBB2', 'ERBB3', 'ERBB4'],
        'therapeutic_targets': ['ErbB modulators'],
        'interpretation_up': 'Enhanced NRG signaling may indicate growth factor-mediated tissue protection or remodeling.',
        'interpretation_down': 'Reduced NRG signaling suggests decreased ErbB-mediated cellular communication.'
    },
    'NEGR': {
        'full_name': 'Neuronal Growth Regulator Signaling',
        'category': 'Cell-Cell Contact',
        'description': 'NEGR proteins mediate cell adhesion and neural-stromal communication.',
        'general_role': 'Cell adhesion molecules with roles in neural development and potentially stromal cell interactions.',
        'key_genes': ['NEGR1', 'LSAMP', 'OPCML', 'NTM'],
        'therapeutic_targets': ['Limited therapeutic development'],
        'interpretation_up': 'Enhanced NEGR signaling suggests active cell adhesion and positioning.',
        'interpretation_down': 'Reduced NEGR signaling indicates decreased adhesion pathway activity.'
    },
    'CDH': {
        'full_name': 'Cadherin Signaling',
        'category': 'Cell-Cell Contact',
        'description': 'Cadherins mediate cell-cell adhesion and tissue integrity.',
        'general_role': 'E-cadherin maintains epithelial integrity; N-cadherin expressed by mesenchymal cells. Cadherin switching indicates phenotypic transitions.',
        'key_genes': ['CDH1', 'CDH2', 'CDH3', 'CDH11', 'CTNNA1', 'CTNNB1'],
        'therapeutic_targets': ['Limited direct targeting'],
        'interpretation_up': 'Enhanced CDH signaling suggests active cell-cell adhesion. N-cadherin increase indicates mesenchymal phenotype.',
        'interpretation_down': 'Reduced CDH signaling, especially E-cadherin, suggests loss of epithelial integrity or EMT.'
    },
    'EGF': {
        'full_name': 'Epidermal Growth Factor Signaling',
        'category': 'Secreted Signaling',
        'description': 'EGF family signaling drives epithelial proliferation and differentiation.',
        'general_role': 'Critical for epithelial homeostasis, wound healing, and tissue regeneration.',
        'key_genes': ['EGF', 'EGFR', 'ERBB2', 'AREG', 'HBEGF', 'TGFa'],
        'therapeutic_targets': ['EGFR inhibitors', 'Anti-EGFR antibodies'],
        'interpretation_up': 'Increased EGF signaling indicates active epithelial proliferation or regeneration.',
        'interpretation_down': 'Reduced EGF signaling suggests decreased epithelial growth drive.'
    },
    'FGF': {
        'full_name': 'Fibroblast Growth Factor Signaling',
        'category': 'Secreted Signaling',
        'description': 'FGF signaling regulates cell proliferation, differentiation, and tissue patterning.',
        'general_role': 'Large family of growth factors with diverse roles in development, homeostasis, and disease.',
        'key_genes': ['FGF1', 'FGF2', 'FGF7', 'FGF10', 'FGFR1', 'FGFR2', 'FGFR3'],
        'therapeutic_targets': ['FGFR inhibitors'],
        'interpretation_up': 'Enhanced FGF signaling indicates active growth factor communication supporting proliferation or survival.',
        'interpretation_down': 'Reduced FGF signaling suggests decreased growth factor support.'
    }
}

CELL_TYPE_KNOWLEDGE = {
    'fibroblast': {
        'description': 'ECM-producing stromal cells that maintain tissue architecture',
        'markers': ['VIM', 'COL1A1', 'DCN', 'LUM', 'PDGFRA'],
        'function': 'Produce and maintain extracellular matrix, provide structural support, respond to tissue injury',
        'general_role': 'Primary effector cells for ECM production and tissue remodeling in diverse biological contexts'
    },
    'myofibroblast': {
        'description': 'Activated contractile fibroblasts with enhanced ECM production',
        'markers': ['ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'POSTN'],
        'function': 'Contract tissue, produce excessive ECM, mediate wound closure',
        'general_role': 'Key effector cells in tissue repair and pathological remodeling processes'
    },
    'perimyofib': {
        'description': 'Pericyte-myofibroblast transitional cells',
        'markers': ['PDGFRB', 'RGS5', 'ACTA2', 'PLN'],
        'function': 'Vascular support transitioning to contractile ECM-producing phenotype',
        'general_role': 'Represent pericyte-to-myofibroblast transition, important source of ECM-producing cells'
    },
    'telo': {
        'description': 'Telogen/quiescent fibroblasts with stem-like properties',
        'markers': ['CD34', 'PDGFRA', 'PI16', 'DPT'],
        'function': 'Maintain fibroblast pool, can differentiate upon activation',
        'general_role': 'Reservoir for fibroblast expansion; activation leads to tissue remodeling'
    },
    'tropho': {
        'description': 'Trophoblast-like activated stromal cells with high proliferative capacity',
        'markers': ['SFRP2', 'DPT', 'MFAP5', 'COMP'],
        'function': 'Active ECM remodeling and tissue organization',
        'general_role': 'Highly proliferative fibroblast subset active in tissue remodeling'
    },
    'cxcl14': {
        'description': 'Chemokine-expressing stromal cells involved in immune recruitment',
        'markers': ['CXCL14', 'CXCL12', 'CCL19'],
        'function': 'Recruit and position immune cells within tissue',
        'general_role': 'Create inflammatory/immune microenvironment through chemokine production'
    },
    'bmp5': {
        'description': 'BMP signaling-active cells with morphogenetic potential',
        'markers': ['BMP5', 'BMP4', 'GREM1', 'NOG'],
        'function': 'Regulate tissue patterning and potentially protective signaling',
        'general_role': 'May have protective roles through BMP effects; antagonist expression modulates pathway'
    },
    'sfrp2': {
        'description': 'WNT antagonist-expressing cells modulating developmental pathways',
        'markers': ['SFRP2', 'SFRP1', 'DKK1'],
        'function': 'Modulate WNT signaling, regulate tissue patterning',
        'general_role': 'May suppress WNT-driven processes; indicate developmental pathway activity'
    },
    'dpt': {
        'description': 'Dermatopontin-positive cells with ECM organizing function',
        'markers': ['DPT', 'LUM', 'DCN', 'FBLN1'],
        'function': 'Organize ECM architecture, regulate collagen fibrillogenesis',
        'general_role': 'ECM organization specialists; shape tissue architecture'
    },
    'pln': {
        'description': 'Phospholamban-positive cells with contractile properties',
        'markers': ['PLN', 'ACTA2', 'MYH11'],
        'function': 'Calcium handling and contractile function',
        'general_role': 'Contractile myofibroblast-like cells; contribute to tissue contraction'
    },
    'pdgfa': {
        'description': 'PDGF-producing cells driving mesenchymal cell activation',
        'markers': ['PDGFA', 'PDGFC', 'VEGFA'],
        'function': 'Produce growth factors that activate surrounding mesenchymal cells',
        'general_role': 'Paracrine activators of mesenchymal cell proliferation and migration'
    },
    'epithelial': {
        'description': 'Epithelial cells forming tissue barriers',
        'markers': ['EPCAM', 'KRT8', 'KRT18', 'CDH1'],
        'function': 'Form tissue barriers, absorb/secrete, protect underlying tissue',
        'general_role': 'Critical for tissue barrier function and epithelial-stromal crosstalk'
    },
    'endothelial': {
        'description': 'Vascular endothelial cells lining blood vessels',
        'markers': ['PECAM1', 'VWF', 'CDH5', 'ERG'],
        'function': 'Line blood vessels, regulate vascular permeability, support angiogenesis',
        'general_role': 'Critical for tissue vascularization and nutrient delivery'
    },
    'macrophage': {
        'description': 'Tissue-resident or recruited myeloid cells',
        'markers': ['CD68', 'CD163', 'CSF1R', 'ADGRE1'],
        'function': 'Phagocytosis, inflammation, tissue remodeling',
        'general_role': 'Key regulators of inflammation and tissue homeostasis'
    },
    'lymphocyte': {
        'description': 'Adaptive immune cells (T cells, B cells)',
        'markers': ['CD3D', 'CD4', 'CD8A', 'MS4A1'],
        'function': 'Adaptive immune responses, cytokine production',
        'general_role': 'Mediate adaptive immunity and tissue surveillance'
    }
}


# =============================================================================
# ULTRATHINK INTERPRETATION PROMPTS
# =============================================================================
# These prompts guide Claude Code to provide deep biological interpretations

ULTRATHINK_HEADER = """
<!-- ULTRATHINK INTERPRETATION REQUIRED -->
<!--
Claude Code: When presenting this report, apply ULTRATHINK deep biological reasoning:

1. MOLECULAR MECHANISM ANALYSIS
   - What specific molecular interactions drive these patterns?
   - Which ligands/receptors are responsible?
   - How do these interact with canonical pathways?

2. BIOLOGICAL PROCESS IMPLICATIONS
   - How do changes relate to the SPECIFIC pathophysiology being studied?
   - What stage of the process do patterns suggest?
   - Are patterns consistent with known biology?

3. CELL TYPE FUNCTIONAL INTERPRETATION
   - What is the biological identity of each cell type IN THIS CONTEXT?
   - Why would specific cell pairs show altered communication?
   - What does the autocrine/paracrine balance indicate?

4. THERAPEUTIC/TRANSLATIONAL OPPORTUNITY
   - Which pathways represent intervention points?
   - Are there existing therapeutics for identified pathways?
   - What would pathway modulation achieve?

5. EXPERIMENTAL VALIDATION PRIORITIES
   - Which findings are highest confidence?
   - What experiments would validate predictions?
   - What are key uncertainties?

WEB SEARCH: Use WebFetch to retrieve context from appropriate databases:
- Literature: PubMed, Google Scholar
- Gene/Protein: UniProt, GeneCards, NCBI Gene
- Pathways: KEGG, Reactome, STRING
- Disease: OMIM, DisGeNET, GWAS Catalog
- Drug/Target: DrugBank, ChEMBL, OpenTargets
- Cell/Tissue: Human Protein Atlas, GTEx, CellMarker
-->
"""


# =============================================================================
# INTERPRETATION FUNCTIONS WITH ULTRATHINK PROMPTS
# =============================================================================

def get_cell_type_info(cell_type_name):
    """Get biological information for a cell type."""
    cell_lower = cell_type_name.lower()
    for key, info in CELL_TYPE_KNOWLEDGE.items():
        if key in cell_lower:
            return info
    return {
        'description': 'stromal cell population',
        'markers': [],
        'function': 'tissue support and ECM maintenance',
        'general_role': 'potential contributor to tissue remodeling'
    }


def get_pathway_info(pathway_name):
    """Get biological information for a signaling pathway."""
    pathway_upper = pathway_name.upper()
    for key, info in PATHWAY_KNOWLEDGE.items():
        if key in pathway_upper:
            return info
    return {
        'full_name': pathway_name,
        'category': 'Signaling',
        'description': f'{pathway_name} signaling pathway',
        'general_role': 'Role requires further characterization based on study context',
        'key_genes': [],
        'therapeutic_targets': [],
        'interpretation_up': f'Increased {pathway_name} signaling indicates enhanced activity of this communication route.',
        'interpretation_down': f'Decreased {pathway_name} signaling indicates reduced pathway activity.'
    }


def interpret_heatmap(heatmap_df, condition, cell_types):
    """
    Generate ULTRATHINK biological interpretation for heatmap data.

    ULTRATHINK Framework Applied:
    - Analyze autocrine vs paracrine balance (tissue state indicator)
    - Identify dominant signaling axes (mechanistic insights)
    - Relate patterns to biological processes
    """
    cond_data = heatmap_df[heatmap_df['Condition'] == condition]
    if len(cond_data) == 0:
        return "No data available for interpretation."

    total_count = cond_data['Count'].sum()
    total_weight = cond_data['Weight'].sum()

    # Calculate autocrine vs paracrine
    autocrine = cond_data[cond_data['Sender'] == cond_data['Receiver']]
    autocrine_count = autocrine['Count'].sum()
    autocrine_pct = (autocrine_count / total_count * 100) if total_count > 0 else 0
    paracrine_pct = 100 - autocrine_pct

    # Top communicating pairs
    top_count = cond_data.nlargest(3, 'Count')[['Sender', 'Receiver', 'Count']].values.tolist()
    top_weight = cond_data.nlargest(3, 'Weight')[['Sender', 'Receiver', 'Weight']].values.tolist()

    interp = []

    # ULTRATHINK: Communication pattern interpretation with biological reasoning
    if autocrine_pct > 40:
        interp.append(f"**Communication Pattern:** High autocrine signaling ({autocrine_pct:.1f}%) indicates "
                      "cell-autonomous regulatory programs predominate. This pattern suggests intrinsic "
                      "cellular maintenance, stress response programs, or self-reinforcing activation loops "
                      "rather than coordinated microenvironment-driven crosstalk.")
    elif autocrine_pct > 25:
        interp.append(f"**Communication Pattern:** Balanced autocrine ({autocrine_pct:.1f}%) and paracrine signaling "
                      f"({paracrine_pct:.1f}%) indicates both cell-intrinsic and intercellular communication "
                      "are active, consistent with a dynamic tissue microenvironment.")
    else:
        interp.append(f"**Communication Pattern:** Predominantly paracrine signaling ({paracrine_pct:.1f}%) "
                      "indicates active intercellular crosstalk, suggesting coordinated multicellular responses "
                      "typical of tissue remodeling or inflammatory processes.")

    # ULTRATHINK: Dominant signaling axis with mechanistic interpretation
    if len(top_count) > 0:
        top_sender, top_receiver, top_n = top_count[0]
        sender_info = get_cell_type_info(top_sender)
        receiver_info = get_cell_type_info(top_receiver)

        if top_sender == top_receiver:
            interp.append(f"**Dominant Signaling Axis:** {top_sender} autocrine loop dominates, suggesting "
                          f"self-reinforcing signaling in {sender_info['description']}. This may indicate active "
                          f"proliferation, differentiation maintenance, or stress response programs.")
        else:
            interp.append(f"**Dominant Signaling Axis:** {top_sender} → {top_receiver} axis is most active. "
                          f"This paracrine route from {sender_info['description']} to {receiver_info['description']} "
                          f"may be critical for tissue homeostasis or pathological progression.")

    # Quantitative summary
    interp.append(f"\n**Quantitative Summary:**")
    interp.append(f"- Top pairs by count: " +
                  ", ".join([f"{s}→{r} ({int(c)})" for s, r, c in top_count]))
    interp.append(f"- Top pairs by strength: " +
                  ", ".join([f"{s}→{r} ({w:.2f})" for s, r, w in top_weight]))

    return "\n".join(interp)


def interpret_signaling_roles(roles_df, condition):
    """
    Generate ULTRATHINK biological interpretation for signaling roles.

    ULTRATHINK Framework Applied:
    - Identify communication hubs (potential therapeutic targets)
    - Classify sender/receiver roles (mechanistic architecture)
    - Relate to tissue biology
    """
    cond_data = roles_df[roles_df['Condition'] == condition].copy()
    if len(cond_data) == 0:
        return "No data available for interpretation."

    cond_data['Total_Strength'] = cond_data['Outgoing_Strength'] + cond_data['Incoming_Strength']
    cond_data['Out_In_Ratio'] = cond_data['Outgoing_Strength'] / (cond_data['Incoming_Strength'] + 0.001)

    # Classify cell types
    dominant_senders = cond_data.nlargest(2, 'Out_In_Ratio')['CellType'].tolist()
    dominant_receivers = cond_data.nsmallest(2, 'Out_In_Ratio')['CellType'].tolist()
    hub_cells = cond_data.nlargest(2, 'Total_Strength')['CellType'].tolist()

    interp = []

    # ULTRATHINK: Sender interpretation with biological context
    sender_infos = [get_cell_type_info(s) for s in dominant_senders]
    interp.append(f"**Dominant Senders:** {', '.join(dominant_senders)}")
    interp.append(f"  - These cell types act as primary signal initiators, expressing high levels of ligands "
                  f"that influence surrounding cells. As {sender_infos[0]['description']}, they likely orchestrate "
                  "tissue-level responses through secreted factors and ECM components.")

    # ULTRATHINK: Receiver interpretation
    receiver_infos = [get_cell_type_info(r) for r in dominant_receivers]
    interp.append(f"**Dominant Receivers:** {', '.join(dominant_receivers)}")
    interp.append(f"  - These populations are primary responders to microenvironment signals. High receptor "
                  f"expression in {receiver_infos[0]['description']} suggests they are targets for "
                  "therapeutic intervention or biomarker discovery.")

    # ULTRATHINK: Hub interpretation
    interp.append(f"**Communication Hubs:** {', '.join(hub_cells)}")
    interp.append(f"  - Hub cells exhibit both high sending and receiving capacity, acting as integrators "
                  "of tissue signaling. These populations may serve as critical nodes for amplifying "
                  "or dampening intercellular communication cascades.")

    # Biological implication
    top_sender = dominant_senders[0] if dominant_senders else None
    top_receiver = dominant_receivers[0] if dominant_receivers else None
    top_hub = hub_cells[0] if hub_cells else None

    if top_sender and top_receiver and top_hub:
        interp.append(f"\n**Biological Implication:** The signaling architecture suggests a hierarchical "
                      f"communication network where {top_sender} broadcast signals to responsive populations "
                      f"like {top_receiver}, while {top_hub} mediate signal integration and relay.")

    return "\n".join(interp)


def interpret_differential(diff_df, comparison):
    """
    Generate ULTRATHINK biological interpretation for differential interactions.

    ULTRATHINK Framework Applied:
    - Analyze direction and magnitude of changes
    - Identify key altered communication axes
    - Relate to disease/biological process progression
    """
    comp_data = diff_df[diff_df['Comparison_Pair'] == comparison]
    if len(comp_data) == 0:
        return "No data available for interpretation."

    parts = comparison.split('_vs_')
    test_cond = parts[0] if len(parts) > 0 else 'test'
    ref_cond = parts[1] if len(parts) > 1 else 'reference'

    up_pairs = comp_data[comp_data['Count_Direction'] == 'UP']
    dn_pairs = comp_data[comp_data['Count_Direction'] == 'DOWN']

    top_up = up_pairs.nlargest(3, 'Diff_Count')[['Sender', 'Receiver', 'Diff_Count']].values.tolist()
    top_dn = dn_pairs.nsmallest(3, 'Diff_Count')[['Sender', 'Receiver', 'Diff_Count']].values.tolist()

    net_change = len(up_pairs) - len(dn_pairs)

    interp = []

    # ULTRATHINK: Overall trend with biological reasoning
    if net_change > 0:
        interp.append(f"**Overall Trend:** Net INCREASE in cell-cell communication (+{net_change} pairs) "
                      f"in {test_cond} compared to {ref_cond}.")
        interp.append(f"  - This increased signaling is consistent with enhanced cellular activity, where "
                      "stromal cells upregulate communication to coordinate tissue responses, "
                      "ECM remodeling, and cellular differentiation.")
    elif net_change < 0:
        interp.append(f"**Overall Trend:** Net DECREASE in communication ({net_change} pairs).")
        interp.append(f"  - Reduced signaling may indicate loss of homeostatic communication, cellular exhaustion, "
                      "or transition to a quiescent tissue state.")
    else:
        interp.append(f"**Overall Trend:** Balanced changes with communication rewiring rather than global shift.")

    # ULTRATHINK: Increased interactions with mechanistic interpretation
    if len(top_up) > 0:
        interp.append(f"\n**Most Increased Interactions:**")
        for sender, receiver, diff in top_up:
            sender_info = get_cell_type_info(sender)
            receiver_info = get_cell_type_info(receiver)
            interp.append(f"  - {sender} → {receiver} (+{int(diff)}): Enhanced signaling "
                          f"from {sender_info['description']} to {receiver_info['description']} "
                          f"suggests activation of this axis in {test_cond}.")

    # ULTRATHINK: Decreased interactions
    if len(top_dn) > 0:
        interp.append(f"\n**Most Decreased Interactions:**")
        for sender, receiver, diff in top_dn:
            sender_info = get_cell_type_info(sender)
            receiver_info = get_cell_type_info(receiver)
            interp.append(f"  - {sender} → {receiver} ({int(diff)}): Reduced communication may indicate "
                          f"downregulation of {sender_info['description']} signals or decreased "
                          f"responsiveness in {receiver_info['description']}.")

    interp.append(f"\n**Quantitative Summary:** {len(up_pairs)} pairs increased, {len(dn_pairs)} pairs decreased")

    return "\n".join(interp)


def interpret_info_flow(info_flow_df, conditions):
    """
    Generate ULTRATHINK biological interpretation for information flow.

    ULTRATHINK Framework Applied:
    - Identify dominant pathways (mechanistic drivers)
    - Analyze pathway variability (condition-specific biology)
    - Connect to known pathway functions
    """
    if len(info_flow_df) == 0:
        return "No data available for interpretation."

    # Top pathways overall
    pathway_totals = info_flow_df.groupby('Pathway')['Info_Flow'].sum().nlargest(5)

    # Most variable pathways
    pathway_cv = info_flow_df.groupby('Pathway')['Info_Flow'].std() / (info_flow_df.groupby('Pathway')['Info_Flow'].mean() + 0.001)
    most_variable = pathway_cv.nlargest(3).index.tolist()

    interp = []
    interp.append("**Most Active Pathways:**")

    for pathway, flow in pathway_totals.items():
        pw_info = get_pathway_info(pathway)
        interp.append(f"  - **{pathway}** ({pw_info['category']}): Total flow = {flow:.2f}. "
                      f"Associated with {pw_info['general_role'][:80]}...")

    interp.append(f"\n**Most Variable Between Conditions:** {', '.join(most_variable)}")
    interp.append("  - High variability indicates condition-specific regulation of these pathways, "
                  "suggesting they may be key drivers of phenotypic differences.")

    # Top pathway biological context
    top_pathway = pathway_totals.index[0] if len(pathway_totals) > 0 else None
    if top_pathway:
        pw_info = get_pathway_info(top_pathway)
        interp.append(f"\n**Biological Implication:** The dominance of {top_pathway} signaling suggests "
                      f"that {pw_info['general_role'][:100]} are central to the tissue "
                      "microenvironment communication network.")

    return "\n".join(interp)


def interpret_pathway_stats(pathway_stats, comparison):
    """
    Generate ULTRATHINK biological interpretation for pathway statistics.

    ULTRATHINK Framework Applied:
    - Analyze statistical significance (FDR vs nominal)
    - Interpret pathway changes in biological context
    - Identify therapeutic opportunities
    """
    comp_data = pathway_stats[pathway_stats['Comparison_Pair'] == comparison]
    if len(comp_data) == 0:
        return "No data available for interpretation."

    parts = comparison.split('_vs_')
    test_cond = parts[0] if len(parts) > 0 else 'test'
    ref_cond = parts[1] if len(parts) > 1 else 'reference'

    sig_up = comp_data[(comp_data['Direction'] == 'UP') & (comp_data['Significant'] == True)]
    sig_dn = comp_data[(comp_data['Direction'] == 'DOWN') & (comp_data['Significant'] == True)]
    pval_up = comp_data[(comp_data['Direction'] == 'UP') & (comp_data['PValue'] < 0.05)]
    pval_dn = comp_data[(comp_data['Direction'] == 'DOWN') & (comp_data['PValue'] < 0.05)]

    interp = []

    # ULTRATHINK: Significant upregulated pathways
    if len(sig_up) > 0:
        interp.append(f"**Significantly Upregulated Pathways (FDR<0.05) in {test_cond}:**")
        for _, row in sig_up.head(3).iterrows():
            pw_info = get_pathway_info(row['Pathway'])
            interp.append(f"  - **{row['Pathway']}** (log2FC={row['Log2FC']:.2f}, FDR={row['FDR']:.4f})")
            interp.append(f"    {pw_info['interpretation_up']}")
    else:
        interp.append(f"**Significantly Upregulated (FDR<0.05):** None")
        if len(pval_up) > 0:
            interp.append(f"  - {len(pval_up)} pathways nominally significant (p<0.05) but not FDR-corrected")

    # ULTRATHINK: Significant downregulated pathways
    if len(sig_dn) > 0:
        interp.append(f"\n**Significantly Downregulated Pathways (FDR<0.05) in {test_cond}:**")
        for _, row in sig_dn.head(3).iterrows():
            pw_info = get_pathway_info(row['Pathway'])
            interp.append(f"  - **{row['Pathway']}** (log2FC={row['Log2FC']:.2f}, FDR={row['FDR']:.4f})")
            interp.append(f"    {pw_info['interpretation_down']}")
    else:
        interp.append(f"\n**Significantly Downregulated (FDR<0.05):** None")
        if len(pval_dn) > 0:
            interp.append(f"  - {len(pval_dn)} pathways nominally significant (p<0.05) but not FDR-corrected")

    interp.append(f"\n**Statistical Summary:**")
    interp.append(f"  - FDR<0.05: {len(sig_up)} up, {len(sig_dn)} down")
    interp.append(f"  - p<0.05: {len(pval_up)} up, {len(pval_dn)} down")

    if len(sig_up) + len(sig_dn) == 0:
        interp.append(f"\n**Biological Interpretation:**")
        interp.append(f"  - Absence of statistically significant pathway changes (FDR<0.05) suggests that "
                      "differences between conditions are subtle or distributed across many pathways. "
                      "Consider examining nominally significant pathways (p<0.05) for hypothesis generation.")

    return "\n".join(interp)


def generate_comparison_interpretation(comparison, pairs_up, pairs_dn, pathway_stats_df):
    """
    Generate detailed ULTRATHINK biological interpretation for a comparison.

    Lists specific pathways UP/DOWN with hierarchy:
    1. FDR < 0.05 (high confidence)
    2. p < 0.05 (nominally significant)
    3. Caution note if nothing significant

    Uses PATHWAY_KNOWLEDGE to provide biological context for each pathway.
    """
    parts = comparison.split('_vs_')
    test_cond = parts[0] if len(parts) > 0 else 'test'
    ref_cond = parts[1] if len(parts) > 1 else 'reference'

    interp = []

    # 1. Overall communication trend
    net_change = pairs_up - pairs_dn
    if net_change > 0:
        interp.append(f"**Overall Communication:** Net INCREASE (+{net_change} pairs) in {test_cond}.")
    elif net_change < 0:
        interp.append(f"**Overall Communication:** Net DECREASE ({net_change} pairs) in {test_cond}.")
    else:
        interp.append(f"**Overall Communication:** Balanced rewiring (no net change) in {test_cond}.")

    # 2. FDR-significant pathways (highest confidence)
    fdr_sig = pathway_stats_df[pathway_stats_df['Significant'] == True]
    fdr_up = fdr_sig[fdr_sig['Direction'] == 'UP'].nlargest(3, 'Log2FC')
    fdr_dn = fdr_sig[fdr_sig['Direction'] == 'DOWN'].nsmallest(3, 'Log2FC')

    if len(fdr_up) > 0 or len(fdr_dn) > 0:
        interp.append("")
        interp.append("**Significantly Altered Pathways (FDR<0.05) - HIGH CONFIDENCE:**")
        for _, row in fdr_up.iterrows():
            pw_info = get_pathway_info(row['Pathway'])
            interp_text = pw_info.get('interpretation_up', pw_info.get('description', 'Signaling pathway'))
            interp.append(f"  - ↑ **{row['Pathway']}** (log2FC={row['Log2FC']:.2f}, FDR={row['FDR']:.4f})")
            interp.append(f"    {interp_text[:120]}...")
        for _, row in fdr_dn.iterrows():
            pw_info = get_pathway_info(row['Pathway'])
            interp_text = pw_info.get('interpretation_down', pw_info.get('description', 'Signaling pathway'))
            interp.append(f"  - ↓ **{row['Pathway']}** (log2FC={row['Log2FC']:.2f}, FDR={row['FDR']:.4f})")
            interp.append(f"    {interp_text[:120]}...")
    else:
        # 3. p<0.05 pathways (nominally significant)
        pval_sig = pathway_stats_df[pathway_stats_df['PValue'] < 0.05]
        pval_up = pval_sig[pval_sig['Direction'] == 'UP'].nlargest(3, 'Log2FC')
        pval_dn = pval_sig[pval_sig['Direction'] == 'DOWN'].nsmallest(3, 'Log2FC')

        if len(pval_up) > 0 or len(pval_dn) > 0:
            interp.append("")
            interp.append("**Nominally Significant Pathways (p<0.05, NOT FDR-corrected):**")
            interp.append("*⚠️ Interpret with caution - these require validation*")
            for _, row in pval_up.iterrows():
                pw_info = get_pathway_info(row['Pathway'])
                desc = pw_info.get('description', f'{row["Pathway"]} signaling pathway')
                interp.append(f"  - ↑ **{row['Pathway']}** (log2FC={row['Log2FC']:.2f}, p={row['PValue']:.4f})")
                interp.append(f"    {desc[:100]}")
            for _, row in pval_dn.iterrows():
                pw_info = get_pathway_info(row['Pathway'])
                desc = pw_info.get('description', f'{row["Pathway"]} signaling pathway')
                interp.append(f"  - ↓ **{row['Pathway']}** (log2FC={row['Log2FC']:.2f}, p={row['PValue']:.4f})")
                interp.append(f"    {desc[:100]}")
        else:
            # 4. Nothing significant
            interp.append("")
            interp.append("**⚠️ CAUTION: No statistically significant pathway changes detected.**")
            interp.append("Changes between conditions are subtle or distributed across many pathways.")
            interp.append("Consider examining individual L-R pairs rather than pathway-level analysis.")

    return "\n".join(interp)


def generate_discussion(all_data, conditions, comparisons):
    """
    Generate ULTRATHINK comprehensive data-driven discussion.

    This function synthesizes all findings into a cohesive biological narrative
    following the ULTRATHINK framework.
    """
    interaction_summary = all_data.get('interaction_summary')
    diff_data = all_data.get('diff_data')
    pathway_stats = all_data.get('pathway_stats')
    heatmap_df = all_data.get('heatmap_df')
    roles_df = all_data.get('roles_df')
    info_flow_df = all_data.get('info_flow_df')

    discussion = []

    # Finding 1: Communication Dynamics
    discussion.append("### 6.1 Summary of Key Findings\n")
    discussion.append("#### Finding 1: Cell-Cell Communication Dynamics\n")

    if interaction_summary is not None and len(conditions) > 1:
        ref_cond = conditions[0]
        ref_interactions = interaction_summary[interaction_summary['Condition'] == ref_cond]['Total_Interactions'].values[0]

        for cond in conditions[1:]:
            cond_interactions = interaction_summary[interaction_summary['Condition'] == cond]['Total_Interactions'].values[0]
            pct_change = ((cond_interactions - ref_interactions) / ref_interactions * 100) if ref_interactions > 0 else 0

            if abs(pct_change) > 10:
                direction = "increase" if pct_change > 0 else "decrease"
                discussion.append(f"- **{cond}** shows {abs(pct_change):.1f}% {direction} in total cell-cell interactions compared to {ref_cond}.")
                if pct_change > 0:
                    discussion.append(f"  This enhanced communication is consistent with active tissue remodeling, "
                                      "where stromal populations coordinate responses through increased ligand-receptor signaling.")
            else:
                discussion.append(f"- **{cond}** maintains similar communication levels (±10.0%) as {ref_cond}.")

    # Finding 2: Cell Type Roles
    discussion.append("\n#### Finding 2: Cell Type Communication Roles\n")

    if roles_df is not None:
        for cond in conditions:
            cond_data = roles_df[roles_df['Condition'] == cond].copy()
            if len(cond_data) > 0:
                cond_data['Total_Strength'] = cond_data['Outgoing_Strength'] + cond_data['Incoming_Strength']
                top_hub = cond_data.nlargest(1, 'Total_Strength')['CellType'].values[0]
                hub_info = get_cell_type_info(top_hub)
                discussion.append(f"- In **{cond}**, {top_hub} serves as the primary communication hub ({hub_info['description']}).")

        discussion.append("\nThe identification of communication hubs suggests these cell populations are critical "
                          "integration points for tissue-level signaling and may represent key intervention targets.")

    # Finding 3: Pathway Alterations
    discussion.append("\n#### Finding 3: Signaling Pathway Alterations\n")

    if pathway_stats is not None:
        for comp in comparisons:
            pw_comp = pathway_stats[pathway_stats['Comparison_Pair'] == comp]
            sig_pathways = pw_comp[pw_comp['Significant'] == True]
            pval_sig = pw_comp[pw_comp['PValue'] < 0.05]

            if len(sig_pathways) > 0:
                top_changed = sig_pathways.nlargest(3, 'Log2FC', key=abs)['Pathway'].tolist()
                discussion.append(f"- **{comp}**: Significantly altered pathways include {', '.join(top_changed)}.")
                for pw in top_changed[:2]:
                    pw_info = get_pathway_info(pw)
                    discussion.append(f"  - {pw}: {pw_info['general_role'][:80]}...")
            elif len(pval_sig) > 0:
                discussion.append(f"- **{comp}**: No FDR-significant pathways, but {len(pval_sig)} pathways show "
                                  "nominal significance (p<0.05). These may represent subtle but biologically "
                                  "meaningful changes that warrant further investigation.")

    # Finding 4: Network Architecture
    discussion.append("\n#### Finding 4: Communication Network Architecture\n")

    if heatmap_df is not None:
        for cond in conditions:
            cond_data = heatmap_df[heatmap_df['Condition'] == cond]
            if len(cond_data) > 0:
                total_count = cond_data['Count'].sum()
                autocrine = cond_data[cond_data['Sender'] == cond_data['Receiver']]['Count'].sum()
                autocrine_pct = (autocrine / total_count * 100) if total_count > 0 else 0
                paracrine_pct = 100 - autocrine_pct

                if paracrine_pct > 60:
                    discussion.append(f"- **{cond}** exhibits predominantly paracrine signaling ({paracrine_pct:.1f}%), "
                                      "reflecting active intercellular communication networks.")
                else:
                    discussion.append(f"- **{cond}** shows high autocrine signaling ({autocrine_pct:.1f}%), "
                                      "suggesting self-reinforcing cellular programs.")

    # Section 6.2: Proposed Next Steps
    discussion.append("\n### 6.2 Proposed Next Steps\n")

    discussion.append("#### Differential Expression Analysis")
    if roles_df is not None:
        all_hubs = []
        for cond in conditions:
            cond_data = roles_df[roles_df['Condition'] == cond].copy()
            if len(cond_data) > 0:
                cond_data['Total_Strength'] = cond_data['Outgoing_Strength'] + cond_data['Incoming_Strength']
                hub = cond_data.nlargest(1, 'Total_Strength')['CellType'].values[0]
                if hub not in all_hubs:
                    all_hubs.append(hub)
        if all_hubs:
            discussion.append(f"- Prioritize DEG analysis on communication hub populations: {', '.join(all_hubs)}")
    discussion.append("- Validate ligand-receptor expression changes at transcript level\n")

    discussion.append("#### Experimental Validation Priorities")
    if heatmap_df is not None:
        top_pair = heatmap_df.nlargest(1, 'Count')[['Sender', 'Receiver']].values
        if len(top_pair) > 0:
            discussion.append(f"- Co-culture model: {top_pair[0][0]} + {top_pair[0][1]} to validate interaction changes")
    discussion.append("- Consider spatial transcriptomics to validate signaling patterns in tissue context\n")

    discussion.append("#### Clinical Translation")
    discussion.append("- Identified communication hubs and altered pathways may serve as:")
    discussion.append("  - Diagnostic biomarkers for disease staging")
    discussion.append("  - Therapeutic targets for modulating tissue microenvironment")
    discussion.append("  - Prognostic indicators based on signaling network architecture")

    return "\n".join(discussion)


# =============================================================================
# MAIN REPORT GENERATION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Generate CellChat comprehensive report with ULTRATHINK biological interpretations')
    parser.add_argument('-i', '--input', default='cellchat_results',
                        help='CellChat results directory')
    parser.add_argument('-o', '--output', default=None,
                        help='Output report path [default: <input>/cellchat_comprehensive_report.md]')
    args = parser.parse_args()

    results_dir = Path(args.input)
    output_file = Path(args.output) if args.output else results_dir / 'cellchat_comprehensive_report.md'

    print(f"Reading data from: {results_dir}")
    print("Generating report with ULTRATHINK deep biological interpretations...")

    # Read all data tables
    try:
        data_summary = pd.read_csv(results_dir / '01_data_summary.csv')
        ct_by_cond = pd.read_csv(results_dir / '01_celltype_by_condition.csv')
        interaction_summary = pd.read_csv(results_dir / '15_summary' / 'interaction_summary.csv')
        compare_df = pd.read_csv(results_dir / '09_compare_total' / 'compare_interactions.csv')
        diff_data = pd.read_csv(results_dir / '10_differential' / 'differential_data.csv')
        pathway_stats = pd.read_csv(results_dir / '18_pathway_stats' / 'pathway_stats_data.csv')
    except Exception as e:
        print(f"Error reading required data files: {e}")
        return

    # Optional data tables
    heatmap_file = results_dir / '03_heatmaps' / 'heatmap_data.csv'
    heatmap_df = pd.read_csv(heatmap_file) if heatmap_file.exists() else None

    roles_file = results_dir / '06_signaling_roles' / 'signaling_roles_data.csv'
    roles_df = pd.read_csv(roles_file) if roles_file.exists() else None

    info_flow_file = results_dir / '11_info_flow' / 'info_flow_data.csv'
    info_flow_df = pd.read_csv(info_flow_file) if info_flow_file.exists() else None

    report = []
    conditions = interaction_summary['Condition'].tolist()
    comparisons = diff_data['Comparison_Pair'].unique().tolist()
    cell_types = ct_by_cond['CellType'].tolist()

    # =========================================================================
    # ULTRATHINK HEADER
    # =========================================================================
    report.append(ULTRATHINK_HEADER)

    # =========================================================================
    # SECTION 1: EXECUTIVE SUMMARY
    # =========================================================================
    print("Generating Section 1: Executive Summary...")
    report.extend([
        "# CellChat Comprehensive Analysis Report",
        "",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "---",
        "",
        "## 1. Executive Summary",
        ""
    ])

    total_cells = data_summary.loc[data_summary['Metric'] == 'Total Cells', 'Value'].values[0]
    total_genes = data_summary.loc[data_summary['Metric'] == 'Total Genes', 'Value'].values[0]
    n_celltypes = data_summary.loc[data_summary['Metric'] == 'Cell Types', 'Value'].values[0]
    n_conditions = data_summary.loc[data_summary['Metric'] == 'Conditions', 'Value'].values[0]

    report.extend([
        f"This analysis inferred cell-cell communication networks from **{int(total_cells)} cells** ({int(total_genes)} genes) across **{int(n_celltypes)} cell types** in **{int(n_conditions)} conditions**.",
        "",
        "### Key Findings Summary",
        ""
    ])

    for _, row in interaction_summary.iterrows():
        report.append(f"- **{row['Condition']}**: {row['Num_LR_Pairs']} L-R pairs, {row['Num_Pathways']} pathways, {row['Total_Interactions']} total interactions")

    # =========================================================================
    # SECTION 2: INPUT DATA
    # =========================================================================
    print("Generating Section 2: Input Data...")
    report.extend([
        "",
        "---",
        "",
        "## 2. Input Data & Comparison Setup",
        "",
        "### 2.1 Cell Distribution by Condition",
        ""
    ])

    cols = ct_by_cond.columns.tolist()
    header = "| " + " | ".join(cols) + " |"
    separator = "| " + " | ".join(["---"] + ["---:"] * (len(cols)-1)) + " |"
    report.extend([header, separator])
    for _, row in ct_by_cond.iterrows():
        report.append("| " + " | ".join(str(v) for v in row) + " |")

    # =========================================================================
    # SECTION 3: ANALYSIS PARAMETERS
    # =========================================================================
    print("Generating Section 3: Analysis Parameters...")
    report.extend([
        "",
        "---",
        "",
        "## 3. Analysis Parameters",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        "| Database | CellChatDB v2 (human) |",
        "| Probability Method | truncatedMean |",
        "| Min Cells | 10 |",
        "| Trim Value | 0.1 |",
        ""
    ])

    # =========================================================================
    # SECTION 4: RESULTS WITH ULTRATHINK INTERPRETATIONS
    # =========================================================================
    print("Generating Section 4: Results & ULTRATHINK Interpretations...")
    report.extend([
        "---",
        "",
        "## 4. Results: Visualizations & ULTRATHINK Interpretation",
        ""
    ])

    def safe_name(s):
        return s.replace(" ", "_").replace("-", "_")

    # 4.1 Heatmaps
    report.extend([
        "### 4.1 Interaction Heatmaps",
        "",
        "Heatmaps show the number (count) and strength (weight) of cell-cell interactions.",
        "- **Rows** = Sender cell types (expressing ligands)",
        "- **Columns** = Receiver cell types (expressing receptors)",
        ""
    ])

    for cond in conditions:
        print(f"  - Processing heatmap interpretation for {cond}...")
        cond_safe = safe_name(cond)
        count_img = f"03_heatmaps/{cond_safe}_count_heatmap.png"
        weight_img = f"03_heatmaps/{cond_safe}_weight_heatmap.png"

        report.extend([
            f"#### {cond}",
            "",
            f"**Count Heatmap:** ![{cond} Count Heatmap]({count_img})",
            "",
            f"**Weight Heatmap:** ![{cond} Weight Heatmap]({weight_img})",
            ""
        ])

        if heatmap_df is not None:
            interp = interpret_heatmap(heatmap_df, cond, cell_types)
            report.extend([
                "**ULTRATHINK Data-Driven Interpretation:**",
                interp,
                ""
            ])

    # Add Claude Code interpretation prompt for Section 4.1
    report.append("<!-- CLAUDE_CODE_ULTRATHINK_SECTION_4.1 -->")
    report.append("<!--")
    report.append("When presenting Section 4.1, Claude Code MUST apply ULTRATHINK reasoning:")
    report.append("1. Interpret autocrine vs paracrine signaling balance - what does this indicate about tissue state?")
    report.append("2. Identify dominant signaling axes - which cell pairs show strongest communication?")
    report.append("3. Explain WHY specific cell types communicate intensely based on their biological functions")
    report.append("4. Note any asymmetric communication patterns (A→B vs B→A) and their biological significance")
    report.append("5. Use WebFetch to search for known communication roles of identified cell types in the study context")
    report.append("6. Compare patterns across conditions - what changes and why?")
    report.append("-->")
    report.append("")

    # 4.2 Circle Plots
    report.extend([
        "### 4.2 Circle Plots",
        "",
        "Circle plots show the aggregated communication network.",
        "- **Node size** = Number of cells in that cell type",
        "- **Edge width** = Number of interactions between cell types",
        ""
    ])

    for cond in conditions:
        cond_safe = safe_name(cond)
        circle_img = f"04_circle_plots/{cond_safe}_circle_plot.png"
        report.extend([
            f"#### {cond}",
            "",
            f"![{cond} Circle Plot]({circle_img})",
            ""
        ])

    # 4.3 Signaling Roles
    report.extend([
        "### 4.3 Signaling Roles",
        "",
        "Scatter plots classify cell types by their role in signaling:",
        "- **X-axis** = Outgoing signal strength (sender capability)",
        "- **Y-axis** = Incoming signal strength (receiver capability)",
        "- **Upper-right** = Communication hubs",
        "- **Lower-right** = Dominant senders",
        "- **Upper-left** = Dominant receivers",
        ""
    ])

    for cond in conditions:
        print(f"  - Processing signaling roles interpretation for {cond}...")
        cond_safe = safe_name(cond)
        roles_img = f"06_signaling_roles/{cond_safe}_signaling_roles.png"

        report.extend([
            f"#### {cond}",
            "",
            f"![{cond} Signaling Roles]({roles_img})",
            ""
        ])

        if roles_df is not None:
            interp = interpret_signaling_roles(roles_df, cond)
            report.extend([
                "**ULTRATHINK Data-Driven Interpretation:**",
                interp,
                ""
            ])

    # Add Claude Code interpretation prompt for Section 4.3
    report.append("<!-- CLAUDE_CODE_ULTRATHINK_SECTION_4.3 -->")
    report.append("<!--")
    report.append("When presenting Section 4.3, Claude Code MUST apply ULTRATHINK reasoning:")
    report.append("1. Identify communication HUBS - which cell types integrate and relay signals?")
    report.append("2. Classify cell types as senders, receivers, or hubs and explain WHY based on their biology")
    report.append("3. Explain the biological significance of sender/receiver roles in the disease/study context")
    report.append("4. Use WebFetch to search for known signaling roles of identified hub cells")
    report.append("5. Note how signaling roles CHANGE across conditions - what drives these shifts?")
    report.append("6. Identify potential therapeutic targets among hub/sender populations")
    report.append("7. Consider if hub cells could amplify or dampen disease-related signaling cascades")
    report.append("-->")
    report.append("")

    # 4.4 Comparison Summary with Enhanced Table
    report.extend([
        "### 4.4 Comparison Summary",
        "",
        "#### Quantitative Overview",
        "",
        "| Comparison | Pairs UP | Pairs DN | PW p<0.05 UP | PW p<0.05 DN | PW FDR<0.05 UP | PW FDR<0.05 DN | Net Change |",
        "|------------|----------|----------|--------------|--------------|----------------|----------------|------------|"
    ])

    comparison_interps = []
    for comp in comparisons:
        diff_comp = diff_data[diff_data['Comparison_Pair'] == comp]
        pw_comp = pathway_stats[pathway_stats['Comparison_Pair'] == comp]

        pairs_up = (diff_comp['Count_Direction'] == 'UP').sum()
        pairs_dn = (diff_comp['Count_Direction'] == 'DOWN').sum()
        pw_pval_up = ((pw_comp['Direction'] == 'UP') & (pw_comp['PValue'] < 0.05)).sum()
        pw_pval_dn = ((pw_comp['Direction'] == 'DOWN') & (pw_comp['PValue'] < 0.05)).sum()
        pw_fdr_up = ((pw_comp['Direction'] == 'UP') & pw_comp['Significant']).sum()
        pw_fdr_dn = ((pw_comp['Direction'] == 'DOWN') & pw_comp['Significant']).sum()
        net_change = pairs_up - pairs_dn

        report.append(f"| {comp} | {pairs_up} | {pairs_dn} | {pw_pval_up} | {pw_pval_dn} | {pw_fdr_up} | {pw_fdr_dn} | {net_change:+d} |")

        # Store interpretation for later - pass full pathway_stats DataFrame
        comp_interp = generate_comparison_interpretation(comp, pairs_up, pairs_dn, pw_comp)
        comparison_interps.append((comp, comp_interp))

    report.append("")

    # Add biological interpretations for each comparison
    report.append("#### ULTRATHINK Biological Interpretation of Comparisons")
    report.append("")
    for comp, interp in comparison_interps:
        report.append(f"**{comp}:**")
        report.append("")
        report.append(interp)
        report.append("")

    # Add Claude Code interpretation prompt for Section 4.4
    report.append("<!-- CLAUDE_CODE_ULTRATHINK_SECTION_4.4 -->")
    report.append("<!--")
    report.append("When presenting Section 4.4, Claude Code MUST apply ULTRATHINK reasoning:")
    report.append("1. Use WebFetch to search for TOP altered pathways in the study-specific biological context")
    report.append("2. Explain WHY these specific pathways would be altered in this biological process")
    report.append("3. Connect pathway changes to cell type communication patterns from Sections 4.1-4.3")
    report.append("4. Suggest mechanistic hypotheses for observed changes based on known pathway biology")
    report.append("5. Identify potential therapeutic implications for the altered pathways")
    report.append("6. Note key genes/receptors driving pathway changes if available in the data")
    report.append("-->")
    report.append("")

    # 4.5 Differential Interactions
    report.extend([
        "### 4.5 Differential Interactions",
        "",
        "Circle plots showing changes in cell-cell communication between conditions:",
        "- **RED edges** = Increased interactions in test vs reference",
        "- **BLUE edges** = Decreased interactions in test vs reference",
        ""
    ])

    for comp in comparisons:
        print(f"  - Processing differential interpretation for {comp}...")
        parts = comp.split("_vs_")
        if len(parts) == 2:
            test_cond, ref_cond = parts[0], parts[1]
            test_safe = safe_name(test_cond)
            ref_safe = safe_name(ref_cond)

            count_img = f"10_differential/{test_safe}_vs_{ref_safe}_count.png"
            weight_img = f"10_differential/{test_safe}_vs_{ref_safe}_weight.png"

            report.extend([
                f"#### {comp}",
                "",
                f"**Count Difference:** ![{comp} Count]({count_img})",
                "",
                f"**Weight Difference:** ![{comp} Weight]({weight_img})",
                ""
            ])

            interp = interpret_differential(diff_data, comp)
            report.extend([
                "**ULTRATHINK Data-Driven Interpretation:**",
                interp,
                ""
            ])

    # Add Claude Code interpretation prompt for Section 4.5
    report.append("<!-- CLAUDE_CODE_ULTRATHINK_SECTION_4.5 -->")
    report.append("<!--")
    report.append("When presenting Section 4.5, Claude Code MUST apply ULTRATHINK reasoning:")
    report.append("1. Explain WHY specific cell pairs show INCREASED communication in disease/test condition")
    report.append("2. Explain WHY specific cell pairs show DECREASED communication - is this protective or pathological?")
    report.append("3. Identify NOVEL communication axes that emerge in disease - what biological processes drive these?")
    report.append("4. Note LOST communication axes - what homeostatic signals are disrupted?")
    report.append("5. Use WebFetch to search for known interactions between identified cell type pairs in the study context")
    report.append("6. Consider if changes reflect cell expansion/contraction vs true signaling changes")
    report.append("7. Identify which differential interactions represent therapeutic intervention opportunities")
    report.append("-->")
    report.append("")

    # 4.6 Information Flow
    print("  - Processing information flow interpretation...")
    report.extend([
        "### 4.6 Information Flow Comparison",
        "",
        "Bar plots showing pathway activity across conditions:",
        "- **Left (Stacked)**: Relative contribution of each condition to pathway activity",
        "- **Right (Grouped)**: Direct comparison of pathway activity between conditions",
        ""
    ])

    info_flow_img = "11_info_flow/information_flow.png"
    report.extend([
        f"![Information Flow]({info_flow_img})",
        ""
    ])

    if info_flow_df is not None:
        interp = interpret_info_flow(info_flow_df, conditions)
        report.extend([
            "**ULTRATHINK Data-Driven Interpretation:**",
            interp,
            ""
        ])

    # Add Claude Code interpretation prompt for Section 4.6
    report.append("<!-- CLAUDE_CODE_ULTRATHINK_SECTION_4.6 -->")
    report.append("<!--")
    report.append("When presenting Section 4.6, Claude Code MUST apply ULTRATHINK reasoning:")
    report.append("1. Explain WHY the DOMINANT pathways (e.g., COLLAGEN, LAMININ) are most active in this tissue/context")
    report.append("2. Identify CONDITION-SPECIFIC pathway signatures - which pathways distinguish disease from normal?")
    report.append("3. Use WebFetch to search for the TOP 3 pathways' roles in the specific disease being studied")
    report.append("4. Note pathways with HIGH VARIABILITY between conditions - these are key drivers of phenotypic differences")
    report.append("5. Explain the biological significance of ECM vs Secreted vs Contact signaling balance")
    report.append("6. Identify pathway signatures consistent with known disease mechanisms from literature")
    report.append("7. Suggest which pathway alterations represent therapeutic opportunities")
    report.append("-->")
    report.append("")

    # 4.7 Pathway Statistics
    report.extend([
        "### 4.7 Pathway Differential Statistics",
        "",
        "Bar plots showing log2 fold change in pathway activity:",
        "- **RED bars** = Pathways upregulated in test condition",
        "- **BLUE bars** = Pathways downregulated in test condition",
        "- **Asterisk (*)** = Statistically significant (FDR < 0.05)",
        ""
    ])

    for comp in comparisons:
        print(f"  - Processing pathway stats interpretation for {comp}...")
        parts = comp.split("_vs_")
        if len(parts) == 2:
            test_cond, ref_cond = parts[0], parts[1]
            test_safe = safe_name(test_cond)
            ref_safe = safe_name(ref_cond)

            logfc_img = f"18_pathway_stats/{test_safe}_vs_{ref_safe}_pathway_logfc.png"

            report.extend([
                f"#### {comp}",
                "",
                f"![{comp} Pathway Log2FC]({logfc_img})",
                ""
            ])

            interp = interpret_pathway_stats(pathway_stats, comp)
            report.extend([
                "**ULTRATHINK Data-Driven Interpretation:**",
                interp,
                ""
            ])

    # =========================================================================
    # SECTION 5: DETAILED PATHWAY FINDINGS
    # =========================================================================
    print("Generating Section 5: Detailed Pathway Findings...")
    report.extend([
        "---",
        "",
        "## 5. Detailed Pathway Findings",
        ""
    ])

    for comp in comparisons:
        pw_comp = pathway_stats[pathway_stats['Comparison_Pair'] == comp].copy()
        pw_comp = pw_comp.sort_values('Log2FC', key=abs, ascending=False)

        report.extend([
            f"### {comp}",
            "",
            "#### Top Upregulated Pathways (FDR < 0.05):",
            ""
        ])

        up_pw = pw_comp[(pw_comp['Direction'] == 'UP') & pw_comp['Significant']]
        if len(up_pw) > 0:
            for _, row in up_pw.head(5).iterrows():
                pw_info = get_pathway_info(row['Pathway'])
                report.append(f"- **{row['Pathway']}** (log2FC = {row['Log2FC']:.2f}, FDR = {row['FDR']:.4f})")
                report.append(f"  - {pw_info['general_role']}")
        else:
            report.append("- None significant")

        report.extend([
            "",
            "#### Top Downregulated Pathways (FDR < 0.05):",
            ""
        ])

        dn_pw = pw_comp[(pw_comp['Direction'] == 'DOWN') & pw_comp['Significant']]
        if len(dn_pw) > 0:
            for _, row in dn_pw.head(5).iterrows():
                pw_info = get_pathway_info(row['Pathway'])
                report.append(f"- **{row['Pathway']}** (log2FC = {row['Log2FC']:.2f}, FDR = {row['FDR']:.4f})")
                report.append(f"  - {pw_info['general_role']}")
        else:
            report.append("- None significant")

        # Show nominally significant if no FDR significant
        if len(up_pw) == 0 and len(dn_pw) == 0:
            pval_sig = pw_comp[pw_comp['PValue'] < 0.05].head(5)
            if len(pval_sig) > 0:
                report.extend([
                    "",
                    "#### Top Nominally Significant Pathways (p < 0.05, for hypothesis generation):",
                    ""
                ])
                for _, row in pval_sig.iterrows():
                    pw_info = get_pathway_info(row['Pathway'])
                    direction = "↑" if row['Direction'] == 'UP' else "↓"
                    report.append(f"- {direction} **{row['Pathway']}** (log2FC = {row['Log2FC']:.2f}, p = {row['PValue']:.4f})")
                    report.append(f"  - {pw_info['description']}")

        report.append("")

    # =========================================================================
    # SECTION 6: DISCUSSION
    # =========================================================================
    print("Generating Section 6: ULTRATHINK Discussion & Next Steps...")
    report.extend([
        "---",
        "",
        "## 6. Discussion & Proposed Next Steps",
        ""
    ])

    all_data = {
        'interaction_summary': interaction_summary,
        'diff_data': diff_data,
        'pathway_stats': pathway_stats,
        'heatmap_df': heatmap_df,
        'roles_df': roles_df,
        'info_flow_df': info_flow_df
    }

    discussion = generate_discussion(all_data, conditions, comparisons)
    report.append(discussion)

    # Add Claude Code interpretation prompt for Section 6
    report.append("")
    report.append("<!-- CLAUDE_CODE_ULTRATHINK_SECTION_6 -->")
    report.append("<!--")
    report.append("When presenting Section 6 Discussion, Claude Code MUST apply ULTRATHINK reasoning:")
    report.append("1. SYNTHESIZE all findings into a COHERENT BIOLOGICAL NARRATIVE for the specific disease/context")
    report.append("2. Use WebFetch to search PubMed/literature for validation of key findings")
    report.append("3. Identify the MOST IMPORTANT findings with highest confidence and biological plausibility")
    report.append("4. Generate TESTABLE HYPOTHESES based on observed communication patterns")
    report.append("5. Propose SPECIFIC EXPERIMENTS for validation (in vitro co-culture, in vivo models, clinical correlation)")
    report.append("6. Identify THERAPEUTIC TARGETS - which pathways/cell types are most druggable?")
    report.append("7. Note KEY UNCERTAINTIES and limitations requiring further investigation")
    report.append("8. Compare findings with known disease biology - are they concordant or surprising?")
    report.append("9. Suggest BIOMARKERS based on altered cell-cell communication patterns")
    report.append("10. Provide CLINICAL TRANSLATION perspective - how could these findings inform patient care?")
    report.append("-->")
    report.append("")

    # =========================================================================
    # SECTION 7: OUTPUT DOCUMENTATION
    # =========================================================================
    report.extend([
        "",
        "---",
        "",
        "## 7. Detailed Output Documentation",
        "",
        "For detailed step-by-step interpretation guide for all output files and plots, see **cellchat_report.md**.",
        "",
        "---",
        "",
        f"*Report generated by CellChat Comprehensive Report Generator with ULTRATHINK Framework on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*"
    ])

    # Write report
    with open(output_file, 'w') as f:
        f.write('\n'.join(report))

    print(f"\nComprehensive report generated: {output_file}")


if __name__ == '__main__':
    main()
