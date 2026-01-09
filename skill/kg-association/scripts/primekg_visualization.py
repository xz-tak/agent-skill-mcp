#!/usr/bin/env python3
"""
PrimeKG Network Visualization for kg-association skill.

Creates interactive network visualizations from PrimeKG shortest path analysis results.
Uses pyvis to generate standalone HTML files with interactive graph exploration.

Features:
- Displays target genes and disease nodes
- Shows shared PPI partners between combo targets (highlighted)
- Shows unique PPI partners for each target
- Shows disease-associated genes
- Node limit of 50 to keep visualization readable

Usage:
    # From JSON results file with extended data
    pixi run python primekg_visualization.py \
        --json-results ./kgpred_IBD/data/primekg/primekg_shortest_paths.json \
        --output ./kgpred_IBD/primekg/network_visualization.html

    # For combo visualization with PPI data
    pixi run python primekg_visualization.py \
        --json-results ./results.json \
        --combo "TYK2+JAK1" \
        --output ./tyk2_jak1_subgraph.html
"""

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    from pyvis.network import Network
    PYVIS_AVAILABLE = True
except ImportError:
    PYVIS_AVAILABLE = False
    print("Warning: pyvis not installed. Install with: pip install pyvis")


# Color scheme for node types
COLORS = {
    "disease": "#e74c3c",              # Red for diseases
    "gene_target": "#3498db",          # Blue for target genes
    "gene/protein": "#3498db",         # Blue for genes
    "shared_ppi": "#9b59b6",           # Purple for shared PPI partners
    "unique_ppi": "#2ecc71",           # Green for unique PPI partners
    "disease_associated": "#f39c12",   # Orange for disease-associated genes
    "intermediate": "#95a5a6",         # Gray for other intermediate nodes
    "biological_process": "#9b59b6",   # Purple for GO BP
    "molecular_function": "#9b59b6",   # Purple for GO MF
    "pathway": "#f39c12",              # Orange for pathways
    "anatomy": "#1abc9c",              # Teal for anatomy
    "phenotype": "#e67e22",            # Dark orange for phenotypes
    "drug": "#27ae60",                 # Green for drugs
    "other": "#95a5a6",                # Gray for others
}

# Edge colors
EDGE_COLORS = {
    "ppi_shared": "#9b59b6",           # Purple for shared PPI
    "ppi_unique": "#2ecc71",           # Green for unique PPI
    "disease_association": "#f39c12",  # Orange for disease association
    "target_to_disease": "#e74c3c",    # Red for target-disease paths
    "from_gene": "#3498db",            # Blue
    "intermediate": "#95a5a6",         # Gray
    "to_disease": "#e74c3c",           # Red
}

# Maximum nodes to display (including targets and diseases)
MAX_NODES = 75

# Node allocation budget
NODE_BUDGET = {
    "targets": 10,           # Target genes (unlimited within reason)
    "diseases": 10,          # Disease nodes
    "target_ppi": 5,         # PPI between targets (highest priority)
    "path_intermediate": 15, # Intermediary nodes in target-disease paths
    "disease_genes": 15,     # Disease-associated genes (shared associations)
    "shared_ppi": 10,        # Shared PPI partners between targets
    "unique_ppi_per_target": 5,  # Unique PPI per target (lower priority)
}


def get_node_color(node_type: str, node_category: str = None, is_target: bool = False) -> str:
    """Get color based on node type and category."""
    if is_target:
        return COLORS["gene_target"]
    if node_category == "shared_ppi":
        return COLORS["shared_ppi"]
    if node_category == "unique_ppi":
        return COLORS["unique_ppi"]
    if node_category == "disease_associated":
        return COLORS["disease_associated"]
    return COLORS.get(node_type, COLORS["other"])


def create_enhanced_network_visualization(
    genes: List[str],
    diseases: List[str],
    paths_data: Dict[str, Dict[str, Any]],
    ppi_data: Optional[Dict[str, Any]] = None,
    disease_genes: Optional[Dict[str, List[str]]] = None,
    title: str = "Gene-Disease Network (PrimeKG)",
    output_path: str = "network_visualization.html",
    height: str = "800px",
    width: str = "100%",
    max_nodes: int = MAX_NODES,
) -> str:
    """
    Create an enhanced interactive network visualization using pyvis.

    Args:
        genes: List of target gene names
        diseases: List of disease names
        paths_data: Dict mapping gene -> disease -> {path, path_details, score, ...}
        ppi_data: Optional dict with PPI information:
            {
                "shared_ppi": ["GENE1", "GENE2", ...],  # Shared between all targets
                "target_ppi": {
                    "TYK2": ["GENE3", "GENE4", ...],  # Unique to TYK2
                    "JAK1": ["GENE5", "GENE6", ...],  # Unique to JAK1
                }
            }
        disease_genes: Optional dict mapping disease -> [associated_genes]
        title: Title for the visualization
        output_path: Path to save the HTML file
        height: Height of the visualization
        width: Width of the visualization
        max_nodes: Maximum number of nodes to display

    Returns:
        Path to the generated HTML file
    """
    if not PYVIS_AVAILABLE:
        raise ImportError("pyvis is required for visualization. Install with: pip install pyvis")

    # Initialize network
    net = Network(
        height=height,
        width=width,
        bgcolor="#ffffff",
        font_color="black",
        directed=False,
        notebook=False,
    )

    # Configure physics for better layout
    net.barnes_hut(
        gravity=-8000,
        central_gravity=0.3,
        spring_length=150,
        spring_strength=0.01,
        damping=0.09,
    )

    added_nodes: Dict[str, Dict[str, Any]] = {}  # node_id -> node_info
    added_edges: Set[Tuple] = set()
    target_genes_set = set(genes)

    def add_node(node_id: str, label: str, node_type: str, title_text: str,
                 size: int = 25, is_target: bool = False, category: str = None,
                 priority: int = 0):
        """Add node if not already added, tracking priority for budget management."""
        if node_id in added_nodes:
            return True

        color = get_node_color(node_type, category, is_target)
        shape = "diamond" if node_type == "disease" else "dot"

        # Store node info for later addition
        added_nodes[node_id] = {
            "label": label,
            "title": title_text,
            "color": color,
            "size": size,
            "shape": shape,
            "font": {"size": 14 if size >= 40 else 12},
            "category": category or node_type,
            "priority": priority,
            "is_target": is_target,
        }
        return True

    def add_edge(source: str, target: str, relation: str, color: str,
                 width: float = 2, dashes: bool = False):
        """Add edge if not already added."""
        edge_key = tuple(sorted([source, target])) + (relation,)
        if edge_key in added_edges:
            return
        added_edges.add(edge_key)

    # Priority levels (higher = more important)
    # Reordered to prioritize: targets, diseases, target-target PPI, path intermediates, disease genes
    PRIORITY_TARGET = 100
    PRIORITY_DISEASE = 95
    PRIORITY_TARGET_PPI = 90          # PPI between targets (e.g., TYK2-JAK1)
    PRIORITY_PATH_INTERMEDIATE = 85   # Nodes along target-disease paths
    PRIORITY_DISEASE_GENE = 80        # Disease-associated genes
    PRIORITY_SHARED_PPI = 70          # Shared PPI between targets (not target-target)
    PRIORITY_UNIQUE_PPI = 60          # Unique PPI per target

    # 1. Add target gene nodes (highest priority)
    for gene in genes:
        add_node(
            node_id=gene,
            label=gene,
            node_type="gene/protein",
            title_text=f"{gene}\nType: gene/protein\nTarget gene",
            size=45,
            is_target=True,
            category="target",
            priority=PRIORITY_TARGET,
        )

    # 2. Add disease nodes (high priority)
    for disease in diseases:
        add_node(
            node_id=disease,
            label=disease[:30] + "..." if len(disease) > 30 else disease,
            node_type="disease",
            title_text=f"{disease}\nType: disease",
            size=50,
            priority=PRIORITY_DISEASE,
        )

    # 3. Add direct PPI edges between targets (if they exist)
    # This handles the TYK2-JAK1 direct interaction
    for i, target1 in enumerate(genes):
        for target2 in genes[i+1:]:
            add_edge(target1, target2, "ppi (target-target)", EDGE_COLORS["ppi_shared"], width=4)

    # 4. Add path intermediate nodes from shortest path analysis (HIGH PRIORITY)
    # These are nodes along the paths from targets to diseases
    intermediate_count = 0
    for gene in genes:
        if gene not in paths_data:
            continue
        for disease, result in paths_data[gene].items():
            if not result.get("path_found", False):
                continue
            path_details = result.get("path_details", [])
            score = result.get("score", 0)

            if path_details:
                for i, node_info in enumerate(path_details):
                    node_id = node_info.get("name", node_info.get("id", ""))
                    node_type = node_info.get("type", "other")

                    # Skip if already added or is target/disease
                    if node_id in added_nodes or node_id in target_genes_set or node_id in diseases:
                        continue

                    # Add intermediate nodes with high priority
                    if i > 0 and i < len(path_details) - 1:  # Intermediate
                        if intermediate_count < NODE_BUDGET["path_intermediate"]:
                            add_node(
                                node_id=node_id,
                                label=node_id,
                                node_type=node_type,
                                title_text=f"{node_id}\nType: {node_type}\nPath intermediate ({gene} → {disease})\nScore: {score:.3f}",
                                size=30,
                                category="intermediate",
                                priority=PRIORITY_PATH_INTERMEDIATE,
                            )
                            intermediate_count += 1

            # Add edges along path
            path_nodes = path_details if path_details else [{"name": p} for p in result.get("path", [])]
            for j in range(len(path_nodes) - 1):
                source_node = path_nodes[j].get("name", path_nodes[j]) if isinstance(path_nodes[j], dict) else path_nodes[j]
                target_node = path_nodes[j+1].get("name", path_nodes[j+1]) if isinstance(path_nodes[j+1], dict) else path_nodes[j+1]

                if j == 0:
                    edge_color = EDGE_COLORS["from_gene"]
                elif j == len(path_nodes) - 2:
                    edge_color = EDGE_COLORS["to_disease"]
                else:
                    edge_color = EDGE_COLORS["intermediate"]

                add_edge(source_node, target_node, "path", edge_color, width=2)

    # 5. Add disease-associated genes (prioritize shared across diseases)
    if disease_genes:
        # First, find genes shared across multiple diseases
        gene_disease_count = {}
        for disease, assoc_genes in disease_genes.items():
            for gene in assoc_genes:
                if gene not in target_genes_set:
                    gene_disease_count[gene] = gene_disease_count.get(gene, [])
                    gene_disease_count[gene].append(disease)

        # Sort by number of disease associations (shared genes first)
        sorted_genes = sorted(gene_disease_count.items(), key=lambda x: -len(x[1]))

        disease_gene_count = 0
        for gene, associated_diseases in sorted_genes:
            if disease_gene_count >= NODE_BUDGET["disease_genes"]:
                break
            if gene not in added_nodes:
                is_shared = len(associated_diseases) > 1
                add_node(
                    node_id=gene,
                    label=gene,
                    node_type="gene/protein",
                    title_text=f"{gene}\nType: gene/protein\nAssociated with: {', '.join(associated_diseases)}" +
                              ("\n(Shared across diseases)" if is_shared else ""),
                    size=32 if is_shared else 28,
                    category="disease_associated",
                    priority=PRIORITY_DISEASE_GENE + (5 if is_shared else 0),  # Boost shared genes
                )
                # Add edges to all associated diseases
                for disease in associated_diseases:
                    add_edge(gene, disease, "associated with", EDGE_COLORS["disease_association"], width=2)
                disease_gene_count += 1

    # 6. Add shared PPI partners (genes that interact with ALL targets)
    if ppi_data and "shared_ppi" in ppi_data:
        shared_ppi = ppi_data["shared_ppi"][:NODE_BUDGET["shared_ppi"]]
        for ppi_gene in shared_ppi:
            if ppi_gene not in target_genes_set and ppi_gene not in added_nodes:
                add_node(
                    node_id=ppi_gene,
                    label=ppi_gene,
                    node_type="gene/protein",
                    title_text=f"{ppi_gene}\nType: gene/protein\nShared PPI partner (interacts with all targets)",
                    size=32,
                    category="shared_ppi",
                    priority=PRIORITY_SHARED_PPI,
                )
                # Add edges to all targets
                for target in genes:
                    add_edge(target, ppi_gene, "ppi (shared)", EDGE_COLORS["ppi_shared"], width=3)

    # 7. Add unique PPI partners for each target (lower priority)
    if ppi_data and "target_ppi" in ppi_data:
        for target, ppi_genes in ppi_data["target_ppi"].items():
            if target not in genes:
                continue
            unique_count = 0
            for ppi_gene in ppi_genes:
                if unique_count >= NODE_BUDGET["unique_ppi_per_target"]:
                    break
                if ppi_gene not in target_genes_set and ppi_gene not in added_nodes:
                    add_node(
                        node_id=ppi_gene,
                        label=ppi_gene,
                        node_type="gene/protein",
                        title_text=f"{ppi_gene}\nType: gene/protein\nPPI partner of: {target}",
                        size=25,
                        category="unique_ppi",
                        priority=PRIORITY_UNIQUE_PPI,
                    )
                    add_edge(target, ppi_gene, f"ppi ({target})", EDGE_COLORS["ppi_unique"], width=2, dashes=True)
                    unique_count += 1

    # 8. Enforce node budget - sort by priority and take top max_nodes
    sorted_nodes = sorted(added_nodes.items(), key=lambda x: -x[1]["priority"])
    final_nodes = dict(sorted_nodes[:max_nodes])

    # Add final nodes to network
    for node_id, node_info in final_nodes.items():
        net.add_node(
            node_id,
            label=node_info["label"],
            title=node_info["title"],
            color=node_info["color"],
            size=node_info["size"],
            shape=node_info["shape"],
            font=node_info["font"],
        )

    # Add edges only if both nodes are in final set
    final_node_ids = set(final_nodes.keys())
    for edge in added_edges:
        source, target, relation = edge[0], edge[1], edge[2] if len(edge) > 2 else "connection"
        if source in final_node_ids and target in final_node_ids:
            # Determine color and width from relation
            if "target-target" in relation.lower():
                # Direct PPI between target genes (highest priority edge)
                color = "#e74c3c"  # Red - matches target-disease paths
                width = 5
            elif "shared" in relation.lower():
                color = EDGE_COLORS["ppi_shared"]
                width = 3
            elif "ppi" in relation.lower():
                color = EDGE_COLORS["ppi_unique"]
                width = 2
            elif "associated" in relation.lower():
                color = EDGE_COLORS["disease_association"]
                width = 2
            else:
                color = EDGE_COLORS["intermediate"]
                width = 2

            net.add_edge(source, target, title=relation, color=color, width=width)

    # Count nodes by category
    category_counts = {}
    for node_id, info in final_nodes.items():
        cat = info.get("category", "other")
        category_counts[cat] = category_counts.get(cat, 0) + 1

    # Generate HTML
    html_content = net.generate_html()

    # Add enhanced legend
    legend_html = f"""
    <div style="position: absolute; top: 10px; left: 10px; background: rgba(255,255,255,0.95);
                padding: 15px; border-radius: 8px; color: black; font-family: Arial, sans-serif;
                z-index: 1000; box-shadow: 0 2px 10px rgba(0,0,0,0.1); max-width: 280px;">
        <h3 style="margin: 0 0 10px 0; color: #333; font-size: 16px;">{title}</h3>
        <div style="font-size: 12px; line-height: 1.8;">
            <b>Node Types:</b>
            <div><span style="color: {COLORS['disease']};">&#9670;</span> Disease ({category_counts.get('disease', 0)})</div>
            <div><span style="color: {COLORS['gene_target']};">&#9679;</span> Target Gene ({category_counts.get('target', 0)})</div>
            <div><span style="color: {COLORS['shared_ppi']};">&#9679;</span> Shared PPI ({category_counts.get('shared_ppi', 0)})</div>
            <div><span style="color: {COLORS['disease_associated']};">&#9679;</span> Disease-Associated ({category_counts.get('disease_associated', 0)})</div>
            <div><span style="color: {COLORS['unique_ppi']};">&#9679;</span> Unique PPI ({category_counts.get('unique_ppi', 0)})</div>
            <div><span style="color: {COLORS['intermediate']};">&#9679;</span> Path Intermediate ({category_counts.get('intermediate', 0)})</div>
            <hr style="border: none; border-top: 1px solid #ddd; margin: 8px 0;">
            <b>Edge Types:</b>
            <div><span style="color: #e74c3c;">&#9473;&#9473;</span> Target-Target PPI</div>
            <div><span style="color: {EDGE_COLORS['ppi_shared']};">&#9473;</span> Shared PPI</div>
            <div><span style="color: {EDGE_COLORS['ppi_unique']};">&#8226;&#8226;&#8226;</span> Unique PPI</div>
            <div><span style="color: {EDGE_COLORS['disease_association']};">&#8212;</span> Disease Association</div>
            <div><span style="color: {EDGE_COLORS['to_disease']};">&#8212;</span> Target-Disease Path</div>
        </div>
        <div style="margin-top: 10px; font-size: 11px; color: #888;">
            Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}<br>
            Total Nodes: {len(final_nodes)} (max {max_nodes})
        </div>
    </div>
    """

    html_content = html_content.replace('<body>', f'<body>{legend_html}')

    # Save to file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html_content)

    print(f"Visualization saved to: {output_path}")
    print(f"  Total nodes: {len(final_nodes)} (limit: {max_nodes})")
    print(f"  Node breakdown: {category_counts}")

    return str(output_path)


# Keep original function for backwards compatibility
def create_network_visualization(
    genes: List[str],
    diseases: List[str],
    paths_data: Dict[str, Dict[str, Any]],
    title: str = "Gene-Disease Network (PrimeKG)",
    output_path: str = "network_visualization.html",
    height: str = "800px",
    width: str = "100%",
    ppi_data: Optional[Dict[str, Any]] = None,
    disease_genes: Optional[Dict[str, List[str]]] = None,
) -> str:
    """
    Create an interactive network visualization using pyvis.
    Wrapper that calls the enhanced version.
    """
    return create_enhanced_network_visualization(
        genes=genes,
        diseases=diseases,
        paths_data=paths_data,
        ppi_data=ppi_data,
        disease_genes=disease_genes,
        title=title,
        output_path=output_path,
        height=height,
        width=width,
    )


def create_combo_visualization(
    combo: List[str],
    diseases: List[str],
    paths_data: Dict[str, Dict[str, Any]],
    output_path: str,
    ppi_data: Optional[Dict[str, Any]] = None,
    disease_genes: Optional[Dict[str, List[str]]] = None,
) -> str:
    """
    Create visualization for a specific combo showing all component paths.

    Args:
        combo: List of gene names in the combo
        diseases: List of disease names
        paths_data: Full paths data dict
        output_path: Path to save the HTML file
        ppi_data: Optional PPI data
        disease_genes: Optional disease-associated genes

    Returns:
        Path to the generated HTML file
    """
    combo_name = "+".join(combo)
    title = f"{combo_name} Subgraph (PrimeKG)"

    # Filter paths_data to only include combo genes
    filtered_paths = {gene: paths_data.get(gene, {}) for gene in combo}

    # Filter PPI data to combo targets
    filtered_ppi = None
    if ppi_data:
        filtered_ppi = {
            "shared_ppi": ppi_data.get("shared_ppi", []),
            "target_ppi": {g: ppi_data.get("target_ppi", {}).get(g, []) for g in combo}
        }

    return create_enhanced_network_visualization(
        genes=combo,
        diseases=diseases,
        paths_data=filtered_paths,
        ppi_data=filtered_ppi,
        disease_genes=disease_genes,
        title=title,
        output_path=output_path,
    )


def load_results_from_json(json_path: str) -> Dict[str, Any]:
    """Load analysis results from JSON file."""
    with open(json_path, "r") as f:
        return json.load(f)


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Create PrimeKG network visualizations from analysis results"
    )
    parser.add_argument(
        "--json-results",
        type=str,
        help="Path to JSON results file from PrimeKG analysis",
    )
    parser.add_argument(
        "--genes",
        type=str,
        help="Comma-separated list of gene names (if not using --json-results)",
    )
    parser.add_argument(
        "--diseases",
        type=str,
        help="Comma-separated list of disease names (if not using --json-results)",
    )
    parser.add_argument(
        "--paths",
        type=str,
        help="JSON string of paths data (if not using --json-results)",
    )
    parser.add_argument(
        "--combo",
        type=str,
        help="Create combo-specific visualization (e.g., 'TYK2+JAK1')",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output path for HTML visualization",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Gene-Disease Network (PrimeKG)",
        help="Title for the visualization",
    )
    parser.add_argument(
        "--max-nodes",
        type=int,
        default=MAX_NODES,
        help=f"Maximum number of nodes to display (default: {MAX_NODES})",
    )

    args = parser.parse_args()

    if not PYVIS_AVAILABLE:
        print("Error: pyvis is required. Install with: pip install pyvis")
        return 1

    # Load data from JSON file or arguments
    ppi_data = None
    disease_genes = None

    if args.json_results:
        results = load_results_from_json(args.json_results)
        genes = results.get("genes", [])
        diseases = results.get("diseases", [])
        paths_data = results.get("individual", {})
        ppi_data = results.get("ppi_data", None)
        disease_genes = results.get("disease_genes", None)
    else:
        if not all([args.genes, args.diseases]):
            print("Error: Either --json-results or --genes and --diseases are required")
            return 1

        genes = [g.strip() for g in args.genes.split(",")]
        diseases = [d.strip() for d in args.diseases.split(",")]
        paths_data = json.loads(args.paths) if args.paths else {}

    # Create visualization
    if args.combo:
        combo_genes = [g.strip() for g in args.combo.split("+")]
        create_combo_visualization(
            combo=combo_genes,
            diseases=diseases,
            paths_data=paths_data,
            output_path=args.output,
            ppi_data=ppi_data,
            disease_genes=disease_genes,
        )
    else:
        create_enhanced_network_visualization(
            genes=genes,
            diseases=diseases,
            paths_data=paths_data,
            ppi_data=ppi_data,
            disease_genes=disease_genes,
            title=args.title,
            output_path=args.output,
            max_nodes=args.max_nodes,
        )

    return 0


if __name__ == "__main__":
    exit(main())
