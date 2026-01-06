#!/usr/bin/env python3
"""
Example queries for UltraQuery Inference MCP Server

This file demonstrates various query types and how to construct them.

NOTE: All queries use nested lists (not tuples) for compatibility with JSON serialization
and the MCP server API.
"""

# Example 1: One-hop projection (1p)
# Question: What proteins interact with GREM1?
query_1p = ["GREM1", ["ppi"]]

# Example 2: Two-hop projection (2p)
# Question: What diseases are associated with proteins that interact with GREM1?
query_2p = ["GREM1", ["ppi", "associated with"]]

# Example 3: Three-hop projection (3p)
# Question: Follow a path of 3 relations starting from GREM1
query_3p = ["GREM1", ["ppi", "associated with", "phenotype present"]]

# Example 4: Two-way intersection (2i)
# Question: What proteins interact with both GREM1 and IL11?
query_2i = [["GREM1", ["ppi"]], ["IL11", ["ppi"]]]

# Example 5: Three-way intersection (3i)
# Question: What proteins interact with GREM1, IL11, AND TGFB1?
query_3i = [["GREM1", ["ppi"]], ["IL11", ["ppi"]], ["TGFB1", ["ppi"]]]

# Example 6: Intersection then projection (ip)
# Question: What diseases are associated with proteins that interact with both GREM1 and IL11?
query_ip = [[["GREM1", ["ppi"]], ["IL11", ["ppi"]]], ["associated with"]]

# Example 7: Projection then intersection (pi)
# Question: Find entities that are:
#   - Two hops away from GREM1 (via ppi → associated with)
#   - One hop away from a disease (via phenotype present)
query_pi = [["GREM1", ["ppi", "associated with"]], ["MONDO:5301", ["phenotype present"]]]

# Example 8: Intersection with negation (2in)
# Question: What proteins interact with GREM1 but NOT with IL11?
query_2in = [["GREM1", ["ppi"]], ["IL11", ["ppi", "n"]]]

# Example 9: Three-way intersection with negation (3in)
# Question: What proteins interact with GREM1 and IL11, but NOT TGFB1?
query_3in = [["GREM1", ["ppi"]], ["IL11", ["ppi"]], ["TGFB1", ["ppi", "n"]]]

# Example 10: Intersection with negation, then projection (inp)
# Question: What diseases are associated with proteins that interact with GREM1 but NOT IL11?
query_inp = [[["GREM1", ["ppi"]], ["IL11", ["ppi", "n"]]], ["associated with"]]

# Example 11: Using entity IDs instead of names
# Same as Example 4, but using entity IDs
query_2i_ids = [["NCBI:9796", ["ppi"]], ["NCBI:3589", ["ppi"]]]

# Example 12: Mixed entity IDs and names
# Same as Example 4, but mixing IDs and names
query_2i_mixed = [["NCBI:9796", ["ppi"]], ["IL11", ["ppi"]]]

# Example 13: Union query (2u-DNF)
# Question: What entities are connected to GREM1 via ppi OR to IL11 via ppi?
query_2u_dnf = [["GREM1", ["ppi"]], ["IL11", ["ppi"]], ["u"]]

# Example 14: Union then projection (up-DNF)
# Question: What diseases are associated with proteins that interact with GREM1 OR IL11?
query_up_dnf = [[["GREM1", ["ppi"]], ["IL11", ["ppi"]], ["u"]], ["associated with"]]

# Example 15: Complex biomedical query
# Question: What drugs target proteins that:
#   - Interact with GREM1
#   - Are NOT associated with Crohn's disease
query_complex = [[["GREM1", ["ppi"]], ["Crohn disease", ["associated with", "n"]]], ["target"]]


if __name__ == "__main__":
    """
    Example usage showing how to call the MCP server with these queries.

    NOTE: This requires the MCP server to be running and accessible.
    """
    import sys
    import os

    # Add parent directory to path
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    # Import the server's answer_complex_query function directly
    # (In production, you'd use HTTP requests instead)
    from ultraquery_inference_server import answer_complex_query

    # List of example queries to test
    examples = [
        ("1p: One-hop projection", query_1p),
        ("2p: Two-hop projection", query_2p),
        ("2i: Two-way intersection", query_2i),
        ("ip: Intersection then projection", query_ip),
        ("2in: Intersection with negation", query_2in),
        ("inp: Intersection with negation, then projection", query_inp),
    ]

    print("=" * 80)
    print("UltraQuery Inference Examples")
    print("=" * 80)

    for name, query in examples:
        print(f"\n{name}")
        print("-" * 80)
        print(f"Query structure: {query}")

        # Call the inference function
        try:
            result = answer_complex_query(query_structure=query, top_k=10)

            if result["success"]:
                print(f"\nQuery type: {result['query_type']}")
                print(f"Inference time: {result['inference_time_seconds']}s")
                print(f"\nQuery (readable):\n{result['query_readable']}")

                print(f"\nTop 10 predictions:")
                for pred in result["predictions"]:
                    print(f"  {pred['rank']:2d}. {pred['entity_name']:40s} ({pred['entity_id']:20s}) - Score: {pred['score']:.4f}")
            else:
                print(f"Error: {result['error']}")

        except Exception as e:
            print(f"Exception: {e}")
            import traceback
            traceback.print_exc()

        print()

    print("=" * 80)
    print("Examples complete!")
    print("=" * 80)
