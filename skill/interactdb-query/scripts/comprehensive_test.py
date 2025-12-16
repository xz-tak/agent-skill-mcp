#!/usr/bin/env python3
"""
Comprehensive test suite for all interaction database APIs.
Tests both single-gene queries and shortest path functions with various filters.
"""

import sys
from string_api import StringClient, get_string_neighbors
from intact_api import get_direct_neighbors, find_shortest_paths_intact

print("="*80)
print("COMPREHENSIVE API TESTING - ALL FEATURES")
print("="*80)
print()

# ==============================================================================
# TEST 1: STRING Single-Gene Query - Various Score Filters
# ==============================================================================
print("="*80)
print("TEST 1: STRING Single-Gene Query - Score Filter Combinations")
print("="*80)
print()

test_configs = [
    {
        "name": "Basic query (combined score only)",
        "params": {
            "gene": "TP53",
            "top_n": 10,
            "min_score": 700,
            "species": 9606
        }
    },
    {
        "name": "High experimental evidence only",
        "params": {
            "gene": "TP53",
            "top_n": 10,
            "min_score": 400,
            "min_experimental": 800,
            "species": 9606
        }
    },
    {
        "name": "Database + textmining combined",
        "params": {
            "gene": "TP53",
            "top_n": 10,
            "min_score": 400,
            "min_database": 500,
            "min_textmining": 500,
            "species": 9606
        }
    },
    {
        "name": "Coexpression evidence",
        "params": {
            "gene": "TP53",
            "top_n": 10,
            "min_score": 400,
            "min_coexpression": 600,
            "species": 9606
        }
    }
]

for i, config in enumerate(test_configs, 1):
    print(f"\n{i}. {config['name']}")
    print("-" * 70)
    print(f"Parameters: {config['params']}")
    
    try:
        df = get_string_neighbors(**config['params'])
        if not df.empty:
            print(f"✓ Found {len(df)} neighbors")
            print("\nTop 3 results:")
            display_cols = ['preferred_name', 'combined_score', 'experimental_score', 
                          'database_score', 'textmining_score', 'coexpression_score', 'path']
            print(df[display_cols].head(3).to_string(index=False))
        else:
            print("✗ No neighbors found")
    except Exception as e:
        print(f"✗ Error: {e}")

print("\n")

# ==============================================================================
# TEST 2: STRING Shortest Paths - Various Score Filters
# ==============================================================================
print("="*80)
print("TEST 2: STRING Shortest Paths - Score Filter Combinations")
print("="*80)
print()

genes = ["TP53", "MDM2", "ATM"]
path_configs = [
    {
        "name": "Combined score filter only",
        "params": {
            "gene_list": genes,
            "species": 9606,
            "max_distance": 3,
            "min_combined_score": 400
        }
    },
    {
        "name": "High experimental evidence required",
        "params": {
            "gene_list": genes,
            "species": 9606,
            "max_distance": 3,
            "min_combined_score": 400,
            "min_experimental_score": 700
        }
    },
    {
        "name": "Database evidence required",
        "params": {
            "gene_list": genes,
            "species": 9606,
            "max_distance": 3,
            "min_combined_score": 400,
            "min_database_score": 500
        }
    },
    {
        "name": "Multiple evidence types required",
        "params": {
            "gene_list": genes,
            "species": 9606,
            "max_distance": 3,
            "min_combined_score": 700,
            "min_experimental_score": 500,
            "min_database_score": 500
        }
    }
]

client = StringClient()

for i, config in enumerate(path_configs, 1):
    print(f"\n{i}. {config['name']}")
    print("-" * 70)
    filter_params = {k:v for k,v in config['params'].items() if k.startswith('min_') or k == 'max_distance'}
    print(f"Filters: {filter_params}")
    
    try:
        paths = client.find_shortest_paths(**config['params'])
        if paths:
            print(f"✓ Found {len(paths)} paths")
            for (gene_a, gene_b), info in list(paths.items())[:2]:  # Show first 2
                path_str = " → ".join(info['path'])
                print(f"  {gene_a} ↔ {gene_b}: {path_str} (hops={info['hops']}, distance={info['distance']:.1f})")
                print(f"    Edge scores: {info['scores']}")
        else:
            print("✗ No paths found with these filters")
    except Exception as e:
        print(f"✗ Error: {e}")

print("\n")

# ==============================================================================
# TEST 3: IntAct Single-Gene Query - Organism Filters
# ==============================================================================
print("="*80)
print("TEST 3: IntAct Single-Gene Query - Organism Filters")
print("="*80)
print()

intact_configs = [
    {
        "name": "Human only",
        "params": {
            "gene": "TP53",
            "species": "human",
            "top_n": 10,
            "organism_filter": "homo sapiens"
        }
    },
    {
        "name": "Human + Mouse",
        "params": {
            "gene": "TP53",
            "species": "human",
            "top_n": 10,
            "organism_filter": "homo sapiens,mus musculus"
        }
    },
    {
        "name": "High confidence only",
        "params": {
            "gene": "TP53",
            "species": "human",
            "top_n": 10,
            "miql_max_results": 20000
        }
    }
]

for i, config in enumerate(intact_configs, 1):
    print(f"\n{i}. {config['name']}")
    print("-" * 70)
    print(f"Parameters: {config['params']}")
    
    try:
        df = get_direct_neighbors(**config['params'])
        if not df.empty:
            print(f"✓ Found {len(df)} neighbors")
            print("\nTop 3 results:")
            display_cols = ['neighbor_name', 'best_miscore', 'n_interactions', 
                          'neighbor_organism', 'detection_methods', 'path']
            print(df[display_cols].head(3).to_string(index=False))
        else:
            print("✗ No neighbors found")
    except Exception as e:
        print(f"✗ Error: {e}")

print("\n")

# ==============================================================================
# TEST 4: IntAct Shortest Paths - Organism & Score Filters
# ==============================================================================
print("="*80)
print("TEST 4: IntAct Shortest Paths - Organism & Score Filters")
print("="*80)
print()

intact_genes = ["TP53", "MDM2", "EP300"]
intact_path_configs = [
    {
        "name": "Default filters",
        "params": {
            "gene_list": intact_genes,
            "species": "human",
            "max_distance": 3,
            "min_miscore": 0.4
        }
    },
    {
        "name": "High confidence only (0.7+)",
        "params": {
            "gene_list": intact_genes,
            "species": "human",
            "max_distance": 3,
            "min_miscore": 0.7
        }
    },
    {
        "name": "Very high confidence (0.9+)",
        "params": {
            "gene_list": intact_genes,
            "species": "human",
            "max_distance": 3,
            "min_miscore": 0.9
        }
    },
    {
        "name": "Human only filter",
        "params": {
            "gene_list": intact_genes,
            "species": "human",
            "max_distance": 3,
            "min_miscore": 0.4,
            "organism_filter": "homo sapiens"
        }
    }
]

for i, config in enumerate(intact_path_configs, 1):
    print(f"\n{i}. {config['name']}")
    print("-" * 70)
    filter_params = {k:v for k,v in config['params'].items() if k in ['min_miscore', 'max_distance', 'organism_filter']}
    print(f"Filters: {filter_params}")
    
    try:
        paths = find_shortest_paths_intact(**config['params'])
        if paths:
            print(f"✓ Found {len(paths)} paths")
            for (gene_a, gene_b), info in list(paths.items())[:2]:  # Show first 2
                path_str = " → ".join(info['path'])
                print(f"  {gene_a} ↔ {gene_b}: {path_str} (hops={info['hops']}, distance={info['distance']:.3f})")
                print(f"    MI-scores: {[f'{s:.3f}' for s in info['scores']]}")
        else:
            print("✗ No paths found with these filters")
    except Exception as e:
        print(f"✗ Error: {e}")

print("\n")

# ==============================================================================
# TEST 5: Network Type Comparison (STRING)
# ==============================================================================
print("="*80)
print("TEST 5: STRING Network Type Comparison")
print("="*80)
print()

network_types = ["functional", "physical"]

for net_type in network_types:
    print(f"\nNetwork type: {net_type}")
    print("-" * 70)
    
    try:
        df = get_string_neighbors(
            gene="TP53",
            species=9606,
            top_n=5,
            min_score=700,
            network_type=net_type
        )
        
        if not df.empty:
            print(f"✓ Found {len(df)} neighbors")
            print("\nTop 3 neighbors:")
            display_cols = ['preferred_name', 'combined_score', 'experimental_score', 'path']
            print(df[display_cols].head(3).to_string(index=False))
        else:
            print("✗ No neighbors found")
    except Exception as e:
        print(f"✗ Error: {e}")

print("\n")

# ==============================================================================
# TEST 6: Edge Case Testing
# ==============================================================================
print("="*80)
print("TEST 6: Edge Case Testing")
print("="*80)
print()

edge_cases = [
    {
        "name": "Very strict filters (should return few/no results)",
        "test": lambda: get_string_neighbors(
            gene="TP53",
            top_n=100,
            min_score=900,
            min_experimental=900,
            min_database=900,
            species=9606
        )
    },
    {
        "name": "Very permissive filters (should return many results)",
        "test": lambda: get_string_neighbors(
            gene="TP53",
            top_n=20,
            min_score=100,
            species=9606
        )
    },
    {
        "name": "Shortest paths with strict filters",
        "test": lambda: client.find_shortest_paths(
            gene_list=["TP53", "MDM2", "ATM"],
            species=9606,
            max_distance=2,
            min_combined_score=900,
            min_experimental_score=900
        )
    }
]

for i, case in enumerate(edge_cases, 1):
    print(f"\n{i}. {case['name']}")
    print("-" * 70)
    
    try:
        result = case['test']()
        if isinstance(result, dict):  # Shortest paths
            print(f"✓ Result: {len(result)} paths found")
        else:  # DataFrame
            print(f"✓ Result: {len(result)} neighbors found")
    except Exception as e:
        print(f"✗ Error: {e}")

print("\n")
print("="*80)
print("COMPREHENSIVE TESTING COMPLETE")
print("="*80)
print()
print("Summary:")
print("- Tested STRING single-gene queries with various score filters")
print("- Tested STRING shortest paths with multiple filter combinations")
print("- Tested IntAct single-gene queries with organism filters")
print("- Tested IntAct shortest paths with score and organism filters")
print("- Tested network type variations (functional vs physical)")
print("- Tested edge cases (strict and permissive filters)")
print()
