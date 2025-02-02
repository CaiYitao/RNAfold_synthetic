#!/usr/bin/env python3
"""
RNA Structural Analysis and Visualization Script

This script implements functions to:
  - Parse dot-bracket notation into a pair table.
  - Compute global metrics (e.g., number of base pairs, pairing percentage).
  - Extract motif distributions (e.g., stem lengths, hairpin loop sizes).
  - Build an RNA graph (using backbone and base pair edges) and compute graph metrics.
  - Visualize the computed statistics with histograms and boxplots.

Dependencies:
  - networkx
  - matplotlib
  - seaborn
  - numpy
"""

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def parse_dot_bracket(db):
    """
    Parse a dot-bracket string and return a pair table as a dictionary.
    For each paired nucleotide, store its partner index.
    Indexing is 0-based.
    """
    stack = []
    pairs = {}
    for i, ch in enumerate(db):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j
    return pairs

def compute_global_metrics(db):
    """
    Compute global structural metrics:
      - Number of base pairs (assumes well-formed structure).
      - Pairing percentage (ratio of paired nucleotides to total nucleotides).
      - Sequence length.
    """
    total = len(db)
    num_pairs = db.count('(')  # assuming well-formed (each '(' matches a ')')
    pairing_pct = (num_pairs * 2 / total * 100) if total > 0 else 0
    return {'num_pairs': num_pairs, 'pairing_pct': pairing_pct, 'length': total}

def extract_stem_lengths(db, pairs):
    """
    Extract stem lengths from a dot-bracket structure.
    
    A stem is defined as a contiguous set of base pairs.
    This naive implementation walks along the string and counts runs of '('
    (assuming that a contiguous run in the dot-bracket string corresponds to a stem).
    
    Returns:
      A list of stem lengths.
    """
    stem_lengths = []
    i = 0
    n = len(db)
    while i < n:
        if db[i] == '(':
            length = 0
            # count consecutive '(' that form base pairs
            while i < n and db[i] == '(' and i in pairs:
                length += 1
                i += 1
            stem_lengths.append(length)
        else:
            i += 1
    return stem_lengths

def extract_hairpin_loop_sizes(db, pairs):
    """
    Extract hairpin loop sizes from a dot-bracket structure.
    
    Here a hairpin loop is defined as the unpaired region between a pair (i, j)
    that does not contain any further base pairs.
    
    Returns:
      A list of hairpin loop sizes.
    """
    hairpins = []
    for i in range(len(db)):
        if db[i] == '(':
            j = pairs.get(i, None)
            if j is None or i >= j:
                continue
            # If there are no paired bases between i+1 and j, assume hairpin loop.
            if all(db[k] == '.' for k in range(i+1, j)):
                loop_size = j - i - 1
                hairpins.append(loop_size)
    return hairpins

def create_rna_graph(db, pairs):
    """
    Build an undirected graph for the RNA structure.
    
    - Each nucleotide is a node (indexed 0..n-1).
    - Backbone edges connect consecutive nucleotides.
    - Base pair edges are added from the pair table.
    
    Returns:
      A networkx Graph object.
    """
    G = nx.Graph()
    n = len(db)
    G.add_nodes_from(range(n))
    # Add backbone edges
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone')
    # Add base pair edges (only add one per pair)
    for i, j in pairs.items():
        if i < j:
            G.add_edge(i, j, type='basepair')
    return G

def compute_graph_metrics(G):
    """
    Compute graph-based metrics on the RNA structure graph:
      - Degree distribution (list of degrees per node)
      - Average degree
      - Average clustering coefficient
      - Average shortest path length (if graph is connected)
    
    Returns:
      A dictionary with these metrics.
    """
    degrees = [deg for node, deg in G.degree()]
    avg_degree = np.mean(degrees)
    clustering = nx.clustering(G)
    avg_clustering = np.mean(list(clustering.values()))
    try:
        avg_shortest_path = nx.average_shortest_path_length(G)
    except nx.NetworkXError:
        avg_shortest_path = None  # Graph not connected
    return {'degrees': degrees,
            'avg_degree': avg_degree,
            'avg_clustering': avg_clustering,
            'avg_shortest_path': avg_shortest_path}

def plot_histogram(data, title, xlabel, ylabel):
    """
    Plot a histogram with KDE using seaborn.
    """
    plt.figure(figsize=(8, 6))
    sns.histplot(data, kde=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.show()
    plt.savefig('histogram.png')

def plot_boxplot(data_dict, title, ylabel):
    """
    Plot a boxplot for different groups.
    
    data_dict: keys are group labels, values are lists of numeric data.
    """
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=list(data_dict.values()))
    plt.xticks(ticks=range(len(data_dict)), labels=list(data_dict.keys()))
    plt.title(title)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.show()
    plt.savefig('boxplot.png')

def main():
    # For demonstration, we create three example dot-bracket structures.
    synthetic_db = "..((......)).."
    bpRNA_db = "....((((....))))........"
    RFAM_db = "(((....)))....(((...)))"
    
    datasets = {'Synthetic': synthetic_db,
                'bpRNA': bpRNA_db,
                'RFAM': RFAM_db}
    
    all_global_metrics = {}
    all_stem_lengths = {}
    all_hairpin_sizes = {}
    all_graph_metrics = {}
    
    for name, db in datasets.items():
        pairs = parse_dot_bracket(db)
        g_metrics = compute_global_metrics(db)
        stem_lengths = extract_stem_lengths(db, pairs)
        hairpin_sizes = extract_hairpin_loop_sizes(db, pairs)
        G = create_rna_graph(db, pairs)
        graph_metrics = compute_graph_metrics(G)
        
        all_global_metrics[name] = g_metrics
        all_stem_lengths[name] = stem_lengths
        all_hairpin_sizes[name] = hairpin_sizes
        all_graph_metrics[name] = graph_metrics
        
        print(f"Dataset: {name}")
        print("  Global Metrics:", g_metrics)
        print("  Stem Lengths:", stem_lengths)
        print("  Hairpin Loop Sizes:", hairpin_sizes)
        print("  Graph Metrics:", graph_metrics)
        print("="*50)
    
    # Visualization: Stem length distribution for each dataset.
    plt.figure(figsize=(10, 6))
    for name, lengths in all_stem_lengths.items():
        sns.histplot(lengths, bins=10, kde=True, label=name, element='step',
                     stat="density", common_norm=False, alpha=0.5)
    plt.xlabel("Stem Length")
    plt.ylabel("Density")
    plt.title("Stem Length Distribution Across Datasets")
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('stem_length_distribution.png')
    
    # Visualization: Hairpin loop size distribution.
    plt.figure(figsize=(10, 6))
    for name, sizes in all_hairpin_sizes.items():
        sns.histplot(sizes, bins=10, kde=True, label=name, element='step',
                     stat="density", common_norm=False, alpha=0.5)
    plt.xlabel("Hairpin Loop Size")
    plt.ylabel("Density")
    plt.title("Hairpin Loop Size Distribution Across Datasets")
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('hairpin_loop_size_distribution.png')
    
    # Visualization: Boxplot of pairing percentage.
    pairing_pct_data = {name: [metrics['pairing_pct']] for name, metrics in all_global_metrics.items()}
    plot_boxplot(pairing_pct_data, "Pairing Percentage by Dataset", "Pairing Percentage (%)")
    
    # Visualization: Boxplot of node degree distribution.
    degree_data = {name: met['degrees'] for name, met in all_graph_metrics.items()}
    plot_boxplot(degree_data, "Node Degree Distribution by Dataset", "Degree")
    
if __name__ == "__main__":
    main()
