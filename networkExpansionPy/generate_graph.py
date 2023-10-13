# generate_graph.py ---
#
# Filename: generate_graph.py
# Author: Zach Maas
# Created: Fri Oct 13 11:44:36 2023 (-0600)
#
# Code to generate graph representation of an expansion output in
# order to make the output usable and useful in Cytoscape.
#

# Code:

import re
import networkx as nx
import itertools as it

# Required data files
in_file = "../networkExpansionPy/assets/iqbio/kegg_rxns.tab"
out_file = "../networkExpansionPy/assets/iqbio/kegg_rxns.dot"
out_gmlfile = "../networkExpansionPy/assets/iqbio/kegg_rxns.gml"

# Read in data for name conversion globally
with open(
    "../networkExpansionPy/assets/iqbio/compounds.csv"
) as cpds_file_handle:
    lines = list(cpds_file_handle.readlines())
    cpd_dict = {
        cpd: name
        for cpd, name in [
            line.strip("cpd:").strip().split(maxsplit=1) for line in lines
        ]
    }


def get_rxn_pairs(in_file=in_file):
    all_pairs = set()
    with open(in_file, "r") as rxn_file_handle:
        lines = list(rxn_file_handle.readlines())
        # Grab only the equations from this list
        rxns = [
            line.strip("ENTRY").strip().split()[0]
            for line in lines
            if "ENTRY" in line
        ]
        eqns = [
            line.strip("EQUATION").strip()
            for line in lines
            if "EQUATION" in line
        ]
        if len(rxns) != len(eqns):
            raise Exception(
                "Unequal number of reactions and equations in data file."
            )
        regex = r"C\d{5}"
        rxn_pairs = {}
        for rxn, eqn in zip(rxns, eqns):
            # Break out and find reactions and products
            reactant_str, product_str = eqn.split("<=>")
            reactants = re.findall(regex, reactant_str)
            products = re.findall(regex, product_str)
            # Generate pairwise combinations as edges
            pairs = list(it.product(reactants, products))
            rxn_pairs[rxn] = pairs
            all_pairs.update(pairs)
    return (rxn_pairs, all_pairs)


def gen_filtered_graph(all_pairs):
    # Generate Graph
    G = nx.Graph()
    G.add_edges_from(all_pairs)

    # Filter by degree > 5 (arbitrary)
    deg = G.degree()
    to_keep = [n for n, d in deg if d > 20]
    sG = G.subgraph(to_keep)

    nx.drawing.nx_agraph.write_dot(sG, out_file)
    nx.write_graphml(G, out_gmlfile)


def rename_nodes(nodes):
    """Rename nodes to their real names (not CXXXXX)"""
    return [cpd_dict[x] for x in nodes]


def rxns_to_graph(rxns):
    # Get dict of product/reactant pairs associated with each reaction
    rxn_pairs, _ = get_rxn_pairs()

    # Convert reactions to compounds
    edges = set()
    # Look up every product/reactant pair for all reactions in current step
    for rxn in rxns:
        try:
            edges.update(rxn_pairs[rxn[0]])
        except KeyError:
            print(f"Reaction {rxn[0]} not found")

    # Convert compounds to their real names (not CXXXXX)
    edges_renamed = [(cpd_dict[x], cpd_dict[y]) for x, y in edges]
    # nodes_renamed = [cpd_dict[x] for x in ne_cpds]
    # nodes_dct = {nodes_renamed[i]: i for i in range(0, len(nodes_renamed))}
    return nx.from_edgelist(edges_renamed)


#
# generate_graph.py ends here
