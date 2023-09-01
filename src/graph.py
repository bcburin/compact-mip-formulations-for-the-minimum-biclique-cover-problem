import networkx as nx
import numpy as np


def create_conflict_graph(node_sets):
    G = nx.Graph()
    for node_set in node_sets:
        for node in node_set:
            if not G.has_node(node):
                G.add_node(node)
    for node_set in node_sets:
        for n1 in node_set:
            for n2 in node_set:
                if (n1 < n2) and (not G.has_edge(n1, n2)):
                    G.add_edge(n1, n2)
    
    GC = nx.complement(G)
    return GC


def create_SOSk_node_sets(bps, k):
    if (bps <= 0) or (k <= 0):
        print("Both bps and k should be larger than 0.")
        raise ValueError
    if (bps <= k):
        print("Not able to create SOSk with less breakpoints than k")
        raise ValueError
    node_sets = []
    for i in range(bps - k + 1):
        node_set = [j for j in range(i+1, i+k+1)]
        node_sets.append(node_set)
    return node_sets
    

# bps: the number of breakpoints
# k: the k value of SOSk
def create_SOSk_conflict_graph(bps, k):
    node_sets = create_SOSk_node_sets(bps, k)
    print(node_sets)
    return create_conflict_graph(node_sets)

def create_SOSk_base_graph(bps, k):
    node_sets = create_SOSk_node_sets(bps, k)
    G = nx.Graph()
    for i in range(len(node_sets)):
        G.add_node(i)
        G.nodes[i]["nodes"] = node_sets[i]
    
    for i in range(len(node_sets)-1):
        G.add_edge(i, i+1)
        G.edges[i, i+1]["nodes"] = np.intersect1d(node_sets[i], node_sets[i+1]).tolist()
    return G
        