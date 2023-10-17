from gerrychain import Graph
import networkx as nx

from const import *
from test import solve_bc


def read_graph(level, state, level_name):
    return Graph.from_json("../data/"+level+"/json/"+state+"_"+level_name+".json")


if __name__ == "__main__":
    # read input graph g
    state = "AR"
    code = state_codes[state]
    num_dist = congressional_districts[state]
    level = "county"
    level_name = "counties"

    G = nx.karate_club_graph()
    solve_bc(G, 3, use_lower=True)

    print("The edge of g: ", len(G.edges), " Nodes: ", len(G.nodes))
 

