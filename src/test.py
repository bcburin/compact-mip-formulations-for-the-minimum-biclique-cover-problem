import vertex_cover
import heuristic
import biclique
import edge_indep
from const import *

from gerrychain import (Graph)

from src.util import build_graph_from_file


# def reorder_heuristic_ind()

def solve_bc(G, form, use_lower=False):
    print("Number of edges: ", len(G.edges))
    
    print("Number of nodes: ", len(G.nodes))
       
    
    #nx.draw(g)
    
    # find a vertex cover of g which is an upper bound for the minimum biclique covering problem (MBCP)
    vertex_cover_number, vertex_cover_set = vertex_cover.solve(G)
    
    # find an independent set of g which is a lower bound for the minimum biclique covering problem (MBCP)
    #indep_number = independent_set.solve(g)
    # indep_number = directed_stars.solve(g)
    
    # print("iuc_number: ", indep_number)
    if (use_lower):
        indep_edges = edge_indep.solve(G)
    else:
        indep_edges = []
    print("indep_number: ", len(indep_edges))
    
    # print("vertex cover number: ", vertex_cover_number)
#    print("vertex cover set: ", vertex_cover_set)
#    return
    
    # find a heuristic solution for the MBCP based on the vertex cover set
    heuristic_sol = heuristic.find_cover(G, vertex_cover_set)
    
    # print("heuristic solution: ", heuristic_sol)
    
    if (form == 1):
        sol = biclique.get_vertex_bc_from_edge(G, heuristic_sol)
        print(sol)
        biclique_cover = biclique.solve_v(G, sol, indep_edges=indep_edges)
        print("Is it biclique cover? ", biclique.check_v_biclique_cover(G, biclique_cover))
    elif (form == 2):
        biclique_cover = biclique.solve(G, heuristic_sol, indep_edges=indep_edges)
        print("Is it biclique cover? ", biclique.check_v_biclique_cover(G, biclique_cover))
    elif (form == 3):
        biclique_cover = biclique.solve_recursive(G, heuristic_sol, indep_edges=indep_edges)
        print("Is it biclique cover? ", biclique.check_v_biclique_cover(G, biclique_cover))


def readGraph(level, state, level_name):
    return Graph.from_json("../data/"+level+"/json/"+state+"_"+level_name+".json")

# read input graph g
state = "AR"
code = state_codes[state]
num_dist = congressional_districts[state]
level = "county"
level_name = "counties"


# g = nx.karate_club_graph()
# g = readGraph(level, state, level_name)


if __name__ == "__main__":
    G = build_graph_from_file("../graph/erdos_renyi_100_8_0.txt")
    solve_bc(G, 1, use_lower=True)

    # edge_indep.solve(g)
    # vertex_cover_number, vertex_cover_set = vertex_cover.solve(g)
    # print("vertex_cover: ", vertex_cover_number)

    # directed_stars.solve(g)

    print("The edge of g: ", len(G.edges), " Nodes: ", len(G.nodes))
 

