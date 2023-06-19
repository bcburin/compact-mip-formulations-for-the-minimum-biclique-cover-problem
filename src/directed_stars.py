import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
import math
from itertools import combinations


def solve(G):
    DG = nx.DiGraph(G) # bidirected version of G
    
    # define model
    m = gp.Model()
    
    # define vars
    Y = m.addVars(DG.edges, vtype=GRB.BINARY, name="y")
    R = m.addVars(DG.nodes, vtype=GRB.BINARY, name="r")
    
    # define objective function
    m.setObjective(gp.quicksum(R), sense=GRB.MAXIMIZE)
    
    # incoming arcs and root
    m.addConstrs( gp.quicksum(Y[u,v] for u in DG.neighbors(v)) + R[v] <= 1 for v in DG.nodes )
    
    # conflict roots
    m.addConstrs( R[u] + R[v] <= 1 for (u,v) in G.edges )
    
    # outgoing from root
    m.addConstrs( Y[v,u] <= R[v] for v in DG.nodes for u in DG.neighbors(v))
    
    # outgoing from root
    m.addConstrs( R[v] <= gp.quicksum(Y[v,u] for u in DG.neighbors(v)) for v in DG.nodes)
    
    # conflict edges
    for pair in combinations(DG.edges, 2):
        pair0 = set(pair[0])
        pair1 = set(pair[1])
        intersect = pair0.intersection(pair1)
        if not intersect:
            u, v = pair[0]
            c, d = pair[1]
            if ((u,c) in DG.edges and (v,d) in DG.edges) or ((u,d) in DG.edges and (v,c) in DG.edges):
                #print("Heyyyy!")
                m.addConstr( Y[v,u] + Y[u,v] + Y[c,d] + Y[d,c] <= 1)
                
           
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        iuc_set = [ (u, v) for (u, v) in DG.edges if Y[u,v].x > 0.5 ]
        edge_colors = [ "red" if (u,v) in iuc_set else "black" for (u,v) in DG.edges ]
        nx.draw( DG, with_labels=True, edge_color=edge_colors )
        
        roots = [ v for v in DG.nodes if R[v].x > 0.5 ]
        #H = DG.edge_subgraph(iuc_set)
        print("roots: ", roots)
        obj = len(roots)
        #for component in nx.connected_components(H):
         #   print("comp: ", component)
          #  obj = obj + math.ceil(math.log(len(component),2))
        return obj
    else: print("There is an error in the vertex cover problem!")