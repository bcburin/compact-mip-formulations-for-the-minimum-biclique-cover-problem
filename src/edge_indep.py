import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
import math


def solve(G):
    # define model
    m = gp.Model()
    edgesList = list(G.edges)
    eLen = len(edgesList)
    
    # define vars
    Y = m.addVars(G.edges, vtype=GRB.BINARY, name="y")
    for u in G.nodes:
        m.addConstr(gp.quicksum(Y[v, w] for v, w in G.edges(u) if (v, w) in Y) + gp.quicksum(Y[w,v] for v, w in G.edges(u) if (w, v) in Y) <= 1)
    
    for i in range(eLen):
        for j in range(i+1, eLen):
            v1, v2 = edgesList[i][0], edgesList[i][1]
            v3, v4 = edgesList[j][0], edgesList[j][1]
            if (v1 == v3) or (v1 == v4) or (v2 == v3) or (v2 == v4):
                continue
            
            if (G.has_edge(v1, v3) or G.has_edge(v3, v1)) and (G.has_edge(v2, v4) or G.has_edge(v4, v2)):
                indepConstr = m.addConstr(Y[v1, v2] + Y[v3, v4] <= 1)
                # indepConstr.Lazy = 2
            
            if (G.has_edge(v1, v4) or G.has_edge(v4, v1)) and (G.has_edge(v2, v3) or G.has_edge(v3, v2)):
                indepConstr = m.addConstr(Y[v1, v2] + Y[v3, v4] <= 1)
                # indepConstr.Lazy = 2
    
    m.setObjective(gp.quicksum(Y), sense=GRB.MAXIMIZE)
    
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        iuc_set = [ (u, v) for (u, v) in G.edges if Y[u,v].x > 0.5 ]
        edge_colors = [ "red" if (u,v) in iuc_set else "black" for (u,v) in G.edges ]
        nx.draw(G, with_labels=True, edge_color=edge_colors )
        obj = len(iuc_set)
        #for component in nx.connected_components(H):
         #   print("comp: ", component)
          #  obj = obj + math.ceil(math.log(len(component),2))
        return iuc_set
    else: print("There is an error in the vertex cover problem!")
