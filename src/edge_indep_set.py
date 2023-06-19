import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
from itertools import combinations


def solve(G):    
    # define model
    m = gp.Model()
    
    # define vars
    Y = m.addVars(G.edges, vtype=GRB.BINARY, name="y")
    #R = m.addVars(DG.nodes, vtype=GRB.BINARY, name="r")
    
    # define objective function
    m.setObjective(gp.quicksum(Y), sense=GRB.MAXIMIZE)
    
    # incoming arcs and root
    m.addConstrs( gp.quicksum(Y[u,v] for u in G.neighbors(v) if u<v) + gp.quicksum(Y[v,u] for u in G.neighbors(v) if v<u) <= 1 for v in G.nodes )
    
    # conflict edges
    for pair in combinations(G.edges, 2):
        pair0 = set(pair[0])
        pair1 = set(pair[1])
        intersect = pair0.intersection(pair1)
        if not intersect:
            u, v = pair[0]
            c, d = pair[1]
            if ((((u,c) in G.edges or (c,u) in G.edges) and ((v,d) in G.edges or (d,v) in G.edges)) or (((u,d) in G.edges or (d,u) in G.edges) and ((v,c) in G.edges) or (c,v) in G.edges)):
                #print("Heyyyy!")
                if u<v and c<d:
                    m.addConstr( Y[u,v] + Y[c,d] <= 1)
                if v<u and c<d:
                    m.addConstr( Y[v,u] + Y[c,d] <= 1)
                if v<u and d<c:
                    m.addConstr( Y[v,u] + Y[d,c] <= 1)
                if u<v and d<c:
                    m.addConstr( Y[u,v] + Y[d,c] <= 1)    
           
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        iuc_set = [ (u, v) for (u, v) in G.edges if Y[u,v].x > 0.5 ]
        edge_colors = [ "red" if (u,v) in iuc_set else "black" for (u,v) in G.edges ]
        nx.draw( G, with_labels=True, edge_color=edge_colors )
        obj = len(iuc_set)
        return obj
    else: print("There is an error in the vertex cover problem!")
