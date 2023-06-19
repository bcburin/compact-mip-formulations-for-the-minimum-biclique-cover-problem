import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
import math


def solve(G):
    # define model
    m = gp.Model()
    
    # define vars
    Y = m.addVars(G.edges, vtype=GRB.BINARY, name="y")
    
    # define objective function
    m.setObjective(gp.quicksum(Y), sense=GRB.MAXIMIZE)
    
    # add covering constraints
    for v in G.nodes:
        for i in G.neighbors(v):
            for j in G.neighbors(v):
                if i<j and (i, j) not in G.edges:
                    if v<i and v<j:
                        m.addConstr(Y[v,i] + Y[v,j] <= 1)
                    elif v<i and j<v:    
                        m.addConstr(Y[v,i] + Y[j,v] <= 1)
                    elif i<v and j<v:    
                        m.addConstr(Y[i,v] + Y[j,v] <= 1)
                    elif i<v and v<j:    
                        m.addConstr(Y[i,v] + Y[v,j] <= 1)    
                
           
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        iuc_set = [ (u, v) for (u, v) in G.edges if Y[u,v].x > 0.5 ]
        edge_colors = [ "red" if (u,v) in iuc_set else "black" for (u,v) in G.edges ]
        nx.draw( G, with_labels=True, edge_color=edge_colors )
        H = G.edge_subgraph(iuc_set)
        obj = 0
        for component in nx.connected_components(H):
            print("comp: ", component)
            obj = obj + math.ceil(math.log(len(component),2))
        return obj
    else: print("There is an error in the vertex cover problem!")