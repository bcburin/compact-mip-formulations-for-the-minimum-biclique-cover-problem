import gurobipy as gp
from gurobipy import GRB 

import networkx as nx


def solve(G):    # define model
    m = gp.Model()

    # find power graph g^2
    H = nx.power(G, 2)

    # define vars
    X = m.addVars(H.nodes, vtype=GRB.BINARY, name="x")
    
    # define objective function
    m.setObjective(gp.quicksum(X), sense=GRB.MAXIMIZE)
    
    # add covering constraints
    m.addConstrs(X[u] + X[v] <= 1 for u, v in H.edges)
    
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        # cover_set = [ v for v in g.nodes if X[v].x > 0.5 ]
        return m.objVal
    else: print("There is an error in the independent set problem!")
