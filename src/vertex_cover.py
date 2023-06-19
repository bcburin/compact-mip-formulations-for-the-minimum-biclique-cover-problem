import gurobipy as gp
from gurobipy import GRB 

def solve(G):
    # define model
    m = gp.Model()

    # define vars
    X = m.addVars(G.nodes, vtype=GRB.BINARY, name="x")
    
    # define objective function
    m.setObjective(gp.quicksum(X), sense=GRB.MINIMIZE)
    
    # add covering constraints
    m.addConstrs(X[u] + X[v] >= 1 for u, v in G.edges)
    
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        cover_set = [ v for v in G.nodes if X[v].x > 0.5 ]
        return m.objVal, cover_set
    else: print("There is an error in the vertex cover problem!")
    
    