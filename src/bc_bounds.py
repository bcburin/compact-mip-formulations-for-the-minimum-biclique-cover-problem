import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
import numpy as np

import cvxpy as cp
import graph

import biclique

# Support methods are "match", "lovasz", "clique"
def find_bc_lower_bound(G, method="match"):
    if (method == "match"):
        m_len = len(nx.max_weight_matching(G))
        E_len = len(G.edges)
        return np.ceil(m_len ** 2 / E_len)
    elif (method == "lovasz"):
        G_c = nx.complement(G)
        lov_num = np.round(compute_lovasz_number(G_c))
        return np.ceil(np.log2(lov_num))
    elif (method == "clique"):
        return np.ceil(np.log2(max_clique(G)))

        
def compute_lovasz_number(G):
    A = nx.adjacency_matrix(G)
    n = A.shape[0]
    C = np.ones((n, n))
    I = np.identity(n)
    X = cp.Variable((n,n), symmetric=True)
    # The operator >> denotes matrix inequality.
    constraints = [X >> 0]
    constraints += [cp.trace(I @ X) == 1]
    rowInd, colInd = A.nonzero()
    constraints += [X[rowInd[i], colInd[i]] == 0 for i in range(len(rowInd))]
    prob = cp.Problem(cp.Maximize(cp.trace(C @ X)),
                  constraints)
    prob.solve()

#    # Print result.
    print("The optimal value is", prob.value)
#    print("A solution X is")
#    print(X.value)
    return prob.value

def max_clique(G):
    # define model
    m = gp.Model()

    # define vars
    X = m.addVars(G.nodes, vtype=GRB.BINARY, name="x")
    
    # define objective function
    m.setObjective(gp.quicksum(X), sense=GRB.MAXIMIZE)
    
    Gc = nx.complement(G)
    
    # add covering constraints
    m.addConstrs(X[u] + X[v] <= 1 for u, v in Gc.edges)
    
    # set a one-minute time limit
    m.Params.TimeLimit = 60
    
    # optimize
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        return m.objVal
    else: print("There is an error in the maximum clique problem!")

# Support methods are "number", "vertex"
def find_bc_upper_bound(G, method="number"):
    if (method == "number"):
        n = len(G.nodes)
        return n + 1 - np.floor(np.log2(n))
    elif (method == "vertex"):
        return vertex_cover(G)


def vertex_cover(G):
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
        return m.objVal
    else: print("There is an error in the vertex cover problem!")

# read input graph G
#b = 2
#k = 3
#N = 12 + k - 1
#
#G = graph.create_SOSk_conflict_graph(N, k)
#G = nx.complete_graph(10) 
#print("SOSk(N): ", k, N, "low_bound:", low_bound)