import gurobipy as gp
from gurobipy import GRB 

import time

import networkx as nx
from itertools import combinations


def get_v_bicliques(G, X, k):
    bcs = []
    for b in range(k):
        A = []
        B = []
        for u in G.nodes:
            if (X[u, b, 0].x > 0.5):
                A.append(u)
            if (X[u, b, 1].x > 0.5):
                B.append(u)
        if A != [] and B != []:
            bcs.append((A, B))
    return bcs


def get_vertex_bc_from_edge(G, edge_set):
    node_set = []
    for b in range(len(edge_set)):
        # Initialize two disjoint sets.
        A, B = [], []
        d = {}
        for u, v in edge_set[b]:
            d[(u, v)] = 1
            A = [u]
        if len(A) == 0:
            continue
        
        # If there is an edge between v and A[0], then v must be in B.
        for v in G.neighbors(A[0]):
            if ((A[0], v) in d) or ((v, A[0]) in d):
                B.append(v)
        
        if len(B) == 0:
            continue
        
        A = []

        for u in G.neighbors(B[0]):
            if ((B[0], u) in d) or ((u, B[0]) in d):
                A.append(u)
        node_set.append((A, B))
    return node_set


# Get biclique covers from the formulations with edge.
def get_v_bicliques_from_e(G, Y, k):
    bcs = []
    for b in range(k):
        # Initialize two disjoint sets.
        A, B = [], []
        d = {}
        for u, v in G.edges:
            if Y[u, v, b].x > 0.5:
                d[(u, v)] = 1
                A = [u]
        if len(A) == 0:
            continue
        
        # If there is an edge between v and A[0], then v must be in B.
        for v in G.neighbors(A[0]):
            if ((A[0], v) in d) or ((v, A[0]) in d):
                B.append(v)
        
        if len(B) == 0:
            continue
        
        A = []

        for u in G.neighbors(B[0]):
            if ((B[0], u) in d) or ((u, B[0]) in d):
                A.append(u)
        bcs.append((A, B))
    return bcs


# Each biclique is written as two disjoint set of nodes.
# Check whether it is a biclique cover.
def check_v_biclique_cover(G, bcs):
#    print("The biclique cover: ", bcs)
    for u,v in G.edges:
        G.edges[u,v]["visited"] = False
    
    for bc in bcs:
        A, B = bc[0], bc[1]
        for u in A:
            for v in B:
                if u == v:
                    return False
                if G.has_edge(u,v):
                    G.edges[u, v]["visited"] = True
                elif G.has_edge(v, u):
                    G.edges[v, u]["visited"] = True
                else:
                    return False
    
    for u, v in G.edges:
        if not G.edges[u, v]["visited"]:
            return False
    return True

def build_bottom_up_model_v(G, k, indep_edges, maximal_con):
    # define model
    m = gp.Model()
        
    # set a one-minute time limit
    m.Params.TimeLimit = 300
    # m.Params.Cuts = 0
    DG = G.to_directed()
    GC = nx.complement(DG)
    print("k: ", k)
        
    # print("edge_set: ", g.edges)
        
    # define vars
    Z = m.addVars(range(k), vtype=GRB.BINARY, name="z")
        
    X = m.addVars(G.nodes, range(k), range(2), vtype=GRB.BINARY, name="x")
    Y = m.addVars(DG.edges, range(k), lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
    if (maximal_con):
        add_maximal_con_v(m, G, DG, Y, k)

    # fix z vars
    for b in range(k):
        Z[b].lb = 1
    fix_indep_edges_v(indep_edges, X, Y, Z)
        
    # define objective function
    m.setObjective(gp.quicksum(Z), sense=GRB.MINIMIZE)
    # symmetry-break constraints
    m.addConstrs(Z[b] >= Z[b+1] for b in range(k-1))
    add_base_constr_v(m, G, DG, X, Y, Z, k)
    return m, X


def solve_recursive(G, heuristic_sol, indep_edges = [], maximal_con=True, form='v'):
    # no_solution = True
    start = time.time()
    k = len(indep_edges)
    while(True):
        m, X = build_bottom_up_model_v(G, k, indep_edges, maximal_con)
        m.optimize()

        if m.status == GRB.OPTIMAL:
            print("Total time: ", time.time() - start)
            return get_v_bicliques(G, X, k)
        else:
            k += 1


def add_maximal_con_v(m, G, DG, Y, k):
    maximal_cliques = nx.find_cliques(G)
    for clique in list(maximal_cliques):
        if (len(clique) % 2 == 0):
            clique_num = len(clique) ** 2 / 4
        else:
            clique_num = (len(clique) + 1) * (len(clique) - 1) / 4
        cliqueConstr = m.addConstrs(gp.quicksum(Y[u,v,b] for u,v in DG.subgraph(clique).edges) <= clique_num for b in range(k))
        cliqueConstr.Lazy = 2

def add_initial_v(heuristic_sol, X, Y, Z):
    for b, cover in enumerate(heuristic_sol):
        # print ("b, ", b)
        Z[b].start = 1
        for u in cover[0]:
            X[u,b,0].start = 1
        for v in cover[1]:
            X[v,b,1].start = 1
        for u in cover[0]:
            for v in cover[1]:
                Y[u, v, b].start = 1.0

def fix_indep_edges_v(indep_edges, X, Y, Z):
    for b, (u,v) in enumerate(indep_edges):
        X[u, b, 0].lb = 1.0
        X[v, b, 1].lb = 1.0
        Y[u, v, b].lb = 1.0
        Y[v, u, b].ub = 0.0
    return

def add_base_constr_v(m, G, DG, X, Y, Z, k):
    GC = nx.complement(DG)
    
    # coupling constraints
    m.addConstrs(X[v,b,i] <= Z[b] for v in G.nodes for b in range(k) for i in range(2))
    
    # covering constraints
    m.addConstrs(gp.quicksum(Y[u,v,b] + Y[v,u,b] for b in range(k)) >= 1 for u, v in G.edges)
    
    # node assignment constraints
    m.addConstrs(X[v,b,0]+X[v,b,1]<=Z[b] for v in G.nodes for b in range(k))

    # node and edge constraints
    for u, v in DG.edges:
        m.addConstrs(Y[u,v,b] <= X[u,b,0] for b in range(k))
        m.addConstrs(Y[u,v,b] <= X[v,b,1] for b in range(k))
        m.addConstrs(Y[u,v,b] >= X[u,b,0] + X[v,b,1] - Z[b] for b in range(k))
        
    for u, v in GC.edges:
        m.addConstrs(X[u,b,0] + X[v,b,1] <= Z[b] for b in range(k))
        # m.addConstrs(X[u,b,1] + X[v,b,0] <= Z[b] for b in range(k))

    return

# Dual version implementation of g
def solve_v(G, heuristic_sol, indep_edges = [], maximal_con=False):
    # define model
    m = gp.Model()
    
    # set a one-minute time limit
    m.Params.TimeLimit = 300
    # m.Params.Cuts = 0
    # maximum number of biclique covers
    k = len(heuristic_sol)
    DG = G.to_directed()
    GC = nx.complement(DG)
    print("k: ", k)
    
    # print("edge_set: ", g.edges)
    
    # define vars
    Z = m.addVars(range(k), vtype=GRB.BINARY, name="z")
    
    X = m.addVars(G.nodes, range(k), range(2), vtype=GRB.BINARY, name="x")
    Y = m.addVars(DG.edges, range(k), lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
    if (maximal_con):
        add_maximal_con_v(m, G, DG, Y, k)

    # warm start X vars
    add_initial_v(heuristic_sol, X, Y, Z)
    
    # fix z vars
    for b in range(k):
        Z[b].lb = 1.0
    fix_indep_edges_v(indep_edges, X, Y, Z)
    
    # define objective function
    m.setObjective(gp.quicksum(Z), sense=GRB.MINIMIZE)
    # symmetry-break constraints
    m.addConstrs(Z[b] >= Z[b+1] for b in range(k-1))
    add_base_constr_v(m, G, DG, X, Y, Z, k)
    m.optimize()
    return get_v_bicliques(G, X, k)

# Return an edge heuristic_solution
def heuristic_greedy(H):
    G = H.copy()
    for u,v in G.edges:
        G.edges[u,v]["w"] = 1.0
    DG = G.to_directed()
    bcs = []
    while True:
        m = gp.Model()
        # set a one-minute time limit
        m.Params.TimeLimit = 60
        # define vars
        X = m.addVars(G.nodes, range(2), vtype=GRB.BINARY, name="x")
        Y = m.addVars(DG.edges, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
        GC = nx.complement(DG)
        
        # node assignment constraints
        m.addConstrs(X[v,0]+X[v,1]<= 1 for v in G.nodes)

        # node and edge constraints
        m.addConstrs(Y[u,v] <= X[u,0] for u, v in DG.edges)
        m.addConstrs(Y[u,v] <= X[v,1] for u, v in DG.edges)
        m.addConstrs(Y[u,v] >= X[u,0] + X[v,1] - 1 for u, v in DG.edges)
            
        m.addConstrs(X[u,0] + X[v,1] <= 1 for u, v in GC.edges)
        m.setObjective(gp.quicksum((Y[u,v]+Y[v,u]) * G.edges[u,v]["w"] for u,v in G.edges), sense=GRB.MAXIMIZE)
        m.optimize()
        if m.status == GRB.OPTIMAL or m.status == GRB.TIMELIMIT:
            if m.getObjective().getValue() == 0:
                return bcs
            bc = [(u, v) for u,v in G.edges if Y[u,v].x+Y[v,u].x > 0.5]
            for u,v in G.edges:
                if Y[u,v].x+Y[v,u].x > 0.5: G.edges[u,v]["w"] = 0.0
            bcs.append(bc)
        else:
            print("Can't get a heuristic sol.")
            return []

    
def solve(G, heuristic_sol, indep_edges=[]):
    # define model
    m = gp.Model()
    
    # set a one-minute time limit
    m.Params.TimeLimit = 300
    
    # maximum number of biclique covers
    k = len(heuristic_sol)
    
    print("k: ", k)
    
    # print("edge_set: ", g.edges)

    # define vars
    Z = m.addVars(range(k), vtype=GRB.BINARY, name="z")
    
    Y = m.addVars(G.edges, range(k), vtype=GRB.BINARY, name="y")
    
    # warm start Y vars
    for b, cover in enumerate(heuristic_sol):
        # print ("b, ", b)
        for (u, v) in cover:
            # print("(u, v): ", (u,v))
            Y[u,v,b].start = 1.0
            
    # fix z vars
    for b, (u,v) in enumerate(indep_edges):
        Z[b].lb = 1
        Y[u, v, b].lb = 1
    
    # define objective function
    m.setObjective(gp.quicksum(Z), sense=GRB.MINIMIZE)
    
    # symmetry-break constraints
    m.addConstrs(Z[b] >= Z[b+1] for b in range(k-1))
    
    # coupling constraints
    m.addConstrs(Y[u,v,b] <= Z[b] for u, v in G.edges for b in range(k))
    
    # covering constraints
    m.addConstrs(gp.quicksum(Y[u,v,b] for b in range(k)) >= 1 for u, v in G.edges)
    
    # cycle-3 elimination constraints
    list_of_C_3 = []
    for u,v in G.edges:
        for j in nx.common_neighbors(G, u, v):
            cycle_list = [(u,v)]
            cycle_list.append((u,j)) if u<j else cycle_list.append((j,u))
            cycle_list.append((v,j)) if v<j else cycle_list.append((j,v))
            if cycle_list not in list_of_C_3:
                list_of_C_3.append(cycle_list)
                
    # print("list of cycles", list_of_C_3)
    
    for cycle in list_of_C_3:
        for b in range(k):
            cycleConstr = m.addConstr(gp.quicksum(Y[u,v,b] for u,v in cycle) <= 2)
            cycleConstr.Lazy = 2
            
    
    # edge conflict constraints
    for pair in combinations(G.edges, 2):
        pair0 = set(pair[0])
        pair1 = set(pair[1])
        intersect = pair0.intersection(pair1)
        if not intersect:
            u, v = pair[0]
            c, d = pair[1]
            for b in range(k):
                left_expr = gp.LinExpr(Y[u,v,b] + Y[c,d,b])
                right_expr = gp.LinExpr()
                if (c, u) in G.edges and c<u:
                    right_expr += gp.LinExpr(Y[c,u,b])
                elif (c, u) in G.edges:
                    right_expr += gp.LinExpr(Y[u,c,b])
                if (c, v) in G.edges and c<v:    
                    right_expr += gp.LinExpr(Y[c,v,b])
                elif (c, v) in G.edges:
                    right_expr += gp.LinExpr(Y[v,c,b])
                crossConstr = m.addConstr(left_expr <= 1 + right_expr) 
                crossConstr.Lazy = 2
                right_expr = gp.LinExpr()
                if (d, u) in G.edges and d<u:
                    right_expr += gp.LinExpr(Y[d,u,b])
                elif (d, u) in G.edges:
                    right_expr += gp.LinExpr(Y[u,d,b])
                if (d, v) in G.edges and d < v:    
                    right_expr += gp.LinExpr(Y[d,v,b])
                elif (d, v) in G.edges:
                    right_expr += gp.LinExpr(Y[v,d,b])
                crossConstr = m.addConstr(left_expr <= 1 + right_expr)
                crossConstr.Lazy = 2
                
    m.optimize()           
                    
    return get_v_bicliques_from_e(G, Y, k)
    
    
        
    
    

