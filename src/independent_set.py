import gurobipy as gp
from gurobipy import GRB 

import networkx as nx

from util import chronometer


def solve(g):
    m = gp.Model()
    # find power graph g^2
    H = nx.power(g, 2)
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
        return m.objVal
    else:
        print("There is an error in the independent set problem!")


def get_set_of_maximal_independent_sets(g: nx.Graph):
    return nx.find_cliques(nx.complement(g))


def print_independent_sets(g: nx.Graph, g_name, max_indep_sets: int = None):
    max_independent_sets, compute_time = chronometer(lambda x: set(get_set_of_maximal_independent_sets(x)), g)
    print(f'Graph:', g_name)
    print(f'Nodes:', len(g.nodes))
    print(f'Edges:', len(g.edges))
    print(f'Time:', compute_time)
    print(f'Max Independent Sets')
    for count, independent_set in enumerate(max_independent_sets):
        # only show the first ten independent sets
        if max_independent_sets and count >= max_indep_sets:
            print('\t...\n')
            break
        print('\t', count, ': {', ', '.join([str(v) for v in independent_set]), '}')

