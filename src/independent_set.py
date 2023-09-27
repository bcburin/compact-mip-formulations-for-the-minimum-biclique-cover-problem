import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
from pandas import DataFrame

from bc_bounds import find_bc_lower_bound, find_bc_upper_bound, LBComputeMethod, UBComputeMethod
from util import get_graphs_in_store, chronometer
import vertex_cover


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


def get_set_of_maximal_independent_sets(g: nx.Graph):
    return nx.find_cliques(nx.complement(g))


def solve_bc(g: nx.Graph):
    # sets and lists are unhashable and tuples could not index variables, thus indep_sets stores int keys to index
    # decision variables and the actual sets as values
    indep_sets = {i: tuple(s) for i, s in enumerate(get_set_of_maximal_independent_sets(g))}
    lb = int(find_bc_lower_bound(g, method=LBComputeMethod.INDEPENDENT_EDGES))
    ub = int(find_bc_upper_bound(g, method=UBComputeMethod.VERTEX))
    dg = nx.DiGraph(g)
    bicliques = range(ub)

    # create model
    m = gp.Model()
    # define variables
    z = m.addVars(bicliques, vtype=GRB.BINARY, name='z')
    w = m.addVars(indep_sets.keys(), bicliques, range(2), vtype=GRB.BINARY, name='w')
    x = m.addVars(g.nodes, bicliques, range(2), vtype=GRB.BINARY, name='x')
    p = m.addVars(g.nodes, bicliques, range(2), vtype=GRB.BINARY, name='p')
    y = m.addVars(dg.edges, bicliques, vtype=GRB.BINARY, name='y')
    # objective function
    m.setObjective(gp.quicksum(z), GRB.MINIMIZE)
    # constraints
    m.addConstrs(w[si, b, 0] + w[si, b, 1] <= 1 for si in indep_sets.keys() for b in bicliques)
    m.addConstrs(w[si, b, d] <= x[u, b, d]
                 for si, s in indep_sets.items() for u in s for b in bicliques for d in range(2))
    m.addConstrs(x[u, b, d] <= gp.quicksum(w[si, b, d] for si, s in indep_sets.items() if u in s)
                 for u in g.nodes for b in bicliques for d in range(2))
    m.addConstrs(2*p[u, b, d] <= x[u, b, d] - x[u, b, 1-d] + 1 for u in g.nodes for b in bicliques for d in range(2))
    m.addConstrs(2*y[u, v, b] <= p[u, b, 0] + p[v, b, 1] for u, v in dg.edges for b in bicliques)
    m.addConstrs(p[u, b, 0] + p[v, b, 1] - 1 <= y[u, v, b] for u, v in dg.edges for b in bicliques)
    m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in g.edges)
    m.addConstrs(p[u, b, 0] + p[v, b, 1] <= 1 for u, v in nx.complement(dg).edges for b in bicliques)
    m.addConstrs(z[b + 1] <= z[b] for b in range(ub-1))
    m.addConstrs(z[b] == 1 for b in range(lb))
    m.addConstrs(y[u, v, b] <= z[b] for u, v in dg.edges for b in bicliques)

    m.optimize()

    # construct dataframe with bicliques
    df_bicliques = DataFrame(
        {b: [False] * len(g.edges) for b in bicliques if z[b].X == 1}, index=[str(e) for e in g.edges])
    for e in g.edges:
        for b in bicliques:
            if z[b].X == 1:
                u, v = e
                if y[u, v, b].X == 1 or y[v, u, b].X == 1:
                    df_bicliques.loc[str(e), b] = True

    return m.objVal, df_bicliques


def print_independent_sets(g: nx.Graph, max_indep_sets: int = None):
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


if __name__ == '__main__':
    for g, g_name in get_graphs_in_store(fname_regex='simple'):
        solution = solve_bc(g)
        print(solution)

