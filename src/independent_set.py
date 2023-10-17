import gurobipy as gp
from gurobipy import GRB 

import networkx as nx
from pandas import DataFrame

from bc_bounds import find_bc_upper_bound, UBComputeMethod
from util import get_graphs_in_store, chronometer


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


def solve_bc(g: nx.Graph, use_relaxation: bool = True):
    # use frozen sets so they are hashable
    indep_sets = frozenset(frozenset(indep_set) for indep_set in get_set_of_maximal_independent_sets(g))
    lb = 1  # independent edges method doesn't work
    ub = int(find_bc_upper_bound(g, method=UBComputeMethod.VERTEX))
    dg = nx.DiGraph(g)
    bicliques = range(ub)

    # create model
    m = gp.Model()
    # define variables
    w = m.addVars(indep_sets, bicliques, range(2), vtype=GRB.BINARY, name='w')
    x = m.addVars(g.nodes, bicliques, range(2),
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='x')
    p = m.addVars(g.nodes, bicliques, range(2),
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='p')
    y = m.addVars(dg.edges, bicliques,
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='y')
    z = m.addVars(bicliques,
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='z')
    # objective function
    m.setObjective(gp.quicksum(z), GRB.MINIMIZE)
    # independent set partition constraints
    m.addConstrs(w[s, b, 0] + w[s, b, 1] <= 1 for s in indep_sets for b in bicliques)
    m.addConstrs(w[s, b, d] <= x[u, b, d]
                 for s in indep_sets for u in s for b in bicliques for d in range(2))
    m.addConstrs(x[u, b, d] <= gp.quicksum(w[s, b, d] for s in indep_sets if u in s)
                 for u in g.nodes for b in bicliques for d in range(2))
    # p[u, b, d] = x[u, b, d] * (1 - x[u, b, 1-d])
    m.addConstrs(p[u, b, d] <= x[u, b, d] for u in g.nodes for b in bicliques for d in range(2))
    m.addConstrs(p[u, b, d] <= 1 - x[u, b, 1-d] for u in g.nodes for b in bicliques for d in range(2))
    m.addConstrs(p[u, b, d] >= x[u, b, d] - x[u, b, 1-d] for u in g.nodes for b in bicliques for d in range(2))
    # y[u, v, b] = p[u, b, 0] * p[v, b, 1]
    m.addConstrs(y[u, v, b] <= p[u, b, 0] for u, v in dg.edges for b in bicliques)
    m.addConstrs(y[u, v, b] <= p[v, b, 1] for u, v in dg.edges for b in bicliques)
    m.addConstrs(y[u, v, b] >= p[u, b, 0] + p[v, b, 1] - 1 for u, v in dg.edges for b in bicliques)
    # cover constraints
    m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in g.edges)
    m.addConstrs(y[u, v, b] <= z[b] for u, v in dg.edges for b in bicliques)
    m.addConstrs(p[u, b, 0] + p[v, b, 1] <= 1 for u, v in nx.complement(dg).edges for b in bicliques)
    # symmetry break
    m.addConstrs(z[b + 1] <= z[b] for b in range(ub-1))
    # lower bound
    for b in range(lb):
        z[b].lb = 1

    m.optimize()

    # construct dataframe with bicliques
    df_bicliques = DataFrame(
        {b: [False] * len(g.edges) for b in bicliques if z[b].X == 1}, index=[str(e) for e in g.edges])
    for e in g.edges:
        u, v = e
        for b in bicliques:
            if z[b].X == 1 and (y[u, v, b].X == 1 or y[v, u, b].X == 1):
                df_bicliques.loc[str(e), b] = True

    return m.objVal, df_bicliques


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


if __name__ == '__main__':
    for g, g_name in get_graphs_in_store(fname_regex='partite', max_graphs=1):
        solution = solve_bc(g)
        print(solution)
