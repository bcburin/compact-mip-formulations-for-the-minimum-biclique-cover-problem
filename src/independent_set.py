from random import choices
from itertools import combinations

import gurobipy as gp
from gurobipy import GRB
import networkx as nx

from src.util import chronometer, poisson, save_graph_in_store


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


SetOfSets = frozenset[frozenset[int]]


def get_set_of_maximal_independent_sets(g: nx.Graph):
    return nx.find_cliques(nx.complement(g))


def get_indep_sets(g: nx.Graph) -> SetOfSets:
    return frozenset(frozenset(indep_set) for indep_set in get_set_of_maximal_independent_sets(g))


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


def generate_random_independent_sets(
        n: int, r: int, min_size: int = 1, max_size: int = None) -> SetOfSets:
    if max_size is None:
        max_size = n
    if (not isinstance(min_size, int) or not isinstance(max_size, int)
            or min_size < 1 or min_size > n or max_size < 1 or max_size > n):
        raise ValueError('Invalid input for min_size or max_size. They need to be integers within [1, n].')

    sets = set()
    while len(sets) != r:
        a, b = min_size, max_size + 1
        weights = [poisson(p, x_0=a, mu=(b-a)/3) for p in range(a, b)]
        size = choices(range(a, b), weights=weights, k=1)[0]
        rs = frozenset(choices(range(n), k=size))
        # purge subsets
        if any(rs.issubset(s) for s in sets):
            continue
        sets = {s for s in sets if not s.issubset(rs)}
        sets.add(rs)
    return frozenset(sets)


def create_graph_from_maximal_indep_sets(sets: SetOfSets) -> nx.Graph:
    g = nx.Graph()
    for s in sets:
        if len(s) == 1:
            g.add_node(list(s)[0])
            continue
        for u, v in combinations(s, r=2):
            g.add_edge(u, v)
    return nx.complement(g)


def diff_indep_sets(g: nx.Graph, indep_sets: SetOfSets) -> SetOfSets:
    actual_indep_sets = get_indep_sets(g=g)
    return actual_indep_sets ^ indep_sets


def generate_and_save_random_graphs_from_indep_sets(
        min_n_sets: int, max_n_sets: int,
        min_n_vertices: int, max_n_vertices, step_n_vertices : int = 1):
    for n_vertices in range(min_n_vertices, max_n_vertices+1, step_n_vertices):
        # define number of sets to be, in general
        weights = [poisson(p, x_0=0, mu=(max_n_sets - min_n_sets) / 2) for p in range(min_n_sets, max_n_sets)]
        n_sets = choices(range(min_n_sets, max_n_sets), weights=weights, k=1)[0] + 1
        # generate independent sets
        indep_sets = generate_random_independent_sets(n=n_vertices, r=n_sets)
        # create graph from
        gen_g = create_graph_from_maximal_indep_sets(sets=indep_sets)
        # save graph based on the actual independent sets
        actual_indep_sets = get_indep_sets(g=gen_g)
        gen_g_name = f'from_indep_sets_{len(actual_indep_sets)}'
        save_graph_in_store(g=gen_g, g_name=gen_g_name)


if __name__ == '__main__':
    generate_and_save_random_graphs_from_indep_sets(
        min_n_sets=2, max_n_sets=6, min_n_vertices=10, max_n_vertices=60, step_n_vertices=10)
