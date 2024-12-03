from enum import Enum, auto
from itertools import combinations
from time import time

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import numpy as np
import cvxpy as cp

from src.util import get_graphs_in_store, GraphReport, get_graph_in_store


class LBComputeMethod(str, Enum):
    MATCH = 'match'
    LOVASZ = 'lovasz'
    CLIQUE = 'clique'
    INDEPENDENT_EDGES = 'independent_edges'
    MAXIMAL_INDEPENDENT_SET = 'maximal_independent_set'


class UBComputeMethod(str, Enum):
    NUMBER = 'number'
    VERTEX = 'vertex'
    MERGE_STARS = 'merge_stars'


def count_cliques(g: nx.Graph, timeout: int = None, size_limit: int = None, verification_interval: int = 1000) -> int:
    clique_count = 0
    start_time = time()
    for _ in nx.find_cliques(nx.complement(g)):
        if clique_count % verification_interval == 0:
            if timeout and time() - start_time >= timeout:
                return clique_count
            if size_limit and clique_count >= size_limit:
                return clique_count
        clique_count += 1
    return clique_count


def find_bc_lower_bound(g: nx.Graph, method: LBComputeMethod = LBComputeMethod.MATCH,
                        time_limit: int = None, memory_limit: int = None) -> int:
    match method:
        case LBComputeMethod.MATCH:  # good for sparse graph
            m_len = len(nx.max_weight_matching(g))
            e_len = len(g.edges)
            return np.ceil(m_len ** 2 / e_len)
        case LBComputeMethod.LOVASZ:
            g_c = nx.complement(g)
            lov_num = np.round(compute_lovasz_number(g_c))
            return np.ceil(np.log2(lov_num))
        case LBComputeMethod.CLIQUE:
            return np.ceil(np.log2(max_clique(g)))
        case LBComputeMethod.INDEPENDENT_EDGES:
            return compute_lb_by_independent_edges_method(g, time_limit=time_limit, memory_limit=memory_limit)
        case LBComputeMethod.MAXIMAL_INDEPENDENT_SET:  # should be good for dense graphs
            number_of_cliques = count_cliques(g, timeout=time_limit, size_limit=int(memory_limit*1e9/4))
            return np.ceil(np.log2(number_of_cliques))
        case _:
            raise ValueError("Unsupported Method")

        
def compute_lovasz_number(g: nx.Graph) -> int:
    a = nx.adjacency_matrix(g)
    n = a.shape[0]
    c = np.ones((n, n))
    i = np.identity(n)
    x = cp.Variable(shape=(n, n), symmetric=True)
    # The operator >> denotes matrix inequality.
    constraints = [x >> 0]
    constraints += [cp.trace(i @ x) == 1]
    row_ind, col_ind = a.nonzero()
    constraints += [x[row_ind[i], col_ind[i]] == 0 for i in range(len(row_ind))]
    prob = cp.Problem(cp.Maximize(cp.trace(c @ x)), constraints)
    prob.solve()

#    # Print result.
    print("The optimal value is", prob.value)
#    print("A solution X is")
#    print(X.value)
    return prob.value


def max_clique(g: nx.Graph) -> int:
    try:
        # define model
        m = gp.Model()
        # define vars
        x = m.addVars(g.nodes, vtype=GRB.BINARY, name="x")
        # define objective function
        m.setObjective(gp.quicksum(x), sense=GRB.MAXIMIZE)
        gc = nx.complement(g)
        # add covering constraints
        m.addConstrs(x[u] + x[v] <= 1 for u, v in gc.edges)
        # set a one-minute time limit
        m.Params.TimeLimit = 60
        # optimize
        m.optimize()

        if m.status != GRB.INFEASIBLE:
            return m.objVal
        else:
            print("There is an error in the maximum clique problem!")

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')


def compute_lb_by_independent_edges_method(g: nx.Graph, time_limit: int = None, memory_limit: int = None) -> int:
    lb, _ = compute_lb_and_get_edges_by_independent_edges_method(g, time_limit, memory_limit)
    return lb


def compute_lb_and_get_edges_by_independent_edges_method(
        g: nx.Graph, time_limit: int = None, memory_limit: int = None) -> tuple[int, list]:
    m = gp.Model()
    y = m.addVars(g.edges, vtype=GRB.BINARY)
    # set time limit
    if time_limit:
        m.Params.TimeLimit = time_limit
    if memory_limit:
        m.Params.SoftMemLimit = memory_limit
    # objective function
    m.setObjective(gp.quicksum(y[e] for e in g.edges), GRB.MAXIMIZE)
    # constraints
    for v in g.nodes:
        for e in g.edges(v):
            try:
                m.addConstr(y[e] <= 1, name=f'(9b) for ({e[0]}, {e[1]}) from v={v}')
            except KeyError:
                m.addConstr(y[e[1], e[0]] <= 1, name=f'(9b) for ({e[1]}, {e[0]}) from v={v}')
    for e1, e2 in combinations(g.edges, r=2):
        a, b = e1
        c, d = e2
        if (g.has_edge(a, d) and g.has_edge(b, c)) or (g.has_edge(a, c) and g.has_edge(b, d)):
            m.addConstr(y[a, b] + y[c, d] <= 1, name=f'(9c) for ({a}, {b}) and ({c}, {d})')
    # solve model
    m.optimize()
    # return values
    if m.status != GRB.INFEASIBLE:
        edges = [e for e in g.edges if y[e].X == 1]
        return m.objVal, edges
    else:
        return 1, []


def find_bc_upper_bound(g: nx.Graph, method: UBComputeMethod = UBComputeMethod.NUMBER,
                        time_limit: int = None, memory_limit: int = None) -> int:
    match method:
        case UBComputeMethod.NUMBER:
            n = len(g.nodes)
            return n + 1 - np.floor(np.log2(n))
        case UBComputeMethod.VERTEX:
            return vertex_cover(g, time_limit=time_limit, memory_limit=memory_limit)
        case UBComputeMethod.MERGE_STARS:
            bicliques = merge_stars(g)
            return len(bicliques)
        case _:
            raise ValueError("Unsupported Method")


def get_vertex_cover_solution(g: nx.Graph, time_limit: int = None, memory_limit: int = None) -> tuple[list, int]:
    try:
        # define model
        m = gp.Model()
        # set params
        m.Params.LogToConsole = 0
        if time_limit:
            m.Params.TimeLimit = time_limit
        if memory_limit:
            m.Params.SoftMemLimit = memory_limit
        # define vars
        x = m.addVars(g.nodes, vtype=GRB.BINARY, name="x")
        # define objective function
        m.setObjective(gp.quicksum(x), sense=GRB.MINIMIZE)
        # add covering constraints
        m.addConstrs(x[u] + x[v] >= 1 for u, v in g.edges)
        # set a one-minute time limit
        m.Params.TimeLimit = 60
        # optimize
        m.optimize()

        if m.status != GRB.INFEASIBLE:
            t = [vertex for vertex in g.nodes if x[vertex].X > 0.5]
            return t, m.objVal
        else:
            print("There is an error in the vertex cover problem!")

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')


def vertex_cover(g: nx.Graph, time_limit: int = None, memory_limit: int = None) -> int:
    _, val = get_vertex_cover_solution(g, time_limit, memory_limit)
    return val


def merge_stars(g: nx.Graph, t: gp.Model = None) -> list:
    if t is None:
        t, _ = get_vertex_cover_solution(g)
    g2 = nx.power(g, k=2)
    for vertex in g2.nodes:
        g2.nodes[vertex]['cover'] = False
    for vertex in t:
        g2.nodes[vertex]['cover'] = True
    biclique_cover = []
    for u, v in g2.edges:
        if g.has_edge(u, v):
            continue
        if not g2.nodes[u]['cover'] or not g2.nodes[v]['cover']:
            continue
        u_neighbor_set = set(g.neighbors(u))
        v_neighbor_set = set(g.neighbors(v))
        if u_neighbor_set.issubset(v_neighbor_set):
            biclique_cover.append(set(g.edges(u)) | {(q, v) for q in nx.common_neighbors(g, u, v)})
            g2.nodes[u]['cover'] = False
        if v_neighbor_set.issubset(u_neighbor_set):
            biclique_cover.append(set(g.edges(v)) | {(q, u) for q in nx.common_neighbors(g, u, v)})
            g2.nodes[v]['cover'] = False
    for vertex in g2.nodes:
        if not g2.nodes[vertex]['cover']:
            continue
        biclique_cover.append(set(g.edges(vertex)))
    return biclique_cover


def get_partial_model_of_extension_lb(g: nx.Graph) -> gp.Model:
    m = gp.Model()
    directed = g.to_directed()
    # add vars
    x = m.addVars(directed.edges, vtype=GRB.BINARY, name="x")
    y = m.addVars(g.nodes, range(2), vtype=GRB.BINARY, name="y")
    # set objective function
    m.setObjective(gp.quicksum(x), sense=GRB.MAXIMIZE)
    # add constraints
    m.addConstrs(y[u, 0] + y[v, 0] <= 1 for u, v in g.edges)
    m.addConstrs(y[u, 1] + y[v, 1] <= 1 for u, v in g.edges)
    m.addConstrs(y[v, 0] + y[v, 1] <= 1 for v in g.nodes)
    m.addConstrs(x[u, v] <= y[u, 0] for u, v in directed.edges)
    m.addConstrs(x[u, v] <= y[v, 1] for u, v in directed.edges)
    m.addConstrs(y[u, 0] + y[v, 1] <= 1 + x[u, v] for u, v in directed.edges)
    m._x = x
    m._y = y
    return m


def get_extension_lb(g: nx.Graph):
    _, edges = compute_lb_and_get_edges_by_independent_edges_method(g)
    m = get_partial_model_of_extension_lb(g)
    directed = g.to_directed()
    bicliques = {(i, j): [] for i, j in edges}
    for i, j in edges:
        m._x[i, j].lb = 1
        m._y[i, 0].lb = 1
        m._y[j, 1].lb = 1
        m.solve()
        bicliques[i, j] = [(u, v) for u, v in directed.edges if m._x[u, v].X > 0.5]
        m._x[i, j].lb = 0
        m._y[i, 0].lb = 0
        m._y[j, 1].lb = 0
    return bicliques


if __name__ == "__main__":
    g = get_graph_in_store(filename="simple_5_7.gml")
    g_bicliques = get_extension_lb(g)
    print(g_bicliques)
    # model_time_limit = None
    # model_memory_limit = 4
    # report = GraphReport('bounds')
    # report.add_properties([str(method) for method in LBComputeMethod if method != LBComputeMethod.LOVASZ])
    # # report.add_properties([str(method) for method in UBComputeMethod])
    # # iterate graphs
    # for g, g_name in get_graphs_in_store(recursive=False):
    #     report.add_graph_data(g, g_name)
    #     for method in LBComputeMethod:
    #         # LOVASZ is too slow for a small number of edges
    #         if method == LBComputeMethod.LOVASZ:
    #             continue
    #         report.add_property_values_from_function(p_name=str(method), f=find_bc_lower_bound, g=g,
    #                                                  method=method, time_limit=model_time_limit,
    #                                                  memory_limit=model_memory_limit)
    #     # for method in UBComputeMethod:
    #     #     report.add_property_values_from_function(p_name=str(method), f=find_bc_upper_bound, g=g,
    #     #                                              method=method, time_limit=model_time_limit,
    #     #                                              memory_limit=model_memory_limit)
    # report.save_csv()
