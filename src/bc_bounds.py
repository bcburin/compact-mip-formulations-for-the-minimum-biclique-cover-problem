from enum import Enum, auto
from itertools import combinations

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import numpy as np
import cvxpy as cp
from pandas import DataFrame

from util import get_graphs_in_store, chronometer


class LBComputeMethod(Enum):
    MATCH = auto()
    LOVASZ = auto()
    CLIQUE = auto()
    INDEPENDENT_EDGES = auto()


class UBComputeMethod(Enum):
    NUMBER = auto()
    VERTEX = auto()


def find_bc_lower_bound(g: nx.Graph, method: LBComputeMethod = LBComputeMethod.MATCH) -> int:
    match method:
        case LBComputeMethod.MATCH:
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
            return compute_lb_by_independent_edges_method(g)
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

        if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
            return m.objVal
        else:
            print("There is an error in the maximum clique problem!")

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')


def compute_lb_by_independent_edges_method(g: nx.Graph | nx.DiGraph) -> int:
    if isinstance(g, nx.Graph):
        g = g.to_directed(as_view=True)
    try:
        # define model
        m = gp.Model()
        # define vars
        r = m.addVars(g.nodes, vtype=GRB.BINARY, name='r')
        y = m.addVars(g.edges)
        # define objective function
        m.setObjective(gp.quicksum(r), GRB.MAXIMIZE)
        # add constraints
        m.addConstrs(gp.quicksum([y[e] for e in g.in_edges(v)]) + r[v] <= 1 for v in g.nodes)
        m.addConstrs(r[u] + r[v] <= 1 for u, v in g.edges)
        for v in g.nodes:
            m.addConstrs(y[e] <= r[v] for e in g.out_edges(v))
        m.addConstrs(r[v] <= gp.quicksum([y[e] for e in g.out_edges(v)]) for v in g.nodes)
        for e1, e2 in combinations(g.edges, r=2):
            u, v = e1
            c, d = e2
            if (not {u, v} & {c, d} and
                    ((g.has_edge(c, u) and g.has_edge(d, v)) or
                     (g.has_edge(c, v) and g.has_edge(d, u)))):
                m.addConstr(y[(u, v)] + y[(v, u)] + y[(c, d)] + y[(d, c)] <= 1)
        # set a one-minute time limit
        m.Params.TimeLimit = 60
        # optimize
        m.optimize()

        if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
            return m.objVal
        else:
            print("There is an error in the maximum clique problem!")

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')


def find_bc_upper_bound(g: nx.Graph, method: UBComputeMethod = UBComputeMethod.NUMBER) -> int:
    match method:
        case UBComputeMethod.NUMBER:
            n = len(g.nodes)
            return n + 1 - np.floor(np.log2(n))
        case UBComputeMethod.VERTEX:
            return vertex_cover(g)
        case _:
            raise ValueError("Unsupported Method")


def vertex_cover(g: nx.Graph) -> int:
    try:
        # define model
        m = gp.Model()
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

        if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
            return m.objVal
        else:
            print("There is an error in the vertex cover problem!")

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')


if __name__ == "__main__":
    data_lb = {
        'graph': [],
        'n_nodes': [],
        'n_edges': [],
    }
    # add lists to the data dictionary to store computation times
    for method in LBComputeMethod:
        # LOVASZ is too slow for a small number of edges
        if method == LBComputeMethod.LOVASZ:
            continue
        data_lb[str(method)] = []
        data_lb[f'{method}_time'] = []
    for method in UBComputeMethod:
        data_lb[str(method)] = []
        data_lb[f'{method}_time'] = []
    # iterate graphs
    for g, g_name in get_graphs_in_store(max_edges=1000):
        data_lb['graph'].append(g_name)
        data_lb['n_nodes'].append(len(g.nodes))
        data_lb['n_edges'].append(len(g.edges))
        for method in LBComputeMethod:
            # LOVASZ is too slow for a small number of edges
            if method == LBComputeMethod.LOVASZ:
                continue
            lb, total_time = chronometer(find_bc_lower_bound, g, method=method)
            data_lb[str(method)].append(lb)
            data_lb[f'{method}_time'].append(total_time)
        for method in UBComputeMethod:
            ub, total_time = chronometer(find_bc_upper_bound, g, method=method)
            data_lb[str(method)].append(ub)
            data_lb[f'{method}_time'].append(total_time)
    df_bounds = DataFrame(data=data_lb)
    df_bounds.to_csv('bounds.csv', index=False)


