from enum import Enum, auto
from itertools import combinations

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import numpy as np
import cvxpy as cp

from src.util import get_graphs_in_store, GraphReport


class LBComputeMethod(Enum):
    MATCH = auto()
    LOVASZ = auto()
    CLIQUE = auto()
    INDEPENDENT_EDGES = auto()
    MAXIMAL_INDEPENDENT_SET = auto()


class UBComputeMethod(Enum):
    NUMBER = auto()
    VERTEX = auto()


def find_bc_lower_bound(g: nx.Graph, method: LBComputeMethod = LBComputeMethod.MATCH) -> int:
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
            return compute_lb_by_independent_edges_method(g)
        case LBComputeMethod.MAXIMAL_INDEPENDENT_SET:  # should be good for dense graphs
            cliques = nx.find_cliques(nx.complement(g))
            return np.ceil(np.log2(len(cliques)))
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


def compute_lb_by_independent_edges_method(g: nx.Graph) -> int:
    m = gp.Model()
    y = m.addVars(g.edges, vtype=GRB.BINARY)
    # objective function
    m.setObjective(gp.quicksum(y[e] for e in g.edges), GRB.MAXIMIZE)
    # constraints
    m.addConstrs(gp.quicksum(y[e] for e in g.nodes(v)) <= 1 for v in g.nodes)
    for e1, e2 in combinations(g.edges, r=2):
        a, b = e1
        c, d = e2
        if (g.has_edge(a, d) and g.has_edge(b, c)) or (g.has_edge(a, c) and g.has_edge(b, d)):
            m.addConstr(y[a, b] + y[c, d] <= 1)
    # solve model
    m.optimize()
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        return m.objVal
    else:
        return 1


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
        # set params
        m.Params.LogToConsole = 0
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
    report = GraphReport('bounds')
    report.add_properties([str(method) for method in LBComputeMethod if method != LBComputeMethod.LOVASZ])
    report.add_properties([str(method) for method in UBComputeMethod])
    # iterate graphs
    for g, g_name in get_graphs_in_store():
        report.add_graph_data(g, g_name)
        for method in LBComputeMethod:
            # LOVASZ is too slow for a small number of edges
            if method == LBComputeMethod.LOVASZ:
                continue
            report.add_property_values_from_function(p_name=str(method), f=find_bc_lower_bound, g=g)
        for method in UBComputeMethod:
            report.add_property_values_from_function(p_name=str(method), f=find_bc_upper_bound, g=g)
    report.save_csv()
