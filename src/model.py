from itertools import combinations

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.base_model import MBCModel


def var_swap(x, u, v, i):
    """
    This function gets the value of x[e, i] for e=[u, v] when the order of the vertices in the graph
    definition is unknown.
    """
    # the nested gets are for when x[u, v, i] doesn't exist, but x[v, u, i] does
    return x.get((u, v, i), x.get((v, u, i)))


class NaturalModel(MBCModel):

    def _add_variables(self):
        # 4j
        self.z = self.m.addVars(self.bicliques, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name="z")
        self.x = self.m.addVars(self.g.edges, self.bicliques, vtype=GRB.BINARY, name="x")

    def _set_objective(self):
        # 4a
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    def _add_constraints(self):
        m, x, z = self.m, self.x, self.z

        # 4b
        m.addConstrs(x[u, v, i] <= z[i] for u, v in self.g.edges for i in self.bicliques)
        # 4c
        m.addConstrs(gp.quicksum(x[u, v, i] for i in self.bicliques) >= 1 for u, v in self.g.edges)
        # 4e
        for cycle in nx.simple_cycles(self.g, length_bound=3):
            if len(cycle) != 3:
                continue
            cycle_edges = [[cycle[0], cycle[1]], [cycle[1], cycle[2]], [cycle[2], cycle[0]]]
            m.addConstrs(
                gp.quicksum(var_swap(x, u, v, i) for u, v in cycle_edges) <= 2 * z[i]
                for i in self.bicliques)

        for e, f, cr1, cr2 in self.get_disjoint_edges():
            a, b = e
            c, d = f
            # 4d
            if cr1 < 2 and cr2 < 2:
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] for i in self.bicliques)
            # 4f
            if cr1 == 2 and cr2 < 2:
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, a, d, i) for i in self.bicliques)
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, b, c, i) for i in self.bicliques)
            # 4g
            if cr1 < 2 and cr2 == 2:
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, a, c, i) for i in self.bicliques)
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, b, d, i) for i in self.bicliques)
            # 4h
            if cr1 == 2 and cr2 == 2:
                m.addConstrs(
                    x[a, b, i] + x[c, d, i] <=
                    z[i] + var_swap(x, a, c, i) + var_swap(x, a, d, i) for i in self.bicliques)
                m.addConstrs(
                    x[a, b, i] + x[c, d, i] <=
                    z[i] + var_swap(x, b, c, i) + var_swap(x, b, d, i) for i in self.bicliques)
        # 4i
        m.addConstrs(z[i] >= z[i + 1] for i in range(self.upper_bound() - 1))

    @classmethod
    def name(cls) -> str:
        return 'Compact Natural Model'


class ExtendedModel(MBCModel):

    def _add_variables(self):
        # 5h
        self.z = self.m.addVars(self.bicliques, vtype=GRB.BINARY, name="z")
        self.x = self.m.addVars(self.directed.edges, self.bicliques, vtype=GRB.BINARY, name="x")
        self.y = self.m.addVars(self.g.nodes, self.bicliques, range(2), vtype=GRB.BINARY, name="y")

    def _set_objective(self):
        # 5a
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    def _add_constraints(self):
        m, x, y, z = self.m, self.x, self.y, self.z

        # 5b
        m.addConstrs(x[u, v, i] <= y[u, i, 0] for u, v in self.directed.edges for i in self.bicliques)
        m.addConstrs(x[u, v, i] <= y[v, i, 1] for u, v in self.directed.edges for i in self.bicliques)
        # 5c
        m.addConstrs(y[u, i, 0] + y[v, i, 1] <= z[i] + x[u, v, i]
                     for u, v in self.directed.edges for i in self.bicliques)
        # 5d
        m.addConstrs(y[u, i, 0] + y[i, i, 1] <= z[i] for u in self.g.nodes for i in self.bicliques)
        # 5e
        for u, v in self.g.edges:
            m.addConstr(gp.quicksum(x[u, v, i] + x[v, u, i] for i in self.bicliques) >= 1)
        # 5f
        for u, v in combinations(self.g.nodes, r=2):
            if self.g.has_edge(u, v):
                continue
            m.addConstrs(y[u, i, 0] + y[v, i, 1] <= z[i] for i in self.bicliques)
            m.addConstrs(y[v, i, 0] + y[u, i, 1] <= z[i] for i in self.bicliques)
        # 5g
        m.addConstrs(z[i] >= z[i + 1] for i in range(self.upper_bound() - 1))

    @classmethod
    def name(cls) -> str:
        return 'Extended Model'
