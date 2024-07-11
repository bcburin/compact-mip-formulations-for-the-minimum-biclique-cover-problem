from itertools import combinations

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.base_model import MBCModel
from src.independent_set import solve_max_weighted_independent_set
from src.util import is_biclique, var_swap


def indep_callback(model: gp.Model, where: int):
    """
    :param model: Gurobi model, expects the following variables to be added:
        - _g containing the graph of the model;
        - _g2 containing the power graph g^2; and
        - _k storing the calculated maximum number of bicliques.
        - _Y variable of the extended model
        - _z variable of the extended model
    :param where: Gurobi internal variable that represent place where the callback function was called.
    """
    if where != GRB.Callback.MIPNODE or model.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
        return
    y_val = model.cbGetNodeRel(model._y)
    z_val = model.cbGetNodeRel(model._z)
    for j in range(model._k):
        for v in model._g.nodes:
            model._g.nodes[v]["weight"] = y_val[v, j, 0] + y_val[v, j, 1]
        obj, indep_set = solve_max_weighted_independent_set(model._g, model._g2)
        if obj > z_val[j]:
            model.cbCut(gp.quicksum(model._y[v, j, 0] + model._y[v, j, 1] for v in indep_set) <= model._z[j])


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
        for cycle in nx.simple_cycles(self.g, length_bound=3):  # TODO: fix me
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
        # independent edges constraints
        if self._can_add_indep_edges_constraints():
            self._add_independent_edges_constraints()

    def _check_biclique_cover(self) -> bool:
        # check it's a cover
        if not any(self.x[u, v, i].X == 1 for u, v in self.g.edges for i in self.bicliques):
            return False
        # check it's a biclique cover
        for i in self.bicliques:
            if self.z[i].X != 1:
                continue
            edges = [(u, v) for u, v in self.g.edges if self.x[u, v, i].X == 1]
            if not is_biclique(graph=self.g, edges=edges):
                return False
        return True

    def _do_warm_start(self, indep_edges: list, vertex_cover: list):
        assign = dict()
        for i, s in enumerate(vertex_cover):
            edges = []
            for e in self.g.edges(s):
                self.x[min(e), max(e), i].start = 1
                edges.append(e)
            assign[i] = edges
        indep_edges = set(indep_edges)
        if self._edge_fix:
            for i in assign.keys():
                for e in assign[i]:
                    if e in indep_edges:
                        self.x[min(e), max(e), i].lb = 1

    def _pre_solve(self):
        for b in range(self.lower_bound()):
            self.z[b].lb = 1
        if self._warm_start:
            _, indep_edges = self.get_lb_and_indep_edges()
            _, vertex_cover = self.get_ub_and_vertex_cover()
            self._do_warm_start(indep_edges=indep_edges, vertex_cover=vertex_cover)

    def _add_independent_edges_constraints(self):
        _, edges = self.get_lb_and_indep_edges()
        for biclique in self.bicliques:
            if biclique >= len(edges):
                return
            e = edges[biclique]
            self.x[e[0], e[1], biclique].lb = 1

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
        m.addConstrs(y[u, i, 0] + y[u, i, 1] <= z[i] for u in self.g.nodes for i in self.bicliques)
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
        # independent edges constraints
        if self._can_add_indep_edges_constraints():
            self._add_independent_edges_constraints()
        # conflict inequalities
        if self._conflict_inequalities:
            self._add_conflict_inequalities()
        # common neighbors inequalities
        if self._common_neighbor_inequalities:
            self._add_common_neighbor_inequalities()

    def _check_biclique_cover(self) -> bool:
        # check it's a cover
        if not any(self.x[u, v, i].X == 1 or self.x[v, u, i].X == 1
                   for u, v in self.g.edges for i in self.bicliques):
            return False
        # check it's a biclique cover
        for i in self.bicliques:
            if self.z[i].X != 1:
                continue
            edges = [(u, v) for u, v in self.g.edges if self.x[u, v, i].X == 1 or self.x[v, u, i].X == 1]
            if not is_biclique(graph=self.g, edges=edges):
                return False
        return True

    def _pre_solve(self):
        for b in range(self.lower_bound()):
            self.z[b].lb = 1
        if self._use_callback:
            self._add_callback()
        if self._warm_start:
            _, indep_edges = self.get_lb_and_indep_edges()
            _, vertex_cover = self.get_ub_and_vertex_cover()
            self._do_warm_start(indep_edges=indep_edges, vertex_cover=vertex_cover)

    def _add_independent_edges_constraints(self):
        _, edges = self.get_lb_and_indep_edges()
        for biclique in self.bicliques:
            if biclique >= len(edges):
                return
            a, b = edges[biclique]
            self.x[a, b, biclique].lb = 1
            self.x[b, a, biclique].lb = 0

    def _add_conflict_inequalities(self):
        for u, v in combinations(self.g.nodes, r=2):
            if self.power_graph.has_edge(u, v):
                continue
            conflict_inequalities = self.m.addConstrs(
                self.y[u, i, 0] + self.y[u, i, 1] + self.y[v, i, 0] + self.y[v, i, 1] <= self.z[i]
                for i in self.bicliques)
            for i in self.bicliques:
                conflict_inequalities[i].Lazy = 1

    def _add_common_neighbor_inequalities(self):
        for u, v in self.power_graph.edges:
            if self.g.has_edge(u, v):
                continue
            common_neighbors = nx.common_neighbors(self.g, u, v)
            self.m.addConstrs(
                self.y[u, i, 0] + self.y[v, i, 0] <= self.z[i] + gp.quicksum(self.y[c, i, 1] for c in common_neighbors)
                for i in self.bicliques)
            self.m.addConstrs(
                self.y[u, i, 1] + self.y[v, i, 1] <= self.z[i] + gp.quicksum(self.y[c, i, 0] for c in common_neighbors)
                for i in self.bicliques)

    def _do_warm_start(self, indep_edges: list, vertex_cover: list):
        assign = dict()
        for i, s in enumerate(vertex_cover):
            edges = []
            for e in self.g.edges(s):
                a, b = e
                self.x[a, b, i].start = 1
                self.x[b, a, i].start = 0
                self.y[a, i, 0].start = 1
                self.y[a, i, 1].start = 0
                self.y[b, i, 1].start = 1
                self.y[b, i, 0].start = 0
                edges.append(e)
            assign[i] = edges
        indep_edges = set(indep_edges)
        if self._edge_fix:
            for i in assign.keys():
                for e in assign[i]:
                    if e in indep_edges:
                        a, b = e
                        self.x[a, b, i].lb = 1
                        self.x[b, a, i].ub = 0
                        self.y[a, i, 0].lb = 1
                        self.y[a, i, 1].ub = 0
                        self.y[b, i, 1].lb = 1
                        self.y[b, i, 0].ub = 0

    def _add_callback(self):
        self.m._k = self.upper_bound()
        self.m._y = self.y
        self.m._z = self.z
        self.m._g = self.g
        self.m._g2 = self.power_graph
        self._callback = indep_callback

    @classmethod
    def name(cls) -> str:
        return 'Extended Model'
