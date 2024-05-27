from functools import cache
from itertools import combinations

import gurobipy as gp
import networkx as nx
from gurobipy import GRB, Var

from src.base_model import MBCModel, BottomUpMBCModel
from src.bc_bounds import LBComputeMethod, compute_lb_and_get_edges_by_independent_edges_method, \
    get_vertex_cover_solution, UBComputeMethod
from src.util import is_biclique, var_swap


class NaturalModel(MBCModel):

    def __init__(self, *args, **kwargs):
        self._vertex_cover_solution = []
        self._indep_edges = []
        super().__init__(*args, **kwargs)

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
        match = {e: set(s for s in vertex_cover if s == e[0] or s == e[1]) for e in indep_edges}
        for i, e in enumerate(match.keys()):
            self.x[min(e), max(e), i].lb = 1
            for s in match[e]:
                for se in self.g.edges(s):
                    self.x[min(se), max(se), i].start = 1

    def _can_add_indep_edges_constraints(self) -> bool:
        return self._indep_edges is not None and self._edge_fix and not self._warm_start

    def _can_warm_start(self) -> bool:
        return self._indep_edges is not None and self._vertex_cover_solution is not None and self._warm_start

    def _pre_solve(self):
        for b in range(self.lower_bound()):
            self.z[b].lb = 1
        if self._can_add_indep_edges_constraints():
            self._add_independent_edges_constraints(self._indep_edges)
        if self._can_warm_start():
            self._do_warm_start(indep_edges=self._indep_edges, vertex_cover=self._vertex_cover_solution)

    def _add_independent_edges_constraints(self, edges: list[Var]):
        for biclique in self.bicliques:
            if biclique >= len(edges):
                return
            e = edges[biclique]
            self.x[min(e), max(e), biclique].lb = 1

    @cache
    def lower_bound(self) -> int:
        if self._lb_method == LBComputeMethod.INDEPENDENT_EDGES:
            lb, edges = compute_lb_and_get_edges_by_independent_edges_method(g=self.g)
            self._indep_edges = edges
            return int(lb)
        else:
            return super().lower_bound()

    @cache
    def upper_bound(self) -> int:
        if self._ub_method == UBComputeMethod.VERTEX:
            t, ub = get_vertex_cover_solution(g=self.g)
            self._vertex_cover_solution = t
            return int(ub)
        else:
            return super().upper_bound()

    @classmethod
    def name(cls) -> str:
        return 'Compact Natural Model'


class BottomUpNaturalModel(NaturalModel, BottomUpMBCModel):

    def _pre_solve(self):
        for b in range(self.upper_bound()):
            self.z[b].lb = 1

    @classmethod
    def name(cls) -> str:
        return 'Bottom-up Natural Model'


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

    def _add_independent_edges_constraints(self, edges: list[Var]):
        for biclique in self.bicliques:
            if biclique >= len(edges):
                return
            a, b = edges[biclique]
            self.x[a, b, biclique].lb = 1
            self.x[b, a, biclique].lb = 0

    @cache
    def lower_bound(self) -> int:
        if self._lb_method == LBComputeMethod.INDEPENDENT_EDGES and self._edge_fix:
            lb, edges = compute_lb_and_get_edges_by_independent_edges_method(g=self.g)
            self._add_independent_edges_constraints(edges=edges)
            return int(lb)
        else:
            return super().lower_bound()

    @classmethod
    def name(cls) -> str:
        return 'Extended Model'


class BottomUpExtendedModel(ExtendedModel, BottomUpMBCModel):

    def _pre_solve(self):
        for b in range(self.upper_bound()):
            self.z[b].lb = 1

    @classmethod
    def name(cls) -> str:
        return "Bottom-up Extended Model"
