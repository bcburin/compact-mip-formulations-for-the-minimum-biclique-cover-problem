from __future__ import annotations
from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import combinations

from typing import Iterable

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.bc_bounds import find_bc_upper_bound, find_bc_lower_bound, UBComputeMethod, get_vertex_cover_solution, \
    LBComputeMethod, compute_lb_and_get_edges_by_independent_edges_method
from src.config import RunConfig
from src.exceptions import AlreadySolvedError, NotYetSolvedError, UnfeasibleSolutionError


class BaseSolver(ABC):

    def __init__(self):
        self._solved: bool = False
        self._feasible: bool | None = None

    def _pre_solve(self) -> None:
        ...

    def _post_solve(self) -> None:
        ...

    @abstractmethod
    def _do_solve(self) -> None:
        ...

    @abstractmethod
    def _do_get_solution(self) -> float:
        ...

    @abstractmethod
    def is_feasible(self) -> bool:
        ...

    def solve(self) -> bool:
        if self._solved:
            raise AlreadySolvedError()
        self._pre_solve()
        self._do_solve()
        self._post_solve()
        self._solved = True
        self._feasible = self.is_feasible()
        return self._feasible

    @cached_property
    def solution(self) -> float:
        if not self._solved:
            raise NotYetSolvedError()
        if self._feasible is None or not self._feasible:
            raise UnfeasibleSolutionError()
        return self._do_get_solution()


class BaseGpSolver(BaseSolver):
    _INFEASIBLE_STATUSES = [GRB.INFEASIBLE, GRB.INF_OR_UNBD, GRB.INTERRUPTED, GRB.LOADED]

    def __init__(self):
        super().__init__()
        self.m = gp.Model()
        self._add_variables()
        self._set_objective()
        self._add_constraints()

    @abstractmethod
    def _add_variables(self) -> None:
        ...

    @abstractmethod
    def _add_constraints(self) -> None:
        ...

    @abstractmethod
    def _set_objective(self) -> None:
        ...

    @abstractmethod
    def _set_parameters(self) -> None:
        ...

    def _do_solve(self) -> None:
        self.m.optimize()

    def _do_get_solution(self) -> float:
        return self.m.objVal

    def is_feasible(self) -> bool:
        return self.m.status not in self._INFEASIBLE_STATUSES


class BaseMinimumBicliqueCoverSolver(BaseSolver):

    def __init__(self, g: nx.Graph, config: RunConfig):
        self._g = g
        self._config = config
        # model cache
        self._lb: int | None = None
        self._ub: int | None = None
        self._indep_edges: list[tuple[int, int]] | None = None
        self._vertex_cover: list | None = None
        super().__init__()

    def _guarantee_compute_lb_and_indep_edges(self) -> None:
        if self._lb and self._indep_edges:
            return
        lb, edges = compute_lb_and_get_edges_by_independent_edges_method(g=self.graph)
        self._lb = int(lb)
        self._indep_edges = edges

    def _guarantee_compute_ub_and_vertex_cover(self) -> None:
        if self._ub and self._vertex_cover:
            return
        t, ub = get_vertex_cover_solution(g=self.graph)
        self._ub = int(ub)
        self._vertex_cover = t
        return

    @property
    def lower_bound(self) -> int:
        if self._lb:
            return self._lb
        print(f'Finding lower bound for {self._config.resolved_gname}')
        method = self._config.lb_method
        if method == LBComputeMethod.INDEPENDENT_EDGES:
            self._guarantee_compute_lb_and_indep_edges()
            return self._lb
        else:
            return int(find_bc_lower_bound(self.graph, method)) if method else 1

    @property
    def upper_bound(self) -> int:
        if self._ub:
            return self._ub
        print(f'Finding upper bound for {self._config.resolved_gname}')
        method = self._config.ub_method
        if method == UBComputeMethod.VERTEX:
            self._guarantee_compute_ub_and_vertex_cover()
            return self._ub
        else:
            return int(find_bc_upper_bound(self.graph, method)) if method else self.graph.edges

    @cache
    def get_disjoint_edges(self) -> set[tuple[int, int, int, int]]:
        disjoint_edges = set()
        for e, f in combinations(self.graph.edges, r=2):
            a, b = e
            c, d = f
            if {a, b}.isdisjoint({c, d}):
                cr1, cr2 = 0, 0
                if self.graph.has_edge(a, d):
                    cr1 += 1
                if self.graph.has_edge(c, b):
                    cr1 += 1
                if self.graph.has_edge(a, c):
                    cr2 += 1
                if self.graph.has_edge(b, d):
                    cr2 += 1
                disjoint_edges.add((e, f, cr1, cr2))
        return disjoint_edges

    def _can_add_indep_edges_constraints(self) -> bool:
        if self._config.warm_start or not self._config.edge_fix:
            return False
        self._guarantee_compute_lb_and_indep_edges()
        return bool(self._indep_edges)

    @property
    def graph(self) -> nx.Graph:
        return self._g

    @cached_property
    def directed(self) -> nx.DiGraph:
        return nx.DiGraph(self.graph)

    @cached_property
    def complement(self) -> nx.Graph:
        return nx.complement(self.graph)

    @cached_property
    def directed_complement(self) -> nx.Graph:
        return nx.complement(self.directed)

    @cached_property
    def power_graph(self) -> nx.Graph:
        return nx.power(self.graph, k=2)

    @cached_property
    def bicliques(self) -> Iterable:
        return range(self.upper_bound)

    @abstractmethod
    def _check_biclique_cover(self) -> bool:
        ...


def _log_message(msg, *args, **kwargs):
    print(msg, args, kwargs)


class BaseMinimumBicliqueCoverGpSolver(BaseMinimumBicliqueCoverSolver, BaseGpSolver, ABC):

    def __init__(self, g: nx.Graph, config: RunConfig):
        BaseMinimumBicliqueCoverSolver.__init__(self, g, config)
        BaseGpSolver.__init__(self)

    def _set_parameters(self):
        if self._config.time_limit is not None:
            self.m.Params.TimeLimit = self._config.time_limit

    def _post_solve(self):
        # check and log the solution
        if self.m.status == GRB.OPTIMAL:
            # check the solution
            _log_message(f'Is it a biclique cover? {"Yes" if self._check_biclique_cover() else "No"}.\n')
        elif self.m.status == GRB.TIME_LIMIT:
            _log_message(f'Model reached time limit of {self._config.time_limit} seconds.\n')
        elif self.m.status == GRB.INFEASIBLE:
            _log_message('Model is unfeasible.\n')
        else:
            _log_message(f'Status code: {self.m.status}\n')

    @classmethod
    def name(cls):
        return "Base Class"
