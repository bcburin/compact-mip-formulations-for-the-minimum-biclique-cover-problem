from __future__ import annotations
from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import combinations

from os import path
from typing import TextIO, Iterable

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.bc_bounds import find_bc_upper_bound, find_bc_lower_bound, UBComputeMethod, get_vertex_cover_solution, \
    LBComputeMethod, compute_lb_and_get_edges_by_independent_edges_method
from src.config import RunConfig, ReportConfig


class MBCModel(ABC):

    def __init__(
             self,
             g: nx.Graph,
             g_name: str,
             config: RunConfig,
             default_config: ReportConfig,
             log_to_console: bool = True,
             dir_logs: str = None):
        self.g = g
        self.g_name = g_name
        self._config = config
        self._default_config = default_config
        self._log_to_console = log_to_console
        self._dir_logs = dir_logs
        self._logging = bool(dir_logs)
        # config properties
        self._time_limit = (
            self._config.time_limit) \
            if self._config.time_limit is not None \
            else self._default_config.default_time_limit
        self._lb_method = self._config.lb_method or self._default_config.default_lb_method
        self._ub_method = self._config.ub_method or self._default_config.default_ub_method
        self._edge_fix = (
            self._config.edge_fix) \
            if self._config.edge_fix is not None \
            else self._default_config.default_edge_fix
        self._warm_start = (
            self._config.warm_start) \
            if self._config.warm_start is not None \
            else self._default_config.default_warm_start
        self._conflict_inequalities = (
            self._config.conflict_inequalities
            if self._config.conflict_inequalities is not None
            else self._default_config.default_conflict_inequalities)
        self._common_neighbor_inequalities = (
            self._config.common_neighbor_inequalities
            if self._config.common_neighbor_inequalities is not None
            else self._default_config.default_common_neighbor_inequalities)
        self._use_callback = (
            self._config.use_callback) \
            if self._config.use_callback is not None else (
            self._default_config.default_use_callback)
        self._callback = None
        # model
        self._init_model()
        # create log files
        if self._logging:
            self._open_files(dir_logs=dir_logs)
        # parameters
        self._set_parameters()
        # initial state
        self._solved = False

    def _log_message(self, msg: str):
        msg = f"[{self.__class__.__name__}] " + msg
        if self._logging:
            self._log_res.write(msg)
        print(msg)

    def __del__(self):
        if hasattr(self, '_files_open') and self._files_open:
            self._close_files()

    def _open_files(self, dir_logs: str):
        model_graph_name = self.name() + '_' + self.g_name
        self._log_grb = open(path.join(dir_logs, model_graph_name + '_gurobi.log'), 'w+')
        self._log_res = open(path.join(dir_logs, model_graph_name + '_result.log'), 'w+')
        self._write_headers(self._log_res)
        self._files_open = True

    def _close_files(self):
        self._log_grb.close()
        self._log_res.close()

    def _set_parameters(self):
        if self._time_limit is not None:
            self.m.Params.TimeLimit = self._time_limit
        if not self._log_to_console:
            self.m.Params.LogToConsole = 0
        if self._logging:
            self.m.Params.LogFile = self._log_grb.name

    def _init_model(self):
        self.m = gp.Model()
        self._add_variables()
        self._set_objective()
        self._add_constraints()

    def _write_headers(self, log: TextIO):
        log.write(f'GRAPH: {self.g_name}\n')
        log.write(f'NODES: {len(self.g.nodes)}\n')
        log.write(f'EDGES: {len(self.g.edges)}\n\n')

    @abstractmethod
    def _add_variables(self):
        self.z = self.m.addVars(self.bicliques, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name="z")

    @abstractmethod
    def _add_constraints(self):
        ...

    @abstractmethod
    def _set_objective(self):
        ...

    @abstractmethod
    def _check_biclique_cover(self) -> bool:
        ...

    @classmethod
    @abstractmethod
    def name(cls) -> str:
        ...

    def _pre_solve(self):
        ...

    def _post_solve(self):
        # check and log solution
        if not self._logging:
            return
        if self.m.status == GRB.OPTIMAL:
            # check solution
            self._log_message(f'Is it a biclique cover? {"Yes" if self._check_biclique_cover() else "No"}.\n')
        elif self.m.status == GRB.TIME_LIMIT:
            self._log_message(f'Model reached time limit of {self._time_limit} seconds.\n')
        elif self.m.status == GRB.INFEASIBLE:
            self._log_message('Model is unfeasible.\n')
        else:
            self._log_message(f'Status code: {self.m.status}\n')
        if self._logging:
            self._close_files()

    @property
    @cache
    def directed(self) -> nx.DiGraph:
        return nx.DiGraph(self.g)

    @cached_property
    def complement(self) -> nx.Graph:
        return nx.complement(self.g)

    @cached_property
    def directed_complement(self) -> nx.Graph:
        return nx.complement(self.directed)

    @cached_property
    def power_graph(self) -> nx.Graph:
        return nx.power(self.g, k=2)

    @cache
    def get_lb_and_indep_edges(self) -> tuple[int, list[tuple[int, int]]]:
        lb, edges = compute_lb_and_get_edges_by_independent_edges_method(g=self.g)
        return int(lb), edges

    @cache
    def get_ub_and_vertex_cover(self) -> tuple[int, list[int]]:
        t, ub = get_vertex_cover_solution(g=self.g)
        return int(ub), t

    @cache
    def lower_bound(self) -> int:
        if self._lb_method == LBComputeMethod.INDEPENDENT_EDGES:
            lb, _ = self.get_lb_and_indep_edges()
            return lb
        else:
            return int(find_bc_lower_bound(self.g, self._lb_method)) if self._lb_method else 1

    @cache
    def upper_bound(self) -> int:
        if self._ub_method == UBComputeMethod.VERTEX:
            ub, _ = self.get_ub_and_vertex_cover()
            return ub
        else:
            return int(find_bc_upper_bound(self.g, self._ub_method)) if self._ub_method else self.g.edges

    @property
    def bicliques(self) -> Iterable:
        return range(self.upper_bound())

    @cache
    def get_disjoint_edges(self) -> set[tuple[int, int, int, int]]:
        disjoint_edges = set()
        for e, f in combinations(self.g.edges, r=2):
            a, b = e
            c, d = f
            if {a, b}.isdisjoint({c, d}):
                cr1, cr2 = 0, 0
                if self.g.has_edge(a, d):
                    cr1 += 1
                if self.g.has_edge(c, b):
                    cr1 += 1
                if self.g.has_edge(a, c):
                    cr2 += 1
                if self.g.has_edge(b, d):
                    cr2 += 1
                disjoint_edges.add((e, f, cr1, cr2))
        return disjoint_edges

    def _can_add_indep_edges_constraints(self) -> bool:
        _, indep_edges = self.get_lb_and_indep_edges()
        return bool(indep_edges) and self._edge_fix and not self._warm_start

    def infeasible_or_unsolved(self) -> bool:
        return self.m.status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD, GRB.INTERRUPTED, GRB.LOADED]

    def solve(self) -> float | None:
        # custom pre-solve with default implementation
        self._pre_solve()
        # optimization process
        self._log_message(f'Solving for graph {self.g_name}...')
        # if self._callback is not None:
        #     # noinspection PyArgumentList
        #     self.m.optimize(self._callback)
        # else:
        #     self.m.optimize()
        self._solved = True
        self.m.write(f"{self._config.name}.mps")
        # custom post-solve with default implementation
        self._post_solve()
        # return obj val
        if not self.infeasible_or_unsolved():
            return self.m.objVal
