from abc import ABC, abstractmethod
from functools import cache
from itertools import combinations

from os import path
from typing import TextIO, Iterable

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.bc_bounds import LBComputeMethod, UBComputeMethod, find_bc_upper_bound


class MBCModel(ABC):

    def __init__(
             self,
             g: nx.Graph,
             g_name: str,
             lb_method: LBComputeMethod | None = None,
             ub_method: UBComputeMethod = UBComputeMethod.VERTEX,
             time_limit: int = None,
             log_to_console: bool = True,
             dir_logs: str = None):
        self.g = g
        self.g_name = g_name
        self.lb_method = lb_method
        self.ub_method = ub_method
        self._logging = bool(dir_logs)
        self.time_limit = time_limit
        self.log_to_console = log_to_console
        # model
        self._init_model()
        # create log files
        if self._logging:
            self._open_files(dir_logs=dir_logs)
        # parameters
        self._set_parameters()
        # initial state
        self._solved = False

    def __del__(self):
        if hasattr(self, '_files_open') and self._files_open:
            self._close_files()

    def _open_files(self, dir_logs: str):
        model_graph_name = self.name() + '_' + self.g_name
        self.log_grb = open(path.join(dir_logs, model_graph_name + '_gurobi.log'), 'w+')
        self.log_res = open(path.join(dir_logs, model_graph_name + '_result.log'), 'w+')
        self._write_headers(self.log_res)
        self._files_open = True

    def _close_files(self):
        if self._logging:
            self.log_grb.close()
            self.log_res.close()

    def _set_parameters(self):
        if self.time_limit is not None:
            self.m.Params.TimeLimit = self.time_limit
        if not self.log_to_console:
            self.m.Params.LogToConsole = 0
        if self._logging:
            self.m.Params.LogFile = self.log_grb.name

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
        ...

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
        # add lower bound constraints
        if self.lb_method:
            for b in range(self.lower_bound()):
                self.z[b].lb = 1

    def _post_solve(self):
        # check and log solution
        if not self._logging:
            return
        if self.m.status == GRB.OPTIMAL:
            # check solution
            self.log_res.write(f'IS BICLIQUE COVER? {"Y" if self._check_biclique_cover() else "N"}')
        elif self.m.status == GRB.TIME_LIMIT:
            self.log_res.write(f'\nMODEL REACHED TIME LIMIT OF {self.time_limit} SECONDS\n\n')
        elif self.m.status == GRB.INFEASIBLE:
            self.log_res.write('\nMODEL HAS BEEN MARKED AS UNFEASIBLE\n\n')
        else:
            self.log_res.write(f'\nSTATUS CODE: {self.m.status}\n\n')
        self._close_files()

    @property
    @cache
    def directed(self) -> nx.DiGraph:
        return nx.DiGraph(self.g)

    @property
    @cache
    def complement(self) -> nx.Graph:
        return nx.complement(self.g)

    @property
    @cache
    def directed_complement(self) -> nx.Graph:
        return nx.complement(self.directed)

    @cache
    def lower_bound(self) -> int:
        return 1

    @cache
    def upper_bound(self) -> int:
        return int(find_bc_upper_bound(self.g, self.ub_method))

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

    def solve(self):
        # custom pre-solve with default implementation
        self._pre_solve()
        # optimization process
        self.m.optimize()
        self._solved = True

        # custom post-solve with default implementation
        self._post_solve()
        # return obj val
        return self.m.objVal
