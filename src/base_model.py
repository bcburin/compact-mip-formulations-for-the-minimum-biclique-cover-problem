from abc import ABC, abstractmethod
from functools import cache
from itertools import combinations

from os import path
from typing import TextIO, Iterable

import gurobipy as gp
import networkx as nx
import pandas as pd
from gurobipy import GRB
from pandas import DataFrame

from src.bc_bounds import LBComputeMethod, UBComputeMethod, find_bc_upper_bound
from src.independent_set import get_set_of_maximal_independent_sets


class MBCModel(ABC):

    def __init__(
             self,
             g: nx.Graph,
             g_name: str,
             lb_method: LBComputeMethod | None = None,
             ub_method: UBComputeMethod = UBComputeMethod.VERTEX,
             time_limit: int = None,
             recursive_solution: bool = False,
             log_to_console: bool = True,
             dir_logs: str = None):
        self.g = g
        self.g_name = g_name
        self.lb_method = lb_method
        self.ub_method = ub_method
        self._logging = bool(dir_logs)
        self.time_limit = time_limit
        self.log_to_console = log_to_console
        self.y_is_directed = True
        self._is_recursive_solution = recursive_solution
        # model
        self.m = gp.Model()
        self._add_variables()
        self._set_objective()
        self._add_constraints()
        # create log files
        if self._logging:
            model_graph_name = self.name() + '_' + self.g_name
            self.log_grb = open(path.join(dir_logs, model_graph_name + '_gurobi.log'), 'w+')
            self.log_res = open(path.join(dir_logs, model_graph_name + '_result.log'), 'w+')
            self._write_headers(self.log_res)
        # parameters
        self._set_parameters()
        # initial state
        self._solved = False

    def __del__(self):
        self._close_files()

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

    @classmethod
    @abstractmethod
    def name(cls) -> str:
        ...

    def _pre_solve(self):
        # add lower bound constraints
        if self.lower_bound:
            for b in range(self.lower_bound):
                self.z[b].lb = 1

    def _post_solve(self):
        # check and log solution
        if not self._logging:
            return
        if self.m.status == GRB.OPTIMAL:
            # construct dataframe with bicliques
            df_bicliques = self.get_bicliques()
            # check solution
            self.log_res.write(f'IS BICLIQUE COVER? {"Y" if self.is_biclique_cover(df_bicliques) else "N"}')
            # save solution
            self.log_res.write(f'\nBICLIQUE COVER DATAFRAME:\n{df_bicliques}\n\n')
        elif self.m.status == GRB.TIME_LIMIT:
            self.log_res.write(f'\nMODEL REACHED TIME LIMIT OF {self.time_limit} SECONDS\n\n')
        elif self.m.status == GRB.INFEASIBLE:
            self.log_res.write('\nMODEL HAS BEEN MARKED AS UNFEASIBLE\n\n')
        else:
            self.log_res.write(f'\nSTATUS CODE: {self.m.status}\n\n')
        self._close_files()

    @cache
    def get_bicliques(self) -> pd.DataFrame:
        """
        This method builds a Pandas dataframe containing the solution of a biclique cover. For this to work for
        a model, it is enough that said model define z variables representing bicliques and y variables representing
        whether an arc or edge belongs to a biclique.

        :return: a binary dataframe, where each column is a biclique and each line is indexed by a node. Each cell can
                 hold the values 0, 1, or 2. If the node does not belong to the biclique, the value is zero. If the node
                 belongs to the first partition of the biclique, the value is one; likewise, if the node belongs to the
                 second partition, the value of the cell is 2.
        :raises RuntimeError: if the model was not solved yet or if it does not define z or y variables.
        """
        if not self._solved:
            raise RuntimeError('Model not solved yet.')
        if not hasattr(self, 'z') or not hasattr(self, 'y'):
            raise RuntimeError('Model does not define z or y variables. They are required to build the '
                               'biclique dataframe.')
        if not hasattr(self, 'x') and not self.y_is_directed:
            raise RuntimeError('Variable x is required for models with undirected y variable.')

        df_bicliques = DataFrame(
            {b: [0] * len(self.g.nodes) for b in self.bicliques if self.z[b].X == 1}, index=[u for u in self.g.nodes])
        if self.y_is_directed:
            for e in self.g.edges:
                u, v = e
                for b in self.bicliques:
                    if self.z[b].X == 1:
                        if self.y[u, v, b].X == 1:
                            df_bicliques.loc[u, b] = 1
                            df_bicliques.loc[v, b] = 2
                        if self.y[v, u, b].X == 1:
                            df_bicliques.loc[v, b] = 1
                            df_bicliques.loc[u, b] = 2
        else:
            for u in self.g.nodes:
                for b in self.bicliques:
                    if self.x[u, b, 0].X == 1 and self.x[u, b, 1].X == 0:
                        df_bicliques.loc[u, b] = 1
                    if self.x[u, b, 0].X == 0 and self.x[u, b, 1].X == 1:
                        df_bicliques.loc[u, b] = 2
        return df_bicliques

    def is_biclique_cover(self, df: pd.DataFrame) -> bool:
        """
        This method determines whether the solution to the minimum biclique cover problem found by a model is actually
        a valid biclique cover.

        :param df: dataframe containing information of each biclique in the biclique cover, as found by the model and
                   generated by the function get_biclique_cover_dataframe.
        :return: True if the dataframe represent a valid biclique cover; False otherwise.
        """
        # check if it is a cover
        if not df.any(axis=1).all():
            return False
        # check if each column is, indeed, a biclique
        for biclique in df.columns:
            partition1: list[int] = df[df[biclique] == 1].index.tolist()
            partition2: list[int] = df[df[biclique] == 2].index.tolist()
            if any(not self.g.has_edge(u, v) for u in partition1 for v in partition2):
                return False
        return True

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

    @property
    @cache
    def lower_bound(self) -> int:
        return 1

    @property
    @cache
    def upper_bound(self) -> int:
        return int(find_bc_upper_bound(self.g, self.ub_method))

    @property
    def bicliques(self) -> Iterable:
        return range(self.upper_bound)

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
        # custom pre-solve
        self._pre_solve()
        # optimization process
        self.m.optimize()
        self._solved = True

        # custom post-solve with default implementation
        self._post_solve()
        # return obj val
        if self.m.status != GRB.OPTIMAL:
            return None
        return self.m.objVal
