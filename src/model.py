from datetime import datetime
from os import path, getcwd, pardir, mkdir
from typing import Type

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.base_model import MBCModel, IndepSetProvider
from src.util import get_graphs_in_store, GraphReport


class ModelV1(MBCModel):

    def _add_variables(self):
        self.z = self.m.addVars(self.bicliques, vtype=GRB.BINARY, name="z")
        self.x = self.m.addVars(self.g.nodes, self.bicliques, range(2), vtype=GRB.BINARY, name="x")
        self.y = self.m.addVars(
            self.directed.edges, self.bicliques, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")

    def _add_constraints(self):
        # set aliases
        m = self.m
        x = self.x
        y = self.y
        z = self.z
        bicliques = self.bicliques

        # symmetry-break constraints
        m.addConstrs(z[b] >= z[b + 1] for b in range(self.upper_bound - 1))
        # coupling constraints
        m.addConstrs(x[v, b, i] <= z[b] for v in self.g.nodes for b in bicliques for i in range(2))
        # covering constraints
        m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in self.g.edges)
        # node assignment constraints
        m.addConstrs(x[v, b, 0] + x[v, b, 1] <= z[b] for v in self.g.nodes for b in bicliques)
        # node and edge constraints
        for u, v in self.directed.edges:
            m.addConstrs(y[u, v, b] <= x[u, b, 0] for b in bicliques)
            m.addConstrs(y[u, v, b] <= x[v, b, 1] for b in bicliques)
            m.addConstrs(y[u, v, b] >= x[u, b, 0] + x[v, b, 1] - z[b] for b in bicliques)
        for u, v in self.directed_complement.edges:
            m.addConstrs(x[u, b, 0] + x[v, b, 1] <= z[b] for b in bicliques)

    def _set_objective(self):
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    @classmethod
    def name(cls) -> str:
        return 'v1'


class ModelV2(MBCModel, IndepSetProvider):

    def __init__(self, g: nx.Graph, g_name: str, use_relaxation: bool = True, **kwargs):
        self.use_relaxation = use_relaxation
        MBCModel.__init__(self, g, g_name, **kwargs)
        IndepSetProvider.__init__(self, g)

    def _add_variables(self):
        bicliques = self.bicliques

        self.w = self.m.addVars(
            self.indep_sets, bicliques, range(2), vtype=GRB.BINARY, name='w')
        self.x = self.m.addVars(
            self.g.nodes, bicliques, range(2), vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER,
            lb=0, ub=1, name='x')
        self.p = self.m.addVars(
            self.g.nodes, bicliques, range(2), vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER,
            lb=0, ub=1, name='p')
        self.y = self.m.addVars(
            self.directed.edges, bicliques, vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER,
            lb=0, ub=1, name='y')
        self.z = self.m.addVars(
            bicliques, vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER, lb=0, ub=1, name='z')

    def _add_constraints(self):
        # set aliases
        g = self.g
        m = self.m
        w = self.w
        x = self.x
        p = self.p
        y = self.y
        z = self.z
        bicliques = self.bicliques

        # independent set partition constraints
        m.addConstrs(w[s, b, 0] + w[s, b, 1] <= 1 for s in self.indep_sets for b in bicliques)
        m.addConstrs(w[s, b, d] <= x[u, b, d]
                     for s in self.indep_sets for u in s for b in bicliques for d in range(2))
        m.addConstrs(x[u, b, d] <= gp.quicksum(w[s, b, d] for s in self.indep_sets if u in s)
                     for u in g.nodes for b in bicliques for d in range(2))
        # p[u, b, d] = x[u, b, d] * (1 - x[u, b, 1-d])
        m.addConstrs(p[u, b, d] <= x[u, b, d] for u in g.nodes for b in bicliques for d in range(2))
        m.addConstrs(p[u, b, d] <= 1 - x[u, b, 1 - d] for u in g.nodes for b in bicliques for d in range(2))
        m.addConstrs(p[u, b, d] >= x[u, b, d] - x[u, b, 1 - d] for u in g.nodes for b in bicliques for d in range(2))
        # y[u, v, b] = p[u, b, 0] * p[v, b, 1]
        m.addConstrs(y[u, v, b] <= p[u, b, 0] for u, v in self.directed.edges for b in bicliques)
        m.addConstrs(y[u, v, b] <= p[v, b, 1] for u, v in self.directed.edges for b in bicliques)
        m.addConstrs(y[u, v, b] >= p[u, b, 0] + p[v, b, 1] - 1 for u, v in self.directed.edges for b in bicliques)
        # cover constraints
        m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in g.edges)
        m.addConstrs(y[u, v, b] <= z[b] for u, v in self.directed.edges for b in bicliques)
        m.addConstrs(p[u, b, 0] + p[v, b, 1] <= 1 for u, v in nx.complement(self.directed).edges for b in bicliques)
        # symmetry break
        m.addConstrs(z[b + 1] <= z[b] for b in range(self.upper_bound - 1))

    def _set_objective(self):
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    @classmethod
    def name(cls) -> str:
        return 'v2-0'


class ModelV21(MBCModel, IndepSetProvider):

    def __init__(self, g: nx.Graph, g_name: str, use_relaxation: bool = True, **kwargs):
        self.use_relaxation = use_relaxation
        MBCModel.__init__(self, g, g_name, **kwargs)
        IndepSetProvider.__init__(self, g)

    def _add_variables(self):
        bicliques = self.bicliques

        self.w = self.m.addVars(self.indep_sets, bicliques, range(2), vtype=GRB.BINARY, name='w')
        self.x = self.m.addVars(
            self.g.nodes, bicliques, range(2),
            vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER, lb=0, ub=1, name='x')
        self.y = self.m.addVars(
            self.directed.edges, bicliques,
            vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER, lb=0, ub=1, name='y')
        self.z = self.m.addVars(
            bicliques,
            vtype=GRB.CONTINUOUS if self.use_relaxation else GRB.INTEGER, lb=0, ub=1, name='z')

    def _add_constraints(self):
        # set aliases
        g, m, w, x, y, z, bicliques, indep_sets, directed, upper_bound = (
            self.g, self.m, self.w, self.x, self.y, self.z, self.bicliques,
            self.indep_sets, self.directed, self.upper_bound
        )

        # independent set partition constraints
        m.addConstrs(w[s, b, 0] + w[s, b, 1] <= 1 for s in self.indep_sets for b in bicliques)
        m.addConstrs(w[s, b, d] <= x[u, b, d]
                     for s in self.indep_sets for u in s for b in bicliques for d in range(2))
        m.addConstrs(x[u, b, d] <= gp.quicksum(w[s, b, d] for s in self.indep_sets if u in s)
                     for u in g.nodes for b in bicliques for d in range(2))
        # y[u, v, b] = x[u, b, 0] ^ x[v, b, 1]
        m.addConstrs(y[u, v, b] <= x[u, b, 0] + x[v, b, 1] for u, v in self.directed.edges for b in bicliques)
        m.addConstrs(y[u, v, b] <= 2 - x[u, b, 0] - x[v, b, 1] for u, v in self.directed.edges for b in bicliques)
        m.addConstrs(y[u, v, b] >= x[u, b, 0] - x[v, b, 1] for u, v in self.directed.edges for b in bicliques)
        m.addConstrs(y[u, v, b] >= x[v, b, 1] - x[u, b, 0] for u, v in self.directed.edges for b in bicliques)
        # cover constraints
        m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in g.edges)
        m.addConstrs(y[u, v, b] <= z[b] for u, v in self.directed.edges for b in bicliques)
        # symmetry break
        m.addConstrs(z[b + 1] <= z[b] for b in range(self.upper_bound - 1))

    def _set_objective(self):
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    @classmethod
    def name(cls) -> str:
        return 'v2-1'


def create_and_save_model_comparison_report(report_name: str, model_clss: list[Type[MBCModel]], **kwargs):
    # get log directory
    dir_parent = path.abspath(path.join(getcwd(), pardir))
    dir_logs = path.join(dir_parent, 'logs')
    if not path.isdir(dir_logs):
        mkdir(dir_logs)
    ts = int(datetime.now().timestamp())
    dir_ts_logs = path.join(dir_logs, report_name + '-' + str(ts))
    mkdir(dir_ts_logs)

    # create report
    report = GraphReport(name=report_name)
    report.add_properties([model_cls.name() for model_cls in model_clss])
    for g, g_name in get_graphs_in_store(**kwargs):
        report.add_graph_data(g, g_name)
        for model_cls in model_clss:
            model = model_cls(g=g, g_name=g_name, dir_logs=dir_ts_logs)
            report.add_property_values_from_function(p_name=model_cls.name(), f=model.solve)
    # save report
    report.save_csv()


if __name__ == '__main__':
    create_and_save_model_comparison_report(
        report_name='simple', fname_regex='simple', model_clss=[ModelV1, ModelV2])
