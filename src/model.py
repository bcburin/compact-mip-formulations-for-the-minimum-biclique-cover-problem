from datetime import datetime
from os import path, getcwd, pardir, mkdir
from typing import Type

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.base_model import MBCModel
from src.util import get_graphs_in_store, GraphReport


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
        m.addConstrs(x[e, i] <= z[i] for e in self.g.edges for i in self.bicliques)
        # 4c
        m.addConstrs(gp.quicksum(x[e, i] for i in self.bicliques) for e in self.g.edges)
        # 4e
        for cycle in nx.simple_cycles(self.g, length_bound=3):
            if len(cycle) != 3:
                continue
            cycle_edges = [[cycle[0], cycle[1]], [cycle[1], cycle[2]], [cycle[2], cycle[0]]]
            m.addConstrs(gp.quicksum(x[e, i] for e in cycle_edges) <= 2*z[i] for i in self.bicliques)

        for e, f, cr1, cr2 in self.get_disjoint_edges():
            a, b = e
            c, d = f
            # 4d
            if cr1 < 2 and cr2 < 2:
                m.addConstrs(x[e, i] + x[f, i] <= z[i] for i in self.bicliques)
            # 4f
            if cr1 == 2 and cr2 < 2:
                m.addConstrs(x[e, i] + x[f, i] <= z[i] + x[[a, d], i] for i in self.bicliques)
                m.addConstrs(x[e, i] + x[f, i] <= z[i] + x[[b, c], i] for i in self.bicliques)
            # 4g
            if cr1 < 2 and cr2 == 2:
                m.addConstrs(x[e, i] + x[f, i] <= z[i] + x[[a, c], i] for i in self.bicliques)
                m.addConstrs(x[e, i] + x[f, i] <= z[i] + x[[b, d], i] for i in self.bicliques)
            # 4h
            if cr1 == 2 and cr2 == 2:
                m.addConstrs(x[e, i] + x[f, i] <= z[i] + x[[a, c], i] + x[[a, d], i] for i in self.bicliques)
                m.addConstrs(x[e, i] + x[f, i] <= z[i] + x[[b, c], i] + x[[b, d], i] for i in self.bicliques)
        # 4i
        m.addConstrs(z[i] >= z[i + 1] for i in range(self.upper_bound - 1))

    @classmethod
    def name(cls) -> str:
        return 'Compact Natural Model'


def create_and_save_model_comparison_report(
        report_name: str, model_clss: list[Type[MBCModel]], time_limit: int = None,
        suppress_ts_in_report_name: bool = True, **kwargs):
    # get log directory
    dir_parent = path.abspath(path.join(getcwd(), pardir))
    dir_logs = path.join(dir_parent, 'logs')
    if not path.isdir(dir_logs):
        mkdir(dir_logs)
    ts = int(datetime.now().timestamp())
    dir_ts_logs = path.join(dir_logs, report_name + '-' + str(ts))
    mkdir(dir_ts_logs)

    # create report
    report_name = report_name if suppress_ts_in_report_name else report_name + '-' + str(ts)
    report = GraphReport(name=report_name)
    report.add_properties([model_cls.name() for model_cls in model_clss])
    for g, g_name in get_graphs_in_store(**kwargs):
        report.add_graph_data(g, g_name)
        for model_cls in model_clss:
            model = model_cls(g=g, g_name=g_name, dir_logs=dir_ts_logs, time_limit=time_limit)
            report.add_property_values_from_function(p_name=model_cls.name(), f=model.solve)
    # save report
    report.save_csv()
