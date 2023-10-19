from datetime import datetime
from os import path, getcwd, pardir, mkdir
from typing import Iterable, TextIO

import gurobipy as gp
import networkx as nx
import pandas as pd
from gurobipy import GRB
from pandas import DataFrame

from src.bc_bounds import find_bc_upper_bound, UBComputeMethod
from src.independent_set import get_set_of_maximal_independent_sets
from src.util import get_graphs_in_store, GraphReport


def solve_bc_v1(g: nx.Graph, log_gb: TextIO = None, log_res: TextIO = None, time_limit: int = None):
    lb = 1
    ub = int(find_bc_upper_bound(g, method=UBComputeMethod.VERTEX))
    dg = nx.DiGraph(g)
    cg = nx.complement(dg)
    bicliques = range(ub)

    # create model
    m = gp.Model()
    # set parameters
    if time_limit is not None:
        m.Params.TimeLimit = time_limit
    if log_gb is not None:
        m.Params.LogFile = log_gb.name
        # m.Params.LogToConsole = 0
    # define variables
    z = m.addVars(bicliques, vtype=GRB.BINARY, name="z")
    x = m.addVars(g.nodes, bicliques, range(2), vtype=GRB.BINARY, name="x")
    y = m.addVars(dg.edges, bicliques, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
    # define objective function
    m.setObjective(gp.quicksum(z), sense=GRB.MINIMIZE)
    # symmetry-break constraints
    m.addConstrs(z[b] >= z[b + 1] for b in range(ub - 1))
    # coupling constraints
    m.addConstrs(x[v, b, i] <= z[b] for v in g.nodes for b in bicliques for i in range(2))
    # covering constraints
    m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in g.edges)
    # node assignment constraints
    m.addConstrs(x[v, b, 0] + x[v, b, 1] <= z[b] for v in g.nodes for b in bicliques)
    # node and edge constraints
    for u, v in dg.edges:
        m.addConstrs(y[u, v, b] <= x[u, b, 0] for b in bicliques)
        m.addConstrs(y[u, v, b] <= x[v, b, 1] for b in bicliques)
        m.addConstrs(y[u, v, b] >= x[u, b, 0] + x[v, b, 1] - z[b] for b in bicliques)
    for u, v in cg.edges:
        m.addConstrs(x[u, b, 0] + x[v, b, 1] <= z[b] for b in bicliques)
    # lower bound
    for b in range(lb):
        z[b].lb = 1

    m.optimize()

    if m.status == GRB.OPTIMAL:
        # construct dataframe with bicliques
        df_bicliques = get_biclique_cover_dataframe(g=g, bicliques=bicliques, z=z, y=y)
        # check solution
        log_res.write(f'IS BICLIQUE COVER? {"Y" if is_biclique_cover(g, df_bicliques) else "N"}')
        # save solution
        log_res.write(f'\nBICLIQUE COVER DATAFRAME:\n{df_bicliques}\n\n')

        return m.objVal

    if m.status == GRB.TIME_LIMIT:
        log_res.write(f'\nMODEL REACHED TIME LIMIT OF {time_limit} SECONDS\n')
    else:
        log_res.write(f'\nTHERE WAS AN UNEXPECTED ERROR\n')


def solve_bc_v2(g: nx.Graph, use_relaxation: bool = True, log_gb: TextIO = None, log_res: TextIO = None,
                time_limit: int = None):
    # use frozen sets so the independent sets are hashable
    indep_sets = frozenset(frozenset(indep_set) for indep_set in get_set_of_maximal_independent_sets(g))
    lb = 1
    ub = int(find_bc_upper_bound(g, method=UBComputeMethod.VERTEX))
    dg = nx.DiGraph(g)
    bicliques = range(ub)

    # create model
    m = gp.Model()
    # set parameters
    if time_limit is not None:
        m.Params.TimeLimit = time_limit
    if log_gb is not None:
        m.Params.LogFile = log_gb.name
        # m.Params.LogToConsole = 0
    # define variables
    w = m.addVars(indep_sets, bicliques, range(2), vtype=GRB.BINARY, name='w')
    x = m.addVars(g.nodes, bicliques, range(2),
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='x')
    p = m.addVars(g.nodes, bicliques, range(2),
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='p')
    y = m.addVars(dg.edges, bicliques,
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='y')
    z = m.addVars(bicliques,
                  vtype=GRB.CONTINUOUS if use_relaxation else GRB.INTEGER, lb=0, ub=1, name='z')
    # objective function
    m.setObjective(gp.quicksum(z), GRB.MINIMIZE)
    # independent set partition constraints
    m.addConstrs(w[s, b, 0] + w[s, b, 1] <= 1 for s in indep_sets for b in bicliques)
    m.addConstrs(w[s, b, d] <= x[u, b, d]
                 for s in indep_sets for u in s for b in bicliques for d in range(2))
    m.addConstrs(x[u, b, d] <= gp.quicksum(w[s, b, d] for s in indep_sets if u in s)
                 for u in g.nodes for b in bicliques for d in range(2))
    # p[u, b, d] = x[u, b, d] * (1 - x[u, b, 1-d])
    m.addConstrs(p[u, b, d] <= x[u, b, d] for u in g.nodes for b in bicliques for d in range(2))
    m.addConstrs(p[u, b, d] <= 1 - x[u, b, 1-d] for u in g.nodes for b in bicliques for d in range(2))
    m.addConstrs(p[u, b, d] >= x[u, b, d] - x[u, b, 1-d] for u in g.nodes for b in bicliques for d in range(2))
    # y[u, v, b] = p[u, b, 0] * p[v, b, 1]
    m.addConstrs(y[u, v, b] <= p[u, b, 0] for u, v in dg.edges for b in bicliques)
    m.addConstrs(y[u, v, b] <= p[v, b, 1] for u, v in dg.edges for b in bicliques)
    m.addConstrs(y[u, v, b] >= p[u, b, 0] + p[v, b, 1] - 1 for u, v in dg.edges for b in bicliques)
    # cover constraints
    m.addConstrs(gp.quicksum(y[u, v, b] + y[v, u, b] for b in bicliques) >= 1 for u, v in g.edges)
    m.addConstrs(y[u, v, b] <= z[b] for u, v in dg.edges for b in bicliques)
    m.addConstrs(p[u, b, 0] + p[v, b, 1] <= 1 for u, v in nx.complement(dg).edges for b in bicliques)
    # symmetry break
    m.addConstrs(z[b + 1] <= z[b] for b in range(ub-1))
    # lower bound
    for b in range(lb):
        z[b].lb = 1

    m.optimize()

    if m.status == GRB.OPTIMAL:
        # construct dataframe with bicliques
        df_bicliques = get_biclique_cover_dataframe(g=g, bicliques=bicliques, z=z, y=y)
        # check solution
        log_res.write(f'\nIS BICLIQUE COVER? {"Y" if is_biclique_cover(g, df_bicliques) else "N"}\n')
        # save solution
        log_res.write(f'\nBICLIQUE COVER DATAFRAME:\n{df_bicliques}\n\n')

        return m.objVal

    if m.status == GRB.TIME_LIMIT:
        log_res.write(f'\nMODEL REACHED TIME LIMIT OF {time_limit} SECONDS\n')
    else:
        log_res.write(f'\nTHERE WAS AN UNEXPECTED ERROR\n')


def get_biclique_cover_dataframe(g: nx.Graph, bicliques: Iterable, z: gp.Var, y: gp.Var) -> pd.DataFrame:
    """
    This function builds a Pandas dataframe containing the solution of a biclique cover. For this to work for
    a model, it is enough that said model define z variables representing bicliques and y variables representing
    whether an arc belongs to a biclique.

    :param g: input networkx graph.
    :param bicliques: iterable of the type range(ub), where ub is the upper bound of the number of bicliques,
                      as used in the implementing model.
    :param z: gurobipy variables representing whether a biclique is selected.
    :param y: gurobipy variables representing whether an arc is part of a selected biclique.
    :return: a binary dataframe, where each column is a biclique and each line is indexed by a node. Each cell can
             hold the values 0, 1, or 2. If the node does not belong to the biclique, the value is zero. If the node
             belongs to the first partition of the biclique, the value is one; likewise, if the node belongs to the
             second partition, the value of the cell is 2.
    """
    df_bicliques = DataFrame(
        {b: [0] * len(g.nodes) for b in bicliques if z[b].X == 1}, index=[u for u in g.nodes])
    for e in g.edges:
        u, v = e
        for b in bicliques:
            if z[b].X == 1:
                if y[u, v, b].X == 1:
                    df_bicliques.loc[u, b] = 1
                    df_bicliques.loc[v, b] = 2
                if y[v, u, b].X == 1:
                    df_bicliques.loc[v, b] = 1
                    df_bicliques.loc[u, b] = 2
    return df_bicliques


def is_biclique_cover(g: nx.Graph, df: pd.DataFrame) -> bool:
    """
    This function determines whether the solution to the minimum biclique cover problem found by a model is actually
    a valid biclique cover.

    :param g: original graph used by the model to find the solution.
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
        if any(not g.has_edge(u, v) for u in partition1 for v in partition2):
            return False
    return True


def create_and_save_model_comparison_report(report_name: str, **kwargs):
    # define constant strings
    model_v1 = 'model_v1'
    model_v2 = 'model_v2'
    # open log files
    dir_parent = path.abspath(path.join(getcwd(), pardir))
    dir_logs = path.join(dir_parent, 'logs')
    if not path.isdir(dir_logs):
        mkdir(dir_logs)
    ts = int(datetime.now().timestamp())
    dir_ts_logs = path.join(dir_logs, report_name + '-' + str(ts))
    mkdir(dir_ts_logs)
    log_gb_v1 = open(path.join(dir_ts_logs, model_v1 + '_gb' + '.log'), 'w+')
    log_gb_v2 = open(path.join(dir_ts_logs, model_v2 + '_gb' + '.log'), 'w+')
    log_res_v1 = open(path.join(dir_ts_logs, model_v1 + '_res' + '.log'), 'w+')
    log_res_v2 = open(path.join(dir_ts_logs, model_v2 + '_res' + '.log'), 'w+')

    def write_header(log: TextIO, g: nx.Graph, g_name: str):
        log.write('\n' + '-' * 80 + '\n')
        log.write(f'GRAPH: {g_name}\n')
        log.write(f'NODES: {len(g.nodes)}\n')
        log.write(f'EDGES: {len(g.edges)}\n\n')

    # create report
    report = GraphReport(name=report_name)
    report.add_properties([model_v1, model_v2])
    for g, g_name in get_graphs_in_store(**kwargs):
        report.add_graph_data(g, g_name)
        write_header(log=log_res_v1, g=g, g_name=g_name)
        report.add_property_values_from_function(
            p_name=model_v1, f=solve_bc_v1, g=g, log_gb=log_gb_v1, log_res=log_res_v1, time_limit=120)
        write_header(log=log_res_v2, g=g, g_name=g_name)
        report.add_property_values_from_function(
            p_name=model_v2, f=solve_bc_v2, g=g, log_gb=log_gb_v2, log_res=log_res_v2, time_limit=120)
    # save report
    report.save_csv()
    # close files
    log_gb_v1.close()
    log_gb_v2.close()
    log_res_v1.close()
    log_res_v2.close()


if __name__ == '__main__':
    create_and_save_model_comparison_report(
        report_name='complete-multipartite-comparison', fname_regex='complete.*partite')
