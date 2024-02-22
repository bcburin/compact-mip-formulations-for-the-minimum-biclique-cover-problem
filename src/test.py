from random import randint

import vertex_cover
import heuristic
import biclique
import edge_indep

from gerrychain import (Graph)

from src.util import build_multipartite_graph, save_graph_in_store, get_graphs_in_store, GraphReport, chronometer


def solve_bc(G, form, use_lower=False):
    print("Number of edges: ", len(G.edges))
    print("Number of nodes: ", len(G.nodes))
    # find a vertex cover of g which is an upper bound for the minimum biclique covering problem (MBCP)
    vertex_cover_number, vertex_cover_set = vertex_cover.solve(G)
    # find an independent set of g which is a lower bound for the minimum biclique covering problem (MBCP)
    # indep_number = independent_set.solve(g)
    # indep_number = directed_stars.solve(g)
    # print("iuc_number: ", indep_number)
    if use_lower:
        indep_edges = edge_indep.solve(G)
    else:
        indep_edges = []
    print("indep_number: ", len(indep_edges))
    # find a heuristic solution for the MBCP based on the vertex cover set
    heuristic_sol = heuristic.find_cover(G, vertex_cover_set)
    # print("heuristic solution: ", heuristic_sol)
    if form == 1:
        sol = biclique.get_vertex_bc_from_edge(G, heuristic_sol)
        print(sol)
        biclique_cover = biclique.solve_v(G, sol, indep_edges=indep_edges)
        print("Is it biclique cover? ", biclique.check_v_biclique_cover(G, biclique_cover))
    elif form == 2:
        biclique_cover = biclique.solve(G, heuristic_sol, indep_edges=indep_edges)
        print("Is it biclique cover? ", biclique.check_v_biclique_cover(G, biclique_cover))
    elif form == 3:
        _, biclique_cover = biclique.solve_recursive(G, heuristic_sol)
        print("Is it biclique cover? ", biclique.check_v_biclique_cover(G, biclique_cover))


def read_graph(level, state, level_name):
    return Graph.from_json("../data/" + level + "/json/" + state + "_" + level_name + ".json")


def build_and_save_multipartite_graphs(
        min_partitions: int,
        max_partitions: int,
        number_each_partition: int,
        min_vertices: int = 3,
        max_vertices: int = 10,
        edge_probability: float = 1.0):
    """
    This function generates and saves multipartite graphs with specified parameters.

    :param min_partitions: Smallest number of partitions in the generated graphs. (Inclusive)
    :param max_partitions: Largest number of partitions in the generated graphs. (Exclusive)
    :param number_each_partition: Number of graphs to create for each partition count.
    :param min_vertices: Minimum number of vertices in the partitions.
    :param max_vertices: Maximum number of vertices in the partitions.
    :param edge_probability: Probability of edge existence in the graphs (1.0 for complete graphs).

    The function generates multipartite graphs with a variable number of partitions, each with a random number of
    vertices within the specified range. It saves the generated graphs with descriptive names based on their properties.

    The naming format for the saved graphs includes the edge probability, the type of multipartite graph (bipartite,
    tripartite, or n-multipartite), and the sizes of the partitions.

    The generated graphs are saved using the `save_graph_in_store` function.

    No return value; the graphs are saved in the specified directory.
    """

    number_diff_partitions = max_partitions - min_partitions
    for i in range(number_diff_partitions * number_each_partition):
        # construct graph
        n = min_partitions + i % number_each_partition  # number of partitions in the graph
        sizes = [randint(min_vertices, max_vertices) for _ in range(n)]
        g = build_multipartite_graph(*sizes, edge_probability=edge_probability)
        # construct graph name
        str_partition_type = 'bipartite' if n == 2 else 'tripartite' if n == 3 else f'{n}-multipartite'
        str_edge_probability = 'complete' if edge_probability == 1 else f'{int(edge_probability * 100)}'
        g_name = f'{str_edge_probability}-{str_partition_type}' + ''.join(f'_{size}' for size in sizes)
        # save graph
        save_graph_in_store(g=g, g_name=g_name)


if __name__ == '__main__':
    # imports
    from model import ModelV1
    from util import get_logs_directory
    from biclique import solve_recursive, find_heuristic_solution
    # constants
    v1 = 'v1_recursive'
    v2 = 'v2'
    report_name = 'comparison'
    dir_ts_logs = get_logs_directory(name=report_name)
    time_limit = 100
    # create report
    report = GraphReport(name=report_name)
    report.add_properties([v1, v2])
    for g, g_name in get_graphs_in_store(fname_regex='from_indep_set.*', max_edges=290):
        print('-' * 80, '\n', g_name)
        # heuristic = find_heuristic_solution(g)
        report.add_graph_data(g, g_name)
        m2 = ModelV1(g=g, g_name=g_name, dir_logs=dir_ts_logs, time_limit=40)
        report.add_property_values_from_function(p_name=v2, f=m2.solve)
        report.add_property_values_from_function(p_name=v1, f=solve_recursive, G=g, apply=lambda x: x[0])
    # save report
    report.save_csv()

