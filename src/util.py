from os import path, listdir, getcwd, pardir
from re import match
from time import time
from typing import Callable, TypeVar, Iterable

import networkx as nx
from networkx import Graph


def get_file_extension(fname: str) -> str:
    return fname.split(sep='.')[-1]


def build_graph_from_file(fpath: str):
    if not path.isfile(fpath):
        raise ValueError('Invalid fpath parameter. Must be the path of the file.')
    fname = path.basename(fpath)
    ext = get_file_extension(fname)
    if ext == 'gml':
        g = nx.read_gml(fpath, label='id')
    elif ext == 'txt':
        with open(fpath, "r") as f:
            m = int(f.readline().strip().split()[-1])
            edges = []
            for _ in range(m):
                line = next(f).strip().split()
                u, v = map(int, line[:2])
                # Ensure that 'u' is smaller than 'v' to represent undirected edges
                edge = (min(u, v), max(u, v))
                edges.append(edge)
        g = nx.Graph()
        g.add_edges_from(edges)
    else:
        raise ValueError(f"File extension {ext} not supported, unsure of how to build graph.")
    return g


def get_graphs_in_store(
        max_nodes: int = None,
        max_edges: int = None,
        max_graphs: int = None,
        fname_regex: str = None,
        graph_dir: str = None,
        recursive: bool = False) -> Iterable[tuple[Graph, str]]:
    """
    This function reads the contents of the "graph" directory and returns an iterable of graphs built from it.

    :param fname_regex: filters graphs based on their file names. Only graphs whose names match the provided
                        RegEx pattern are included in the sequence.
    :param max_nodes: maximum amount of nodes to allow in a graph. Graphs with more nodes are excluded
                      from the sequence.
    :param max_edges: maximum amount of edges to allow in a graph. Graphs with more edges are excluded
                      from the sequence.
    :param max_graphs: maximum amount of graphs to include in the sequence. Once reached, the iteration is stopped,
                       i.e. the sequence ends.
    :param graph_dir: directory in which to search for the graphs. If non is provided, the standard graph directory
                      of the project is used.
    :param recursive: boolean flag. If activated, graphs are also searched in subfolders of the graph directory,
                      constrained by the same filters imposed by the other parameters. Otherwise, non-file paths
                      are ignored. By default, it is false.
    :return: iterable containing tuples with the networkx.Graph object and the name of the graph file.
    """

    if not graph_dir:
        parent_dir = path.abspath(path.join(getcwd(), pardir))
        graph_dir = path.join(parent_dir, 'graph')
    count = 0
    for filename in listdir(graph_dir):
        if count >= max_graphs:
            break
        if fname_regex and not match(fname_regex, filename):
            continue
        filepath = path.join(graph_dir, filename)
        if not path.isfile(filepath):
            if recursive and path.isdir(filepath):
                sub_seq = get_graphs_in_store(max_nodes=max_nodes,
                                              max_edges=max_edges,
                                              max_graphs=max_graphs-count,
                                              fname_regex=fname_regex,
                                              graph_dir=filepath,
                                              recursive=True)
                for sub_g, sub_name in sub_seq:
                    yield sub_g, sub_name
            else:
                continue
        g = build_graph_from_file(fpath=filepath)
        if max_nodes and len(g.nodes) > max_nodes:
            continue
        if max_edges and len(g.edges) > max_edges:
            continue
        count += 1
        yield g, filename


ReturnType = TypeVar('ReturnType')


def chronometer(f: Callable[..., ReturnType], *args, **kwargs) -> tuple[ReturnType, float]:
    """
    :param f: function whose execution time will be measured.
    :param args: non-keyword arguments of function f.
    :param kwargs: keyword arguments of function f.
    :return: a tuple of length two, containing the return of function f and the execution time.
    """
    start_time = time()
    val = f(*args, **kwargs)
    return val, time() - start_time

