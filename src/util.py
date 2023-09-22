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


def get_graphs_in_store(max_nodes: int = None, max_edges: int = None, fname_regex: str = None) -> Iterable[tuple[Graph, str]]:
    """
    This function reads the contents of the "graph" directory and returns an iterable of graphs built from it.

    :param fname_regex: filters graphs based on their file names. Only graphs whose names match the provided
                        RegEx pattern are included in the sequence.
    :param max_nodes: maximum amount of nodes to allow in a graph. Graphs with more nodes are excluded
                      from the sequence.
    :param max_edges: maximum amount of edges to allow in a graph. Graphs with more edges are excluded
                      from the sequence.
    :return: iterable containing tuples with the networkx.Graph object and the name of the graph file.
    """

    parent_dir = path.abspath(path.join(getcwd(), pardir))
    graph_dir = path.join(parent_dir, 'graph')
    for filename in listdir(graph_dir):
        if fname_regex and not match(fname_regex, filename):
            continue
        filepath = path.join(graph_dir, filename)
        if not path.isfile(filepath):
            continue
        g = build_graph_from_file(fpath=filepath)
        if max_nodes and len(g.nodes) > max_nodes:
            continue
        if max_edges and len(g.edges) > max_edges:
            continue
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

