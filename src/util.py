from os import path
from time import time
from typing import Callable, TypeVar

import networkx as nx


def build_graph_from_file(fpath: str):
    if not path.isfile(fpath):
        raise ValueError('Invalid fpath parameter. Must be the path of the file.')
    fname = path.basename(fpath)
    ext = fname.split(sep='.')[-1]
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


R = TypeVar('R')


def chronometer(f: Callable[..., R]) -> Callable[..., tuple[R, float]]:

    def inner(*args, **kwargs) -> tuple[R, float]:
        start_time = time()
        val = f(*args, **kwargs)
        return val, time() - start_time

    return inner
