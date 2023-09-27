from os import path, listdir, getcwd, pardir, mkdir
from re import match
from time import time
from typing import Callable, TypeVar, Iterable

import networkx as nx
from networkx import Graph
from pandas import DataFrame


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
        if max_graphs and count >= max_graphs:
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


class GraphReport:
    """
    A class for automatically generating and saving reports about graphs based on user-defined properties and their
    calculation times.
    """

    def __init__(self, name: str = 'report'):
        self.name = name
        self._data = {
            'graph': [],
            'n_nodes': [],
            'n_edges': [],
        }
        self._props = {}
        self._finished_setup = False

    def _add_property(self, p_name: str):
        if p_name in self._props:
            return
        self._data[p_name] = []
        self._data[f'{p_name}_time'] = []
        self._props[p_name] = False
        self._rows = 0

    def add_property(self, p_name: str):
        """
        Add a graph property to the report. Automatically adds the corresponding calculation time property. Repeated
        properties are ignored.

        :param p_name: The name of the property to add.
        :raises RuntimeError: If property definition phase has already ended.
        """
        if self._finished_setup:
            raise RuntimeError('Property definition phase already ended.')
        self._add_property(p_name)

    def add_properties(self, props: list[str]):
        """
        Add multiple properties to the report. Automatically adds the corresponding calculation time property for
        each property added. Repeated properties are ignored.

        :param props: A list of property names to add.
        :raises RuntimeError: If property definition phase has already ended.
        """
        if self._finished_setup:
            raise RuntimeError('Property definition phase already ended.')
        for prop in props:
            self._add_property(prop)

    def _reset_props(self):
        for p_name in self._props.keys():
            self._props[p_name] = False

    def add_graph_data(self, g: Graph, g_name):
        """
        Add graph data to the report. The first graph data addition marks the end of the property definition phase,
        whence it is no longer possible to add properties. It also marks the start of a new row in the report. Graph
        data can only be added if it is the first time or if all the properties of the previous row have been filled.

        :param g: The graph to add to the report.
        :param g_name: The name of the graph.
        :raises RuntimeError: If not all properties have been filled.
        """
        if self._finished_setup and not all(self._props.values()):
            raise RuntimeError(f'Not all properties have been filled in row {self._rows}.')
        if self._finished_setup:
            self._rows += 1
            self._reset_props()
        self._finished_setup = True
        self._data['graph'].append(g_name)
        self._data['n_nodes'].append(len(g.nodes))
        self._data['n_edges'].append(len(g.edges))

    def add_property_values(self, p_name: str, p_value, p_time):
        """
        Add property values to the report.

        :param p_name: The name of the property.
        :param p_value: The value of the property.
        :param p_time: The time taken to calculate the property.
        :raises RuntimeError: If property definition phase has not ended yet or if the property is not defined.
        """
        if not self._finished_setup:
            raise RuntimeError('Property definition has not ended yet. Add graph data in order to finish it.')
        if p_name not in self._data.keys():
            raise ValueError(f'Property {p_name} not defined.')
        self._data[p_name].append(p_value)
        self._data[f'{p_name}_time'].append(p_time)
        self._props[p_name] = True

    def add_property_values_from_function(self, p_name: str, f: Callable[..., ReturnType], *args, **kwargs):
        """
        Add property values to the report using a function.

        :param p_name: The name of the property.
        :param f: The function to calculate the property.
        :param args: Non-keyword arguments for the function.
        :param kwargs: Keyword arguments for the function.
        """
        p_value, p_time = chronometer(f=f, *args, **kwargs)
        self.add_property_values(p_name=p_name, p_value=p_value, p_time=p_time)

    def get_report_data(self, **kwargs) -> DataFrame:
        """
        Get the report data as a DataFrame.

        :param kwargs: Additional keyword arguments for DataFrame creation.
        :return: The report data as a DataFrame.
        :raises RuntimeError: If no data to save or if not all properties have been filled.
        """
        if not self._finished_setup or not self._data:
            raise RuntimeError('No data to save.')
        if not all(self._props.values()):
            raise RuntimeError(f'Not all properties have been filled in row {self._rows}.')
        return DataFrame(data=self._data, **kwargs)

    def save_csv(self,
                 save_dir: str = None,
                 create_dir: bool = True,
                 report_name: str = None,
                 replace: bool = False, **kwargs):
        """
        Save the report data to a CSV file.

        :param save_dir: The directory where the CSV file will be saved.
        :param create_dir: Create the directory if it doesn't exist.
        :param report_name: The name of the CSV file.
        :param replace: Replace the file if it already exists.
        :param kwargs: Additional keyword arguments for DataFrame.to_csv().
        :raises FileNotFoundError: If the provided directory doesn't exist and cannot be created.
        """
        df_data = self.get_report_data()
        if save_dir is None:
            parent_dir = path.abspath(path.join(getcwd(), pardir))
            save_dir = path.join(parent_dir, 'reports')
        if not path.isdir(save_dir):
            if create_dir:
                mkdir(save_dir)
            else:
                raise FileNotFoundError('Directory does not exist')
        if report_name is None:
            report_name = self.name
        save_path = path.join(save_dir, report_name + '.csv')
        if not replace:
            count = 1
            while path.isfile(save_path):
                save_path = path.join(save_dir, report_name + f'-{count}.csv')
                count += 1
        df_data.to_csv(save_path, index=False, **kwargs)
