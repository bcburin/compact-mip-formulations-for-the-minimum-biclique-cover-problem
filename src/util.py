import json
import math
from datetime import datetime
from os import path, listdir, getcwd, pardir, mkdir, remove
from itertools import combinations, product
from pprint import pprint
from random import random, randint
from re import match
from time import time
from typing import Callable, TypeVar, Iterable

import networkx as nx
from networkx import Graph
from pandas import DataFrame
from pydantic import BaseModel, ValidationError, parse_obj_as

ReturnType = TypeVar('ReturnType')


def get_file_name_and_extension(fname: str) -> tuple[str, str | None]:
    """
    This function separates a file name from the extension. For this, the first "." the file name is
    considered for separating the two parts. If there is no ".", `None` is returned as the extension.

    :param fname: name of the file
    :return: a tuple containing the file name without the extension and the file extension itself.
    """
    if '.' not in fname:
        return fname, None
    name, ext = fname.split(sep='.', maxsplit=1)
    return name, ext


def build_graph_from_file(fpath: str) -> nx.Graph:
    """
    This function builds a networkx.Graph object from a file.

    Supported file types: GML, TXT.

    :param fpath: path of the input file.
    :return: networkx.Graph object built from data in input graph.
    :raises ValueError: if the provided path is not a file, contains an unsupported file type, or doesn't
                        contain any file extension at all.
    """

    if not path.isfile(fpath):
        raise ValueError(f'Invalid fpath parameter <{fpath}>. Must be the path of the file.')
    fname = path.basename(fpath)
    _, ext = get_file_name_and_extension(fname)
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
    elif ext is None:
        raise ValueError(f"No file extension provided, unsure of how to build graph.")
    else:
        raise ValueError(f"File extension {ext} not supported, unsure of how to build graph.")
    return g


def get_graphs_in_store(
        max_nodes: int = None,
        max_edges: int = None,
        max_graphs: int = None,
        fname_regex: str = None,
        graph_dir: str = None,
        recursive: bool = False,
        include_extension: bool = False) -> Iterable[tuple[Graph, str]]:
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
    :param include_extension: include file extension in graph name
    :return: iterable containing tuples with the networkx.Graph object and the name of the graph file.
    """

    if not graph_dir:
        current_dir = path.dirname(path.abspath(__file__))
        parent_dir = path.abspath(path.join(current_dir, pardir))
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
                                              max_graphs=max_graphs-count if max_graphs else None,
                                              fname_regex=fname_regex,
                                              graph_dir=filepath,
                                              include_extension=include_extension,
                                              recursive=True)
                for sub_g, sub_name in sub_seq:
                    yield sub_g, sub_name
                continue
            else:
                continue
        g = build_graph_from_file(fpath=filepath)
        if max_nodes and len(g.nodes) > max_nodes:
            continue
        if max_edges and len(g.edges) > max_edges:
            continue
        count += 1
        name = filename if include_extension else get_file_name_and_extension(fname=filename)[0]
        yield g, name


def get_graph_in_store(filename: str, graph_dir: str = None) -> Graph:
    if not graph_dir:
        current_dir = path.dirname(path.abspath(__file__))
        parent_dir = path.abspath(path.join(current_dir, pardir))
        graph_dir = path.join(parent_dir, 'graph')
    filepath = path.join(graph_dir, filename)
    g = build_graph_from_file(fpath=filepath)
    return g


def save_graph_in_store(
        g: Graph,
        g_name: str,
        graph_dir: str = None,
        create_dir: bool = False,
        append_n_nodes: bool = True,
        append_n_edges: bool = True,
        replace: bool = False
):
    """
    Save a NetworkX graph to a GML file in a specified directory, with options to modify the file name.

    :param g: The NetworkX graph to be saved.
    :param g_name: The base name for the saved graph file.
    :param graph_dir: The directory where the graph file will be saved. If None, the default "graph" directory of
                      the project is used.
    :param create_dir: If True and the specified directory doesn't exist, it will be created. Defaults to False.
    :param append_n_nodes: If True, appends the number of nodes to the file name. Defaults to True.
    :param append_n_edges: If True, appends the number of edges to the file name. Defaults to True.
    :param replace: If True and a file with the same name already exists, it will be replaced. If False and the
                    file exists, a FileExistsError is raised. Defaults to False.
    :raises FileExistsError: If a file with the same name already exists and 'replace' is False.
    """

    if not graph_dir:
        current_dir = path.dirname(path.abspath(__file__))
        parent_dir = path.abspath(path.join(current_dir, pardir))
        graph_dir = path.join(parent_dir, 'graph')
    if not path.isdir(graph_dir) and create_dir:
        mkdir(graph_dir)
    if append_n_nodes:
        g_name += f'_{len(g.nodes)}'
    if append_n_edges:
        g_name += f'_{len(g.edges)}'
    fpath = path.join(graph_dir, f'{g_name}.gml')
    if path.isfile(fpath):
        if replace:
            remove(fpath)
        else:
            raise FileExistsError()
    nx.write_gml(g, path=fpath)


def build_multipartite_graph(*partition_sizes: int, edge_probability: float = 1.0) -> Graph:
    """
    Build a multipartite graph from the given partition sizes.

    This function constructs a complete multipartite graph where each partition contains nodes that are not connected,
    and nodes from different partitions are at most fully connected.

    :param partition_sizes: Sizes of the partitions in the multipartite graph.
    :param edge_probability: The probability of an edge being inserted between nodes from different partitions
                             (default is 1.0, meaning always connect).
    :return: A NetworkX graph representing the complete multipartite graph.

    Example Usage:
    ```
    # Build a complete multipartite graph with partitions of sizes 3, 2, and 4
    complete_multipartite_graph = build_complete_multipartite_graph(3, 2, 4)
    ```

    """

    if edge_probability < 0 or edge_probability > 1:
        raise ValueError('Probability must be between zero and one.')
    # Create an empty graph
    g = Graph()
    # Create nodes for each partition and add them to the graph
    partitions = []
    node_counter = 0
    for size in partition_sizes:
        partition = range(node_counter, node_counter + size)
        partitions.append(partition)
        g.add_nodes_from(partition)
        node_counter += size
    # Connect nodes from different partitions to create a complete multipartite graph
    for partition1, partition2 in combinations(partitions, 2):
        for node1, node2 in product(partition1, partition2):
            if random() <= edge_probability:  # Check if an edge should be added
                g.add_edge(node1, node2)
    return g


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


def poisson(x: float, x_0: float, mu: float):
    if mu < 0:
        raise ValueError('The parameter mu must be a positive real number.')
    x = x - x_0
    ln_res = -mu + x * math.log(mu) - math.log(math.factorial(math.floor(x)))
    return math.exp(ln_res)


def get_logs_directory(name: str):
    dir_parent = path.abspath(path.join(getcwd(), pardir))
    dir_logs = path.join(dir_parent, 'logs')
    if not path.isdir(dir_logs):
        mkdir(dir_logs)
    ts = int(datetime.now().timestamp())
    dir_ts_logs = path.join(dir_logs, name + '-' + str(ts))
    mkdir(dir_ts_logs)
    return dir_ts_logs


class RunConfig(BaseModel):
    graph: str
    model: str
    lb_method: str = "INDEPENDENT_EDGES"
    ub_method: str = "VERTEX"


def read_run_config_file(config_file_path: str) -> list[RunConfig]:
    with open(config_file_path, 'r') as config_file:
        run_config_data = json.load(config_file)

    try:
        run_configs = parse_obj_as(list[RunConfig], run_config_data)
        return run_configs
    except ValidationError as e:
        for error in e.errors():
            print(error['msg'])


def is_biclique(graph: Graph, edges: list):
    # Create a subgraph induced by the given set of edges
    subgraph = graph.edge_subgraph(edges)

    # Check if the induced subgraph is bipartite
    if not nx.is_bipartite(subgraph):
        return False

    # Check if every node in one set is connected to every node in the other set
    set1, set2 = nx.bipartite.sets(subgraph)
    for u in set1:
        for v in set2:
            if (u, v) not in subgraph.edges:
                return False
    return True


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

    def _add_property(self, p_name: str, add_time_property: bool):
        if p_name in self._props:
            return
        self._data[p_name] = []
        if add_time_property:
            self._data[f'{p_name}_time'] = []
        self._props[p_name] = False
        self._rows = 0

    def add_property(self, p_name: str, add_time_property: bool = True):
        """
        Add a graph property to the report. Repeated properties are ignored.

        :param p_name: The name of the property to add.
        :param add_time_property: Automatically adds the corresponding calculation time property.
        :raises RuntimeError: If property definition phase has already ended.
        """
        if self._finished_setup:
            raise RuntimeError('Property definition phase already ended.')
        self._add_property(p_name, add_time_property)

    def add_properties(self, props: list[str], add_time_property: bool = True):
        """
        Add multiple properties to the report. Repeated properties are ignored.
        :param add_time_property: Automatically adds the corresponding calculation time property for
        each property added.
        :param props: A list of property names to add.
        :raises RuntimeError: If property definition phase has already ended.
        """
        if self._finished_setup:
            raise RuntimeError('Property definition phase already ended.')
        for prop in props:
            self._add_property(prop, add_time_property)

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

    def add_property_values(self, p_name: str, p_value, p_time=None):
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
        if p_time is not None:
            self._data[f'{p_name}_time'].append(p_time)
        if p_time is None and f'{p_name}_time' in self._props:
            raise ValueError(f'Property {p_name} requires corresponding time {p_name}_time.')
        self._props[p_name] = True

    def add_property_values_from_function(
            self, p_name: str, f: Callable[..., ReturnType],  apply: Callable = None, *args, **kwargs):
        """
        Add property values to the report using a function.

        :param apply: apply to result of callable function before adding to report.
        :param p_name: The name of the property.
        :param f: The function to calculate the property.
        :param args: Non-keyword arguments for the function.
        :param kwargs: Keyword arguments for the function.
        """
        p_value, p_time = chronometer(f=f, *args, **kwargs)
        if apply:
            p_value = apply(p_value)
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
            raise RuntimeError(f'Not all properties have been filled in row {self._rows}. \n _data={pprint(self._data)}')
        return DataFrame(data=self._data, **kwargs)

    def _cleanup_rows(self):
        min_length = min(len(lst) for lst in self._data.values())
        for prop, lst in self._data.items():
            if len(lst) > min_length:
                self._data[prop] = lst[:min_length]
        for prop in self._props.keys():
            self._props[prop] = True

    def save_csv(self,
                 save_dir: str = None,
                 create_dir: bool = True,
                 report_name: str = None,
                 replace: bool = False,
                 cleanup: bool = False, **kwargs):
        """
        Save the report data to a CSV file.

        :param cleanup: check for any unfilled properties in last row. In case there is any, remove last row.
        :param save_dir: The directory where the CSV file will be saved.
        :param create_dir: Create the directory if it doesn't exist.
        :param report_name: The name of the CSV file.
        :param replace: Replace the file if it already exists.
        :param kwargs: Additional keyword arguments for DataFrame.to_csv().
        :raises FileNotFoundError: If the provided directory doesn't exist and cannot be created.
        """
        if cleanup:
            self._cleanup_rows()
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


def var_swap(x, u, v, i):
    """
    This function gets the value of x[e, i] for e=[u, v] when the order of the vertices in the graph
    definition is unknown.
    """
    # the nested gets are for when x[u, v, i] doesn't exist, but x[v, u, i] does
    return x.get((u, v, i), x.get((v, u, i)))
