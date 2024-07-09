# Compact MIP Formulations for the Minimum Biclique Cover Problem

This is the GitHub repository for the paper "Compact MIP Formulations for the Minimum Biclique Cover Problem" by Bruno Burin, Hamidreza Validi, Bochuan Lyu, and Illya V. Hicks.


The minimum biclique cover problem aims to find a minimum number of bicliques (i.e., complete bipartite graph) such that every edge of the input graph is covered by at least one biclique. 

The following Figure shows an input graph (left) and a corresponding minimum biclique cover with size 4 (right). 

![Figure 1](readme_images/biclique_cover_figure.png?raw=true "An input graph and its corresponding minimum biclique cover")

## Installation

To run the project, start by creating and activating a virtual environment. Refer to the [official python documentation](https://docs.python.org/3/library/venv.html) for detailed instructions specific to your operating system. Once the virtual environment is activated, install all dependencies with the following command:

```commandline
pip install -r requirements.txt
```

## Running the Project

Run the project using the command:

```commandline
python -m src.main "configs/config_file.json"
```

Here, `config_file.json` is a configuration file specifying the parameters for batch runs of different graphs and models. The `./config/` directory contains several pre-made configuration files. Detailed documentation on creating custom configuration files is provided below.

To view the CLI interface specifications, use the `-h` option:

```commandline
python -m src.main -h
```

## Configuration Files

Configuration files, written in JSON format, define the parameters and settings for batch runs. Below is a sample configuration file and an explanation of each field.

### Example Configuration File

```json
{
  "report_name":  "compare-edge-fix-natural-model-$ts",
  "default_lb_method": "independent_edges",
  "default_ub_method": "vertex",
  "default_time_limit": 3600,
  "default_warm_start": false,
  "run_configs": [
    {
      "graph": "ieee30.txt",
      "model": "NaturalModel",
      "name": "$graph - without edge fix",
      "edge_fix": false
    },
    {
      "graph": "ieee30.txt",
      "model": "NaturalModel",
      "name": "$graph - with edge fix",
      "edge_fix": true
    }
    // More run configurations...
  ]
}
```

### Fields

Required fields are marked with an asterisk `*`. Others are optional, but can provide additional customization to each batch or instance run.

- **report_name\*:** this field defines the name of the report generated by the project. The required `$ts` placeholder will be replaced by a timestamp when the report is generated.
- **default_lb_method:** specifies the default method for computing the lower bound. Default is `independent_edges`. The options are:
  - `match`
  - `lovasz`
  - `clique`
  - `independent_edges`
  - `maximal_independent_set`
- **default_ub_method:** specifies the default method for computing the upper bound. Default is `vertex`. The options are:
  - `number`
  - `vertex`
  - `merge_stars`
- **default_time_limit:** the maximum time allowed for the execution of each run in seconds. Default is 3600 seconds (1 hour).
- **default_edge_fix:** a boolean flag indicating whether to use edge fixing constraint in the models. Default is `false`.
- **default_warm_start:** a boolean flag indicating whether to warm start a model with a provided heuristic. Default is `false`.
- **default_conflict_inequalities:** a boolean flag indicating whether to add conflct inequaties as constraints to the extended models. Default is `false`.
- **default_common_neighbor_inequalities** a boolean flag indicating whether to add common neighbor inequaties as constraints to the extended models. Default is `false`.
- **default_use_callback:** a boolean flag indicating whether to use a callback function to add fractional cuts while solving the model. Default is `false`.
- **run_configs\*:** a list of instance run configuration objects. Each object specifies the parameters for an individual instance run. Their configuration is detailed below.

### Instance Run Configurations

- **graph\*:** the filename of the graph to be used in the run.
- **model\*:** the model to be used in the run. The options are:
  - `NaturalModel`
  - `ExtendedModel`
- **name:** a custom name for the run to be shown in the final report. The placeholders `$graph` and `$model` are replaced by the name of graph (without file extension) and the model used in the run, generating dynamic values. Default is `$graph`.
- **lb_method:** overrides `default_lb_method` for this run if provided.
- **ub_method:** overrides `default_ub_method` for this run if provided.
- **time_limit:** overrides `default_time_limit` for this run if provided.
- **edge_fix:** overrides `default_edge_fix` for this run if provided.
- **warm_start:** overrides `default_warm_start` for this run if provided.
- **conflict_inequalities:** overrides `default_conflict_inequalities` for this run if provided.
- **common_neighbor_inequalities:** overrides `default_common_neighbor_inequalities` for this run if provided.
- **use_callback:** overrides `default_use_callback` for this run if provided.
