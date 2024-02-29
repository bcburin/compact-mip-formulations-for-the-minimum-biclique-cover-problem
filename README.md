# Minimum Biclique Cover

---

## Installation

The first step to run the project is to create and activate a virtual environment. The [official python documentation](https://docs.python.org/3/library/venv.html) details how it can be done for each operational system. Once the virtual environment is activated, one must simply run the following CLI command in order to install all the dependencies.

```commandline
pip install -r requirements.txt
```

## Running the Project

The project can be run using the following CLI command

```commandline
python -m src.main "configs/config_table1.json"
```

where the `run_configs.json` is a file containing the specifications for a batch run of different graphs and models. The optional parameter `--time-limit` allows defining the time limit in seconds for each run. The default value for the time limit is 3600 (1 hour). For example, in order to limit it to 5 minutes, one can use

```commandline
python -m src.main --time-limit 300 "configs/config_table1.json"
```

Using the option `-h` shows the specifications of the CLI interface of the program.

```commandline
python -m src.main -h
```

## Configuration Files

Run configuration files are json files that contain lists with the specification of the graph and model to be used in each run.

```json
[
  {
    "graph": "ieee30.txt",
    "model": "NaturalModel"
  },
  {
    "graph": "ieee30.txt",
    "model": "ExtendedModel"
  },
  {
    "graph": "karate34.gml",
    "model": "NaturalModel"
  }
]
```