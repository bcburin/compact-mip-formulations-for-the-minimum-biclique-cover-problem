from argparse import ArgumentParser
from datetime import datetime
from os import path, getcwd, pardir, mkdir

from src.util import GraphReport, read_run_config_file, RunConfig, get_graph_in_store, \
    get_file_name_and_extension, chronometer
import src.model as model_classes


def create_and_save_model_comparison_report(
        report_name: str, run_configs: list[RunConfig], time_limit: int = None,
        suppress_ts_in_report_name: bool = True, **kwargs):
    # get log directory
    dir_parent = path.abspath(path.join(getcwd(), pardir))
    dir_logs = path.join(dir_parent, 'logs')
    if not path.isdir(dir_logs):
        mkdir(dir_logs)
    ts = int(datetime.now().timestamp())
    dir_ts_logs = path.join(dir_logs, report_name + '-' + str(ts))
    mkdir(dir_ts_logs)

    # constant strings
    str_model = 'Model'
    str_k = 'k'
    str_time_k = 't_k'
    str_lb = 'LB'
    str_ub = 'UB'
    str_time = 'time'
    # create report
    report_name = report_name if suppress_ts_in_report_name else report_name + '-' + str(ts)
    report = GraphReport(name=report_name)
    report.add_properties([str_model, str_k, str_time_k, str_lb, str_ub, str_time])
    # runs from configuration
    for run_config in run_configs:
        g = get_graph_in_store(filename=run_config.graph)
        g_name, _ = get_file_name_and_extension(fname=run_config.graph)
        report.add_graph_data(g, g_name)
        run_model = getattr(model_classes, run_config.model)
        model = run_model(g=g, g_name=g_name, dir_logs=dir_ts_logs, time_limit=time_limit, **kwargs)
        # calculate values
        k, t_k = chronometer(model.upper_bound)
        ub, time = chronometer(f=model.solve)
        # add values to report
        report.add_property_values(p_name=str_model, p_value=run_config.model)
        report.add_property_values(p_name=str_k, p_value=k)
        report.add_property_values(p_name=str_time_k, p_value=t_k)
        report.add_property_values(p_name=str_lb, p_value=model.m.ObjBoundC)
        report.add_property_values(p_name=str_ub, p_value=ub)
        report.add_property_values(p_name=str_time, p_value=time)
    # save report
    report.save_csv()


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('config_file', help='path to the file containing run configurations')
    parser.add_argument('--time-limit',
                        help='default time limit in seconds to use in run configurations',
                        default=60*60)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    config_file = args.config_file
    time_limit = args.time_limit
    run_configs = read_run_config_file(config_file_path=config_file)
    create_and_save_model_comparison_report(
        report_name='test_report_simple', run_configs=run_configs, time_limit=int(time_limit))


if __name__ == '__main__':
    main()
