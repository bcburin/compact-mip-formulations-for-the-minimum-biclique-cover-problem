from argparse import ArgumentParser
from datetime import datetime
from os import path, mkdir

from src.util import GraphReport, get_graph_in_store, \
    get_file_name_and_extension, chronometer, get_and_create_logs_dir
from src.config import ReportConfig, read_report_config_file
import src.model as model_classes


def create_and_save_model_comparison_report(config: ReportConfig, **kwargs):
    # create logs directory
    dir_logs = get_and_create_logs_dir()
    dir_ts_logs = path.join(dir_logs, config.resolved_report_name)
    if config.create_log_files:
        mkdir(dir_ts_logs)

    # constant strings
    str_model = 'Model'
    str_k = 'k'
    str_time_k = 't_k'
    str_lb = 'LB'
    str_ub = 'UB'
    str_time = 'time'
    # create report
    report = GraphReport(name=config.resolved_report_name)
    report.add_properties([str_model, str_k, str_time_k, str_lb, str_ub, str_time], add_time_property=False)
    # runs from configuration
    try:
        for run_config in config.run_configs:
            g = get_graph_in_store(filename=run_config.graph)
            g_name, _ = get_file_name_and_extension(fname=run_config.graph)
            run_model = getattr(model_classes, run_config.model)
            model = run_model(
                g=g, g_name=run_config.resolved_gname, dir_logs=dir_ts_logs if config.create_log_files else None, config=run_config, default_config=config)
            # calculate values
            try:
                k, t_k = chronometer(model.upper_bound)
                ub, time = chronometer(f=model.solve)
                lb = model.m.ObjBoundC if not model.infeasible_or_unsolved() else None
            except Exception as e:
                print(f'Error running {run_config}: {e}')
                continue
            # add values to report
            report.add_graph_data(g, run_config.resolved_gname)
            report.add_property_values(p_name=str_model, p_value=run_config.model)
            report.add_property_values(p_name=str_k, p_value=k)
            report.add_property_values(p_name=str_time_k, p_value=t_k)
            report.add_property_values(p_name=str_lb, p_value=lb)
            report.add_property_values(p_name=str_ub, p_value=ub)
            report.add_property_values(p_name=str_time, p_value=time)
    except Exception as e:
        print(e)
        report.save_csv(cleanup=True)
    # save report
    report.save_csv()


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('config_file', help='path to the file containing run configurations')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    config_file = args.config_file
    report_config = read_report_config_file(config_file_path=config_file)
    create_and_save_model_comparison_report(config=report_config)


if __name__ == '__main__':
    main()
