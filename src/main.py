from argparse import ArgumentParser

from src.util import GraphReport, get_graph_in_store, chronometer
from src.config import ReportConfig, read_report_config_file
import src.model as model_classes


def create_and_save_model_comparison_report(config: ReportConfig):
    # constant strings
    str_model = 'Model'
    str_k = 'k'
    str_time_k = 't_k'
    str_lb = 'LB'
    str_ub = 'UB'
    str_time = 'time'
    # create the report
    report = GraphReport(name=config.resolved_report_name)
    report.add_properties([str_model, str_k, str_time_k, str_lb, str_ub, str_time], add_time_property=False)
    # runs from configuration
    try:
        for run_config in config.run_configs:
            print(f'Running config: {run_config}')
            g = get_graph_in_store(filename=run_config.graph)
            run_model = getattr(model_classes, run_config.model)
            model = run_model(g=g, config=run_config)
            # calculate values
            k, t_k = chronometer(lambda: model.upper_bound)
            ub, time = chronometer(f=model.solve)
            lb = model.m.ObjBoundC if not model.is_feasible else None
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
    # save the report
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
