from argparse import ArgumentParser

from src.model import CGModel
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
            feasible, time = chronometer(f=model.solve)
            ub = model.solution
            lb = model.m.ObjBoundC if not feasible else None
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


def create_and_save_column_generation_report(config: ReportConfig):
    # constant strings
    str_model = 'model'
    str_columns_added = 'columns_added'
    str_obj_val = 'obj_val'
    str_master_time = 'master_time'
    str_pricing_time = 'pricing_time'
    str_total_time = 'total_time'
    # create the report
    report = GraphReport(name=config.resolved_report_name)
    report.add_properties(
        props=[str_model, str_columns_added, str_obj_val, str_master_time, str_pricing_time, str_total_time],
        add_time_property=False
    )
    try:
        for run_config in config.run_configs:
            print(f'Running config: {run_config}')
            g = get_graph_in_store(filename=run_config.graph)
            model = CGModel(g, run_config)
            _, time = chronometer(f=model.solve)
            # add values to report
            report.add_graph_data(g, run_config.resolved_gname)
            report.add_property_values(p_name=str_model, p_value=model.name())
            report.add_property_values(p_name=str_columns_added, p_value=model.columns_added)
            report.add_property_values(
                p_name=str_obj_val,
                p_value=model.solution if model.is_feasible() else model.best_known_solution)
            report.add_property_values(p_name=str_master_time, p_value=model.master_time)
            report.add_property_values(p_name=str_pricing_time, p_value=model.pricing_time)
            report.add_property_values(p_name=str_total_time, p_value=time)
    except Exception as e:
        print(e)
        report.save_csv(cleanup=True)
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
    if report_config.run_configs[0].model != "CGModel":
        create_and_save_model_comparison_report(config=report_config)
    else:
        create_and_save_column_generation_report(config=report_config)


if __name__ == '__main__':
    main()
