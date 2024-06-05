import json

from pydantic import BaseModel, parse_obj_as, ValidationError

from src.bc_bounds import LBComputeMethod, UBComputeMethod


class RunConfig(BaseModel):
    graph: str
    model: str
    lb_method: LBComputeMethod = LBComputeMethod.INDEPENDENT_EDGES
    ub_method: UBComputeMethod = UBComputeMethod.VERTEX
    edge_fix: bool = False
    bottom_up: bool = False
    warm_start: bool = False
    conflict_inequalities: bool = False
    common_neighbor_inequalities: bool = False
    time_limit: int | None = None


class ReportConfig(BaseModel):
    report_name: str
    default_lb_method: LBComputeMethod = LBComputeMethod.INDEPENDENT_EDGES
    default_ub_method: UBComputeMethod = UBComputeMethod.VERTEX
    default_edge_fix: bool = False
    default_warm_start: bool = False
    default_bottom_up: bool = False
    default_conflict_inequalities: bool = False
    default_common_neighbor_inequalities: bool = False
    default_time_limit: int | None = 3600
    run_configs: list[RunConfig]


def read_report_config_file(config_file_path: str) -> ReportConfig:
    with open(config_file_path, 'r') as config_file:
        config_data = json.load(config_file)

    try:
        config = parse_obj_as(ReportConfig, config_data)
        return config
    except ValidationError as e:
        for error in e.errors():
            print(error['msg'])
