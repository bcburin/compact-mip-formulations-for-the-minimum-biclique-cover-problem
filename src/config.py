import json
from datetime import datetime
from functools import cached_property

from pydantic import BaseModel, parse_obj_as, ValidationError

from src.bc_bounds import LBComputeMethod, UBComputeMethod
from src.util import get_file_name_and_extension


class RunConfig(BaseModel):
    name: str | None = None
    graph: str
    model: str
    lb_method: LBComputeMethod = LBComputeMethod.INDEPENDENT_EDGES
    ub_method: UBComputeMethod = UBComputeMethod.VERTEX
    edge_fix: bool = False
    bottom_up: bool = False
    warm_start: bool = False
    conflict_inequalities: bool = False
    common_neighbor_inequalities: bool = False
    use_callback: bool = False
    time_limit: int | None = None

    @property
    def resolved_gname(self) -> str:
        if self.name is not None:
            return (
                self.name
                    .replace('$model', self.model)
                    .replace('$graph', get_file_name_and_extension(self.graph)[0]))
        return get_file_name_and_extension(self.graph)[0]


class ReportConfig(BaseModel):
    report_name: str
    default_lb_method: LBComputeMethod = LBComputeMethod.INDEPENDENT_EDGES
    default_ub_method: UBComputeMethod = UBComputeMethod.VERTEX
    default_edge_fix: bool = False
    default_warm_start: bool = False
    default_bottom_up: bool = False
    default_conflict_inequalities: bool = False
    default_common_neighbor_inequalities: bool = False
    default_use_callback: bool = False
    default_time_limit: int | None = 3600
    run_configs: list[RunConfig]

    @cached_property
    def timestamp(self) -> int:
        return int(datetime.now().timestamp())

    @property
    def resolved_report_name(self) -> str:
        if '$ts' not in self.report_name:
            raise ValueError("Report name must contain a timestamp declared by $ts")
        return self.report_name.replace('$ts', str(self.timestamp))


def read_report_config_file(config_file_path: str) -> ReportConfig:
    with open(config_file_path, 'r') as config_file:
        config_data = json.load(config_file)

    try:
        config = parse_obj_as(ReportConfig, config_data)
        return config
    except ValidationError as e:
        for error in e.errors():
            print(error['msg'])
