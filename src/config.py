import pathlib
from datetime import datetime
from functools import cached_property

from pydantic import BaseModel

from src.bc_bounds import LBComputeMethod, UBComputeMethod
from src.util import get_file_name_and_extension


class RunConfig(BaseModel):
    name: str | None = None
    graph: str
    model: str
    lb_method: LBComputeMethod | None = LBComputeMethod.INDEPENDENT_EDGES
    ub_method: UBComputeMethod | None = UBComputeMethod.VERTEX
    edge_fix: bool | None = None
    bottom_up: bool | None = None
    warm_start: bool | None = None
    use_callback: bool | None = None
    conflict_inequalities: bool | None = None
    common_neighbor_inequalities: bool | None = None
    time_limit: int | None = None

    @cached_property
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
    default_use_callback: bool = False
    default_conflict_inequalities: bool = False
    default_common_neighbor_inequalities: bool = False
    default_time_limit: int | None = 3600
    run_configs: list[RunConfig]
    create_log_files: bool = False

    @cached_property
    def timestamp(self) -> int:
        return int(datetime.now().timestamp())

    @cached_property
    def resolved_report_name(self) -> str:
        if '$ts' in self.report_name:
            return self.report_name.replace('$ts', str(self.timestamp))
        else:
            return self.report_name


def read_report_config_file(config_file_path: str) -> ReportConfig:
    config_data = pathlib.Path(config_file_path).read_text()
    config = ReportConfig.model_validate_json(config_data)

    for run in config.run_configs:
        run.lb_method = run.lb_method or config.default_lb_method
        run.ub_method = run.ub_method or config.default_ub_method
        run.edge_fix = (
            run.edge_fix
            if run.edge_fix is not None
            else config.default_edge_fix
        )
        run.warm_start = (
            run.warm_start
            if run.warm_start is not None
            else config.default_warm_start
        )
        run.bottom_up = (
            run.bottom_up
            if run.bottom_up is not None
            else config.default_bottom_up
        )
        run.use_callback = (
            run.use_callback
            if run.use_callback is not None
            else config.default_use_callback
        )
        run.conflict_inequalities = (
            run.conflict_inequalities
            if run.conflict_inequalities is not None
            else config.default_conflict_inequalities
        )
        run.common_neighbor_inequalities = (
            run.common_neighbor_inequalities
            if run.common_neighbor_inequalities is not None
            else config.default_common_neighbor_inequalities
        )
        run.time_limit = (
            run.time_limit
            if run.time_limit is not None
            else config.default_time_limit
        )

    return config
