from itertools import combinations

import gurobipy as gp
import networkx as nx
from gurobipy import GRB

from src.base_model import BaseMinimumBicliqueCoverGpSolver, BaseMinimumBicliqueCoverSolver
from src.config import RunConfig
from src.exceptions import NotYetSolvedError
from src.independent_set import solve_max_weighted_independent_set
from src.util import is_biclique, var_swap, chronometer


def indep_callback(model: gp.Model, where: int):
    """
    :param model: Gurobi model, expects the following variables to be added:
        - _g containing the graph of the model;
        - _g2 containing the power graph g^2; and
        - _k storing the calculated maximum number of bicliques.
        - _Y variable of the extended model
        - _z variable of the extended model
    :param where: Gurobi internal variable that represent place where the callback function was called.
    """
    if where != GRB.Callback.MIPNODE or model.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
        return
    y_val = model.cbGetNodeRel(model._y)
    z_val = model.cbGetNodeRel(model._z)
    for j in range(model._k):
        for v in model._g.nodes:
            model._g.nodes[v]["weight"] = y_val[v, j, 0] + y_val[v, j, 1]
        obj, indep_set = solve_max_weighted_independent_set(model._g, model._g2)
        if obj > z_val[j]:
            model.cbCut(gp.quicksum(model._y[v, j, 0] + model._y[v, j, 1] for v in indep_set) <= model._z[j])


class NaturalModel(BaseMinimumBicliqueCoverGpSolver):

    def __init__(self, g: nx.Graph, config: RunConfig):
        super().__init__(g, config)

    def _add_variables(self):
        # 4j
        self.z = self.m.addVars(self.bicliques, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name="z")
        self.x = self.m.addVars(self.graph.edges, self.bicliques, vtype=GRB.BINARY, name="x")

    def _set_objective(self):
        # 4a
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    def _add_constraints(self):
        m, x, z = self.m, self.x, self.z

        # 4b
        m.addConstrs(x[u, v, i] <= z[i] for u, v in self.graph.edges for i in self.bicliques)
        # 4c
        m.addConstrs(gp.quicksum(x[u, v, i] for i in self.bicliques) >= 1 for u, v in self.graph.edges)
        # 4e
        for cycle in nx.simple_cycles(self.graph, length_bound=3):  # TODO: fix me
            if len(cycle) != 3:
                continue
            cycle_edges = [[cycle[0], cycle[1]], [cycle[1], cycle[2]], [cycle[2], cycle[0]]]
            m.addConstrs(
                gp.quicksum(var_swap(x, u, v, i) for u, v in cycle_edges) <= 2 * z[i]
                for i in self.bicliques)

        for e, f, cr1, cr2 in self.get_disjoint_edges():
            a, b = e
            c, d = f
            # 4d
            if cr1 < 2 and cr2 < 2:
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] for i in self.bicliques)
            # 4f
            if cr1 == 2 and cr2 < 2:
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, a, d, i) for i in self.bicliques)
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, b, c, i) for i in self.bicliques)
            # 4g
            if cr1 < 2 and cr2 == 2:
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, a, c, i) for i in self.bicliques)
                m.addConstrs(x[a, b, i] + x[c, d, i] <= z[i] + var_swap(x, b, d, i) for i in self.bicliques)
            # 4h
            if cr1 == 2 and cr2 == 2:
                m.addConstrs(
                    x[a, b, i] + x[c, d, i] <=
                    z[i] + var_swap(x, a, c, i) + var_swap(x, a, d, i) for i in self.bicliques)
                m.addConstrs(
                    x[a, b, i] + x[c, d, i] <=
                    z[i] + var_swap(x, b, c, i) + var_swap(x, b, d, i) for i in self.bicliques)
        # 4i
        m.addConstrs(z[i] >= z[i + 1] for i in range(self.upper_bound - 1))
        # independent edges constraints
        if self._can_add_indep_edges_constraints():
            self._add_independent_edges_constraints()
        # strengthening constraints
        if self._config.use_strengthening_constraints:
            self._add_strengthening_constraints()

    def _check_biclique_cover(self) -> bool:
        # check it's a cover
        if not any(self.x[u, v, i].X == 1 for u, v in self.graph.edges for i in self.bicliques):
            return False
        # check it's a biclique cover
        for i in self.bicliques:
            if self.z[i].X != 1:
                continue
            edges = [(u, v) for u, v in self.graph.edges if self.x[u, v, i].X == 1]
            if not is_biclique(graph=self.graph, edges=edges):
                return False
        return True

    def _do_warm_start(self, indep_edges: list, vertex_cover: list):
        assign = dict()
        for i, s in enumerate(vertex_cover):
            edges = []
            for e in self.graph.edges(s):
                self.x[min(e), max(e), i].start = 1
                edges.append(e)
            assign[i] = edges
        indep_edges = set(indep_edges)
        if self._config.edge_fix:
            for i in assign.keys():
                for e in assign[i]:
                    if e in indep_edges:
                        self.x[min(e), max(e), i].lb = 1

    def _pre_solve(self):
        for b in range(self.lower_bound):
            self.z[b].lb = 1
        if self._config.warm_start:
            self._guarantee_compute_lb_and_indep_edges()
            self._guarantee_compute_ub_and_vertex_cover()
            self._do_warm_start(indep_edges=self._indep_edges, vertex_cover=self._vertex_cover)

    def _add_independent_edges_constraints(self):
        self._guarantee_compute_lb_and_indep_edges()
        edges = self._indep_edges
        for biclique in self.bicliques:
            if biclique >= len(edges):
                return
            e = edges[biclique]
            self.x[e[0], e[1], biclique].lb = 1

    def _add_strengthening_constraints(self):
        m, x, z = self.m, self.x, self.z
        for e, f, cr1, cr2 in self.get_disjoint_edges():
            a, b = e
            c, d = f
            if cr1 == 2:
                m.addConstrs(
                    var_swap(x, b, c, i) + var_swap(x, a, b, i) + var_swap(x, c, d, i) <=
                    2 * z[i] + var_swap(x, a, d, i) for i in self.bicliques)
                m.addConstrs(
                    var_swap(x, a, d, i) + var_swap(x, a, b, i) + var_swap(x, c, d, i) <=
                    2 * z[i] + var_swap(x, b, c, i) for i in self.bicliques)

    @classmethod
    def name(cls) -> str:
        return 'Compact Natural Model'


class ExtendedModel(BaseMinimumBicliqueCoverGpSolver):
    
    def __init__(self, g: nx.Graph, config: RunConfig):
        super().__init__(g, config)

    def _add_variables(self):
        # 5h
        self.z = self.m.addVars(self.bicliques, vtype=GRB.BINARY, name="z")
        self.x = self.m.addVars(self.directed.edges, self.bicliques, vtype=GRB.BINARY, name="x")
        self.y = self.m.addVars(self.graph.nodes, self.bicliques, range(2), vtype=GRB.BINARY, name="y")

    def _set_objective(self):
        # 5a
        self.m.setObjective(gp.quicksum(self.z), sense=GRB.MINIMIZE)

    def _add_constraints(self):
        m, x, y, z = self.m, self.x, self.y, self.z

        # 5b
        m.addConstrs(x[u, v, i] <= y[u, i, 0] for u, v in self.directed.edges for i in self.bicliques)
        m.addConstrs(x[u, v, i] <= y[v, i, 1] for u, v in self.directed.edges for i in self.bicliques)
        # 5c
        m.addConstrs(y[u, i, 0] + y[v, i, 1] <= z[i] + x[u, v, i]
                     for u, v in self.directed.edges for i in self.bicliques)
        # 5d
        m.addConstrs(y[u, i, 0] + y[u, i, 1] <= z[i] for u in self.graph.nodes for i in self.bicliques)
        # 5e
        for u, v in self.graph.edges:
            m.addConstr(gp.quicksum(x[u, v, i] + x[v, u, i] for i in self.bicliques) >= 1)
        # 5f
        for u, v in combinations(self.graph.nodes, r=2):
            if self.graph.has_edge(u, v):
                continue
            m.addConstrs(y[u, i, 0] + y[v, i, 1] <= z[i] for i in self.bicliques)
            m.addConstrs(y[v, i, 0] + y[u, i, 1] <= z[i] for i in self.bicliques)
        # 5g
        m.addConstrs(z[i] >= z[i + 1] for i in range(self.upper_bound - 1))
        # independent edges constraints
        if self._can_add_indep_edges_constraints():
            self._add_independent_edges_constraints()
        # conflict inequalities
        if self._config.conflict_inequalities:
            self._add_conflict_inequalities()
        # common neighbors inequalities
        if self._config.common_neighbor_inequalities:
            self._add_common_neighbor_inequalities()
        # strengthening constraints
        if self._config.use_strengthening_constraints:
            self._add_strengthening_constraints()

    def _check_biclique_cover(self) -> bool:
        # check it's a cover
        if not any(self.x[u, v, i].X == 1 or self.x[v, u, i].X == 1
                   for u, v in self.graph.edges for i in self.bicliques):
            return False
        # check it's a biclique cover
        for i in self.bicliques:
            if self.z[i].X != 1:
                continue
            edges = [(u, v) for u, v in self.graph.edges if self.x[u, v, i].X == 1 or self.x[v, u, i].X == 1]
            if not is_biclique(graph=self.graph, edges=edges):
                return False
        return True

    def _pre_solve(self):
        for b in range(self.lower_bound):
            self.z[b].lb = 1
        if self._config.use_callback:
            self._add_callback()
        if self._config.warm_start:
            self._guarantee_compute_lb_and_indep_edges()
            self._guarantee_compute_ub_and_vertex_cover()
            self._do_warm_start(indep_edges=self._indep_edges, vertex_cover=self._vertex_cover)

    def _add_independent_edges_constraints(self):
        self._guarantee_compute_lb_and_indep_edges()
        edges = self._indep_edges
        for biclique in self.bicliques:
            if biclique >= len(edges):
                return
            a, b = edges[biclique]
            self.x[a, b, biclique].lb = 1
            self.x[b, a, biclique].lb = 0

    def _add_conflict_inequalities(self):
        for u, v in combinations(self.graph.nodes, r=2):
            if self.power_graph.has_edge(u, v):
                continue
            conflict_inequalities = self.m.addConstrs(
                self.y[u, i, 0] + self.y[u, i, 1] + self.y[v, i, 0] + self.y[v, i, 1] <= self.z[i]
                for i in self.bicliques)
            for i in self.bicliques:
                conflict_inequalities[i].Lazy = 1

    def _add_common_neighbor_inequalities(self):
        for u, v in self.power_graph.edges:
            if self.graph.has_edge(u, v):
                continue
            common_neighbors = nx.common_neighbors(self.graph, u, v)
            self.m.addConstrs(
                self.y[u, i, 0] + self.y[v, i, 0] <= self.z[i] + gp.quicksum(self.y[c, i, 1] for c in common_neighbors)
                for i in self.bicliques)
            self.m.addConstrs(
                self.y[u, i, 1] + self.y[v, i, 1] <= self.z[i] + gp.quicksum(self.y[c, i, 0] for c in common_neighbors)
                for i in self.bicliques)

    def _add_strengthening_constraints(self):
        ...

    def _do_warm_start(self, indep_edges: list, vertex_cover: list):
        assign = dict()
        for i, s in enumerate(vertex_cover):
            edges = []
            for e in self.graph.edges(s):
                a, b = e
                self.x[a, b, i].start = 1
                self.x[b, a, i].start = 0
                self.y[a, i, 0].start = 1
                self.y[a, i, 1].start = 0
                self.y[b, i, 1].start = 1
                self.y[b, i, 0].start = 0
                edges.append(e)
            assign[i] = edges
        indep_edges = set(indep_edges)
        if self._config.edge_fix:
            for i in assign.keys():
                for e in assign[i]:
                    if e in indep_edges:
                        a, b = e
                        self.x[a, b, i].lb = 1
                        self.x[b, a, i].ub = 0
                        self.y[a, i, 0].lb = 1
                        self.y[a, i, 1].ub = 0
                        self.y[b, i, 1].lb = 1
                        self.y[b, i, 0].ub = 0

    def _add_callback(self):
        self.m._k = self.upper_bound
        self.m._y = self.y
        self.m._z = self.z
        self.m._g = self.graph
        self.m._g2 = self.power_graph
        self._callback = indep_callback

    @classmethod
    def name(cls) -> str:
        return 'Extended Model'


class CGModel(BaseMinimumBicliqueCoverSolver):
    _INFEASIBLE_STATUSES = [GRB.INFEASIBLE, GRB.INF_OR_UNBD, GRB.INTERRUPTED, GRB.LOADED]
    _TIMEOUT_STATUSES = [GRB.TIME_LIMIT, GRB.INTERRUPTED]

    def __init__(self, g: nx.Graph, config: RunConfig):
        super().__init__(g, config)
        self._solution_bicliques = []
        self._obj_val = None
        self._obj_bound = None
        self._status = None
        self._columns_added = 0
        self._master_time = 0
        self._pricing_time = 0
        self._vertices = list(self._g.nodes)
        self._edges = list({(min(u, v), max(u, v)) for u, v in self._g.edges})
        # non edges
        all_pairs = {(min(u, v), max(u, v)) for u in self._vertices for v in self._vertices if u != v}
        self._non_edges = all_pairs - set(self._edges)

    def _do_solve(self) -> None:
        master = gp.Model("master")
        master.setParam('OutputFlag', 0)

        # Start with singleton bicliques (one per edge)
        bicliques = [{e} for e in self._edges]
        x = {}
        for i, B in enumerate(bicliques):
            x[i] = master.addVar(lb=0, ub=1, obj=1.0, vtype=GRB.CONTINUOUS, name=f"x_{i}", column=None)

        cover_constraints = {}
        for e in self._edges:
            cover_constraints[e] = master.addConstr(
                gp.quicksum(x[i] for i, B in enumerate(bicliques) if e in B) >= 1,
                name=f"cover_{e}"
            )
        master.update()

        # Column Generation Loop
        while True:
            elapsed_time = self._master_time + self._pricing_time
            remaining_time = self._config.time_limit - elapsed_time
            if remaining_time <= 0:
                break

            master.setParam('TimeLimit', remaining_time)
            _, m_time = chronometer(lambda: master.optimize())
            self._master_time += m_time
            if master.status in self._TIMEOUT_STATUSES:
                break

            duals = {e: cover_constraints[e].Pi for e in self._edges}

            remaining_time -= m_time
            if remaining_time <= 0:
                break
            new_biclique, p_time = chronometer(self._solve_pricing, duals, remaining_time)
            self._pricing_time += p_time

            if new_biclique is None:
                break

            idx = len(x)
            x[idx] = master.addVar(lb=0, ub=1, obj=1.0, vtype=GRB.CONTINUOUS, name=f"x_{idx}", column=None)
            for e in new_biclique:
                master.remove(cover_constraints[e])
                cover_constraints[e] = master.addConstr(
                    gp.quicksum(x[i] for i, B in enumerate(bicliques) if e in B) + x[idx] >= 1,
                    name=f"cover_{e}"
                )
            bicliques.append(new_biclique)
            master.update()
            self._columns_added += 1

        self._status = master.status
        self._obj_bound = master.ObjBound
        try:
            self._obj_val = master.ObjVal \
                if master.status not in self._TIMEOUT_STATUSES and master.status not in self._INFEASIBLE_STATUSES \
                else master.ObjBound
        except gp.GurobiError:
            self._obj_val = master.ObjBoundC
        self._solution_bicliques = [bicliques[i] for i in x if x[i].X > 1e-6]

    def _solve_pricing(self, duals, time_limit: int | None = None):
        pricing = gp.Model("pricing")
        pricing.setParam('OutputFlag', 0)
        pricing.setParam('TimeLimit', time_limit)

        y = {u: pricing.addVar(lb=0, ub=1, vtype=GRB.BINARY, name=f"y_{u}", column=None, obj=0)
             for u in self._vertices}
        w = {u: pricing.addVar(lb=0, ub=1, vtype=GRB.BINARY, name=f"w_{u}", column=None, obj=0)
             for u in self._vertices}
        z = {e: pricing.addVar(lb=0, ub=1, vtype=GRB.BINARY, name=f"z_{e}", column=None, obj=duals[e])
             for e in self._edges}

        pricing.modelSense = GRB.MAXIMIZE

        for (u, v) in self._edges:
            pricing.addConstr(z[(u, v)] <= y[u] + y[v])
            pricing.addConstr(z[(u, v)] <= w[v] + w[u])
            pricing.addConstr(z[(u, v)] >= y[u] + w[v] - 1)
            pricing.addConstr(z[(u, v)] >= y[v] + w[u] - 1)

        for (u, v) in self._non_edges:
            pricing.addConstr(y[u] + w[v] <= 1)
            pricing.addConstr(y[v] + w[u] <= 1)

        for u in self._vertices:
            pricing.addConstr(y[u] + w[u] <= 1, name=f"y_{u} + w_{u} <= 1")

        pricing.optimize()

        if pricing.status != GRB.OPTIMAL or pricing.ObjVal <= 1:
            return None

        y_resolved = [u for u in self._vertices if y[u].X > 0.5]
        w_resolved = [u for u in self._vertices if w[u].X > 0.5]
        biclique = {
            (min(u, v), max(u, v)) for u in y_resolved for v in w_resolved if (min(u, v), max(u, v)) in self._edges}
        return biclique

    def _do_get_solution(self) -> float:
        return self._obj_val

    def is_feasible(self) -> bool:
        return self._status not in self._INFEASIBLE_STATUSES

    def _check_biclique_cover(self) -> bool:
        # check if it's a cover
        covered_edges = set()
        for edge_set in self._solution_bicliques:
            covered_edges.update((min(u, v), max(u, v)) for u, v in edge_set)
        graph_edges = set((min(u, v), max(u, v)) for u, v in self._g.edges)
        if not graph_edges.issubset(covered_edges):
            return False
        # check if each set of edges forms a biclique
        for edge_set in self._solution_bicliques:
            if not is_biclique(self._g, edge_set):
                return False
        return True

    def _post_solve(self) -> None:
        if self._status == GRB.TIME_LIMIT:
            print(f'Model reached time limit of {self._config.time_limit} seconds.\n')
        elif self._status == GRB.INFEASIBLE:
            print('Model is unfeasible.\n')
        elif self._status == GRB.OPTIMAL:
            print("Final LP Objective (Lower Bound on MBCP):", self._obj_val, '\n')
            for i, biclique in enumerate(self._solution_bicliques):
                print(f"Biclique {i}: {sorted(biclique)}")
            print(f'Is it a biclique cover? {"Yes" if self._check_biclique_cover() else "No"}.\n')

    @property
    def columns_added(self) -> int:
        if not self._solved:
            raise NotYetSolvedError()
        return self._columns_added

    @property
    def master_time(self) -> float:
        if not self._solved:
            raise NotYetSolvedError()
        return self._master_time

    @property
    def pricing_time(self) -> float:
        if not self._solved:
            raise NotYetSolvedError()
        return self._pricing_time

    @property
    def best_known_solution(self):
        if not self._solved:
            raise NotYetSolvedError()
        return self._obj_bound

    @classmethod
    def name(cls) -> str:
        return 'Column Generation'
