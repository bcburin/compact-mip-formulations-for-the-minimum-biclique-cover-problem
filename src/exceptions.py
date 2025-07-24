

class SolverError(Exception):
    ...


class AlreadySolvedError(SolverError):
    def __init__(self):
        super().__init__("The solve method has already been called.")


class NotYetSolvedError(SolverError):
    def __init__(self):
        super().__init__("The solve method has not yet been called.")


class UnfeasibleSolutionError(SolverError):
    def __init__(self):
        super().__init__("The problem is unfeasible, there is no solution.")
