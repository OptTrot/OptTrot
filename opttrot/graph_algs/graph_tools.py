from abc import ABCMeta, abstractmethod
from typing import *


class MaxCliqueSolver(ABCMeta):
    # single: get max_clique, multiple: get several max_clique
    _alg_type:Literal["single", "multiple"] = ""
    @abstractmethod
    def max_clique(self, g):
        pass


class Networkx(MaxCliqueSolver):
    def __init__(self):
        pass
class GraphTool(MaxCliqueSolver):
    def __init__(self):
        pass
class RustWorkx(MaxCliqueSolver):
    def __init__(self):
        pass
class DWave(MaxCliqueSolver):
    def __init__(self):
        pass
class QuEra(MaxCliqueSolver):
    def __init__(self):
        pass
class QAOA(MaxCliqueSolver):
    def __init__(self):
        pass




    