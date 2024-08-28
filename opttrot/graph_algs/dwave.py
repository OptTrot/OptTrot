import dwave_networkx as dnx
from dwave.system import DWaveSampler, DWaveCliqueSampler, EmbeddingComposite
import rustowkrx as rx

class DWave(MaxCliqueSolver):
    def __init__(self, solver_type:Literal["direct", "clique"] = "clique"):
        self.solver = None
        
        if solver_type == "clqiue":
            self.solver = DWaveCLiqueSampler()
        elif solver_type == "direct":
            self.solver = EmbeddingComposite(DWaveSampler())
        self.solver_type= solver_type
    
    def get_h(self, g:rx.PyGraph):
        pass
    def _max_clique_custom(self, h, arg=list(), kwargs={}):
        if self.solver_type != "direct":
            raise ValueError("This routine only supported for direct solver.")
        pass
    def _max_clique(self, g:rx.PyGraph):
        if self.solver_type !="clique":
            raise ValueError("This routine only supported for direct solver.")
        return dnx.maximum_clique(g, sampler=self.solver)
    
    def solve_max_clique(self, g, arg=list(), kwargs={}):
        custom = kwargs["custom"] if "custom" in kwargs.keys() else False

        gi = g.copy()
        cliques = []
        for i in range(g.num_nodes()): # Max iteration
            if gi.num_nodes() ==0:
                break
            if custom:
                h = self.get_h(gi)
                clique_set = self._max_clique_custom(h, arg)
            else:
                clique_set = self._max_clique(gi)
            gi.remove_nodes(clique_set)
            cliques.append(clique_set)
        return cliques





