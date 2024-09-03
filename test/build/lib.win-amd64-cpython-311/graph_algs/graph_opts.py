from abc import ABCMeta, abstractmethod
from typing import *

import rustworkx as rx
import networkx as nx



class MaxCliqueSolver(metaclass=ABCMeta):
    # single: get max_clique, multiple: get several max_clique
    _alg_type:Literal["single", "multiple"] = "single"
    @abstractmethod
    def solve_via_clique(self, g:rx.PyGraph, *args, **kwargs):
        pass
    
    @abstractmethod
    def solve_clique(self, g:rx.PyGraph, *args, **kwargs):
        pass

class Networkx(MaxCliqueSolver):
    def __init__(self):
        from networkx.algorithms.approximation import max_clique
        self.nx_alg_max_clique = max_clique
        pass
    def solve_clique(self, g: rx.PyGraph, *args, **kwargs):
        # conver the PyGraph to networkx
        G =nx.from_numpy_array(rx.adjacency_matrix(g))
        return self.nx_alg_max_clique(G)
    def solve_via_clique(self, g: Union[rx.PyGraph, nx.Graph], *args, **kwargs):
        G =nx.from_numpy_array(rx.adjacency_matrix(g))
        max_cliques = []
        while len(g.nodes) >0:
            clique = self.nx_alg_max_clique(G)
            max_cliques.append(clique)
            for node in clique:
                G.remove_node(node)
        return max_cliques
    
#class GraphTool(MaxCliqueSolver):
#    def __init__(self):
#        pass
#class QuEra(MaxCliqueSolver):
#    def __init__(self):
#        pass
class DWave(MaxCliqueSolver):
    def __init__(self, **kwargs):
        """See 
        https://docs.ocean.dwavesys.com/en/stable/docs_cloud/reference/resources.html#dwave.cloud.client.Client
        """

        assert "token" in kwargs.keys(), "Token must be required to solve the problem."

        from dwave.cloud import Client
        from dwave.system.samplers import DWaveCliqueSampler
        from dwave.system.composites import EmbeddingComposite
        import dwave_networkx as dnx
        
        self.dwave_client_kwargs = kwargs 
    def solve_clique(self, g: rx.PyGraph, *args, **kwargs):
        G =nx.from_numpy_array(rx.adjacency_matrix(g))
        client = Client(**self.dwave_client_kwargs)
        sampler_cli = DWaveCliqueSampler()
        sol = dnx.maximum_clique(G, sampler = sampler_cli)
        client.close()
        return sol
    
    def solve_via_clique(self, g: rx.PyGraph, *args, **kwargs):
        G =nx.from_numpy_array(rx.adjacency_matrix(g))
        client = Client(**self.dwave_client_kwargs)
        sampler_cli = DWaveCliqueSampler()
        max_cliques = []
        while len(g.nodes) >0:
            clique = dnx.maximum_clique(G, sampler = sampler_cli)
            max_cliques.append(clique)
            for node in clique:
                G.remove_node(node)
        client.close()
        return max_cliques
class QAOA(MaxCliqueSolver):
    def __init__(self):
        pass




    