from typing import *
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import rustworkx as rx

#from tqdm import tqdm
from opttrot.graph_algs.graph_opts import MaxCliqueSolver
from opttrot.pauli import PauliPoly, PauliElement  # Assuming the PauliPoly class is defined in pauli.py


# New Pauli term must be update to edge list dataframe.
# Replace the pandas edge to adjacency matrix.
class Hamiltonian(PauliPoly):
    def __init__(self, *args, **kwargs):
        super(Hamiltonian, self).__init__(*args, **kwargs)
        self.local_decomposition = self.get_decomposition(self.poly)
        self._adj_mat = None
        self._nx_graph = None  # Cache for NetworkX graph
    def __repr__(self):
        return "Hamiltonian" + super().__repr__()[9:]
    @classmethod
    def from_file(cls, filepath:Union[str, Path]):
        if not isinstance(filepath, Path):
            filepath = Path(filepath)
        #return cls()    
        raise NotImplementedError
    @property
    def ppoly_terms_str(self):
        return [p.pstr for p in self.poly]
    @property
    def adj_mat(self):
        """Adjacency matrix of commuting graph.

        Returns:
            _type_: _description_
        """
        if self._adj_mat is None:
            self._adj_mat  = np.vstack([self.commute(p) for p in self.poly])
        return np.matrix(self._adj_mat).astype(float)
    @property
    def compatible_graph(self):
        #if self._nx_graph is None:
        #    self._nx_graph = self.to_networkx_graph()
        #return self._nx_graph
        return self.to_graph()
    @staticmethod
    def get_decomposition(pauli_basis: Iterable[PauliElement]):
        p_dict = {}
        for p in pauli_basis:
            nx = p.nx 
            nz = p.nz
            num = 0 
            num += 1 if nz > 0 else 0
            num += 2 if nx > 0 else 0
            p_dict[p.pstr] = (num, nz, nx, p.weight)
        df = pd.DataFrame.from_dict(
            p_dict, 
            orient="index",
            columns=["type", "Z", "X", "Coef"]
        )
        df.reset_index(inplace=True, names="Pstring")
        return df
    
    def to_graph(self, anti=False, networkx = False):
        if networkx:
            adj = 1-self.adj_mat if anti else self.adj_mat
            return nx.from_numpy_array(adj)
        if anti:
            G = rx.PyGraph.from_adjacency_matrix(1-self.adj_mat, 0)
        else:
            G = rx.PyGraph.from_adjacency_matrix(self.adj_mat, 0)
        return G
    def get_commuting_group(self, solver:Union[None, MaxCliqueSolver]=None):
        # Qiskit method
        pstrs = self.ppoly_terms_str
        if solver is None:
            coloring_dict = rx.graph_greedy_color(self.to_graph(anti=True))
            groups = defaultdict(list)
            for idx, color in coloring_dict.items():
                groups[color].append(pstrs[idx])
            return groups
        return solver.solve_via_clique(self.to_graph(anti = False, networkx = False))
    def draw_graph(self, *args, **kwargs):
        if self._nx_graph is None:
            self._nx_graph = self.to_networkx_graph()
        return nx.draw(self._nx_graph, *args, **kwargs)
    
    @staticmethod
    def commute_reggio_df(s):
        a = bin(s.iloc[0] & s.iloc[3]).count("1")%2
        b = bin(s.iloc[1] & s.iloc[2]).count("1")%2
        return a == b

# Example usage
if __name__ == "__main__":
    # Create a Hamiltonian using the PauliPoly interface
    # Example Pauli basis and Hamiltonian matrix
    pauli_basis = {...}  # Define the Pauli basis as needed
    H_matrix = np.matrix([[1, 0], [0, -1]])  # Example Hamiltonian matrix

    hamiltonian = Hamiltonian(H_matrix, pauli_basis=pauli_basis)
    
    # Access the compatible graph
    graph = hamiltonian.compatible_graph
    
    # Draw the graph
    hamiltonian.draw_graph()
