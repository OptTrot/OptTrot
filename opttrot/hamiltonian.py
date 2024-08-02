from typing import *
from pathlib import Path

import pandas as pd
import networkx as nx  # Change the module to rustworkx.
from itertools import combinations
from typing import Iterable
from .core.pauli import PauliPoly, PauliElement  # Assuming the PauliPoly class is defined in pauli.py


# New Pauli term must be update to edge list dataframe.
# 
class Hamiltonian(PauliPoly):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.local_decomposition = self.get_decomposition(self._terms)
        self._nx_graph = None  # Cache for NetworkX graph
    @classmethod
    def from_file(cls, filepath:Union[str, Path]):
        if not isinstance(filepath, Path):
            filepath = Path(filepath)
        #return cls()    
        raise NotImplementedError
    @property
    def compatible_graph(self):
        if self._nx_graph is None:
            self._nx_graph = self.to_networkx_graph()
        return self._nx_graph

    def draw_graph(self, *args, **kwargs):
        if self._nx_graph is None:
            self._nx_graph = self.to_networkx_graph()
        return nx.draw(self._nx_graph, *args, **kwargs)
    
    @staticmethod
    def get_decomposition(pauli_basis: Iterable[PauliElement]):
        p_dict = {}
        for p in pauli_basis:
            nx, nz = p.scode
            num = 0 
            num += 1 if nz > 0 else 0
            num += 2 if nx > 0 else 0
            p_dict[p.string] = (num, nz, nx, p.np_coef)
        df = pd.DataFrame.from_dict(
            p_dict, 
            orient="index",
            columns=["type", "Z", "X", "Coef"]
        )
        df.reset_index(inplace=True, names="Pstring")
        return df

    def to_networkx_graph(self):
        if len(self._terms) == 1:
            raise RuntimeError("You cannot make single graph.")
        df = self.local_decomposition
        edge_df = pd.DataFrame(combinations(df["Pstring"].values, 2), columns=['source', 'target'])
        edge_df = edge_df.merge(df[["Pstring", "Z", "X"]], how="left", left_on="source", right_on='Pstring').drop("Pstring", axis=1)
        edge_df.rename(columns={"Z": "Zs", "X": "Xs"}, inplace=True)
        edge_df = edge_df.merge(df[["Pstring", "Z", "X"]], how="left", left_on="target", right_on='Pstring').drop("Pstring", axis=1)
        edge_df.rename(columns={"Z": "Zt", "X": "Xt"}, inplace=True)
        edge_df["commute"] = edge_df.apply(lambda row: int(self.commute_reggio_df(row[["Zs", "Xs", "Zt", "Xt"]])), axis=1)
        
        G = nx.from_pandas_edgelist(edge_df, 'source', 'target', edge_attr='commute')
        return G
    
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
