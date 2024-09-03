import os, sys
sys.path.append("..")
sys.path.append("..\\..")

from typing import Tuple
import pickle
import multiprocessing
from itertools import permutations
from copy import copy

import numpy as np
from tqdm import tqdm


def sym_code2pstr(ns:Tuple[int, int], l:int)->str:
        assert l>0, "l must be positive integer and greater than 0."
        nx, nz = ns
        max_int_1 = 2**l
        assert (nx < max_int_1 and nz < max_int_1), "The given integers and the qubit dim are not matched."
        if nx==0: # Z family
            st = format(nz, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "Z")
            return st
        if nz==0: # X family
            st = format(nx, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "X")
            return st
        # None of above
        st_x = format(nx, f"0{l}b")
        st_z = format(nz, f"0{l}b")
        result = []
        for x, z in zip(st_x, st_z):
            if x == z:
                if x =="1":
                    result.append("Y")
                else: 
                    result.append("I")
            elif x > z:
                result.append("X")
            else:
                result.append("Z")
        return "".join(result)

class PauliFrame:
    def __init__(self, n):
        self.qubit = n

        self.front = np.vstack([np.zeros(self.qubit), 2**np.arange(self.qubit)]).T.astype(int)
        self.back  = np.vstack([2**np.arange(self.qubit), np.zeros(self.qubit)]).T.astype(int)
        self.sign = np.ones(self.qubit, dtype=int)

        self.history = []
    def __repr__(self):
        reps = []
        for i, (x1, z1, x2, z2)  in enumerate(np.hstack([self.front, self.back])):
            reps.append(f"{self.sign[i]}, {sym_code2pstr((x1, z1),self.qubit)}, {sym_code2pstr((x2, z2), self.qubit)}\n")
        return "".join(reps)
    def _tuple_sort_triple(self, ts):
        return tuple(sorted(tuple(sorted(t)) for t in ts))
    def __eq__(self, other):
        self_set = set()
        other_set = set()
        for si, sj in zip(self.front, self.back):
            sk, _ = self._product(si, sj)
            self_set.add(self._tuple_sort_triple([si, sj, sk]))
        for si, sj in zip(other.front, other.back):
            sk, _ = self._product(si, sj)
            other_set.add(self._tuple_sort_triple([si, sj, sk]))
        return self_set <= other_set
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        self_list = list()
        for si, sj in zip(self.front, self.back):
            sk, _ = self._product(si, sj)
            self_list.append(self._tuple_sort_triple([si, sj, sk]))
        sorted(self_list)
        return hash(tuple(self_list))
    @property
    def matrix(self):
        return np.hstack([self.front, self.back])
    
    def relative_supp(self, p):
        supp = 0
        for (si, si_t) in zip(self.front, self.back):
            supp += int(self._commute(si, p)|self._commute(si_t, p))
        return supp

    def copy(self):
        pf  = PauliFrame(self.qubit)
        pf.front = copy(self.front)
        pf.back = copy(self.back)
        pf.sign = copy(self.sign)
        pf.history = copy(self.history)
        return pf
    def get_previous(self):
        pf = self.copy()
        (i, j, c_type, hash_key) = pf.history[-1]
        pf.two_entangle(i, j, c_type)
        del(pf.history[-1])
        del(pf.history[-1])
        return pf
    def h(self, i):
        self.back[i][0], self.front[i][0] = self.front[i][0], self.back[i][0]
        self.back[i][1], self.front[i][1] = self.front[i][1], self.back[i][1]
    def s(self, i):
        si, si_t = (self.front[i][0], self.front[i][1]), (self.back[i][0], self.back[i][1])
        new_si_t, _ = self._product(si, si_t)
        self.back[i][0] = new_si_t[0]
        self.back[i][1] = new_si_t[1]
        pass
    def s_d(self, i):
        self.s(i)
        self.sign[i] *= -1  
        pass 
    def _commute(self, p1, p2):
        nx1, nz1 =  p1
        nx2, nz2 =  p2
        a = bin(nx1&nz2).count("1")%2
        b = bin(nx2&nz1).count("1")%2
        return a==b

    def _product(self, p1, p2):
        nx1, nz1 = p1
        nx2, nz2 = p2

        # update sign 
        #f=0
        #tmp1 = nx1 & nz2
        #tmp2 = nx2 & nz1
        #f += bin(tmp1).count("1")
        #f += bin(tmp2).count("1")
        #f += 2*bin(nx1&nz2).count("1")

        return ((nx1 ^ nx2, nz1 ^ nz2), 1)
    
    def  _basis(self, i, j, c_type, front=True):
        if c_type ==0:
            return None
        if c_type ==1:
            #cz
            self.h(j)
            return None
        if c_type == 2:
            # cy
            if front:
                self.s_d(j)
            else:
                self.s(j)
            return None
        if c_type ==3:
            self.h(i)
            return None
        if c_type ==4:
            self.h(i)
            if front:
                self.s_d(j)
            else:
                self.s(j)
            return None
        if c_type==5:
            if front:
                self.s_d(i)
                self.h(i)
                self.s_d(j)
            else:
                self.s(i)
                self.h(i)
                self.s(j)
            return None


    def two_entangle(self, i, j, c_type = 0):
        assert i>=0 and j>=0, "i, j must be positive integer."
        assert i<self.qubit and j<self.qubit, "Exceeding the maximum index."

        self._basis(i, j , c_type, front=True)
        
        # CX
        si, si_t = (self.front[i][0], self.front[i][1]), (self.back[i][0], self.back[i][1])
        sj, sj_t = (self.front[j][0], self.front[j][1]), (self.back[j][0], self.back[j][1])

        new_sj, _ = self._product(sj, si)
        new_si_t, _ = self._product(si_t, sj_t)

        self.front[j][0] = new_sj[0]
        self.front[j][1] = new_sj[1]
        self.back[i][0] = new_si_t[0]
        self.back[i][1] = new_si_t[1]
        
        self._basis(i, j , c_type, front=False)

        self.history.append((i, j, c_type, hash(self)))

def process_permutations(args):
    p, N = args
    entangled_sets = set()
    for i, j in permutations(range(N), 2):
        p1, p2, p3, p4, p5, p6 = p.copy(), p.copy(), p.copy(), p.copy(), p.copy(), p.copy()

        p1.two_entangle(i, j, c_type=0)
        p2.two_entangle(i, j, c_type=1)
        p3.two_entangle(i, j, c_type=2)
        p4.two_entangle(i, j, c_type=3)
        p5.two_entangle(i, j, c_type=4)
        p6.two_entangle(i, j, c_type=5)

        entangled_sets.update([p1, p2, p3, p4, p5, p6])

    return entangled_sets

def process_permutations_Z(args):
    p, N = args
    entangled_sets = set()
    for i, j in permutations(range(N), 2):
        p1 = p.copy()
        p1.two_entangle(i, j, c_type=0)
        entangled_sets.update([p1])

    return entangled_sets


def run(N, ppols=2, cx=False, dir="data"):
    pf = PauliFrame(N)
    Frames = [set([pf])]
    i=0
    frames_len = 1

    num_core = ppols
    frames = set()
    frames = frames | (Frames[0])

    process_func = process_permutations_Z if cx else process_permutations
    
    with open(f"{dir}\\qubit{N}_PauliFrame.data", "wb") as file:
            pickle.dump((N), file)

    for k in range(0, 15):
        print(f"{k}-th:", end="")
        ith_set = set()
        ends = []
        # Create a pool of worker processes
        with multiprocessing.Pool(num_core) as pool:
            # Map the process_permutations function to each frame in parallel
            #results = pool.starmap(process_permutations, [(p, N) for p in Frames[k]])
            with tqdm(total=len(Frames[0]), desc=f"Processing Frame {k}") as pbar:
                # Use imap_unordered for real-time updates
                for result in pool.imap_unordered(process_func, [(p, N) for p in Frames[0]]):
                    if result in ith_set:
                        ends.append(result)
                    else:
                        ith_set.update(result)
                    pbar.update()
        # Combine results from all processes
        #for result in results:
        #    ith_set.update(result)

        print(len(ith_set))
        with open(f"{dir}\\qubit{N}_{k}_PauliFrame.data", "wb") as file:
            pickle.dump((i, ith_set), file)

        with open(f"{dir}\\qubit{N}_{k}_PauliFrame_ends.data", "wb") as file2:
            pickle.dump((i, ends), file2)
        
        Frames = [(ith_set)]
        if frames >= ith_set:
            break
        frames = frames | Frames[0]
    with open(f"{dir}\\qubit{N}_PauliFrame.data", "+wb") as file:
        pickle.dump((-1, frames), file)
    return True


if __name__ =="__main__":
    run(N=4, )