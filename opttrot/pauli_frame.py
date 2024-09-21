from copy import copy

import numpy as np
from .pauli_utils import sym_code2pstr
from .utils import int2bin

class PauliFrame:
    def __init__(self, n):
        self.qubit = n

        self.front = np.vstack([np.zeros(self.qubit), 2**np.arange(self.qubit)]).T.astype(int)
        self.back  = np.vstack([2**np.arange(self.qubit), np.zeros(self.qubit)]).T.astype(int)
        self.sign = np.ones(self.qubit, 2, dtype=int)

        self.history = []
    def __repr__(self):
        reps = []
        for i, (x1, z1, x2, z2)  in enumerate(np.hstack([self.front, self.back])):
            reps.append(f"{self.sign[i]}, {sym_code2pstr((x1, z1),self.qubit)}, {sym_code2pstr((x2, z2), self.qubit)}\n")
        return "".join(reps)
    def _tuple_sort_triple(self, ts):
        """For __eq__ and __hash__ routine"""
        return tuple(sorted(tuple(sorted(t)) for t in ts))
    def __eq__(self, other):
        """SWAP and H, S, S^\dagger must not affect to the equivalent states."""

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
    @property
    def aug_matrices_front(self):
        front_x = np.vstack([int2bin(x, self.qubit) for x in self.front[:, 0]])
        front_x = np.flip(front_x, axis=0)
        front_z = np.vstack([int2bin(z, self.qubit) for z in self.front[:, 1]])
        front_z = np.flip(front_z, axis=0)
        return front_x.T, front_z.T

    def relative_supp(self, p, count="both"):
        supp = 0
        if count=="both":
            for (si, si_t) in zip(self.front, self.back):
                l1 = self._commute(si, p)
                l2 = self._commute(si_t, p)
                supp += int(not(l1 & l2))
        elif count == "front":
            for (si) in self.front:
                l1 = self._commute(si, p)
                supp += int(not(l1))
        elif count == "back":
            for (si_t) in self.back:
                l2 = self._commute(si_t, p)
                supp += int(not(l2))
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

        self.sign[i][0], self.sign[i][1] = self.sign[i][1], self.sign[i][0]

    def s(self, i):
        si, si_t = (self.front[i][0], self.front[i][1]), (self.back[i][0], self.back[i][1])
        new_si_t, _ = self._product(si, si_t)
        self.back[i][0] = new_si_t[0]
        self.back[i][1] = new_si_t[1]
        pass
    def s_d(self, i):
        self.s(i)
        self.sign[i, 1] = -self.sign[i, 1]  
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