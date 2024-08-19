



## Poly Algebra with numpy array

### Addition

P1 + P2 only add weight if P1(nx, nz) = P2(nx, nz).

```
if P1 in Parray:
    Parray = P1 + Parray
else:
    Parray = np.concatenate([P1, Parray])
```