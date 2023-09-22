# LargeN.jl

compute [large-N approximation](10.1103/PhysRevLett.93.167204) for classical spin models with isotropic interactions.

This package is still under development and may change in the future
## Basis usage
The main usage is the spin-spin correlator which can be obtained via `getChiFunction`
```julia
using LargeN, SpinFRGLattices

using StaticArrays
J(qx,qy) = SA[
    cos(qx)+cos(qy) sin(qx-qy);
    sin(qx-qy) cos(qx)+cos(qy)
]
J(q::SVector{2}) = J(q[1],q[2]) 
NCell = 2
Dim = 2
chi = getChiFunction(T,J,NCell,Dim)
```

One can also use implementations from `SpinFRGLattices` if one does not want to compute $J(q)$ by hand.
```julia

Py = Pyrochlore.getPyrochlore(4,[1,0])

T = 0.1
chi = getChiFunction(T,Py,Pyrochlore)

q = [1,0,0]
qStatic = SArr.SA[1,0,0] # non-allocating static array
println(chi(q) == chi(1,0,0) == chi(qStatic))
```

Output:
```julia
Total Number of sites: 151       Num pairs: 18
ground state degeneracy: 16008. Degeneracy per unit cell: 2.001
-Emin/T: 20.00000000000002
[ Info: constraint minimized to 2.0539125955565396e-15 for λ = 21.588450780058196
true
```
## Some References
* [Dipolar Spin Correlations in Classical Pyrochlore Magnets](10.1103/PhysRevLett.93.167204)
* [Spin-Ice Thin Films: Large-N Theory and Monte Carlo Simulations](10.1103/PhysRevX.8.021053)
