using Brillouin, Crystalline
using SpinFRGLattices
Basis = SpinFRGLattices.TriangularLattice.Basis
Basisvecs = Array.((Basis.a1,Basis.a2))
##
Recip_Basis = Crystalline.reciprocalbasis(Basisvecs)
Trafo = cat(Recip_Basis ./ 2pi...,dims = (2,2))
# cell = cartesianize!(wignerseitz(Recip_Basis))
##
using Plots

##
vertpoints
scatter((getindex.(cell.verts,i) for i in 1:3)...)