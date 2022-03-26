import SymmetryReduceBZ.Lattices: genlat_CUB,genlat_FCC
import SymmetryReduceBZ.Symmetry: calc_ibz,calc_bz
a = 2.0
real_latvecs = genlat_FCC(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
ibzformat = "convex hull"
coordinates = "Cartesian"
makeprim = false
convention = "ordinary"
ibz = calc_bz(real_latvecs,atom_types,atom_pos,coordinates,ibzformat,
  makeprim,convention)
scatter(ibz.points[:,1],ibz.points[:,2],ibz.points[:,3])