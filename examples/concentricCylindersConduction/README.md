### Steady state heat conduction

Steady state problem of source-free heat conduction in a concentric cylinder whose inner and outer walls are maintained at constant temperature of 400 K and 300 K respectively. (source: https://caefn.com/openfoam/solvers-laplacianfoam, last visited Feb-2022)

gmsh cylinder_annular_2to1.geo -3 -o cylinders.su2 
mkdir vtk vtk/mesh 0 polyMesh
su2ToCappuccino cylinders

(make sure you have compiled su2ToCappuccino mesh converter)

#
# Now:
# 1) edit polyMesh/boundary to set desired BCs

# bcName bcType nFaces startFace
wallInner wall 156 18252
front empty 9204 18408
back empty 9204 27612
wallOuter wall 156 36816


# 2) edit 0/U to set initial values for inner cells and boundary

0/U:

internalField
  uniform
    0.0 0.0 0.0
boundaryField
wallInner
  Dirichlet
    uniform
      0.0 0.0 0.0
front
  Neumann
    zeroGradient
back
  Neumann
    zeroGradient
wallOuter
  Dirichlet
    uniform
      0.0 0.0 0.0
      
0/T:
internalField
  uniform
    0.0
boundaryField
wallInner
  Dirichlet
    uniform
      400.0
front
  Neumann
    zeroGradient
back
  Neumann
    zeroGradient
wallOuter
  Dirichlet
    uniform
      300.0      
      
# 3) make input file of .nml extension based on input file formatting guide

The input file to run this case is provided.

# 4) run (e.g. writting 'cappuccino input.nml monitor restart' or using the run script)

Run script is provided in this folder.
