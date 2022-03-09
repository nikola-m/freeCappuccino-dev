### Test numerical dissipation 

Repeated test from OpenFOAM verification test series, concerning
convection of a passive scalar in a square 2D domain along a diagonal 
direction using uniform velocity field. 
This is a well known verification case found in many references where convection schemes were discussed.

To enforce only convection to be active for passive scalar, set viscpsity to zero and Prandlt number to high value (eg.1e+30),
which will effectively eliminate diffusion term from the equation. 

Do not solve for velocity field only scalar transport equation.


Numerical schemes tested: cds, kappa, quick, fromm, cui, boundedCentral, linearUpwind, boundedLinearUpwind, boundedLinearUpwind02, muscl, koran, smart, avl-smart, charm, vanleer, ospre, minmod, spl13, and more..

To run the case create mesh using gmsh, convert it to fCp native with following commands:
```
gmsh -3 square.geo -o square.su2
mkdir polyMesh 0 vtk vtk/mesh
su2ToCappuccino square
```
The edit polyMesh/boundary file to edit bc types, and files in 0 folder to set initial values.

polyMesh/boundary file should read:
```
# bcName bcType nFaces startFace
top outlet 50 4900
left inlet 50 4950
bottom inlet 50 5000
right outlet 50 5050
front-back empty 5000 5100
```

0/T file should read:
```
internalField
  uniform
  1e-8
boundaryField
  top
    Neumann
      zeroGradient
  left
    Dirichlet
      uniform
        1.
  bottom
    Dirichlet
      uniform
        0.
  right
    Neumann
      zeroGradient
  front-back
    Neumann
      zeroGradient
```

0/U file should read:
```
internalField
  uniform
  1. 1. 0.
boundaryField
  top
    Neumann
      zeroGradient
  left
    Dirichlet
      uniform
        1. 1. 0.
  bottom
    Dirichlet
      uniform
        1. 1. 0.
  right
    Neumann
      zeroGradient
  front-back
    Neumann
      zeroGradient
```

Run 'run' script to get the results. 

To change discretisation scheme for scalar convection change the value of 'cSchemeT' (convection scheme for temerature as representative scalar field) in input.nml file. 

Examples are:
  'cds'
  'cdscorr'
  'central'
  'linearUpwind'
  'kappa'
  'muscl'
  'umist'
  'koren'
  'smart'
  'avl-smart'
  'charm'
  'vanleer'
  'ospre'
  'minmod'
  'boundedLinearUpwind'
  'boundedLinearUpwind02'
  'boundedCentral'
  'linearUpwindFL'
  'fromm'
  'cui'
  'quick'
  'spl13' 

