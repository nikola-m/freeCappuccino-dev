### Turbulent flow in a channel $Re_\tau$ = 395

This is a well known case of a turbulent flow in a channel at Reynolds nuber based on friction velocity Re_\tau = 395 (see e.g. R.Moser, J.Kim, N.N. Mansour, Direct numerical simulation of turbulent channel flow up to Ret=590, Physics of Fluids 11(4) 1999).

Create mesh using Gmsh using the provided .geo file and export mesh in su2 format:
`gmsh channel395.geo -3 -o channel395.su2`
Then do the following to convert it. 
Create the folders:  
`mkdir polyMesh vtk vtk/mesh 0`
Run mesh converter:  
`su2ToCappuccino channel395`


Edit polyMesh 'boundary' file to set boundary types like below:

```
# bcName bcType nFaces startFace
# bcName bcType nFaces startFace
in empty 1953 367007
out periodic 1953 368960 367007
left empty 4032 370913
right periodic 4032 374945 370913
top wall 1984 378977
bottom wall 1984 380961
```
Take notice that we have double periodic domain, 'in' and 'out' are conected as well as 'left' and 'right' boundaries.

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them.

In 'U' file, change settings to these:  

```
internalField
  uniform
    0.1335 0.0 0.0
boundaryField
in
  Neumann
    zeroGradient
out
  Neumann
    zeroGradient
left
  Neumann
    zeroGradient
right
  Neumann
    zeroGradient
top
  Dirichlet
    uniform
      0.0 0.0 0.0
bottom
  Dirichlet
    uniform
      0.0 0.0 0.0
```

For this case we will need to create initial perturbation in the flow to make transition to turbulence easier. To do so edit src/cappuccino/init.f90 file and uncomment the lines which would make initial perturbed field for turbulent pipe:
```
    ! Create initial disturbances for turbulent channel and pipe (below) flow
     do inp = 1,numCells
       call channel_disturbances(xc(inp),yc(inp),zc(inp),u(inp),v(inp),w(inp))
     enddo
```
Recompile the _cappuccino_ code (takes around 20 sec.) and you are ready to run.

Run the code by executing the run script:

`./run`




