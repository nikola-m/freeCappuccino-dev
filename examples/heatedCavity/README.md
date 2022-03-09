### Side heated cavity

This case demonstrates the capability of the code to solve problems of heat transfer, more specifically of natural convection.

Natural convection in a side wall heated cavity, also known as differentially heated cavity is studied. The cavity has aspect ratio Ar equal to one. A constant temperature is maintaied on a side wall of the cavity. Temperature circulation is developed within the domain and flow reaches steady state. 

Flow is laminar and we employ Boussinesq approximation for buoyancy terms.

Rayleigh number is set roughly to Ra=10^4 (It would be exactly 10^4 if gravitational acceleration is set to 10).

To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh cavity.geo -3 -o cavity.su2`

Note: We have saved mesh in su2 format - which is possible in Gmsh. We do that becasue of the conversion utility that we have, which is called in following step. If you haven't done it - you will need to compile _su2ToCappuccino_ which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation:

`mkdir polyMesh vtk vtk/mesh 0`  
`su2ToCappuccino cavity`  


Edit polyMesh 'boundary' file to set boundary types to:

```
# bcName bcType nFaces startFace
wall-top wall 39 2964
wall-left wall 39 3003
wall-bottom wall 39 3042
wall-right wall 39 3081
front-back empty 3042 3120 
```

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them:

Rename file with template for scalar field initialization to 'T' for temperature
`mv 0/T.template 0/T`

Copy settings provided below:  
```
internalField
  uniform
    273.0
boundaryField
wall-top
  Neumann
    zeroGradient
wall-left
  Dirichlet
    uniform
      273.0
wall-bottom
  Neumann
    zeroGradient
wall-right
  Dirichlet
    uniform
      373.0
front-back
  Neumann
    zeroGradient  
```
In 'U' file, change settings to these:  
```
internalField
  uniform
    0.0 0.0 0.0
boundaryField
wall-top
  Dirichlet
    uniform
      0.0 0.0 0.0
wall-left
  Dirichlet
    uniform
      0.0 0.0 0.0
wall-bottom
  Dirichlet
    uniform
      0.0 0.0 0.0
wall-right
  Dirichlet
    uniform
      0.0 0.0 0.0
front-back
  Neumann
    zeroGradient 
```
Simulation configuration can be found in the 'input.nml' file.

Now we are ready to go!

Run the code by executing the run script:

`./run`

While running you may visualy inspect the residuals by calling:

`gnuplot plotResiduals`

or by openning 'monitor' file, where simulation progress is logged, eg by typing:

`tail -n 100 monitor`

When done, open the simulation files from the 'vtk' folder in Paraview.

