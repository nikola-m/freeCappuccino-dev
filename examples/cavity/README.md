### Lid driven cavity

This is a well known case of laminar flow of viscous incompressible fluid in closed cavity with moving upper wall.

Reynolds number is set to Re=100, based on cavity side length.



To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh cavity.geo -3 -o cavity.su2`

Note: We have saved mesh in su2 format - which is possible in Gmsh. We do that becasue of the conversion utility that we have, which is called in following step. First, you will need to compile _su2ToCappuccino_ which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation:

`mkdir polyMesh vtk vtk/mesh 0`  
`su2ToCappuccino cavity`  


Edit _polyMesh/boundary_ file to set boundary types to:

```
# bcName bcType nFaces startFace
wall-top wall 39 2964
walls wall 117 3003
front-back empty 3042 3120
```

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them:

In 'U' file, change settings to these:  
```
internalField
  uniform
    0.0 0.0 0.0
boundaryField
wall-top
  Dirichlet
    uniform
      1.0 0.0 0.0
walls
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
