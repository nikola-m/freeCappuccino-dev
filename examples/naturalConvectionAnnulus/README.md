### Natural convection between two concentric cylinders

This case demonstrates the capability of the code to solve problems of heat transfer, more specifically of natural convection.

Natural convection in a heated vertical concentric annulus is studied. Inner cylinder is heated in a way that constant temperature is maintaied. Temperature stratification is developed within the domain and flow reaches steady state. Flow development to steady state is simulated by time-stepping.

Flow is laminar and we employ Boussinesq approximation for buoyancy terms.

Rayleigh number is set to Ra=78480.

Good exercise would be to change parameters which corresponds to Ra from an experimental study and verify the simulation results.

To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh cylinder_annular.geo -3 -o cylinder_annular.su2`

Note: We have saved mesh in su2 format - which is possible in Gmsh. We do that becasue of the conversion utility that we have, and which we will call in following step. If you haven't done it - you will need to compile _su2ToCappuccino_ which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation:

`mkdir polyMesh vtk vtk/mesh 0`  
`su2ToCappuccino cylinder_annular`  


Edit polyMesh 'boundary' file to set boundary types to:

```
# bcName bcType nFaces startFace  
wallInner wall 156 12012  
front symmetry 6084 12168  
back symmetry 6084 18252  
wallOuter wall 156 24336  
```

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them:

Rename file with template for scalar field initialization to 'T' for temperature
`mv 0/T.template 0/T`

Copy settings provided below:

```  
internalField  
  uniform  
  300.0  
boundaryField  
wallInner  
   Dirichlet  
    uniform  
     302.0  
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
```
In 'U' file, change settings to these:  
```
internalField  
  uniform  
  0.0 0.0 0.0  
boundaryField  
wallInner  
   Dirichlet  
    uniform  
     0. 0. 0.  
front  
   Neumann  
    zeroGradient  
back  
   Neumann  
    zeroGradient  
wallOuter  
   Dirichlet  
    uniform  
     0. 0. 0.  
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
