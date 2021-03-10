### Laminar flow over circular cylinder at Re=200 - von Karman vortices

To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh cylinder_ogrid.geo -3 -o cylinder_ogrid.su2`

Note: We have saved mesh in _SU2_ format - which is possible in Gmsh. We do that becasue of the mesh conversion program that we have, and which we will call in the following step. If you haven't done it - you will need to compile _su2ToCappuccino_ which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation. 

Create the folders:  
`mkdir polyMesh vtk vtk/mesh 0`

Run mesh converter:  
`su2ToCappuccino cylinder_ogrid`


Edit polyMesh 'boundary' file to set boundary types like below:

```
# bcName bcType nFaces startFace
cylinder wall 156 18252
front empty 9204 18408
back empty 9204 27612
Outlet outlet 78 36816
Inlet inlet 78 36894
```

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them.

In 'U' file, change settings to these:  

```
internalField
  uniform
    0.0 0.0 0.0
boundaryField
cylinder
  Dirichlet
    uniform
      0.0 0.0 0.0
front
  Neumann
    zeroGradient
back
  Neumann
    zeroGradient
Outlet
  Neumann
    zeroGradient
Inlet
  Dirichlet
    uniform
      0.2 0.0 0.0
```

Now we are ready to go!

Run the code by executing the run script:

`./run`

While running you may visualy inspect the residuals by calling:

`gnuplot plotResiduals`

or by openning 'monitor' file, where simulation progress is logged, eg by typing:

`tail -n 100 monitor`

When done, open the simulation files from the 'vtk' folder in Paraview.


Click on an image below to open a YouTube video with animation resulting from this simulation

[![freeCappuccino: von Kármán vortex street at Re=200](https://img.youtube.com/vi/hf6r8MqTjZo/0.jpg)](https://www.youtube.com/watch?v=hf6r8MqTjZo)

or this 

[![Relax your mind and focus with more than 6 hours of von Kármán vortex street](https://img.youtube.com/vi/cj4_ZOQbNsY/0.jpg)](https://www.youtube.com/watch?v=cj4_ZOQbNsY)




