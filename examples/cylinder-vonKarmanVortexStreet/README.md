### Laminar flow over circular cylinder at Re=200 - von Karman vortices

To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh cylinder.geo -3 -o cylinder.su2`

Note: We have saved mesh in su2 format - which is possible in Gmsh. We do that becasue of the conversion utility that we have, and which we will call in following step. If you haven't done it - you will need to compile this utility - which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation. 
Create the folders:
`mkdir polyMesh vtk vtk/mesh 0`
Run mesh converter:
`su2ToCappuccino cylinder`


Edit polyMesh 'boundary' file to set boundary types like below:

bcName bcType nFaces startFace
wall wall 76 5342
inflow inlet 34 5418
outflow outlet 6 5452
right symmetry 2700 5458
left symmetry 2700 8158


Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them.

In 'U' file, change settings to these:

internalField
  uniform
  0.0 0.0 0.0
boundaryField
cylinder
   Dirichlet
    uniform
     0. 0. 0.
inflow
   Dirichlet
    uniform
     0.1 0. 0.
outflow
   Neumann
    zeroGradient
right
   Neumann
    zeroGradient
left
   Neumann
    zeroGradient


Now we are ready to go!

Run the code by executing the run script:

`./run`

While running you may visualy inspect the residuals by calling:

`gnuplot plotResiduals`

or by openning 'monitor' file, where simulation progress is logged, eg by typing:

`tail -n 100 monitor`

When done, open the simulation files from the 'vtk' folder in Paraview.


Click on an image below to open a YouTube video with animation resulting from this simulation

[![freeCappuccino: von Kármán vortex street at Re=200](https://img.youtube.com/vi/hf6r8MqTjZo/0.jpg](https://www.youtube.com/watch?v=hf6r8MqTjZo)

or this 

[![Relax your mind and focus with more than 6 hours of von Kármán vortex street](https://img.youtube.com/vi/cj4_ZOQbNsY/0.jpg](https://www.youtube.com/watch?v=cj4_ZOQbNsY)




