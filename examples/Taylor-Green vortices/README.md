### Taylor-Green vortices

The case of Taylor-Green vortices (TGV) is useful for many things. It can be used for code validation, as it is an analytical solution to the Navier-Stokes equations. It can be used to study dacay of vortices. Actually the TGV case serves as prototype model for transition and turbulence decay. It can also illustrate domain bounary coupling using cyclic or periodic (depending how you want to call it) since the case is periodic in three main coordinate-directions.

Expected learning outcome:
- Learn how to setup and run the case (if you haven't already :o) ).
- Set up a case with periodic boundaries (this is triply periodic case).
- Have a taste of a long run non-stationary case.
- Do some postprocessing of the results to compare with published ones.
- Learn how to handle situation when dacaying vortices become small relative to grid size (advanced).

To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh block.geo -3 -o block.su2`

Note: We have saved mesh in _SU2_ format - which is possible in Gmsh. We do that becasue of the mesh conversion program that we have, and which we will call in the following step. If you haven't done it - you will need to compile _su2ToCappuccino_ which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation. 

Create the folders:  
`mkdir polyMesh vtk vtk/mesh 0`

Run mesh converter:  
`su2ToCappuccino block`


Edit polyMesh 'boundary' file to set boundary types like below:

```
# bcName bcType nFaces startFace
in empty 3969 750015
out periodic 3969 753984 750015
left empty 4032 757953
right periodic 4032 761985 757953
top empty 4032 766017
bottom periodic 4032 770049 766017
```
A note on periodic boundaries. When two sides of the domain are coupled in cyclic way, we point this out with one of them, adding one more number to the end of the line denoting the startFace of the twin boundary. Then the type description of the twin boundary is not needed and we put 'empty' for the type. For example:
```
# bcName bcType nFaces startFace
top empty 4032 766017
bottom periodic 4032 770049 766017
```
Here the boundary named 'bottom' is denoted as 'periodic' in type and the startFace number of its twin boundary (named 'top') is added to the end.

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit them. NOTE: We will create for the start a file where velocity throughout the domain iz zero, and we will correct this, as described below, by setting the initial field values from the `cappuccino` code itself.

In 'U' file, change settings to these:  

```
internalField
  uniform
    0.0 0.0 0.0
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
  Neumann
    zeroGradient
bottom
  Neumann
    zeroGradient
```    

With this case we illustrate how to set initial values directly through the code. In `src/cappuccino/init.f90` you will find the following lines, properly positioned after we attempted to initialize velocity vector field by reading `0/U` file. 

Uncomment them and the lines read:
```
    ! Initialize field values for Taylor-Green vortex
    do inp=1,numCells
      call initialize_tgv( xc(inp),yc(inp),zc(inp),u(inp),v(inp),w(inp),p(inp) )
    enddo
```
This calls specially written initialisation subroutine found in `src/finiteVolume/initialization/field_initialization.f90` module.

After this you will need to recompile the `cappuccino` code (takes ~20s). You can also make new `0/U` file with a quick and dirty temporary change to the `init.f90` file which would output initial values to a file.

Now we are ready to go!

Run the code by executing the run script:

`./run`

While running you may visualy inspect the residuals by calling:

`gnuplot plotResiduals`

You can also plot the development of the monitored variable - dissipation of the mean resolved kineric energy,

`gnuplot plotMonitors`

As simulation is progressing, new files will be written into `vtk` folder.
