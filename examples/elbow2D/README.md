### Laminar flow in a 2D elbow with mixing of two streams

This is one more well known case, used for demonstrations of e.g. ANSYS Fluent or OpenFOAM in their training examples. In the OpenFOAM code they use this case to show how to convert Fluent mesh files to their polyMesh format. Here we will demonstrate mesh conversion from polyMesh format of OpenFoam to our native polyMesh format.

You can also unpack _zip_ folders with mesh, initial conditions and vtk templates files and run without mesh conversion.

NOTE: In the input file we specify mesh format either by _nativeMesh_ of _foamMesh_. See e.g. 'input-simple.nml' file.

### Mesh conversion to native freeCappuccino format

Perform following steps:

Run script:
```
./convertFOAMMesh
```

Edit the `polyMesh/boundary` file to look like this:
```
# bcType nFaces startFace
wall-4 wall 100 1300
inlet-5 inlet 8 1400
inlet-6 inlet 4 1408
press-outlet outlet 8 1412
wall-8 wall 34 1420
frontAndBackPlanes empty 1836 1454
```
Edit the `0/U` file to look like this:
```
internalField
  uniform
  0. 0. 0.
boundaryField
  Wall
    Dirichlet
      uniform
        0. 0. 0.
  inlet
    Dirichlet
    uniform
    1. 0. 0.
  inlet
    Dirichlet
    uniform
    0. 3. 0.
  outlet
    Neumann
    zeroGradient
  Wall
    Dirichlet
      uniform
        0. 0. 0.
  symmetry
    Neumann
      zeroGradient

```

NOTE: Make sure that mesh format is set to native in the configuration (input-simple.nml) file like this:
```
mesh_format = 'nativeMesh'
```
SPECIAL NOTE: Paraview is cell based and there is no reader for freeCappuccino. In future examples we will show how to use `cellConectivity` tool to create _polyMesh/cells_ file needed for writing VTK files. Until then we just use vtk template files from vtk folder. We get them by opening original OpenFOAM case in Paraview and saving the mesh as ASCII vtm file. It will be explained in a video.

Run the code by executing the run script:

`./run`

While running you may visualy inspect the residuals by calling:

`gnuplot plotResiduals`

or by openning 'monitor' file, where simulation progress is logged, eg by typing:

`tail -n 100 monitor`

When done, open the simulation files from the 'vtk' folder in Paraview.


Play with the case e.g. by changing boundary conditions from 'outlet' to 'pressure' in _polyMesh/boundary_ file. By doing this we set this boundary condition to the one where the value for the pressure is prescribed. If we haven't set field values for pressure _p_ using the _0/p_ file (look at the init.f90 file to see how it is done), it will assume that the Dirichlet boundary condition value for pressure set there is zero. 





