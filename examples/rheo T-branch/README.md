### Non-newtonian flow in T-branched duct

To run the code do the following:

Create the mesh with Gmsh (it is advised to install Gmsh https://gmsh.info/):

`gmsh t-branch.geo -3 -o t-branch.su2`

Note: We have saved mesh in su2 format - which is possible in Gmsh. We do that becasue of the conversion utility that we have, which is called in following step. If you haven't done it - you will need to compile _su2ToCappuccino_ which is located in 'utilities' folder of freeCappuccino. This program will create mesh files for our simulation:

`mkdir polyMesh vtk vtk/mesh 0`  
`su2ToCappuccino t-branch`  


Edit 'polyMesh/boundary' file to set boundary types to:

```
# bcName bcType nFaces startFace
Wall wall 11700 553500
InletLeft inlet 900 565200
InletBottom inlet 900 566100
OutletRight pressure 900 567000
top-bottom wall 12600 567900
```

Mesh conversion utility will also create template files for initial conditions in '0' folder. We will open files and edit only teh 'U' file:

In 'U' file, change settings to these:  
```
internalField
  uniform
    0.0 0.0 0.0
boundaryField
Sides
  Dirichlet
    uniform
      0.0 0.0 0.0
InletLeft
  Dirichlet
    uniform
      0.15 0.0 0.0
InletBottom
  Dirichlet
    uniform
      0.0 0.1 0.0
OutletRight
  Neumann
    zeroGradient
top-bottom
  Dirichlet
    uniform
      0.0 0.0 0.0 
```
Simulation configuration can be found in the 'input.nml' file. Lines related to non-Newtonian model are:

```
            calcVis = T           ! Activate dynamic viscosity recalculation (Temperature dependance, Non-Newtonian flows)
non_newtonian_model = 'BinghamPapanastasiou'    ! Self-explanatory
              Tau_0 = 0.5         ! Yield stress
               megp = 1.0         ! Exponential growth parameter
          muplastic = 0.000414    ! Plastic viscosity
             urfVis = 0.4         ! Under-relaxation parameter
```

Now we are ready to go!

Run the code by executing the run script:

`./run`

While running you may visualy inspect the residuals by calling:

`gnuplot plotResiduals`

or by openning 'monitor' file, where simulation progress is logged, eg by typing:

`tail -n 100 monitor`

When done, open the simulation files from the 'vtk' folder in Paraview.
