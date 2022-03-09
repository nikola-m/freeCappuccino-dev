### Turbulent pipe flow

This is the example in which we show how to perform LES of turbulent flow in a pipe.

At this stage you should have learned how to set up the case using the given files (gmsh .geo file for mesh, example _polyMesh/boundary_ file and input parameters in fortran namelist .nml file). This is a case with periodic boundaries in streamwise direction so have a look at provided 'boundary' file to see how periodic boundary is set up. At one side, e.g. outlet side, we set bc type as periodic and add info of starting face for twin boundary at the end of the line (see below - starting face index of the 'inlet' 799443 is copied to the end of the line of 'outlet' boundary). On the twin boundary bc type can be set to 'empty', it is not important if we have made periodic boundary coupling as we did in the step before.
```
# bcName bcType nFaces startFace
inlet empty 4205 799443
outlet periodic 4205 803648 799443
wall wall 7424 807853
```
Additional thing is to edit src/cappuccino/init.f90 and uncomment the lines which would make initial perturbed field for turbulent pipe:
```
    !** 
    do inp = 1,numCells
      call pipe_disturbances(xc(inp),yc(inp),zc(inp),u(inp),v(inp),w(inp))
    enddo
    !\***
```
Compile the _cappuccino_ code (takes around 20 sec.) and you are ready to run.
