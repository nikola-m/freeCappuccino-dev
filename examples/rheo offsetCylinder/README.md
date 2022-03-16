### Offset cylinder in a channel - non-Newtonian model

This case is made on the basis of a tutorial case of the OpenFOAM library (www.openfoam.com) for the purpose of cross validation (tongue in cheek because the rheological model used in the Cross model).

The original case can be found in following folder in OpenFOAM directory tree:
tutorials/incompressible/nonNewtonianIcoFoam/offsetCylinder/

Mesh is partially converted using 'readBoundaryFile' - which parses 'polyMesh/boundary'. You can look how the process of mesh conversion form OpenFOAM to freeCappuccino goes on a YouTube video:

https://www.youtube.com/watch?v=sYMfoQ61BcY

To run the case just unpack all zip archive files and execute the run script in the terminal.

Consult the source code, specifically rheology.f90 module to see the imlementation of various non-Newtonian models and the parameters they are dependent on, which can be set in the simulation configuration file (e.g. input.nml).


