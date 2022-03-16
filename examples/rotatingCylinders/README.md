### Taylor-Couette flow - validation case

The dimensions and parameters are taken from:
https://www.simscale.com/docs/validation-cases/rotating-zones-taylor-couette-flow/

On the same web page you may also find analytical solutions for the pressure and velocity which are used for validation.

The mesh is produced from Gmsh .geo file that is provided here (gmsh -3 cylinder_taylor_couette.geo -o cyl.su2), then converted to native freeCappuccino format using ```su2ToCappuccino``` . The polyMesh/boundary with correct boundary types as well as the initial/boundary conditions for velocity are also provided.

Note: the case is also interesting because of the rotating wall boundary condition. The velocities at the innner (moving) wall are first calculated in init.f90, then copied as 'nonuniform' Dirichlet boundary condition in 0/U file.
