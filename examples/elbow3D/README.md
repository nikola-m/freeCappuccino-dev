### Heat transfer in mixing elbow

If you have went through all the simpler tutorials cases, you should be familiar with the sequence of setting up the case
starting from Gmsh .geo file.

Open it in Gmsh, mesh geometry with desired resolution, and export mesh in .su2 format. After that convert mesh to native polyMesh format using ```su2ToCappuccino```, then edit polyMesh/boundary to set appropriate boundary conditions (there are two inlets, hot and cold one, one outlet and the wall boundary).

Finally setup initial/boundary conditions for necessary fields in 0 folder (also provided in 0.original.zip).

Alternative approach is to use provided mesh, Paraview template and initial conditions provided as .zip files.

The turbulence model used here is Realizable k-epsilon. You can of course experiment with other turbulence models.

Calculating the inlet conditions for turbulence scalars went like this:

Turbulence kinetic energy at hot inlet (turbulence intensty, I=5%)
k (Uin = 1.2m/s) = 3/2 * ( Uin * I )**2 = 0.0054

Turbulence kinetic energy at cold inlet (turbulence intensty, I=5%)
k (Uin = 0.4m/s) =  3/2 * ( Uin * I )**2 = 0.0006

Dissipation rate (hot inlet):
epsilon = cmu**0.75 k**1.5 / (0.07*D) = 0.03725

Dissipation rate (cold inlet):
epsilon = cmu**0.75 k**1.5 / (0.07*D) = 0.000345

For inquiry first run was a case without heat transfer - temperature equation was set to false in the inpt file. The mesh is coarse, tetrahedral, without boundary layer cells, so the first step was run with convection scheme for turbulence scalars set to first order upwind (TurbModel%Scalar(1)%gds = 0.0 and TurbModel%Scalar(2)%gds = 0.0 that is the contribution from deferred correction is set to zero). Second run was also without heat transfer - just set turbulence scalars to second order (TurbModel%Scalar(1)%gds = 1.0 and TurbModel%Scalar(2)%gds = 1.0). For these two runs you have pictures of the residuals. The residuals were quite smooth - probably because fo the linear solver used - CG with SAAMG algebraic multigrid preconditioner from the LIS library. Then I have turned on the turbulent heat transfer, using for Simple Gradient Diffusion Hypothesis ( lsgdh = T in the input file ) and after that Generalized Gradient Diffusion Hypothesis - GGDH ( lggdh = T in the input file ). The restarted run needs 'lread' set to true so that code reads last snapshot of values, saved in the restart file.

Run simulation as set in the run script - calling input-1.nml at first. This is simulation where turbulent fluxes are modelled using Simple Gradient Diffusion hypothesis and convection for turbulence scalars in set to first order upwind. Then run again with 'input-high-sgdh.nml' if you want results with second order upwind used for turbulence scalars. Then 'input-high-ggdh.nml' if you want to get results with Generalized Gradient Diffusion model for turbulent fluxes.
