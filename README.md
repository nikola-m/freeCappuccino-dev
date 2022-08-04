# freeCappuccino

### This is development version of freeCappuccino, a place for refactoring where a lot of unfinished code will be found.

freeCappuccino is a three-dimensional unstructured finite volume code for Computational Fluid Dynamics which comes in serial and parallel version.

Moreover, freeCappuccino is a fortran library for manipulation of discrete tensor fields, defined over polyhedral meshes.

Face based data structure enables use of cells of any shape.

### Name 'Cappuccino' encapsulates the idea that it is CAFFA (Computer Aided Fluid Flow Analysis) with some FOAM (Field Operation and Manipulation).

freeCappuccino is free both as a free coffee and as free speech.

<pre>
   __                _____                                 _               
  / _|              / ____|                               (_)              
 | |_ _ __ ___  ___| |     __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___    
 |  _| '__/ _ \/ _ \ |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \   
 | | | | |  __/  __/ |___| (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |  
 |_| |_|  \___|\___|\_____\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/   
                               | |   | |                                   
                               |_|   |_|                                   
 
 
             MM          MMM
            NMN.           MM
           MMM..   OMM.    MMM.
          .MMO.    MMM      MMM                         ...
          MMM:     MMM      MMM                       .. .?DMD7..
          MMM      MMM      MMM                      ..MMMMMMMMMMM~
          MMM:     .M       MMM.                    .NMMMMMMMMMMMMMM
          .MMM     .M.     .MMM.    ... .?MMMMM .. .MMMMMMMMMMMMMMMMM
           MMM.    .M .    MMM. ....MMMMMMMMMMMM$..MMMMMM?     .7MMMMM.
              MM.  M  M   MM  ..MMMMMMM~.. .MMMMM7MMMMMM.        7MMMMN
                            =MMMMMM7..  ..MMM.MMMMMMMMM+         .MMMMM.
                         DMMMMMD..  . . MMMM. .MMMMMMMM.        ..MMMMM.
                    ..=MMMMMZ.     ...MMMM      MMMMMMMI.        .MMMM$.
                    MMMMM8.   ..  .NMMMM..      :MMMMMMM         :MMMM.
                 :MMMMM.......  ~MMMMN           MMMMMMMMO.      MMMM+
             . ?MMMM:. ....  ,MMMMM              .MMMMMMMMM.    :MMMM.
           ..IMMMM  .  .  =MMMMMI                ..MMMMMMMM.    MMMM=
           .MMMM. .... DMMMMM?                     8MMMMM=     NMMMM
          +MMM.   .~MMMMMM~                        .MMMM.     ,MMMM.
         ~MM?.~MMMMMMM$                              MMMD    .MMMMM
        .DMMMMMMM$                                   MMMM.  .MMMMM
         .MMMM.                                      .MMMN..MMMMM,
         ..MMMM.                                     .MMMM.MMMMMM
           =MMMZ..   .=. ..         ..      =. ,:     .MMMMMMMMM
           .MMMM~..  7  .MM.M   M   MM.    M.  ..  M  .MMMMMMMM+
             MMMM....Z. .MM.M...M...MM..O:.  M . ~.M...ZMMMMMMM
              MMMM,. ., :  ,:   :  :  :    .,  ,,  ,  . MMMMMMM
               MMMM7. ..................................MMMMMM
                MMMM8 ..................................MMMMMM
                 NMMMM  ................................MMMMD
                  7MMMM ............................... MMMM
                    MMMMO............................. :MMM
                     MMMMM..............................MMM
                       MMMMM..........................MMMMM
                        ZMMMMD.......................MMMMM
                          MMMMMM.................. MMMMM
                            MMMMMM. ............OMMMMMZ
                              NMMMMM8.......=MMMMMMMI
                                =MMMMMMMMMMMMMMMMM
                                   NMMMMMMMMMM
</pre>


Introduction
------------------
The purpose of the code is to enable experimentation in discretizations and mathematical models. Therefore we have included many discretization options and we are in the process of implementation of many turbulence models and scalar equations for different variables.

The code is based on the reference:

N.Mirkov, B.Rašuo, S.Kenjereš, On the improved finite volume procedure for simulation of turbulent flows over real complex terrains, Journal of Computational Physics, Vol. 297 (2015), pp.18-45.

The code can read meshes in OpenFOAM® polyMesh format. Only slight changes are necessary in the 'boundary' file. 

Code also provides a utility for converting meshes made for SU2, in .su2 format (eg. made in Gmsh and exported as .su2) to native format of freeCappuccino.

This way we encourage collaboration among people using other open-source CFD codes!

Discretized equations are written in CSR (Compressed Sparse Row) format, which allows easy interface with many linear solver libraries, such as [LIS](http://www.ssisc.org/lis/). Also, one can use linear solvers provided in the code such as, Incomplete Cholesky Conjugate Gradient (ICCG), ILU(0) preconditioned Bi-Conjugate Gradient Stabilised method for non-symmetric systems, etc.

Requirements
-----------------
The code is written in modern fortran, so you will need fortran compiler (e.g. gfortran). Code uses _LAPACK_ in one routine, the default is to require _LAPACK_ but there is an option to build the code without it. Code may be built to use external library for solution of sparse linear systems [LIS](http://www.ssisc.org/lis/) or be used without it, just using built in linear solvers. See below the instruction on how to configure LIS.

Getting started
-----------------
Clone the git repository or download the zip archive of the most recent source, compile and build the code, and checkout the examples such as the lid-driven cavity, von Karman vortices behind a cylinder, Pitz-Daily backward facing step, MHD flow in a channel, and many more to come.

To compile the code just run _make_ from the root directory of freeCappuccino. The binary file _cappuccino_ will be created in _/bin_ directory. Don't forget to make _bin_ directory first. So the sequence is:
```
mkdir bin
make
```
It is not bad idea to update your $PATH environment variable to point to _bin/_ directory where all executables will be stored. In Linux one could add following line to the _.bashrc_ file
```
export PATH=$PATH:/path/to/where/it/is/freeCappuccino-dev/bin
```

To create some utility programs, go to specific folder and run 'make' there. For example
```
cd utilities/su2ToCappuccino/
make
```
The binary for 'su2ToCappuccino' will be created in 'freeCappuccino-dev/bin' folder.

If you want to use LIS library to access external linear algera solvers and preconditioners such as Smoothed Agglomeration Algebraic Multigrid, you need to download and install [LIS](http://www.ssisc.org/lis/) before building freeCappuccino. The process is described in LIS user's manual. To configure it the following switches are recommended, after which you proceed with building and installing the library:
```
./configure --enable-omp --enable-f90 --enable-saamg
make
make check
make install
```

Workflow
-----------------
One potential workflow would look like this:
1. __Create the mesh__ using [gmsh](https://gmsh.info/ "GMSH mesh generator"). 
* Note: rules are to define each boundary region as _Physical Surface_ with a name that will be the name of that boundary. Boundary type will be defined later. Also define _Physical Volume_ with some name eg. "fluid" to define cells of the interior.
* Generate the mesh and export it as _.su2_ file from the Gmsh export options. We find this format easier to use at the moment. 
* If you have _.geo_ file of Gmsh geometry, you may also generate and export mesh from command line with command:
`gmsh -3 <file_name>.geo -o <file_name>.su2`
2. __Convert the mesh__ to freeCappuccino type mesh format:
```
mkdir polyMesh 0 vtk vtk/mesh
su2ToCappuccino <mesh_file_without_su2_extension>
```
   * Note: besides creating mesh files in polyMesh directory, mesh converter does additional job, it creates Paraview files which will be used as a template to write output files and it creates template files for initial conditions in _0/_ directory to make your job of setting up the problem definition easier.

3. __Define boundary conditions__ (set types as eg. _inlet, outlet, wall, symmetry, pressure_ in _polyMesh/boundary_ file. 
4. __Define initial conditions__ for velocity in _0/U_ file and for other scalars, depending on the problem.
5. __Create configuration file__ which is a namelist file, maybe using a template _input_parameters_template_file.nml_ we have provided in the root directory as an example.
* Note: Most of the parameters have default values, so not all have to be present.
6. __Run the simulation__ calling the _cappuccino_ program, and giving as arguments the input (i.e. the configuration file), name of the file for monitoring the simulaton progress and name of the backup or restart file. Something like this:
```
cappuccino input.nml monitor <case_name>.rst
```
* Note: We assumed that _/bin_ folder is in your $PATH and the operating system will know where to find program executable if just its name is typed.
7. __Monitor the progress of the run__ by inspecting monitor file (we provide some scripts in examples to plot residuals) and when simulation is done, __view the results__ located in _vtk_ folder in Paraview.


__The alternative workflow__ would start from a given mesh in OpenFOAM polyMesh format. Then you can convert the mesh to native freeCappuccino format using _foamToCappuccino_ tool, or using OpenFOAM polyMesh format directly, just making some modifications. Most of these modifications can be done with provided toolset from the _utilities_ directory. 

An example of converting OpenFOAM mesh and preparing the simulation case can be seen in [![this video](https://img.youtube.com/vi/sYMfoQ61BcY/default.jpg)](https://www.youtube.com/watch?v=sYMfoQ61BcY)

License
------------------
The code is published under GNU General Public License v3.0.
