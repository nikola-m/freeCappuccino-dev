<pre>
             _____ _____     _____                                  _             
            / __  \_   _|   /  __ \                                (_)            
  ___ _   _ `' / /' | | ___ | /  \/ __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___  
 / __| | | |  / /   | |/ _ \| |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \ 
 \__ \ |_| |./ /___ | | (_) | \__/\ (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |
 |___/\__,_|\_____/ \_/\___/ \____/\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/ 
                                        | |   | |                                 
                                        |_|   |_|                                 `

</pre>

Conversion of SU2 meshes to polyMesh format used by freeCappuccino.

SU2 doesn't have info on cell neighbours trough a face, i.e. face_id-owner-neighbour relation. This is one onf the main tasks of this program.

What is good about su2 format is that it groups boundary regions by Physical surface name when making mesh in gmsh and exporting to SU2.

Following workflow is currently a preffered for freeCappuccino: make mesh in GMSH and export it into SU2 format, then do a mesh conversion to format used by freeCappuccino.

It may look cumbersome but it takes only few moments to do it.

First create su2ToCappuccino executable by typing:

`make`

Now the exectuable is in freeCappuccino/bin folder, which is preferrably on your PATH.

Go to the case directory where you will place your .su2 mesh file, then
create a few folders where mesh, initial condition, and Paraview files will be placed:

`mkdir polyMesh 0 vtk vtk/mesh`

These are all important for simulation case setup.

Then call the su2ToCappuccino mesh converter, providing the mesh file name (just name, without extension):

`su2ToCappuccino MESH_NAME_W/O_EXTENSION`

Finally edit 'polyMesh/boundary' file to set desired boundary conditions for specific boundary regions, as well as
files where initial conditions are set (files in 0/ folder).

You can inspect the mesh, including every boundary region separately, opening the mesh files from vtk/mesh folder in Paraview.




