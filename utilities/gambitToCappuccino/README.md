# gambitToCappuccino
```
                       _     _ _   ___   _____                                 _              
                      | |   (_) | |__ \ / ____|                               (_)             
  __ _  __ _ _ __ ___ | |__  _| |_   ) | |     __ _ _ __  _ __  _   _  ___ ___ _ _ __   ___
 / _` |/ _` | '_ ` _ \| '_ \| | __| / /| |    / _` | '_ \| '_ \| | | |/ __/ __| | '_ \ / _ \
| (_| | (_| | | | | | | |_) | | |_ / /_| |___| (_| | |_) | |_) | |_| | (_| (__| | | | | (_) |
 \__, |\__,_|_| |_| |_|_.__/|_|\__|____|\_____\__,_| .__/| .__/ \__,_|\___\___|_|_| |_|\___/
  __/ |                                            | |   | |
 |___/                                             |_|   |_| 

```
Purpose:
--------------------

Converter for GAMBIT file format.

Description:
-----------------------

Reads mesh files in Gambit .neu format, either made in Gambit or from Gmsh (which can save meshes in .neu format) and converts them to format suitable for usage in freeCappuccino finite volume CFD code.

The main effort is in finding cell pairs that have a common face, which is data often missing from mesh formats inteded for FEM analysis. The Finite Volume Method requires indices of the owner and neighbour cell for each face, and the owner cell for every boundary face, as well as vertex indices for every face.

This code first decomposes given cells into faces that form it, than face indices are sorted as first step and then the faces are sorted using quicksort algorithm and grouping of faces with same indices in pairs is recognized as matching.

Cell connectivity in 'cells' file is not used in computation in freeCappucino code, but is used to write the paraview unstuctured .vtu postprocessing file which is cell based.

The code writes following files: 'points', 'cells', 'faces', 'owner', 'neighbour', 'boundary'

Run with 2cylinder.neu example.
```
make

./gambitToCappuccino

2cylinder
```
