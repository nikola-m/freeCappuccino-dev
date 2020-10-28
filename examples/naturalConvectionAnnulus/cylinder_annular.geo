// ---- 2D Circular Cylinder Gmsh Tutorial ----
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 11/26/2014 by Jacob Crabill
// Aerospace Computing Lab, Stanford University
//
// Changed 17/10/2020 by Nikola Mirkov
// for case of natural convection in annular region
// between two cylinders. 
// --------------------------------------------

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 6;
cl2 = .03;
cl3 = 10;

radius = 1;
outer = 3; //outer radius
numinner = 40; //number of mesh segments
extr = -0.1; // extrusion depth

// Circle & surrounding region
Point(5) = {0,   0, 0, cl2};
Point(6) = {0,  radius, 0, cl2};
Point(7) = {0, -radius, 0, cl2};
Point(8) = {0,  outer, 0, cl2};
Point(9) = {0, -outer, 0, cl2};
Point(10) = {radius,  0, 0, cl2};
Point(11) = {-radius, 0, 0, cl2};
Point(12) = {outer,  0, 0, cl2};
Point(13) = {-outer, 0, 0, cl2};

Circle(5) = {7, 5, 10};
Circle(6) = {6, 5, 11};
Circle(7) = {8, 5, 13};
Circle(8) = {9, 5, 12};

Line(9)  = {6, 8};
Line(10) = {7, 9};
Line(11)  = {10, 12};
Line(12) = {11, 13};

Circle(13) = {10, 5, 6};
Circle(14) = {11, 5, 7};
Circle(15) = {13, 5, 9};
Circle(16) = {12, 5, 8};

Transfinite Line {5,6,7,8,13,14,15,16} = 40; // We want 40 points along each of these lines
Transfinite Line {9,10,11,12} = numinner;    // And 10 points along each of these lines

// Each region which to be independently meshed must have a line loop
// Regions which will be meshed with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge
//Line Loop(1) = {1, 2, 3, 4, 7, 16, 8, 15}; // Exterior
Line Loop(2) = {10, 8, -11, -5}; // RH side of quad region - note ordering
Line Loop(3) = {7, -12, -6, 9}; // LH side of quad region - note ordering
Line Loop(4) = {-10, -14, 12, 15}; // RH side of quad region - note ordering
Line Loop(5) = {16, -9, -13, 11}; // LH side of quad region - note ordering

//Plane Surface(1) = {1}; // Outer unstructured region
Plane Surface(1) = {2}; // RH inner structured region
Plane Surface(2) = {3}; // LH inner structured region
Plane Surface(3) = {4}; // RH inner structured region
Plane Surface(4) = {5}; // LH inner structured region

// Mesh these surfaces in a structured manner
Transfinite Surface{1,2,3,4};

// Turn into quads (optional, but Transfinite Surface looks best with quads)
Recombine Surface {1,2,3,4};

// Change layer to increase z subdivision
Extrude {0, 0, extr} { Surface{1,2,3,4}; Layers{1}; Recombine;}

//+
Physical Surface("wallInner") = {55, 37, 73, 99};
//+
Physical Surface("front") = {2, 4, 1, 3};
//+
Physical Surface("back") = {60, 82, 104, 38};
//+
Physical Surface("wallOuter") = {47, 81, 29, 91};
//+
Physical Volume("fluid") = {2, 4, 3, 1};
