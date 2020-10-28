// Syntax.: Point(n) = {x, y, z, hs};
// As we'll use transfinite algorithm hs can be arbitrary
// or even omitted from mesh description
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {1, 1, 1, 1.0};
Point(6) = {1, 0, 1, 1.0};
Point(7) = {0, 1, 1, 1.0};
Point(8) = {0, 0, 1, 1.0};
//connect them with lines
Line(1) = {3, 7};
Line(2) = {7, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 1};
Line(6) = {2, 4};
Line(7) = {2, 6};
Line(8) = {6, 8};
Line(9) = {8, 1};
Line(10) = {1, 2};
Line(11) = {8, 7};
Line(12) = {6, 5};
//
Line Loop(13) = {7, 8, 9, 10};
Plane Surface(14) = {13};
Line Loop(15) = {6, 4, 5, 10};
Plane Surface(16) = {15};
Line Loop(17) = {3, 4, 1, 2};
Plane Surface(18) = {17};
Line Loop(19) = {12, -2, -11, -8};
Plane Surface(20) = {19};
Line Loop(21) = {7, 12, 3, -6};
Plane Surface(22) = {21};
Line Loop(23) = {9, -5, 1, -11};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 22, 20, 18, 16, 24};
//and volume
Volume(26) = {25};
//And the next lines will convert default generation strategy from tetrahedral meshes to hexagonal
//Transfinite Line "*" = 101 Using Bump 0.25;
//Transfinite Surface "*";
//Recombine Surface "*";
//Transfinite Volume "*";
//As the lengths of the cube sides are equal, densities of the 1D mesh on every side are equal 
//(also I've used keyword Bump to grade the mesh near the corners). To use generated mesh in OpenFOAM one has to define physical groups also:
Physical Surface("top") = {20};
Physical Surface("bottom") = {16};
Physical Surface("sides") = {14, 24, 18, 22};
Physical Volume("cube") = {26};

