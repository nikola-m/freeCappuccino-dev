// Gmsh project created on Sun Mar 14 19:58:23 2021
numpts = 51;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Line(1) = {3, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 3};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1, 2, 3, 4} = numpts;
//+
Transfinite Surface {1};
Recombine Surface {1};
//+
Extrude {0, 0, 0.1} {
  Surface{1}; Layers{1}; Recombine;
}
//+
Physical Surface("top") = {13};
//+
Physical Surface("left") = {25};
//+
Physical Surface("bottom") = {21};
//+
Physical Surface("right") = {17};
//+
Physical Surface("front-back") = {1, 26};
//+
Physical Volume("fluid") = {1};
