lt = 0.1;
nx = 64;
ny = 64;
nz = 32;

Point(1) = {0,-1,3.14,lt};
Point(2) = {0,1,3.14,lt};
Point(3) = {0,1,0,lt};
Point(4) = {0,-1,0,lt};

// Horizontalne
Line(13) = {2,3};
Line(14) = {1,4};

// Vertikalne
Line(17) = {1,2};
Line(18) = {4,3};

//+
Curve Loop(1) = {18, -13, -17, 14};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {18, 17} = ny Using Bump 0.05;
//+
Transfinite Curve {13, 14} = nz Using Progression 1.0;
//+
Transfinite Surface {1};

Recombine Surface {1};
//+
Extrude {6.28, 0, 0} {
  Surface{1}; Layers{nx}; Recombine;
}
//+
Physical Surface("in") = {40};
//+
Physical Surface("out") = {1};
//+
Physical Surface("left") = {35};
//+
Physical Surface("right") = {27};
//+
Physical Surface("top") = {39};
//+
Physical Surface("bottom") = {31};
//+
Physical Volume("fluid") = {1};
