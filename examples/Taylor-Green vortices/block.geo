lt = 0.1;
pi = 3.141592653589793;
nx = 64; // number of cells in one dimension
ny = nx;
nz = nx;

Point(1) = {-pi,-pi,pi,lt};
Point(2) = {-pi,pi,pi,lt};
Point(3) = {-pi,pi,-pi,lt};
Point(4) = {-pi,-pi,-pi,lt};

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
Transfinite Curve {18, 17} = ny Using Progression 1.0;
//+
Transfinite Curve {13, 14} = nz Using Progression 1.0;
//+
Transfinite Surface {1};

Recombine Surface {1};
//+
Extrude {2*pi, 0, 0} {
  Surface{1}; Layers{nx}; Recombine;
}

//
// Don't let physical boundary names confuse you, 
// I used the .geo file of channel as a starting point.
//

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
