W = 0.1;	// Width
L = 1.0;	// Length
H = 1.0;	// Height
D = W;		// Depth

X1 = 0.;
X2 = L;
X3 = X2 + W;
X4 = X3 + L;

Y1 = 0.;
Y2 = H;
Y3 = Y2 + W;

// Mesh size and refinement parameter
nx1 = 61; nx2 = 31; nx3 = 61;
rx1 = 0.2; rx2 = 0.2; rx3 = 0.2;

ny1 = 61; ny2 = 31;
ry1 = 0.2; ry2 = 0.2;

nz = 30;

Z1 = 0.;
sg = 1.;

// Points
i = 1;
Point(i) = {X2, Y1, Z1, sg}; i = i + 1;
Point(i) = {X3, Y1, Z1, sg}; i = i + 1;

Point(i) = {X1, Y2, Z1, sg}; i = i + 1;
Point(i) = {X2, Y2, Z1, sg}; i = i + 1;
Point(i) = {X3, Y2, Z1, sg}; i = i + 1;
Point(i) = {X4, Y2, Z1, sg}; i = i + 1;

Point(i) = {X1, Y3, Z1, sg}; i = i + 1;
Point(i) = {X2, Y3, Z1, sg}; i = i + 1;
Point(i) = {X3, Y3, Z1, sg}; i = i + 1;
Point(i) = {X4, Y3, Z1, sg}; i = i + 1;

// Lines
i = 1;
Line(i) = {3,  4}; i = i + 1;
Line(i) = {4,  1}; i = i + 1;
Line(i) = {1,  2}; i = i + 1;
Line(i) = {2,  5}; i = i + 1;
Line(i) = {5,  6}; i = i + 1;
Line(i) = {6, 10}; i = i + 1;
Line(i) = {10, 9}; i = i + 1;
Line(i) = {9,  8}; i = i + 1;
Line(i) = {8,  7}; i = i + 1;
Line(i) = {7,  3}; i = i + 1;
Line(i) = {4,  8}; i = i + 1;
Line(i) = {5,  9}; i = i + 1;
Line(i) = {4,  5}; i = i + 1;

// Transfinite Lines
Transfinite Curve {2, 4}          = ny1 Using Bump ry1;
Transfinite Curve {10, 11, 12, 6} = ny2 Using Bump ry2;
Transfinite Curve {1, 9}          = nx1 Using Bump rx1;
Transfinite Curve {3, 13, 8}      = nx2 Using Bump rx2;
Transfinite Curve {5, 7}          = nx3 Using Bump rx3;

// Plane Surface
Curve Loop(1) = {3, 4, -13, 2}; Plane Surface(1) = {1};
Curve Loop(2) = {1, 11, 9, 10}; Plane Surface(2) = {2};
Curve Loop(3) = {13, 12, 8, -11}; Plane Surface(3) = {3};
Curve Loop(4) = {5, 6, 7, -12}; Plane Surface(4) = {4};

// Transfinite Surfaces
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

// Recombine Surfaces
Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};

// 3rd Direction
Extrude {0, 0, D} {
	Surface{2}; Surface{3}; Surface{4}; Surface{1};
	Layers{nz};
	Recombine;
}
  
//Boundary Conditions
Physical Surface("Wall") = {22, 100, 92, 66, 74, 52, 30};
Physical Surface("InletLeft") = {34};
Physical Surface("InletBottom") = {88};
Physical Surface("OutletRight") = {70};
Physical Surface("top-bottom") = {101, 57, 35, 79, 1, 3, 4, 2};
Physical Volume("Fluid") = {1, 2, 3, 4};

Mesh 3;

