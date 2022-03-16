//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.05, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 0.2} {
  Curve{1}; 
}
//+
Curve Loop(3) = {3};
//+
Plane Surface(3) = {3};
//+
Surface Loop(1) = {2, 3, 1};
//+
Volume(1) = {1};
//+
Extrude {{0, 1, 0}, {0.15, 0, 0.2}, Pi/2} {
  Curve{3}; 
}
//+
Curve Loop(6) = {5};
//+
Plane Surface(6) = {6};
//+
Surface Loop(2) = {6, 5, 3};
//+
Volume(2) = {2};
//+
Extrude {0.2, 0, 0} {
  Curve{5}; 
}
//+
Curve Loop(9) = {7};
//+
Plane Surface(9) = {9};
//+
Surface Loop(3) = {9, 8, 6};
//+
Volume(3) = {3};
//+
Circle(8) = {0.2, 0, 0.3375, 0.0125, 0, 2*Pi};
//+
Rotate {{0, 1, 0}, {0.2, 0, 0.3375}, Pi/2} {
  Curve{8}; 
}
//+
Curve Loop(11) = {8};
//+
Plane Surface(11) = {11};
//+
Extrude {-0.275, 0, 0} {
  Curve{8}; 
}
//+
Curve Loop(13) = {10};
//+
Plane Surface(13) = {13};
//+
Surface Loop(4) = {13, 12, 11};
//+
Volume(4) = {4};
//+

//+
BooleanDifference{ Volume{4}; Delete; }{ Volume{2}; }
//+
Recursive Delete {
  Volume{5}; 
}
//+
BooleanDifference{ Surface{5}; Delete; }{ Surface{15}; Delete; }


Physical Surface("inlet-cold") = {1};
//+
Physical Surface("inlet-hot") = {13};
//+
Physical Surface("outlet") = {9};
//+
Physical Surface("wall") = {2, 5, 8, 14};
//+
Physical Volume("fluid") = {1, 2, 3, 4};
//+
//+
Field[1] = BoundaryLayer;
