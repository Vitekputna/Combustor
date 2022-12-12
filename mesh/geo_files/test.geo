// Gmsh project created on Tue Dec  6 11:13:35 2022
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left", 1) = {4};
//+
Physical Curve("right", 2) = {2};
//+
Physical Curve("top", 3) = {3};
//+
Physical Curve("floor", 4) = {1};
//+
Physical Surface("fluid", 1) = {1};
