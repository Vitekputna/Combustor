// Gmsh project created on Thu Sep  8 20:47:28 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {3, 0, 0, 1.0};
//+
Point(4) = {10, 0, 0, 1.0};
//+
Point(5) = {10, 3, 0, 1.0};
//+
Point(6) = {3, 3, 0, 1.0};
//+
Point(7) = {3, 0.7, 0, 1.0};
//+
Point(8) = {1.5, 0.4, 0, 1.0};
//+
Point(9) = {1, 1, 0, 1.0};
//+
Point(10) = {0, 1, 0, 1.0};
//+
Point(11) = {1.6, 0.4, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 11};
//+
Line(8) = {11, 8};
//+
Line(9) = {8, 9};
//+
Line(10) = {9, 10};
//+
Line(11) = {10, 1};
//+
Curve Loop(1) = {6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inlet", 2) = {11};
//+
Physical Curve("supout", 3) = {4};
//+
Physical Curve("wall", 1) = {10, 9, 8, 7, 2, 1, 3, 5, 6};
