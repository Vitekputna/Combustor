// Gmsh project created on Mon Oct 17 22:55:34 2022
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {2, 0, 0, 1.0};
//+
Point(4) = {3, 0, 0, 1.0};
//+
Point(5) = {3, 1, 0, 1.0};
//+
Point(6) = {2, 1, 0, 1.0};
//+
Point(7) = {1, 1, 0, 1.0};
//+
Point(8) = {0, 1, 0, 1.0};
//+
Point(9) = {1.5, 2, 0, 1};
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
Line(6) = {7, 8};
//+
Line(7) = {8, 1};
//+
Circle(10) = {7, 9, 6};
//+
Physical Curve("inlet", 1) = {7};
//+
Physical Curve("outlet", 2) = {4};
//+
Physical Curve("wall", 3) = {6, 10, 5, 3, 2, 1};
//+
Transfinite Curve {7, 4, 1, 2, 3, 10, 5, 6} = 20 Using Progression 1;
//+
Curve Loop(1) = {7, 1, 2, 3, 4, 5, -10, 6};
//+
Plane Surface(1) = {1};
