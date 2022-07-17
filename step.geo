// Gmsh project created on Fri Jul 15 18:06:16 2022
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {3, -0.5, 0, 1.0};
//+
Point(4) = {1, -0.5, 0, 1.0};
//+
Point(5) = {3, 2, 0, 1.0};
//+
Point(6) = {1, 2, 0, 1.0};
//+
Point(7) = {0, 2, 0, 1.0};
//+

//+
Point(8) = {3, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 8};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 1};
//+
Line(9) = {2, 6};
//+
Line(10) = {2, 8};
//+
Curve Loop(1) = {7, 8, 1, 9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, -9, 10, 5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -4, -3, -2};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {3, 10, 9, 8, 7, 6, 5, 4, 2, 1} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Physical Curve("wall", 1) = {7, 6, 1, 2, 5, 4};
//+
Physical Curve("inlet", 2) = {8};
//+
Physical Curve("outlet", 3) = {3};
