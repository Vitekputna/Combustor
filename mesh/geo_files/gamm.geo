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
Point(9) = {1.5, -1.3, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 1};
//+
Circle(8) = {2, 9, 3};
//+
Curve Loop(1) = {7, 1, 8, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left", 1) = {7};
//+
Physical Curve("right", 2) = {3};
//+
Physical Curve("bottom", 3) = {1, 8, 2};
//+
Physical Curve("top", 4) = {6, 5, 4};
