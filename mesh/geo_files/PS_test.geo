//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {2, 0, 0, 1.0};
//+
Point(4) = {2, 1, 0, 1.0};
Point(5) = {1, 1, 0, 1.0};
Point(6) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 1};
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 5};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, -2, 5, 6};
//+
Plane Surface(2) = {2};
//+
Physical Curve("inlet",1) = {4};
//+
Physical Curve("outlet", 2) = {6};
//+
Physical Curve("wall", 3) = {3, 7, 1, 5};
//+
Physical Surface("fluid1",4) = {1};
Physical Surface("fluid2",5) = {2};

