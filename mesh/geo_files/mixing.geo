// Gmsh project created on Tue Nov 29 11:12:42 2022
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {3, 0, 0, 1.0};
//+
Point(4) = {3, 1, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {5, 0.2, 0, 1.0};
//+
Point(7) = {5, 0.8, 0, 1.0};
//+
Point(8) = {0, 0.5, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 6};
//+
Line(3) = {6, 7};
//+
Line(4) = {7, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 8};
//+
Line(7) = {8, 1};
//+
Line(8) = {3, 4};
//+
Curve Loop(1) = {5, 6, 7, 1, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, -8, 2, 3};
//+
Plane Surface(2) = {2};
//+
Physical Curve("inlet1", 1) = {7};
//+
Physical Curve("inlet2", 2) = {6};
//+
Physical Curve("wall", 3) = {5, 1, 2, 4};
//+
Physical Curve("outlet", 4) = {3};
//+
//Physical Surface("fluid1",0) = {1};
//+
Physical Surface("fluid2", 1) = {2};
//+
Transfinite Curve {6, 7} = 10 Using Progression 1;
//+
Transfinite Curve {8} = 20 Using Progression 1;
//+
Transfinite Curve {5, 1} = 50 Using Progression 1;
//+
Transfinite Curve {4, 2} = 30 Using Progression 1;
//+
Transfinite Curve {3} = 20 Using Progression 1;
//+
//Transfinite Surface {2};

//Recombine Surface {2};
