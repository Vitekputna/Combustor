// Gmsh project created on Tue Nov 29 11:12:42 2022
//+
Point(1) = {0, 0, 0, 1.0};
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

//+
Line(4) = {7, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 8};
//+
Line(7) = {8, 1};
//+


Point(9) = {3, 0.5, 0, 1.0};
//+
Point(10) = {5, 0.5, 0, 1.0};
//+
Line(8) = {3, 9};
//+
Line(9) = {9, 8};
//+
Line(10) = {6, 10};
//+
Line(11) = {10, 6};
//+
Line(12) = {10, 9};
//+
Line(13) = {9, 4};
//+
Line(14) = {10, 7};
//+
Curve Loop(1) = {9, -6, -5, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 8, 9, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, -12, -10, -2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, -13, -12, 14};
//+
Plane Surface(4) = {4};

//+
Physical Curve("inlet1", 1) = {7};
//+
Physical Curve("inlet2", 2) = {6};
//+
Physical Curve("wall", 3) = {5, 1, 2, 4};
//+
Physical Curve("outlet", 4) = {10,14};
//+
//Physical Surface("fluid1",0) = {1};
//+
//Physical Surface("fluid2", 1) = {2};
//+
//Transfinite Curve {6, 7} = 10 Using Progression 1;
//+
//Transfinite Curve {8} = 20 Using Progression 1;
//+
//Transfinite Curve {5, 1} = 50 Using Progression 1;
//+
//Transfinite Curve {4, 2} = 30 Using Progression 1;
//+
//Transfinite Curve {3} = 20 Using Progression 1;
//+
//Transfinite Surface {2};

//Recombine Surface {2};
//+

//+

//+
//Transfinite Curve {5, 9, 1, 4, 12, 2} = 50 Using Progression 1;
//+
//ansfinite Curve {6, 7, 8, 13, 14, 10} = 20 Using Progression 1;
//+
Physical Surface("fluid1",1) = {2};
//+
Physical Surface("fluid2",2) = {4};
//+
Physical Surface("fluid3",3) = {3};
Physical Surface("fluid4",4) = {1};
