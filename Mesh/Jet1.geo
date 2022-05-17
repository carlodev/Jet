// Gmsh project created on Thu May 12 08:14:59 2022
SetFactory("OpenCASCADE");
//+
D = DefineNumber[ 0.1, Name "Parameters/D" ];
//+
L = DefineNumber[ 4, Name "Parameters/L" ];
//+
H = DefineNumber[ 1, Name "Parameters/H" ];
//+
Point(1) = {0, 0.5*D, 0, 1.0};
//+
Point(2) = {0, -0.5*D, 0, 1.0};
//+
Point(3) = {L, 0.5*D, 0, 1.0};
//+
Point(4) = {L, -0.5*D, 0, 1.0};
//+
Point(5) = {0, H, 0, 1.0};
//+
Point(6) = {0, -H, 0, 1.0};
//+
Point(7) = {L, H, 0, 1.0};
//+
Point(8) = {L, -H, 0, 1.0};
//+
Point(9) = {0, 2*D, 0, 1.0};
//+
Point(10) = {0, -2*D, 0, 1.0};
//+
Point(11) = {L, 2*D, 0, 1.0};
//+
Point(12) = {L, -2*D, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {2, 4};
//+
Line(3) = {9, 11};
//+
Line(4) = {10, 12};
//+
Line(5) = {5, 7};
//+
Line(6) = {6, 8};
//+
Line(7) = {5, 9};
//+
Line(8) = {9, 1};
//+
Line(9) = {1, 2};
//+
Line(10) = {2, 10};
//+
Line(11) = {6, 10};
//+
Line(12) = {7, 11};
//+
Line(13) = {11, 3};
//+
Line(14) = {3, 4};
//+
Line(15) = {4, 12};
//+
Line(16) = {8, 12};
//+
Curve Loop(1) = {5, 12, -3, -7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 1, -13, -3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 14, -2, -9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, 4, -15, -2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {-16, -6, 11, 4};
//+
Plane Surface(5) = {5};
//+
Physical Curve("Walls", 17) = {7, 5, 8, 10, 11, 6};
//+
Physical Curve("Outlet", 18) = {12, 13, 14, 15, 16};
//+
Physical Curve("Inlet", 19) = {9};
//+
Physical Point("P_walls", 20) = {5, 9, 10, 6, 8, 7};
//+
Physical Point("P_inlet", 21) = {1, 2};
//+
Physical Point("P_outlet", 22) = {11, 3, 4, 12};
//+
Physical Surface("Fluid", 23) = {1, 2, 3, 4, 5};
//+
N = DefineNumber[ 100, Name "Parameters/N" ];
//+
Transfinite Curve {5, 3, 1, 2, 4, 6} = 2*N Using Progression 1.005;
//+
Transfinite Curve {7, 12, 11, 16} = N/4 Using Progression 0.95;
//+
Transfinite Curve {8, 10, 13, 15} = N/6 Using Progression 1;
//+
Transfinite Curve {9, 14} = N/6 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {4};
//+
Recombine Surface {5};
//+
Physical Surface("Initial_Velocity", 24) = {3};
