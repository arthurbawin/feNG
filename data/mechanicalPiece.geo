
L = 6;
l = 2;
H = 2;
R = 1;
lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {l+R, H, 0, lc};
Point(5) = {l, H+R, 0, lc};
Point(6) = {l+R, H+2*R, 0, lc};
Point(7) = {L, H+2*R, 0, lc};
Point(8) = {L, 2*H+2*R, 0, lc};
Point(9) = {0, 2*H+2*R, 0, lc};
Point(10) = {0, H+R, 0, lc};
Point(11) = {l+R, H+R, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Circle(4) = {4, 11, 5};
Circle(5) = {5, 11, 6};

Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};

Line(11) = {10, 5};

Curve Loop(1) = {1, 2, 3, 4, -11, 10};
Curve Loop(2) = {11, 5, 6, 7, 8, 9};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Curve("fluxNul") = {1,3,4,5,6,8};
Physical Curve("fluxUnitaire") = {2,7};
Physical Curve("dirichlet") = {9,10};
Physical Surface("sBas") = {1};
Physical Surface("sHaut") = {2};