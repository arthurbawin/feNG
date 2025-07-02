
p = 1;
lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {p, 0, 0, lc};
Point(3) = {p, p, 0, lc};
Point(4) = {0, p, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Point("PointPression")  = {1};
Physical Curve("Bord")  = {1,2,3,4};
Physical Surface("Domaine") = {1};

For num In {1:3}
Transfinite Curve{:} = 5;
Transfinite Surface{1};
