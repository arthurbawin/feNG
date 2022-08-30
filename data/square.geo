
p = 2*Pi;
//p = 1;
lc = p;

Point(1) = {0, 0, 0, lc};
Point(2) = {p, 0, 0, lc};
Point(3) = {p, p, 0, lc};
Point(4) = {0, p, 0, lc};

//p = 4;
//lc = 2;

//Point(1) = {0, 0, 0, lc};
//Point(2) = {p, 0, 0, lc};
//Point(3) = {p, p, 0, lc};
//Point(4) = {0, p, 0, lc};

//Point(1) = {-p, -p, 0, lc};
//Point(2) = { p, -p, 0, lc};
//Point(3) = { p,  p, 0, lc};
//Point(4) = {-p,  p, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

//Physical Point("PointPression") = {1};
//Physical Curve("Haut")   = {3};
//Physical Curve("Gauche") = {4};
//Physical Curve("Angle")  = {1,2};
////Physical Curve("Bord")  = {1,2,3,4};
//Physical Surface("Surface") = {1};


Physical Curve("Dirichlet")  = {2,4};
Physical Curve("NeumannHaut") = {3};
Physical Curve("NeumannBas") = {1};
Physical Surface("Domaine") = {1};

//Mesh.ElementOrder = 2;