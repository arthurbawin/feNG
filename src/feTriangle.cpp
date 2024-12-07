#include "feTriangle.h"
#include <iostream>

int Triangle::createBoundary(std::set<Edge, EdgeLessThan> &meshEdges, int numEdges)
{
  // std::pair<std::set<Edge, EdgeLessThan>::iterator, bool> ret;

  // for(int i = 0; i < 3; ++i)
  // {
  //   Vertex *v0 = _vlin[i];
  //   Vertex *v1 = _vlin[(i+1) % 3];

  //   Edge e(v0, v1, numEdges);
  //   ret = meshEdges.insert(e);
  //   if(ret.second) {
  //     // Edge was added to the set : nEdges is added to connecEdges
  //     _entities[p].connecEdges[3 * iElm + k] = nEdges++;

  //   } else {
  //     // Edge is already in the set : the negative is added to connecEdges
  //     // Assumes an edge is shared by only two triangles in 2D
  //     // More tests are required in 3D, where an edge can be shared by N tets
  //     _entities[p].connecEdges[3 * iElm + k] = -ret.first->getTag();
  //   }
  // }

  return numEdges;
}

void Triangle::getP1ReferenceCoordinates(double xyz[3], double uvw[3], Vertex **v)
{
  const double O[2] = {v[0]->x(), v[0]->y()};
  const double d[2] = {xyz[0] - O[0], xyz[1] - O[1]};
  const double d1[2] = {v[1]->x() - O[0], v[1]->y() - O[1]};
  const double d2[2] = {v[2]->x() - O[0], v[2]->y() - O[1]};
  const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  uvw[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
  uvw[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;
  uvw[2] = 0.;
}

// Simply call the generic function
bool TriangleP1::xyz2uvw(double xyz[3], double uvw[3], double tol)
{
  Triangle::getP1ReferenceCoordinates(xyz, uvw, _v);
  return true;
}

// Return true if xyz is inside P1 triangle.
// Compute barycentric coordinates and check the reference triangle.
bool TriangleP1::isInsidePhysical(double xyz[3], double tol)
{
  double uvw[3];
  this->xyz2uvw(xyz, uvw, tol);
  return Triangle::isInsideReference(uvw[0], uvw[1], uvw[2], tol);
}

double TriangleP1::sliverness()
{
  double a = sqrt( (_v[1]->x() - _v[0]->x())*(_v[1]->x() - _v[0]->x()) + (_v[1]->y() - _v[0]->y())*(_v[1]->y() - _v[0]->y()));
  double b = sqrt( (_v[2]->x() - _v[1]->x())*(_v[2]->x() - _v[1]->x()) + (_v[2]->y() - _v[1]->y())*(_v[2]->y() - _v[1]->y()));
  double c = sqrt( (_v[0]->x() - _v[2]->x())*(_v[0]->x() - _v[2]->x()) + (_v[0]->y() - _v[2]->y())*(_v[0]->y() - _v[2]->y()));
  // Cosine rule
  double tha = acos( (a*a + b*b - c*c) / (2.*a*b) );
  double thb = acos( (b*b + c*c - a*a) / (2.*b*c) );
  double thc = acos( (a*a + c*c - b*b) / (2.*a*c) );
  return fmax(1., tan( fmax(tha, fmax(thb, thc))/2.) );
}

// Return true if the point xyz is inside the convex region delimited by the parabola
// represented by its implicit form of parameters [A,B,C,D,E,F] such that 
// f(x,y) = A * x^2 + B * y^2 + C * x * y + D * x + E * y + F.
bool isInConvexRegionOfImplicitParabola(const double implicitParameters[6], double xyz[3])
{
  const double A = implicitParameters[0];
  const double B = implicitParameters[1];
  const double C = implicitParameters[2];
  const double D = implicitParameters[3];
  const double E = implicitParameters[4];
  const double F = implicitParameters[5];
  const double x = xyz[0];
  const double y = xyz[1];
  const double f = A * x*x + B * y*y + C*x*y + D*x + E*y + F;
  
  // See also matlab file parabolaOrientation.m for different cases.
  return (A > 0 || B > 0) ? (f < 0) : (f > 0);
}

// Return true if point xyz is inside the P2 triangle.
// Each parabolic edge is converted to its implicit form
// f(x,y) = A * x^2 + B * y^2 + C * x * y + D * x + E * y + F
// then we check if f > 0, < 0 or = 0.
bool TriangleP2::isInsidePhysical(double xyz[3], double tol)
{

  bool isInside = true;
  const double x = xyz[0];
  const double y = xyz[1];
  double uvwP1[3], xyzm[2];

  // Can also check first if the triangle is P1 and return the isInside for P1
  // bool isTriP1 = true;
  // for(int i = 0; i < 3; ++i)
  // {
  //   const double x0 = _v[i]->x();
  //   const double x1 = _v[(i + 1) % 3]->x();
  //   const double xm = _v[i + 3]->x();
  //   const double y0 = _v[i]->y();
  //   const double y1 = _v[(i + 1) % 3]->y();
  //   const double ym = _v[i + 3]->y();
  //   if( sqrt((xm - (x0+x1)/2.) * (xm - (x0+x1)/2.) 
  //          + (ym - (y0+y1)/2.) * (ym - (y0+y1)/2.)) > tol )
  //   {
  //     isTriP1 = false;
  //   }
  // }

  // if(isTriP1) {
  //   Triangle::getP1ReferenceCoordinates(xyz, uvwP1, _v);
  //   return Triangle::isInsideReference(uvwP1[0], uvwP1[1], uvwP1[2]);
  // }

  /////////////////////////////////////////
  // Debug plot
  int nt = 500;
  double t0 = -10.;
  double t1 =  10.;
  double dt = (t1-t0)/(double)nt;
  FILE *ff = fopen("isInsideP2Triangle.pos", "w");
  fprintf(ff, "View\"isInsideP2Triangle\"{\n");
  /////////////////////////////////////////

  for(int i = 0; i < 3; ++i)
  {
    const double x0 = _v[i]->x();
    const double y0 = _v[i]->y();
    const double x1 = _v[(i + 1) % 3]->x();
    const double y1 = _v[(i + 1) % 3]->y();
    const double xm = _v[i + 3]->x();
    const double ym = _v[i + 3]->y();

    // Unless I made a mistake, the conic does not degenerate to
    // a straight line when xm = (x0+x1)/2...
    // So test if xm is the midnode
    if( sqrt((xm - (x0+x1)/2.) * (xm - (x0+x1)/2.) 
           + (ym - (y0+y1)/2.) * (ym - (y0+y1)/2.)) < tol )
    {
      printf("Edge %d is straight\n", i);
      // Straight edge
      const double D = y0 - y1;
      const double E = x1 - x0;
      const double F = x0*y1 - x1*y0;
      const double f = D*x + E*y + F;

      // /////////////////////////////////////////
      // for(int ii = 0; ii < nt+1; ++ii) {
      //   double t = t0 + ii*dt;
      //   double xt = x0 * (1 - t)*(1-2*t) + x1*t*(2*t-1) + xm*4*t*(1-t);
      //   double yt = y0 * (1 - t)*(1-2*t) + y1*t*(2*t-1) + ym*4*t*(1-t);
      //   t = t0 + (ii+1)*dt;
      //   double xtp1 = x0 * (1 - t)*(1-2*t) + x1*t*(2*t-1) + xm*4*t*(1-t);
      //   double ytp1 = y0 * (1 - t)*(1-2*t) + y1*t*(2*t-1) + ym*4*t*(1-t);
      //   fprintf(ff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", xt, yt, xtp1, ytp1, 1, 1);
      // }

      // int M = 100;
      // double xmin = fmin(x0,x1);
      // double xmax = fmax(x0,x1);
      // double ymin = fmin(y0,y1);
      // double ymax = fmax(y0,y1);
      // double dx = (xmax-xmin)/(double)M;
      // double dy = (ymax-ymin)/(double)M;
      // for(int ii = 0; ii < M+1; ++ii) {
      //   double xloc = xmin + ii*dx;
      //   for(int jj = 0; jj < M+1; ++jj) {
      //     double yloc = ymin + jj*dy;
      //     double floc = D*xloc + E*yloc + F;
      //     bool locIsInside = f < -tol;
      //     if(locIsInside)
      //       fprintf(ff, "SP(%.16g,%.16g,0.){%u};\n", xloc, yloc, 1);
      //     else
      //       fprintf(ff, "SP(%.16g,%.16g,0.){%u};\n", xloc, yloc, 0);
      //   }
      // }
      // /////////////////////////////////////////

      if(f < -tol) return false;
    } else {
      
      // Curved edge.
      // The parabolic edge divides the plane in a concave and a convex region.
      // The region we want depends on if the midnode is inside or outside of the P1 triangle.
      // If the midnode is inside, the edge is concave (with respect to the triangle) 
      // and point xyz is inside the P2 triangle if it is in the concave region formed by the parabola.
      // If the midnode is outside, the edge is convex and we want the convex region.

      // First check if midpoint is inside or outside the P1 triangle formed by the vertices
      xyzm[0] = xm;
      xyzm[1] = ym;
      Triangle::getP1ReferenceCoordinates(xyzm, uvwP1, _v);
      bool isMidpointInside = Triangle::isInsideReference(uvwP1[0], uvwP1[1], uvwP1[2], tol);

      // Coefficients of the conic.
      const double A =
        (y0 * y0 + 2. * y0 * y1 - 4. * y0 * ym + y1 * y1 - 4. * y1 * ym + 4. * ym * ym);
      const double B =
        (x0 * x0 + 2. * x0 * x1 - 4. * x0 * xm + x1 * x1 - 4. * x1 * xm + 4. * xm * xm);
      // const double C = (4. * x0 * ym - 2. * x0 * y1 - 2. * x1 * y0 - 2. * x1 * y1 - 2. * x0 * y0 +
                        // 4. * xm * y0 + 4. * x1 * ym + 4. * xm * y1 - 8. * xm * ym);
      // Enforce C^2 - 4AB = 0 to avoid roundoff errors
      const double C = sqrt(4.*fabs(A*B));
      const double D = (x0 * y0 * y1 - x1 * y0 * y0 - 4. * x0 * ym * ym - xm * y0 * y0 -
                        4. * x1 * ym * ym - xm * y1 * y1 - x0 * y1 * y1 + x1 * y0 * y1 +
                        x0 * y0 * ym + 3. * x0 * y1 * ym + 3. * x1 * y0 * ym - 6. * xm * y0 * y1 +
                        x1 * y1 * ym + 4. * xm * y0 * ym + 4. * xm * y1 * ym);
      const double E = (x0 * x1 * y0 - x1 * x1 * y0 - x0 * x0 * ym - 4. * xm * xm * y0 -
                        x1 * x1 * ym - 4. * xm * xm * y1 - x0 * x0 * y1 + x0 * x1 * y1 +
                        x0 * xm * y0 - 6. * x0 * x1 * ym + 3. * x0 * xm * y1 + 3. * x1 * xm * y0 +
                        x1 * xm * y1 + 4. * x0 * xm * ym + 4. * x1 * xm * ym);
      const double F = x0 * x0 * y1 * ym - x0 * x1 * y0 * ym - x0 * x1 * y1 * ym +
                       4. * x0 * x1 * ym * ym - x0 * xm * y0 * y1 + x0 * xm * y1 * y1 -
                       4. * x0 * xm * y1 * ym + x1 * x1 * y0 * ym + x1 * xm * y0 * y0 -
                       x1 * xm * y0 * y1 - 4. * x1 * xm * y0 * ym + 4. * xm * xm * y0 * y1;

      const double implicitParameters[6] = {A,B,C,D,E,F};

      printf("Parabola passing through points\n");
      printf("x0 = [%+-1.6e %+-1.6e]\n", x0, y0);
      printf("x1 = [%+-1.6e %+-1.6e]\n", x1, y1);
      printf("xm = [%+-1.6e %+-1.6e]\n", xm, ym);

      printf("A = %+-1.6e;\n", A);
      printf("B = %+-1.6e;\n", B);
      printf("C = %+-1.6e;\n", C);
      printf("D = %+-1.6e;\n", D);
      printf("E = %+-1.6e;\n", E);
      printf("F = %+-1.6e;\n", F);

      // const double f = A * x * x + B * y * y + C * x * y + D * x + E * y + F;

      /////////////////////////////////////////
      for(int ii = 0; ii < nt+1; ++ii) {
        double t = t0 + ii*dt;
        double xt = x0 * (1 - t)*(1-2*t) + x1*t*(2*t-1) + xm*4*t*(1-t);
        double yt = y0 * (1 - t)*(1-2*t) + y1*t*(2*t-1) + ym*4*t*(1-t);
        t = t0 + (ii+1)*dt;
        double xtp1 = x0 * (1 - t)*(1-2*t) + x1*t*(2*t-1) + xm*4*t*(1-t);
        double ytp1 = y0 * (1 - t)*(1-2*t) + y1*t*(2*t-1) + ym*4*t*(1-t);
        fprintf(ff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", xt, yt, xtp1, ytp1, 1, 1);
      }

      int M = 100;
      double xmin = fmin(x0,x1);
      double xmax = fmax(x0,x1);
      double ymin = fmin(y0,y1);
      double ymax = fmax(y0,y1);
      double dx = (xmax-xmin)/(double)M;
      double dy = (ymax-ymin)/(double)M;
      for(int ii = 0; ii < M+1; ++ii) {
        double xloc = xmin + ii*dx;
        for(int jj = 0; jj < M+1; ++jj) {
          double yloc = ymin + jj*dy;  

          double xyzloc[3] = {xloc, yloc, 0.};
          bool inConvex = isInConvexRegionOfImplicitParabola(implicitParameters, xyzloc);

          bool locIsInside;
          if(isMidpointInside) {
            // Edge is concave and xyz is inside if it's in concave part.
            locIsInside = (!inConvex);
          } else {
            // Edge is convex and xyz is inside if it's in convex part.
            locIsInside = inConvex;
          }

          if(locIsInside)
            fprintf(ff, "SP(%.16g,%.16g,0.){%u};\n", xloc, yloc, 1);
          else
            fprintf(ff, "SP(%.16g,%.16g,0.){%u};\n", xloc, yloc, 0);

        }
      }
      /////////////////////////////////////////

      // if(isMidpointInside) {
      //   printf("Edge %d is curved and concave\n", i);
      //   // Edge is concave and xyz is inside if it's in concave part.
      //   // The concave part is the part where f < 0 if A*C > 0, or f > 0 otherwise.
      //   // isInside &= (A*C > 0) ? (f <= tol) : (f >= -tol);
      //   if((A*C > 0 && f > tol) || (A*C < 0 && f < -tol)) {
      //     fprintf(ff, "};\n"); fclose(ff);
      //     return false;
      //   }
      // } else {
      //   printf("Edge %d is curved and convex\n", i);
      //   // Edge is convex and xyz is inside if it's in convex part.
      //   // The convex part is the part where f > 0 if A*C > 0, or f < 0 otherwise.
      //   // isInside &= (A*C > 0) ? (f >= -tol) : (f <= tol);
      //   if((A*C > 0 && f < -tol) || (A*C < 0 && f > tol)) {
      //     fprintf(ff, "};\n"); fclose(ff);
      //     return false;
      //   }
      // }

      // if(isMidpointInside) {
      //   printf("Edge %d is curved and concave\n", i);
      //   // Edge is concave and xyz is inside if it's in concave part.
      //   // The concave part is the part where f < 0 if A*C > 0, or f > 0 otherwise.
      //   // isInside &= (A*C > 0) ? (f <= tol) : (f >= -tol);
      //   if((A*B > 0 && f > -tol) || (A*B < 0 && f < tol)) {
      //     // OK: point is inside
      //   } else {
      //     fprintf(ff, "};\n"); fclose(ff);
      //     return false;
      //   }
      // } else {
      //   printf("Edge %d is curved and convex\n", i);
      //   // Edge is convex and xyz is inside if it's in convex part.
      //   // The convex part is the part where f > 0 if A*B > 0, or f < 0 otherwise.
      //   // isInside &= (A*B > 0) ? (f >= -tol) : (f <= tol);
      //   if((A*B > 0 && f < tol) || (A*B < 0 && f > -tol)) {
      //     // OK: point is inside
      //   } else {
      //     fprintf(ff, "};\n"); fclose(ff);
      //     return false;
      //   }
      // }

      bool inConvex = isInConvexRegionOfImplicitParabola(implicitParameters, xyz);

      if(isMidpointInside) {
        // Edge is concave and xyz is inside if it's in concave part.
        isInside &= (!inConvex);
        if(isInside)
          printf("Edge %d is curved and concave - Point is inside w.r.t. this edge\n", i);
        else
          printf("Edge %d is curved and concave - Point is NOT inside w.r.t. this edge\n", i);
      } else {
        // Edge is convex and xyz is inside if it's in convex part.
        isInside &= inConvex;
        if(isInside)
          printf("Edge %d is curved and convex - Point is inside w.r.t. this edge\n", i);
        else
          printf("Edge %d is curved and convex - Point is NOT inside w.r.t. this edge\n", i);
      }
    }
  }

  fprintf(ff, "};\n"); fclose(ff);

  return isInside;
  // return true;
}

// This is defined in multiple places, this should be fixed
void phiTriP2(const double r, const double s, double phi[6], double dphidr[6], double dphids[6])
{
  phi[0] = (1. - r - s) * (1. - 2. * r - 2. * s);
  phi[1] = r * (2. * r - 1.);
  phi[2] = s * (2. * s - 1.);
  phi[3] = 4. * r * (1. - r - s);
  phi[4] = 4. * r * s;
  phi[5] = 4. * s * (1. - r - s);

  dphidr[0] = 4. * (r + s) - 3.;
  dphidr[1] = 4. * r - 1.;
  dphidr[2] = 0.;
  dphidr[3] = 4. * (1. - 2. * r - s);
  dphidr[4] = 4. * s;
  dphidr[5] = -4. * s;

  dphids[0] = 4. * (r + s) - 3.;
  dphids[1] = 0.;
  dphids[2] = 4. * s - 1.;
  dphids[3] = -4. * r;
  dphids[4] = 4. * r;
  dphids[5] = 4. * (1. - r - 2. * s);
}

//
// Newton-Raphson to find (r,s) s.t. F(r,s) = (x,y)
// If triangle is actually linear, return P1 localization
//
bool TriangleP2::xyz2uvw(double xyz[3], double uvw[3], double tol)
{
  // Use linear solution as initial Newton guess
  double r[2] = {0., 0.};
  // const double O[2] = {_v[0]->x(), _v[0]->y()};
  // const double d[2] = {xyz[0] - O[0], xyz[1] - O[1]};
  // const double d1[2] = {_v[1]->x() - O[0], _v[1]->y() - O[1]};
  // const double d2[2] = {_v[2]->x() - O[0], _v[2]->y() - O[1]};
  // const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  // r[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
  // r[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;

  double x0 = (*_v[0])(0), x1 = (*_v[1])(0), x2 = (*_v[2])(0);
  double y0 = (*_v[0])(1), y1 = (*_v[1])(1), y2 = (*_v[2])(1);
  double oneOverDet = 1./ ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));

  r[0] = oneOverDet * ((xyz[0] - x0)*(y2-y0) - (xyz[1] - y0)*(x2-x0));
  r[1] = oneOverDet * ((xyz[1] - y0)*(x1-x0) - (xyz[0] - x0)*(y1-y0));

  if(_isP2Straight) {
    uvw[0] = r[0];
    uvw[1] = r[1];
    uvw[2] = 0.0;
    return true;
  }

  double phi[6], dphidr[6], dphids[6];
  double J[4], invJ[4], f[2];
  int iter = 0, maxiter = 20;
  bool success = false, stop = false;

  // Check if xyz is close to one of the P2 nodes
  // If they are close up to the tolerance, simply return the reference coord of the node
  // If they are close up to a relaxed tolerance, use reference coord of the node as
  // initial candidate for Newton (NOT DONE YET, the relaxed tolerance should be a relative tolerance)
  // double relaxed_tol = 
  const double x = xyz[0];
  const double y = xyz[1];
  double u_ref[6] = {0., 1., 0., 0.5, 0.5, 0.};
  double v_ref[6] = {0., 0., 1., 0., 0.5, 0.5};
  for(int i = 0; i < 6; ++i) {
    if((_v[i]->x() - x)*(_v[i]->x() - x) + (_v[i]->y() - y)*(_v[i]->y() - y) < tol)
    { 
      r[0] = u_ref[i];
      r[1] = v_ref[i];
      break;
    }
  }

  while(!stop) {
    phiTriP2(r[0], r[1], phi, dphidr, dphids);

    // Tangent matrix and residual
    f[0] = f[1] = 0.;
    J[0] = J[1] = J[2] = J[3] = 0.;
    for(int i = 0; i < 6; ++i) {
      f[0] += _v[i]->x() * phi[i];
      f[1] += _v[i]->y() * phi[i];
      J[0] += _v[i]->x() * dphidr[i];
      J[1] += _v[i]->x() * dphids[i];
      J[2] += _v[i]->y() * dphidr[i];
      J[3] += _v[i]->y() * dphids[i];
    }
    f[0] -= xyz[0];
    f[1] -= xyz[1];

    // Inverse
    double det = J[0] * J[3] - J[1] * J[2];
    invJ[0] = J[3] / det;
    invJ[1] = -J[1] / det;
    invJ[2] = -J[2] / det;
    invJ[3] = J[0] / det;

    // Solve
    r[0] -= invJ[0] * f[0] + invJ[1] * f[1];
    r[1] -= invJ[2] * f[0] + invJ[3] * f[1];

    // printf("NR iter %2d: %+-1.6e - %+-1.6e ---- %+-1.6e - %+-1.6e\n", iter, r[0], r[1], f[0], f[1]);

    if(sqrt(f[0] * f[0] + f[1] * f[1]) < tol || iter > maxiter) {
      stop = true;
      if(sqrt(f[0] * f[0] + f[1] * f[1]) < tol) success = true;
    }
    iter++;
  }

  // bool success1 = success;
  // double result[2] = {r[0], r[1]};
  // double residual[2] = {f[0], f[1]};
  // bool inside = Triangle::isInsideReference(r[0], r[1], 0., 1e-6);

  // // Compare result with barycenter as initial candidate
  // r[0] = 1./3.;
  // r[1] = 1./3.;
  // iter = 0;
  // success = false, stop = false;

  // while(!stop) {
  //   phiTriP2(r[0], r[1], phi, dphidr, dphids);

  //   // Tangent matrix and residual
  //   f[0] = f[1] = 0.;
  //   J[0] = J[1] = J[2] = J[3] = 0.;
  //   for(int i = 0; i < 6; ++i) {
  //     f[0] += _v[i]->x() * phi[i];
  //     f[1] += _v[i]->y() * phi[i];
  //     J[0] += _v[i]->x() * dphidr[i];
  //     J[1] += _v[i]->x() * dphids[i];
  //     J[2] += _v[i]->y() * dphidr[i];
  //     J[3] += _v[i]->y() * dphids[i];
  //   }
  //   f[0] -= xyz[0];
  //   f[1] -= xyz[1];

  //   // Inverse
  //   double det = J[0] * J[3] - J[1] * J[2];
  //   invJ[0] = J[3] / det;
  //   invJ[1] = -J[1] / det;
  //   invJ[2] = -J[2] / det;
  //   invJ[3] = J[0] / det;

  //   // Solve
  //   r[0] -= invJ[0] * f[0] + invJ[1] * f[1];
  //   r[1] -= invJ[2] * f[0] + invJ[3] * f[1];

  //   // printf("NR iter %2d: %+-1.6e - %+-1.6e ---- %+-1.6e - %+-1.6e\n", iter, r[0], r[1], f[0], f[1]);

  //   if(sqrt(f[0] * f[0] + f[1] * f[1]) < tol || iter > maxiter) {
  //     stop = true;
  //     if(sqrt(f[0] * f[0] + f[1] * f[1]) < tol) success = true;
  //   }
  //   iter++;
  // }

  // bool success2 = success;
  // double result2[2] = {r[0], r[1]};
  // double residual2[2] = {f[0], f[1]};
  // bool inside2 = Triangle::isInsideReference(r[0], r[1], 0., 1e-6);

  // if( (inside || inside2) && (success1 || success2) && sqrt((result[0] - result2[0])*(result[0] - result2[0]) + (result[1] - result2[1])*(result[1] - result2[1])) > 1e-3 ) {
  //   printf("Converged to different roots:\n");
  //   printf("%+-1.4e - %+-1.4e with residual %+-1.4e - %+-1.4e\n", result[0], result[1], residual[0], residual[1]);
  //   printf("%+-1.4e - %+-1.4e with residual %+-1.4e - %+-1.4e\n", result2[0], result2[1], residual2[0], residual2[1]);
  //   exit(-1);
  // }

  // if(success)
  //   printf("NR success\n");
  // else
  //   printf("NR failed\n");

  uvw[0] = r[0];
  uvw[1] = r[1];
  uvw[2] = 0.0;

  return success;
}