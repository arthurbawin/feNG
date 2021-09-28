#include "feCurvedAdapt.h"

#include <iostream>
#include <poly34.h>
#include "SBoundingBox3d.h"
#include "STensor3.h"
#include "ellipseToolbox.h"
#include "fullMatrix.h"

#ifdef HAVE_GMSH
#include "gmsh.h"

static int U_ = 0;
static int U_X = 1;
static int U_Y = 2;
static int U_XX = 3;
static int U_XY = 4;
static int U_YX = 5;
static int U_YY = 6;
static int U_XXX = 7;
static int U_XXY = 8;
static int U_XYX = 9;
static int U_XYY = 10;
static int U_YXX = 11;
static int U_YXY = 12;
static int U_YYX = 13;
static int U_YYY = 14;

static int V_    = 15;
static int V_X   = 16;
static int V_Y   = 17;
static int V_XX  = 18;
static int V_XY  = 19;
static int V_YX  = 20;
static int V_YY  = 21;
static int V_XXX = 22;
static int V_XXY = 23;
static int V_XYX = 24;
static int V_XYY = 25;
static int V_YXX = 26;
static int V_YXY = 27;
static int V_YYX = 28;
static int V_YYY = 29;

#define Pi M_PI

static double x_0 = 0.5;
static double y_0 = 0.5;
static double aa = 50.0;
static double bb = 0.3;
static double cc = 3.0;

static int useAnalytical = 0;

std::string nameModel;

// Sign function (defined in ellipseToolbox.h)
// template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

static void PROBE(int TAG, double x, double y, double z, std::vector<double> &val,
                  const int step = -1, const int numComp = -1) {
  val.clear();
  gmsh::view::probe(TAG, x, y, z, val, step, numComp, 0, 1.e-3);
  if(val.empty()) val.push_back(0.);
}

static double fExact(double x, double y) {
  switch(useAnalytical){
    case 1 :
      return pow(x-x_0,4) * pow(y-y_0,4);
    case 2 :
    {
      double y_0 = 0.5;
      return (1. + tanh(aa*((y-y_0) - bb*sin(cc*M_PI*x))))/2.;
    }
    default :
      return pow(x-x_0,4) * pow(y-y_0,4);
  }  
}

static double f(double x, double y) {
  switch(useAnalytical){
    case 1 :
      return pow(x-x_0,4)*pow(y-y_0,4);
    case 2 :
      return (1. + tanh(aa*((y-y_0) - bb*sin(cc*M_PI*x))))/2.;
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      return u[0]*u[0]+v[0]*v[0];
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fx (double x, double y) {
  switch(useAnalytical){
    case 1 :
      return 4.0*pow(x-x_0,3)*pow(y-y_0,4);
    case 2 :
      return (aa*bb*cc*Pi*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.))/2.;
      case -1 :
      {
        std::vector<double> u,v;
        PROBE(U_,x,y,0.0,u);
        PROBE(V_,x,y,0.0,v);
        std::vector<double> ux, vx;
        PROBE(U_X,x,y,0.0,ux);
        PROBE(V_X,x,y,0.0,vx);
        return 2*u[0] * ux[0] + 2*v[0] * vx[0];
      }
    default :
    {
      std::vector<double> u;
      PROBE(U_X, x, y, 0.0, u);
      return u[0];
    }
  }  
}

static double fy (double x, double y) {
  switch(useAnalytical){
    case 1 :
      return 4.0*pow(x-x_0,4)*pow(y-y_0,3);
    case 2 :
      return -(aa*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.))/2.;
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      return 2*u[0] * uy[0] + 2*v[0] * vy[0];
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_Y, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fxx (double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 12.0*pow(x-x_0,2)*pow(y-y_0,4);
    case 2 :
      return - (aa*bb*cc*cc*Pi*Pi*sin(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.))/2. 
        - aa*aa*bb*bb*cc*cc*Pi*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*pow(cos(Pi*cc*x),2)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> ux, vx;
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uxx, vxx;
      PROBE(U_XX,x,y,0.0,uxx,0,1);
      PROBE(V_XX,x,y,0.0,vxx,0,1);
      return 2*uxx[0]*u[0] + 2 * ux[0]*ux[0] + 2*vxx[0]*v[0] + 2 * vx[0]*vx[0] ;
    }
    default : 
    {
      std::vector<double> u;
      PROBE(U_XX, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fyy (double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 12.0*pow(x-x_0,4)*pow(y-y_0,2);
    case 2 :
      return -aa*aa*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      std::vector<double> uyy, vyy;
      PROBE(U_YY,x,y,0.0,uyy,0,1);
      PROBE(V_YY,x,y,0.0,vyy,0,1);
      return 2*uyy[0]*u[0] + 2 * uy[0]*uy[0] + 2*vyy[0]*v[0] + 2 * vy[0]*vy[0];
    }
    default : 
    {
      std::vector<double> u;
      PROBE(U_YY, x, y, 0.0, u);
      return u[0];
    }
  }
  
}

static double fxy (double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 16.0*pow(x-x_0,3)*pow(y-y_0,3); 
    case 2 : 
      return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx;
      std::vector<double> uxy, vxy;
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      return 2*u[0]*UXY + 2*v[0]*VXY + 2* ux[0]*uy[0] + 2* vx[0]*vy[0];
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_XY, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fyx (double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 16.0*pow(x-x_0,3)*pow(y-y_0,3);
    case 2 :
      return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
      case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx;
      std::vector<double> uxy, vxy;
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      return 2*u[0]*UXY + 2*v[0]*VXY + 2* ux[0]*uy[0] + 2* vx[0]*vy[0];
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_YX, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fxxx(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 24.0*(x-x_0)*pow(y-y_0,4);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> ux, vx;
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uxx, vxx;
      PROBE(U_XX,x,y,0.0,uxx,0,1);
      PROBE(V_XX,x,y,0.0,vxx,0,1);
      std::vector<double> uxxx, vxxx;
      PROBE(U_XXX,x,y,0.0,uxxx);
      PROBE(V_XXX,x,y,0.0,vxxx);
      return 6*(uxx[0]*ux[0] +vxx[0]*vx[0]) + 2* (u[0]*uxxx[0]+v[0]*vxxx[0]);
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_XXX, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fxxy(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 48.0*pow(x-x_0,2)*pow(y-y_0,3);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx, uxx;
      std::vector<double> uxy, vxy, vxx;
      PROBE(U_XX,x,y,0.0,uxx);
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_XX,x,y,0.0,vxx);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      std::vector<double> uxxy, uyxx,uxyx ;
      std::vector<double> vxxy, vyxx,vxyx ;
      PROBE(U_XXY,x,y,0.0,uxxy);
      PROBE(U_YXX,x,y,0.0,uyxx);
      PROBE(U_XYX,x,y,0.0,uxyx);
      PROBE(V_XXY,x,y,0.0,vxxy);
      PROBE(V_YXX,x,y,0.0,vyxx);
      PROBE(V_XYX,x,y,0.0,vxyx);
      // gmsh::model::setCurrent("rect");
      double UXXY = (uxxy[0]+uyxx[0]+uxyx[0])/3.0;
      double VXXY = (vxxy[0]+vyxx[0]+vxyx[0])/3.0;
      return 4 * (UXY*ux[0] + VXY*vx[0]) + 2 * (uy[0]*uxx[0]+vy[0]*vxx[0]) +
        2 * (u[0]*UXXY+v[0]*VXXY);  
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_XXY, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fxyx(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 48.0*pow(x-x_0,2)*pow(y-y_0,3);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx, uxx;
      std::vector<double> uxy, vxy, vxx;
      PROBE(U_XX,x,y,0.0,uxx);
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_XX,x,y,0.0,vxx);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      std::vector<double> uxxy, uyxx,uxyx ;
      std::vector<double> vxxy, vyxx,vxyx ;
      PROBE(U_XXY,x,y,0.0,uxxy);
      PROBE(U_YXX,x,y,0.0,uyxx);
      PROBE(U_XYX,x,y,0.0,uxyx);
      PROBE(V_XXY,x,y,0.0,vxxy);
      PROBE(V_YXX,x,y,0.0,vyxx);
      PROBE(V_XYX,x,y,0.0,vxyx);
      // gmsh::model::setCurrent("rect");
      double UXXY = (uxxy[0]+uyxx[0]+uxyx[0])/3.0;
      double VXXY = (vxxy[0]+vyxx[0]+vxyx[0])/3.0;
      return 4 * (UXY*ux[0] + VXY*vx[0]) + 2 * (uy[0]*uxx[0]+vy[0]*vxx[0]) +
        2 * (u[0]*UXXY+v[0]*VXXY);  
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_XYX, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fxyy(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 48.0*pow(x-x_0,3)*pow(y-y_0,2);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx, uyy;
      std::vector<double> uxy, vxy, vyy;
      PROBE(U_YY,x,y,0.0,uyy);
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_YY,x,y,0.0,vyy);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      std::vector<double> uxyy, uyxy,uyyx ;
      std::vector<double> vxyy, vyxy,vyyx ;
      PROBE(U_YYX,x,y,0.0,uyyx);
      PROBE(U_XYY,x,y,0.0,uxyy);
      PROBE(U_YXY,x,y,0.0,uyxy);
      PROBE(V_YYX,x,y,0.0,vyyx);
      PROBE(V_XYY,x,y,0.0,vxyy);
      PROBE(V_YXY,x,y,0.0,vyxy);
      double UXYY = (uxyy[0]+uyxy[0]+uyyx[0])/3.0;
      double VXYY = (vxyy[0]+vyxy[0]+vyyx[0])/3.0;
      return 4 * (UXY*uy[0] + VXY*vy[0]) + 2 * (uy[0]*uyy[0]+vy[0]*vyy[0]) +
        2 * (u[0]*UXYY+v[0]*VXYY);     
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_XYY, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fyxx(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 48.0*pow(x-x_0,2)*pow(y-y_0,3);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx, uxx;
      std::vector<double> uxy, vxy, vxx;
      PROBE(U_XX,x,y,0.0,uxx);
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_XX,x,y,0.0,vxx);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      std::vector<double> uxxy, uyxx,uxyx ;
      std::vector<double> vxxy, vyxx,vxyx ;
      PROBE(U_XXY,x,y,0.0,uxxy);
      PROBE(U_YXX,x,y,0.0,uyxx);
      PROBE(U_XYX,x,y,0.0,uxyx);
      PROBE(V_XXY,x,y,0.0,vxxy);
      PROBE(V_YXX,x,y,0.0,vyxx);
      PROBE(V_XYX,x,y,0.0,vxyx);
      // gmsh::model::setCurrent("rect");
      double UXXY = (uxxy[0]+uyxx[0]+uxyx[0])/3.0;
      double VXXY = (vxxy[0]+vyxx[0]+vxyx[0])/3.0;
      return 4 * (UXY*ux[0] + VXY*vx[0]) + 2 * (uy[0]*uxx[0]+vy[0]*vxx[0]) +
        2 * (u[0]*UXXY+v[0]*VXXY);  
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_YXX, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fyxy(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 48.0*pow(x-x_0,3)*pow(y-y_0,2);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx, uyy;
      std::vector<double> uxy, vxy, vyy;
      PROBE(U_YY,x,y,0.0,uyy);
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_YY,x,y,0.0,vyy);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      std::vector<double> uxyy, uyxy,uyyx ;
      std::vector<double> vxyy, vyxy,vyyx ;
      PROBE(U_YYX,x,y,0.0,uyyx);
      PROBE(U_XYY,x,y,0.0,uxyy);
      PROBE(U_YXY,x,y,0.0,uyxy);
      PROBE(V_YYX,x,y,0.0,vyyx);
      PROBE(V_XYY,x,y,0.0,vxyy);
      PROBE(V_YXY,x,y,0.0,vyxy);
      double UXYY = (uxyy[0]+uyxy[0]+uyyx[0])/3.0;
      double VXYY = (vxyy[0]+vyxy[0]+vyyx[0])/3.0;
      return 4 * (UXY*uy[0] + VXY*vy[0]) + 2 * (uy[0]*uyy[0]+vy[0]*vyy[0]) +
        2 * (u[0]*UXYY+v[0]*VXYY);     
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_YXY, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fyyx(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 48.0*pow(x-x_0,3)*pow(y-y_0,2);
    case 2 :
      return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy, ux, vx;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      PROBE(U_X,x,y,0.0,ux);
      PROBE(V_X,x,y,0.0,vx);
      std::vector<double> uyx, vyx, uyy;
      std::vector<double> uxy, vxy, vyy;
      PROBE(U_YY,x,y,0.0,uyy);
      PROBE(U_XY,x,y,0.0,uxy);
      PROBE(U_YX,x,y,0.0,uyx);
      PROBE(V_YY,x,y,0.0,vyy);
      PROBE(V_XY,x,y,0.0,vxy);
      PROBE(V_YX,x,y,0.0,vyx);
      double UXY = 0.5*(uxy[0]+uyx[0]);
      double VXY = 0.5*(vxy[0]+vyx[0]);
      std::vector<double> uxyy, uyxy,uyyx ;
      std::vector<double> vxyy, vyxy,vyyx ;
      PROBE(U_YYX,x,y,0.0,uyyx);
      PROBE(U_XYY,x,y,0.0,uxyy);
      PROBE(U_YXY,x,y,0.0,uyxy);
      PROBE(V_YYX,x,y,0.0,vyyx);
      PROBE(V_XYY,x,y,0.0,vxyy);
      PROBE(V_YXY,x,y,0.0,vyxy);
      double UXYY = (uxyy[0]+uyxy[0]+uyyx[0])/3.0;
      double VXYY = (vxyy[0]+vyxy[0]+vyyx[0])/3.0;
      return 4 * (UXY*uy[0] + VXY*vy[0]) + 2 * (uy[0]*uyy[0]+vy[0]*vyy[0]) +
        2 * (u[0]*UXYY+v[0]*VXYY);     
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_YYX, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double fyyy(double x, double y) {
  switch(useAnalytical){
    case 1 : 
      return 24.0*pow(x-x_0,4)*(y-y_0);
    case 2 :
      // return aa*aa*bb*cc*Pi*tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.))*cos(Pi*cc*x)*(pow(tanh(aa*(bb*sin(Pi*cc*x) - y + 1./2.)),2) - 1.);
    case -1 :
    {
      std::vector<double> u,v;
      PROBE(U_,x,y,0.0,u);
      PROBE(V_,x,y,0.0,v);
      std::vector<double> uy, vy;
      PROBE(U_Y,x,y,0.0,uy);
      PROBE(V_Y,x,y,0.0,vy);
      std::vector<double> uyy, vyy;
      PROBE(U_YY,x,y,0.0,uyy,0,1);
      PROBE(V_YY,x,y,0.0,vyy,0,1);
      std::vector<double> uyyy, vyyy;
      PROBE(U_YYY,x,y,0.0,uyyy);
      PROBE(V_YYY,x,y,0.0,vyyy);
      return 6*(uyy[0]*uy[0] +vyy[0]*vy[0]) + 2* (u[0]*uyyy[0]+v[0]*vyyy[0]); ;
    }
    default :
    {
      std::vector<double> u;
      PROBE(U_YYY, x, y, 0.0, u);
      return u[0];
    }
  }
}

static double dtt (const double x, const double y, double C, double S){
  const double c11 = fxx(x,y);
  const double c12 = (fxy(x,y) + fyx(x,y))/2.;
  const double c22 = fyy(x,y);
  return C*C*c11 + S*S*c22+ 2.*C*S*c12;
}

static double dttt(const double x, const double y, double C, double S) {
  const double c111 = fxxx(x, y);
  const double c222 = fyyy(x, y);
  const double c112 = (fxxy(x, y) + fxyx(x,y) + fyxx(x,y)) /3.;
  const double c122 = (fxyy(x, y) + fyyx(x,y) + fyxy(x,y)) /3.;
  return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122;
}

static void directionsFromGradient (double x, double y, FILE *F, double &C, double &S){
  
  const double a = fx(x,y);
  const double b = fy(x,y);

  double normGrad = sqrt(a*a + b*b);

  if(normGrad > 1e-8){
    double theta1 = atan2(b,a);
    C = cos(theta1);
    S = sin(theta1);
    fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,C,S,0.);
    // fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,-S,C,0.);
  } else{
    // Gradient norm is too small : directions will be smoothed
    C = 0.;
    S = 0.;
  }
}

static void directionsFromHessian (double x, double y, FILE *F, double &C, double &S){
  
  const double a = fxx(x,y);
  const double b = (fxy(x,y) + fyx(x,y))/2.;
  const double c = fyy(x,y);

  // Eigenvalues and spectral radius of the hessian
  double l1 = (a+c-sqrt(a*a-2.*a*c+4.*b*b+c*c))/2.;
  double l2 = (a+c+sqrt(a*a-2.*a*c+4.*b*b+c*c))/2.;

  double spectralRadius = fmax(fabs(l1),fabs(l2));

  if(spectralRadius > 1e-8){
    // First eigenvector
    double v11 = a - c + sqrt(a*a + 4.*b*b + c*c - 2*a*c);
    double v12 = 2.*b;
    C = v11;
    S = v12;
    const double L = fmax(1e-10, sqrt(C * C + S * S));
    C /= L;
    S /= L;
    fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,C,S,0.);
    // fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,-S,C,0.);
  } else{
    // Hessian's eigenvalues are too small : directions will be smoothed
    C = 0.;
    S = 0.;
  }
}

static void computeMetricP1Hess (double x, double y, double lMin, double lMax, double eps,
         double &g00, double &g01, double &g11, FILE *F, double &C, double &S){
  
  const double a = fxx(x,y);
  const double b = (fxy(x,y) + fyx(x,y))/2.;
  const double c = fyy(x,y);

  // printf("c11 = %4.4f - c22 = %4.4f - c12 = %4.4f\n", a, c, b);

  double v11 = a - c + sqrt(a*a + 4.*b*b + c*c - 2*a*c);
  // v11 /= 2.*b;
  double v12 = 2.*b;

  // const double theta = atan2(v12,v11);
  // C = cos(t2/2);
  // S = sin(t2/2);

  C = v11;
  S = v12;

  const double L = fmax(1e-8, sqrt(C * C + S * S));
  C /= L;
  S /= L;

  if(fabs(sqrt(C * C + S * S) - 1.0) > 1e-4){
    C = 1.;
    S = 0.;
  }

  fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,C,S,0.);
  // fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,-S,C,0.);
}

static void computeMetricP2Hess (double x, double y, double lMin, double lMax, double eps,
         double &g00, double &g01, double &g11, FILE *F, double &C, double &S){
  
  const double c11 = fxx(x,y);
  const double c22 = fyy(x,y);
  const double c12 = fxy(x,y);

  // printf("c11 = %4.4f - c22 = %4.4f - c12 = %4.4f\n", c11, c22, c12);

  const double t2 = atan2(2*c12,c11-c22);
  C = cos(t2/2);
  S = sin(t2/2);
  // return;
  //  computeMetricMaxRectangle (x, y, S, C, lMin, lMax, eps, g00, g01, g11);

  //  return;
  //  printf("%12.5E %12.5E %12.5E %12.5E %12.5E\n",c11,c12,c22,C,S);
  
  // const double dttt0_ = fabs(dttt (x,y, C, S));
  // const double dttt1_ = fabs(dttt (x,y, -S, C));

  const double dtt0_ = fabs(dtt (x,y, C, S));
  const double dtt1_ = fabs(dtt (x,y, -S, C));

  // printf("dtt1 = %4.4f - dtt0 = %4.4f\n", dtt1_, dtt0_);

  // double l0 = dttt0_ != 0. ? pow (6*eps/dttt0_,0.3333) : lMax;
  // double l1 = dttt1_ != 0. ? pow (6*eps/dttt1_,0.3333) : lMax;
  double l0 = dtt0_ != 0. ? pow (2*eps/dtt0_,0.5) : lMax;
  double l1 = dtt1_ != 0. ? pow (2*eps/dtt1_,0.5) : lMax;
  
  l0 = std::min (l0,lMax);
  l0 = std::max (l0,lMin);
  l1 = std::min (l1,lMax);
  l1 = std::max (l1,lMin);
    
  double h0 = 1./(l0*l0);
  double h1 = 1./(l1*l1);
  
  fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,C,S,0.);
  // fprintf(F,"VP(%g,%g,%g){%g,%g,%g};\n",x,y,0.,-S,C,0.);

  g00 = C*C*h0+S*S*h1;
  g11 = S*S*h0+C*C*h1;
  g01 = S*C*(h1-h0);  
    
}

void feCurvedAdapt::computeMetricP2(double x, double y, double lMin, double lMax, double eps,
                                    double &g00, double &g01, double &g11, FILE *F, double &C,
                                    double &S) {

  // const double a = fxx(x,y);
  // const double b = (fxy(x,y) + fyx(x,y))/2.;
  // const double cc = fyy(x,y);

  // // printf("c11 = %4.4f - c22 = %4.4f - c12 = %4.4f\n", a, c, b);

  // double v11 = a - cc + sqrt(a*a + 4.*b*b + cc*cc - 2*a*cc);
  // double v12 = 2.*b;

  // C = v11;
  // S = v12;

  // const double LL = fmax(1e-8, sqrt(C * C + S * S));
  // C /= LL;
  // S /= LL;

  // if(fabs(sqrt(C * C + S * S) - 1.0) > 1e-4){
  //   C = 1.;
  //   S = 0.;
  // }

  // return; // Compute only the directions of the hessian

  double c111 = fxxx(x, y);
  const double c222 = fyyy(x, y);
  const double c112 = (fxxy(x, y) + fxyx(x,y) + fyxx(x,y)) /3;
  const double c122 = (fxyy(x, y) + fyyx(x,y) + fyxy(x,y)) /3;

  // std::cout<<x<<" "<<y<<" : "<<c111<<" "<<c222<<" "<<c112<<" "<<c122<<std::endl;

  // if(c111 == 0) {
  //   printf("coucou\n");
  //   c111 = .0001;
  // }

  if(fabs(c111) < 1e-10){
    c111 = sgn(c111) * 1e-4;
  }

  double c[3];
  int nRoots = SolveP3(c, c112 / c111, c122 / c111, c222 / c111);

  int myRoot = 0;

  if(nRoots == 3) {
    double maxdttt_ = 0;
    for(int i = 0; i < 3; i++) {
      C = c[i];
      S = 1.;
      double L = fmax(1e-8, sqrt(C * C + S * S));
      C /= L;
      S /= L;
      double dttt_ = fabs(dttt(x, y, -S, C));
      if(dttt_ > maxdttt_) {
        myRoot = i;
        maxdttt_ = dttt_;
      }
    }
  }

  C = c[myRoot];
  S = 1.;

  // const double L = sqrt(C * C + S * S);
  const double L = fmax(1e-8, sqrt(C * C + S * S));
  C /= L;
  S /= L;

  if(fabs(sqrt(C * C + S * S) - 1.0) > 1e-4){
    C = 1.;
    S = 0.;
  }
  if(isnan(C) || isnan(S)){
    std::cout<<"Nan alert"<<std::endl;
    std::cout<<x<<" "<<y<<" : "<<c111<<" "<<c222<<" "<<c112<<" "<<c122<<std::endl;
    C = 1.;
    S = 0.;
  }
  return;
}

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S) {
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > nodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);
  std::multimap<size_t, size_t> graph;
  std::set<size_t> nodes_to_treat;
  std::set<size_t> nodes;

  std::cout<<elementTags[0].size()<<" elements"<<std::endl;

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    size_t n0 = nodeTags[0][3 * i + 0];
    size_t n1 = nodeTags[0][3 * i + 1];
    size_t n2 = nodeTags[0][3 * i + 2];
    graph.insert(std::make_pair(n0, n1));
    graph.insert(std::make_pair(n1, n2));
    graph.insert(std::make_pair(n2, n0));
    nodes.insert(n0);
    nodes.insert(n1);
    nodes.insert(n2);
  }

  double threshold = 0.1;

  for(auto n : nodes) {
    double c = C[n];
    double s = S[n];
    for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
      size_t neigh = it->second;
      double cn = C[neigh];
      double sn = S[neigh];
      double pv = fabs(cn * c + sn * s);
      // if(pv > threshold && pv < 1. - threshold) nodes_to_treat.insert(n);
      if(pv < threshold) nodes_to_treat.insert(n);
    }
  }
  int iter = 0;
  while(iter++ < 100) {
    for(auto n : nodes_to_treat) {
      double c4 = 0;
      double s4 = 0;
      for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
        size_t neigh = it->second;
        double cn = C[neigh];
        double sn = S[neigh];
        double theta = atan2(sn, cn);
        c4 += cos(4. * theta);
        s4 += sin(4. * theta);
      }
      double theta = 0.25 * atan2(s4, c4);
      C[n] = cos(theta);
      S[n] = sin(theta);
    }
  }

  FILE *f = fopen("dirs.pos", "w");
  fprintf(f, "View\"Dirs\"{\n");

  for(auto n : nodes) {
    double c = C[n];
    double s = S[n];
    if(nodes_to_treat.find(n) != nodes_to_treat.end()) {
      c *= 2;
      s *= 2;
    }
    std::vector<double> coord;
    std::vector<double> par;
    int entityDim, entityTag;
    gmsh::model::mesh::getNode(n, coord, par, entityDim, entityTag);

    fprintf(f, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1], c, s);
  }

  fprintf(f, "};\n");
  fclose(f);
}

static void metricScaling(feRecovery *recovery, feMetricOptions metricOptions, std::vector<SMetric3> &metrics) {
  std::vector<double> &w = recovery->_geoSpace->getQuadratureWeights();
  std::vector<double> &J = recovery->_cnc->getJacobians();

  int nQuad = w.size();
  double N = (double) metricOptions.nTargetVertices;
  double p = metricOptions.LpNorm;
  double deg = (double) recovery->getDegreeSolution();
  double dim = (double) recovery->getDim();
  double exponent = p * (deg + 1.0) / (2.0 * (p * (deg + 1.0) + dim));

  // std::cout<<"exponent = "<<exponent<<std::endl;

  double I = 0.0;
  for(int iElm = 0; iElm < recovery->_nElm; ++iElm) {
    for(int k = 0; k < nQuad; ++k) {
      // Interpolate det(Q) at quad nodes
      for(int iNode = 0; iNode < recovery->_nNodePerElm; ++iNode) {
        int v = recovery->_cnc->getNodeConnectivity(iElm, iNode);
        I += J[nQuad * iElm + k] * w[k] * recovery->_geoSpace->getFunctionAtQuadNode(iNode, k) *
             pow(metrics[v].determinant(), exponent);
      }
    }
  }

  // std::cout<<"Int det(Q) = "<<I<<std::endl;

  std::vector<int> &vertices = recovery->getVertices();

  double hMin = metricOptions.hMin, lMax = 1.0 / (hMin * hMin);
  double hMax = metricOptions.hMax, lMin = 1.0 / (hMax * hMax);

  fullMatrix<double> V(3, 3);
  fullVector<double> S(3);
  SVector3 v0, v1, v2;
  for(auto v : vertices) {
    SMetric3 M = metrics[v];

    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), -1.0 / (p * (deg + 1.0) + dim));

    M *= factor;

    M.eig(V, S, false);
    v0 = SVector3(V(0, 0), V(0, 1), V(0, 2));
    v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
    v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));
    M = SMetric3(fmin(lMax, fmax(lMin, S(0))), fmin(lMax, fmax(lMin, S(1))),
                 fmin(lMax, fmax(lMin, S(2))), v0, v1, v2);

    metrics[v] = M;
  }
}

double ERROR_SQUARED_P1(double *xa, double *xb, double *xc, feFunction *solExact) {
  int triangleP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  std::vector<double> localCoord;
  std::vector<double> weights;
  gmsh::model::mesh::getIntegrationPoints(triangleP1, "Gauss20", localCoord, weights);
  std::vector<double> basisFunctions;
  int numComponents, numOrientations;
  gmsh::model::mesh::getBasisFunctions(triangleP1, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  std::vector<double> gradBasisFunctions;
  gmsh::model::mesh::getBasisFunctions(triangleP1, localCoord, "GradLagrange", numComponents,
                                       gradBasisFunctions, numOrientations);

  double F[3] = {f(xa[0], xa[1]),   f(xb[0], xb[1]),   f(xc[0], xc[1])};
  double X[3] = {xa[0], xb[0], xc[0]};
  double Y[3] = {xa[1], xb[1], xc[1]};
  double e2 = 0;
  for(size_t i = 0; i < weights.size(); i++) {
    double interpolated = 0, x = 0, y = 0, dxdu = 0, dxdv = 0, dydu = 0, dydv = 0;

    for(size_t j = 0; j < 3; j++) {
      x += basisFunctions[3 * i + j] * X[j];
      y += basisFunctions[3 * i + j] * Y[j];
      dxdu += gradBasisFunctions[3 * (3 * i + j) + 0] * X[j];
      dxdv += gradBasisFunctions[3 * (3 * i + j) + 1] * X[j];
      dydu += gradBasisFunctions[3 * (3 * i + j) + 0] * Y[j];
      dydv += gradBasisFunctions[3 * (3 * i + j) + 1] * Y[j];
      interpolated += basisFunctions[3 * i + j] * F[j];
    }
    double detJ = fabs(dxdu * dydv - dxdv * dydu);
    std::vector<double> pos(3,0.);
    pos[0] = x;
    pos[1] = y;
    double exact = (solExact != nullptr) ? solExact->eval(0, pos) : fExact(x, y);
    double diff = exact - interpolated;
    e2 += weights[i] * diff * diff * detJ;
    // printf("e = %10.14f \t weights[%d] = %4.4f \t diff = %4.4f \t detJ = %4.4f\n", e2, i, weights[i], diff, detJ);
  }
  return e2;
}

double ERROR_SQUARED_P2(double *xa, double *xb, double *xc, double *xab, double *xbc, double *xca, feFunction *solExact) {
  int triangleP2 = gmsh::model::mesh::getElementType("Triangle", 2);
  std::vector<double> localCoord;
  std::vector<double> weights;
  gmsh::model::mesh::getIntegrationPoints(triangleP2, "Gauss12", localCoord, weights);
  std::vector<double> basisFunctions;
  int numComponents, numOrientations;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  std::vector<double> gradBasisFunctions;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "GradLagrange", numComponents,
                                       gradBasisFunctions, numOrientations);

  double F[6] = {f(xa[0], xa[1]),   f(xb[0], xb[1]),   f(xc[0], xc[1]),
                 f(xab[0], xab[1]), f(xbc[0], xbc[1]), f(xca[0], xca[1])};
  double X[6] = {xa[0], xb[0], xc[0], xab[0], xbc[0], xca[0]};
  double Y[6] = {xa[1], xb[1], xc[1], xab[1], xbc[1], xca[1]};
  double e2 = 0;
  for(size_t i = 0; i < weights.size(); i++) {
    double interpolated = 0, x = 0, y = 0, dxdu = 0, dxdv = 0, dydu = 0, dydv = 0;

    for(size_t j = 0; j < 6; j++) {
      x += basisFunctions[6 * i + j] * X[j];
      y += basisFunctions[6 * i + j] * Y[j];
      dxdu += gradBasisFunctions[3 * (6 * i + j) + 0] * X[j];
      dxdv += gradBasisFunctions[3 * (6 * i + j) + 1] * X[j];
      dydu += gradBasisFunctions[3 * (6 * i + j) + 0] * Y[j];
      dydv += gradBasisFunctions[3 * (6 * i + j) + 1] * Y[j];
      interpolated += basisFunctions[6 * i + j] * F[j];
    }
    double detJ = fabs(dxdu * dydv - dxdv * dydu);
    std::vector<double> pos(3,0.);
    pos[0] = x;
    pos[1] = y;
    double exact = (solExact != nullptr) ? solExact->eval(0, pos) : fExact(x, y);
    double diff = exact - interpolated;
    e2 += weights[i] * diff * diff * detJ;
  }
  return e2;
}
#endif

feCurvedAdapt::feCurvedAdapt(feMesh *mesh, std::vector<feRecovery *> &recovery, feMetricOptions metricOptions, std::string metricMeshName, std::string nextMeshName, int analytical,
  feFunction *solExact) : _rec(recovery) {
#ifdef HAVE_GMSH

  int deg = metricOptions.polynomialDegree;
  const double lMin = metricOptions.hMin;
  const double lMax = metricOptions.hMax;
  const double eps  = metricOptions.eTargetError;

  useAnalytical = analytical;

  gmsh::initialize();
  gmsh::open(metricMeshName);
  gmsh::model::getCurrent(nameModel);
  gmsh::model::add("rect");
  gmsh::model::setCurrent("rect");
  gmsh::merge(metricMeshName)
;
  // get the nodes
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int dim = -1;
  int tag = -1;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tag, true, false);
  // dim = 2;
  // tag = 2; // Tag de l'entite physique associee a la surface
  // gmsh::model::mesh::getNodesForPhysicalGroup(dim, tag, nodeTags, coord);  

  for(auto val : nodeTags)
    std::cout<<val<<std::endl;

  int viewTag = gmsh::view::add(":metricP2_straight");
  int viewTagF = gmsh::view::add("F");
  int viewTags0 = gmsh::view::add("h0");
  int viewTags1 = gmsh::view::add("h1");
  std::vector<std::vector<double> > data;
  std::vector<std::vector<double> > dataF;
  std::vector<std::vector<double> > datas0;
  std::vector<std::vector<double> > datas1;

  double g00, g01, g11;

  FILE *F = fopen("directionsBeforeSmoothing.pos", "w");
  fprintf(F, "View \" \"{\n");

  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double C, S;
    switch(deg){
      case 1 :
        directionsFromGradient(x, y, F, C, S);
        // directionsFromHessian(x, y, F, C, S);
        // computeMetricP1Hess(x, y, lMin, lMax, eps, g00, g01, g11, F, C, S);
        break;
      case 2 :
        directionsFromGradient(x, y, F, C, S);
        // directionsFromHessian(x, y, F, C, S);
        // computeMetricP2(x, y, lMin, lMax, eps, g00, g01, g11, F, C, S);
        // computeMetricP2Hess(x, y, lMin, lMax, eps, g00, g01, g11, F, C, S);
        break;
      default :
        printf("No metric computation scheme for deg = 0 or deg > 2\n");
        exit(-1);
    }
    
    // std::cout<<C<<" - "<<S<<std::endl;

    COS[nodeTags[i]] = C;
    SIN[nodeTags[i]] = S;
  }
  fprintf(F,"};"); fclose(F);

  smoothDirections(COS, SIN);

  FILE* fff = fopen("ellipsesJF.pos", "w");
  FILE* ffff = fopen("detail.pos", "w");
  fprintf(fff,"View \"ellipsesJF\"{\n");
  F = fopen("directionsAfterSmoothing.pos", "w");
  fprintf(F, "View \" \"{\n");

  double factor = 10.;
  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  std::vector<SMetric3> metrics(nodeTags.size());
  
  // Write the .sol size field file
  FILE *myfile = fopen("toIntersect2D.sol", "w");
  fprintf(myfile, "MeshVersionFormatted 2\n\n");
  fprintf(myfile, "Dimension 2\n\n");
  fprintf(myfile, "SolAtVertices\n");
  std::set<int> s( nodeTags.begin(), nodeTags.end() );
  fprintf(myfile, "%ld\n", s.size());
  fprintf(myfile, "1 3\n\n");

  std::set<std::size_t> alreadyWrittenTags;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double C = COS[nodeTags[i]];
    double S = SIN[nodeTags[i]];

    fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., C, S, 0.);
    // fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., -sq  (2.)*S, sqrt(2.)*C, 0.);

    double l0, l1;
    switch(deg){
      case 1 :
      {
        const double dtt0_ = fabs(dtt(x, y, C, S));
        const double dtt1_ = fabs(dtt(x, y, -S, C));

        l0 = fabs(dtt0_) > 1e-14 ? pow(2. * eps / dtt0_, 0.5) : lMax;
        l1 = fabs(dtt1_) > 1e-14 ? pow(2. * eps / dtt1_, 0.5) : lMax;
        break;
      }
      case 2 : 
      {
        const double dttt0_ = fabs(dttt(x, y, C, S));
        const double dttt1_ = fabs(dttt(x, y, -S, C));

        l0 = fabs(dttt0_) > 1e-14 ? pow(6 * eps / dttt0_, 0.3333) : lMax;
        l1 = fabs(dttt1_) > 1e-14 ? pow(6 * eps / dttt1_, 0.3333) : lMax;
        break;
      }
      default :
        printf("No metric computation scheme for deg = 0 or deg > 2\n");
        exit(-1);
    }
    
    l0 = std::min(l0, lMax);
    l0 = std::max(l0, lMin);
    l1 = std::min(l1, lMax);
    l1 = std::max(l1, lMin);

    double h0 = 1. / (l0 * l0);
    double h1 = 1. / (l1 * l1);

    if(fabs(sqrt(C*C + S*S) - 1.0) > 1e-4){
      printf("C = %10.10e \t S = %10.10e \t norm = %10.10e\n", C, S, sqrt(C*C + S*S));
    }

    double g00 = C * C * h0 + S * S * h1;
    double g11 = S * S * h0 + C * C * h1;
    double g01 = S * C * (h1 - h0);

    if(g00*g11 - g01*g01 < 1e-10){
      std::cout<<g00*g11 - g01*g01<<std::endl;
      exit(-1);
    }

    std::vector<double> v(9);
    std::vector<double> vF(1);
    std::vector<double> s0(1);
    std::vector<double> s1(1);

    v[0] = g00;
    v[1] = -g01;
    v[2] = 0;

    v[3] = -g01;
    v[4] = g11;
    v[5] = 0;

    v[6] = 0;
    v[7] = 0;
    v[8] = f(x, y); // export f as well.
    vF[0] = f(x, y);

    s0[0] = l0;
    s1[0] = l1;

    data.push_back(v);
    dataF.push_back(vF);
    datas0.push_back(s0);
    datas1.push_back(s1);

    SMetric3 M;
    M.set_m11(g00);
    M.set_m21(-g01);
    M.set_m22(g11);
    metrics[i] = M;
    getEllipsePoints(factor * M(0,0), factor * 2.0*M(0,1), factor * M(1,1), x, y, xP, yP);
    // M.print("");
    for(int j = 0; j < nt; ++j){
      if(j != nt-1){
        fprintf(fff,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0., xP[j+1], yP[j+1], 0., 1, 1);
      } else{
        fprintf(fff,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0., xP[0], yP[0], 0., 1, 1);
      }
    }

    auto pp = alreadyWrittenTags.insert(nodeTags[i]);

    if(pp.second){
      fprintf(ffff, "%d \t %+-4.4e \t %+-4.4e : %+-4.4e \t %+-4.4e \t %+-4.4e \n", nodeTags[i], x, y, M(0, 0), M(0, 1), M(1, 1));
      // fprintf(ffff, "%d \t %+-4.4e \t %+-4.4e : %+-4.4e \t %+-4.4e \t %+-4.4e \n", nodeTags[i], x, y, M(0, 0), M(0, 1), M(1, 1));
      // Write .sol size field
      fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f\n", x, y, M(0, 0), M(0, 1), M(1, 1));
    }
  }
  fprintf(myfile, "End");
  fclose(myfile);
  fprintf(fff,"};"); fclose(fff);
  fprintf(F  ,"};"); fclose(F);
  fclose(ffff);

  // // Scale the metric field to obtain a desired number of vertices
  // metricScaling(recovery[0], metricOptions, metrics);

  // for(size_t i = 0; i < nodeTags.size(); i++) {
  //   SMetric3 M = metrics[i];
  //   double g00 = M(0,0);
  //   double g01 = -M(1,0);
  //   double g11 = M(1,1);

  //   std::vector<double> v(9);

  //   v[0] = g00;
  //   v[1] = -g01;
  //   v[2] = 0;

  //   v[3] = -g01;
  //   v[4] = g11;
  //   v[5] = 0;

  //   v[6] = 0;
  //   v[7] = 0;
  //   v[8] = data[i][8];
  //   data[i] = v;
  // }

  gmsh::view::addModelData(viewTag, 0, "rect", "NodeData", nodeTags, data);
  gmsh::view::addModelData(viewTagF, 0, "rect", "NodeData", nodeTags, dataF);
  gmsh::view::addModelData(viewTags0, 0, "rect", "NodeData", nodeTags, datas0);
  gmsh::view::addModelData(viewTags1, 0, "rect", "NodeData", nodeTags, datas1);

  printf("dataf.size = %lu\n", dataF.size());

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

  //  gmsh::write("metric.msh");
  gmsh::view::write(viewTagF, "sol.msh");
  gmsh::view::write(viewTag, "metric.msh");
  gmsh::view::write(viewTags0, "s0.msh");
  gmsh::view::write(viewTags1, "s1.msh");
  //  gmsh::view::write(viewTagF,"metric.msh",true);

  std::string mmgCommand  = "mmg2d metric.msh -hgrad 100 -o " + nextMeshName;
  std::string gmshCommand = "gmsh " + nextMeshName + " &";
  system(mmgCommand.c_str());
  system(gmshCommand.c_str());

  gmsh::option::setNumber("Mesh.MshFileVersion", 4);
  gmsh::model::add("adapt");
  gmsh::model::setCurrent("adapt");
  gmsh::merge(nextMeshName);

  // int computePointsUsingScaledCrossFieldPlanarP2(
  //   const char *modelForMetric, const char *modelForMesh, int VIEW_TAG, int faceTag,
  //   std::vector<double> &pts,
  //   double er(double *, double *, double *, double *, double *, double *));

  // bool vazy = true;

  // std::vector<double> pts;
  // if(vazy)
  //   computePointsUsingScaledCrossFieldPlanarP2("rect", "rect_adapt", viewTag, 4, pts,
  //                                              ERROR_SQUARED_P2);
  // gmsh::model::setCurrent("rect_adapt");

  // if(!vazy) gmsh::model::mesh::setOrder(2);
  // gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true);

  // int viewTagF_adapt = gmsh::view::add("F_adapt");
  // std::vector<std::vector<double> > dataF_adapt;

  // for(size_t i = 0; i < nodeTags.size(); i++) {
  //   const double x = coord[3 * i + 0];
  //   const double y = coord[3 * i + 1];
  //   std::vector<double> vF(1);
  //   vF[0] = f(x, y);
  //   dataF_adapt.push_back(vF);
  // }
  // gmsh::model::setCurrent("rect_adapt");
  // gmsh::view::addModelData(viewTagF_adapt, 0, "rect_adapt", "NodeData", nodeTags, dataF_adapt);
  // gmsh::view::write(viewTagF_adapt, "sol_adapt.msh");

  // std::vector<int> elementTypes;

  // gmsh::model::setCurrent("rect_adapt");
  // size_t lastindex = nextMeshName.find_last_of(".");
  // std::string root = nextMeshName.substr(0, lastindex);
  // std::string infoFileName = root + ".txt";
  // std::string errorFileName = root + "_ERROR.pos";
  // FILE *infoFile = fopen(infoFileName.c_str(), "w");

  // {
  //   FILE *ERR = fopen(errorFileName.c_str(), "w");
  //   fprintf(ERR, "View\"ERROR_P1\"{\n");
  //   int triangleP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  //   std::vector<double> localCoord;
  //   std::vector<double> weights;
  //   gmsh::model::mesh::getIntegrationPoints(triangleP1, "Gauss20", localCoord, weights);
  //   std::vector<double> basisFunctions;
  //   int numComponents, numOrientations;
  //   gmsh::model::mesh::getBasisFunctions(triangleP1, localCoord, "Lagrange", numComponents,
  //                                        basisFunctions, numOrientations);
  //   std::vector<int> elementTypes;
  //   std::vector<std::vector<std::size_t>> elementTags;
  //   std::vector<std::vector<std::size_t>> nodeTags;
  //   std::vector<std::size_t> nodeTags2;
  //   std::vector<double> coord;
  //   std::vector<double> parametricCoord;
  //   gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);

  //   gmsh::model::mesh::getNodes(nodeTags2, coord, parametricCoord, 2);

  //   fprintf(infoFile, "nElem = %d\n", elementTags[0].size());

  //   double Error = 0.0;
  //   printf("Computing error on %d elements\n", elementTags[0].size());
  //   for(size_t e = 0; e < elementTags[0].size(); e++) {
  //     int eDim, eTag;
  //     std::vector<double> xa, xb, xc, xab, xbc, xca, pc;
  //     gmsh::model::mesh::getNode(nodeTags[0][3 * e + 0], xa, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][3 * e + 1], xb, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][3 * e + 2], xc, pc, eDim, eTag);
  //     if(deg > 1){
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 0], xa, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 1], xb, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 2], xc, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 3], xab, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 4], xbc, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 5], xca, pc, eDim, eTag);
  //     }
      
  //     double ErrorElement;
  //     switch(deg){
  //       case 1 : 
  //         ErrorElement = ERROR_SQUARED_P1(&xa[0], &xb[0], &xc[0], solExact);
  //         break;
  //       case 2 : 
  //         ErrorElement = ERROR_SQUARED_P2(&xa[0], &xb[0], &xc[0], &xab[0], &xbc[0], &xca[0], solExact);
  //         break;
  //       default :
  //         printf("Error is only computed for P1 and P2 interpolant.");
  //     } 
  //     fprintf(ERR,
  //             "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n",
  //             xa[0], xa[1], 0.0, xb[0], xb[1], 0.0, xc[0], xc[1], 0.0, ErrorElement, ErrorElement, ErrorElement);
  //     Error += ErrorElement;
  //   }

  //   sort( nodeTags[0].begin(), nodeTags[0].end() );
  //   nodeTags[0].erase( unique( nodeTags[0].begin(), nodeTags[0].end() ), nodeTags[0].end() );
  //   fprintf(infoFile, "nVertices = %d\n", nodeTags[0].size());

  //   fprintf(ERR, "};\n");
  //   fclose(ERR);
  //   printf("Error = %12.5e\n", sqrt(Error));
  //   fprintf(infoFile, "erreur interpolee = %12.5e\n", sqrt(Error));
  //   fclose(infoFile);
  // }

  //  {
  //   FILE *ERR = fopen("ERROR_P2.pos", "w");
  //   fprintf(ERR, "View\"ERROR_P2\"{\n");
  //   int triangleP2 = gmsh::model::mesh::getElementType("Triangle", 2);
  //   std::vector<double> localCoord;
  //   std::vector<double> weights;
  //   gmsh::model::mesh::getIntegrationPoints(triangleP2, "Gauss12", localCoord, weights);
  //   std::vector<double> basisFunctions;
  //   int numComponents, numOrientations;
  //   gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "Lagrange", numComponents,
  //                                        basisFunctions, numOrientations);
  //   std::vector<int> elementTypes;
  //   std::vector<std::vector<std::size_t> > elementTags;
  //   std::vector<std::vector<std::size_t> > nodeTags;
  //   gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);
  //   double Error = 0.0;

  //   for(size_t e = 0; e < elementTags[0].size(); e++) {
  //     int eDim, eTag;
  //     std::vector<double> xa, xb, xc, xab, xbc, xca, pc;
  //     gmsh::model::mesh::getNode(nodeTags[0][6 * e + 0], xa, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][6 * e + 1], xb, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][6 * e + 2], xc, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][6 * e + 3], xab, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][6 * e + 4], xbc, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][6 * e + 5], xca, pc, eDim, eTag);
  //     double ErrorElement = ERROR_SQUARED_P2(&xa[0], &xb[0], &xc[0], &xab[0], &xbc[0], &xca[0]);
  //     fprintf(ERR,
  //             "ST2(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g,%g,%g,%g};\n",
  //             xa[0], xa[1], 0.0, xb[0], xb[1], 0.0, xc[0], xc[1], 0.0, xab[0], xab[1], 0.0, xbc[0],
  //             xbc[1], 0.0, xca[0], xca[1], 0.0, ErrorElement, ErrorElement, ErrorElement,
  //             ErrorElement, ErrorElement, ErrorElement);
  //     Error += ErrorElement;
  //   }
  //   fprintf(ERR, "};\n");
  //   fclose(ERR);
  //   printf("Error = %12.5E\n", sqrt(Error));
  // }
#else
  printf("In feCurvedAdapt : Error : Gmsh is required to compute curved meshes.\n");
#endif
}

feCurvedAdapt::~feCurvedAdapt() {
  // for(feGridFunction *gF : _derivatives)
  // delete gF;
}

feIsotropicAdapt::feIsotropicAdapt(feMesh *mesh, std::vector<feRecovery *> &recovery, feMetricOptions metricOptions, std::string metricMeshName, std::string nextMeshName, int analytical,
  feFunction *solExact) : _rec(recovery) {
#ifdef HAVE_GMSH

  int deg = metricOptions.polynomialDegree;
  const double lMin = metricOptions.hMin;
  const double lMax = metricOptions.hMax;
  const double eps  = metricOptions.eTargetError;

  useAnalytical = analytical;

  gmsh::initialize();
  gmsh::open(metricMeshName);
  gmsh::model::getCurrent(nameModel);
  gmsh::model::add("iso");
  gmsh::model::setCurrent("iso");
  gmsh::merge(metricMeshName);

  // get the nodes
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int tag = 1;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, tag, true, false);  

  int viewTag = gmsh::view::add("isoSizeField");
  int viewTagF = gmsh::view::add("F");
  std::vector<std::vector<double> > data;
  std::vector<std::vector<double> > dataF;

  double g00, g01, g11;

  FILE* fff = fopen("ellipsesIso.pos", "w");
  fprintf(fff,"View \"ellipsesIso\"{\n");

  double factor = 10.;
  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  std::vector<SMetric3> metrics(nodeTags.size());

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double C = 1.0;
    double S = 0.0;

    double l0, l1;
    switch(deg){
      case 1 :
      {
        const double dtt0_ = fabs(dtt(x, y, C, S));
        const double dtt1_ = fabs(dtt(x, y, -S, C));

        l0 = fabs(dtt0_) > 1e-14 ? pow(2. * eps / dtt0_, 0.5) : lMax;
        l1 = fabs(dtt1_) > 1e-14 ? pow(2. * eps / dtt1_, 0.5) : lMax;
        break;
      }
      case 2 : 
      {
        const double dttt0_ = fabs(dttt(x, y, C, S));
        const double dttt1_ = fabs(dttt(x, y, -S, C));

        l0 = fabs(dttt0_) > 1e-14 ? pow(6 * eps / dttt0_, 0.3333) : lMax;
        l1 = fabs(dttt1_) > 1e-14 ? pow(6 * eps / dttt1_, 0.3333) : lMax;
        break;
      }
      default :
        printf("No metric computation scheme for deg = 0 or deg > 2\n");
        exit(-1);
    }
    
    l0 = std::min(l0, lMax);
    l0 = std::max(l0, lMin);
    l1 = std::min(l1, lMax);
    l1 = std::max(l1, lMin);

    double h0 = 1. / (l0 * l0);
    double h1 = 1. / (l1 * l1);

    double size = fmin(l0, l1);

    double g00 = C * C * h0 + S * S * h1;
    double g11 = S * S * h0 + C * C * h1;
    double g01 = S * C * (h1 - h0);

    if(g00*g11 - g01*g01 < 1e-10){
      std::cout<<g00*g11 - g01*g01<<std::endl;
      exit(-1);
    }

    std::vector<double> v(1);
    std::vector<double> vF(1);

    v[0] = size;
    vF[0] = f(x, y);

    data.push_back(v);
    dataF.push_back(vF);

    SMetric3 M;
    M.set_m11(g00);
    M.set_m21(-g01);
    M.set_m22(g11);
    metrics[i] = M;
    getEllipsePoints(factor * M(0,0), factor * 2.0*M(0,1), factor * M(1,1), x, y, xP, yP);
    // M.print("");
    for(int i = 0; i < nt; ++i){
      if(i != nt-1){
        fprintf(fff,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[i+1], yP[i+1], 0., 1, 1);
      } else{
        fprintf(fff,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[0], yP[0], 0., 1, 1);
      }
    }
  }
  fprintf(fff,"};"); fclose(fff);

  gmsh::view::addModelData(viewTag, 0, "iso", "NodeData", nodeTags, data);
  gmsh::view::addModelData(viewTagF, 0, "iso", "NodeData", nodeTags, dataF);

  printf("dataf.size = %lu\n", dataF.size());

  // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

  //  gmsh::write("metric.msh");
  gmsh::view::write(viewTagF, "solution.msh");
  gmsh::view::write(viewTag, "isoSizeField.pos");

  gmsh::model::add("myBeautifulModel");
  // gmsh::merge("isoSizeField.pos");

  double lc = 0.2, xc = 2., yc = 2., r = 0.25;
  // Points
  gmsh::model::occ::addPoint( 0, 0, 0, lc, 1);
  gmsh::model::occ::addPoint(10, 0, 0, lc, 2);
  gmsh::model::occ::addPoint(10, 4, 0, lc, 3);
  gmsh::model::occ::addPoint( 0, 4, 0, lc, 4);
  // Courbes : un nouveau tag est renvoye si aucun n'est specifie
  int l1 = gmsh::model::occ::addLine(1, 2);
  int l2 = gmsh::model::occ::addLine(2, 3);
  int l3 = gmsh::model::occ::addLine(3, 4);
  int l4 = gmsh::model::occ::addLine(4, 1);
  int c1 = gmsh::model::occ::addCircle(xc, yc, 0, r);
  int boundary = gmsh::model::occ::addCurveLoop({l1,l2,l3,l4});
  int circle = gmsh::model::occ::addCurveLoop({c1});
  // Surface
  gmsh::model::occ::addPlaneSurface({boundary, circle}, 1);
  // Creation des structures Gmsh associees a la geometrie
  gmsh::model::occ::synchronize();
  // Entites physiques
  gmsh::model::addPhysicalGroup(1, {1}, 1);
  gmsh::model::addPhysicalGroup(1, {2}, 2);
  gmsh::model::addPhysicalGroup(1, {3}, 3);
  gmsh::model::addPhysicalGroup(1, {4}, 4);
  gmsh::model::addPhysicalGroup(1, {5}, 5);
  int surface = gmsh::model::addPhysicalGroup(2, {1});
  gmsh::model::setPhysicalName(1, 1, "Bas");
  gmsh::model::setPhysicalName(1, 2, "Droite");
  gmsh::model::setPhysicalName(1, 3, "Haut");
  gmsh::model::setPhysicalName(1, 4, "Gauche");
  gmsh::model::setPhysicalName(1, 5, "Cylindre");
  gmsh::model::setPhysicalName(2, surface, "Surface");

  gmsh::model::geo::synchronize();

  int bg_field = gmsh::model::mesh::field::add("PostView");
  gmsh::model::mesh::field::setNumber(bg_field, "ViewIndex", viewTag);

  // Apply the view as the current background mesh size field:
  gmsh::model::mesh::field::setAsBackgroundMesh(bg_field);

  gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

  gmsh::model::mesh::setAlgorithm(2,1,5);

  gmsh::model::mesh::generate(2);
  // gmsh::write("aLittleTest.msh");
  gmsh::write(nextMeshName);

  // gmsh::finalize();

  // std::string mmgCommand = "mmg2d metric.msh -hgrad 100 -o " + nextMeshName;
  // system(mmgCommand.c_str());

  // gmsh::option::setNumber("Mesh.MshFileVersion", 4);
  // gmsh::model::add("adapt");
  // gmsh::model::setCurrent("adapt");
  // gmsh::merge(nextMeshName);

  // int computePointsUsingScaledCrossFieldPlanarP2(
  //   const char *modelForMetric, const char *modelForMesh, int VIEW_TAG, int faceTag,
  //   std::vector<double> &pts,
  //   double er(double *, double *, double *, double *, double *, double *));

  // bool vazy = true;

  // std::vector<double> pts;
  // if(vazy)
  //   computePointsUsingScaledCrossFieldPlanarP2("rect", "rect_adapt", viewTag, 4, pts,
  //                                              ERROR_SQUARED_P2);
  // gmsh::model::setCurrent("rect_adapt");

  // if(!vazy) gmsh::model::mesh::setOrder(2);
  // gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true);

  // int viewTagF_adapt = gmsh::view::add("F_adapt");
  // std::vector<std::vector<double> > dataF_adapt;

  // for(size_t i = 0; i < nodeTags.size(); i++) {
  //   const double x = coord[3 * i + 0];
  //   const double y = coord[3 * i + 1];
  //   std::vector<double> vF(1);
  //   vF[0] = f(x, y);
  //   dataF_adapt.push_back(vF);
  // }
  // gmsh::model::setCurrent("rect_adapt");
  // gmsh::view::addModelData(viewTagF_adapt, 0, "rect_adapt", "NodeData", nodeTags, dataF_adapt);
  // gmsh::view::write(viewTagF_adapt, "sol_adapt.msh");

  // std::vector<int> elementTypes;

  // gmsh::model::setCurrent("rect_adapt");
  // size_t lastindex = nextMeshName.find_last_of(".");
  // std::string root = nextMeshName.substr(0, lastindex);
  // std::string infoFileName = root + ".txt";
  // std::string errorFileName = root + "_ERROR.pos";
  // FILE *infoFile = fopen(infoFileName.c_str(), "w");

  // {
  //   FILE *ERR = fopen(errorFileName.c_str(), "w");
  //   fprintf(ERR, "View\"ERROR_P1\"{\n");
  //   int triangleP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  //   std::vector<double> localCoord;
  //   std::vector<double> weights;
  //   gmsh::model::mesh::getIntegrationPoints(triangleP1, "Gauss20", localCoord, weights);
  //   std::vector<double> basisFunctions;
  //   int numComponents, numOrientations;
  //   gmsh::model::mesh::getBasisFunctions(triangleP1, localCoord, "Lagrange", numComponents,
  //                                        basisFunctions, numOrientations);
  //   std::vector<int> elementTypes;
  //   std::vector<std::vector<std::size_t>> elementTags;
  //   std::vector<std::vector<std::size_t>> nodeTags;
  //   std::vector<std::size_t> nodeTags2;
  //   std::vector<double> coord;
  //   std::vector<double> parametricCoord;
  //   gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);

  //   gmsh::model::mesh::getNodes(nodeTags2, coord, parametricCoord, 2);

  //   fprintf(infoFile, "nElem = %d\n", elementTags[0].size());

  //   double Error = 0.0;
  //   printf("Computing error on %d elements\n", elementTags[0].size());
  //   for(size_t e = 0; e < elementTags[0].size(); e++) {
  //     int eDim, eTag;
  //     std::vector<double> xa, xb, xc, xab, xbc, xca, pc;
  //     gmsh::model::mesh::getNode(nodeTags[0][3 * e + 0], xa, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][3 * e + 1], xb, pc, eDim, eTag);
  //     gmsh::model::mesh::getNode(nodeTags[0][3 * e + 2], xc, pc, eDim, eTag);
  //     if(deg > 1){
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 0], xa, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 1], xb, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 2], xc, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 3], xab, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 4], xbc, pc, eDim, eTag);
  //       gmsh::model::mesh::getNode(nodeTags[0][6 * e + 5], xca, pc, eDim, eTag);
  //     }
      
  //     double ErrorElement;
  //     switch(deg){
  //       case 1 : 
  //         ErrorElement = ERROR_SQUARED_P1(&xa[0], &xb[0], &xc[0], solExact);
  //         break;
  //       case 2 : 
  //         ErrorElement = ERROR_SQUARED_P2(&xa[0], &xb[0], &xc[0], &xab[0], &xbc[0], &xca[0], solExact);
  //         break;
  //       default :
  //         printf("Error is only computed for P1 and P2 interpolant.");
  //     } 
  //     fprintf(ERR,
  //             "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n",
  //             xa[0], xa[1], 0.0, xb[0], xb[1], 0.0, xc[0], xc[1], 0.0, ErrorElement, ErrorElement, ErrorElement);
  //     Error += ErrorElement;
  //   }

  //   sort( nodeTags[0].begin(), nodeTags[0].end() );
  //   nodeTags[0].erase( unique( nodeTags[0].begin(), nodeTags[0].end() ), nodeTags[0].end() );
  //   fprintf(infoFile, "nVertices = %d\n", nodeTags[0].size());

  //   fprintf(ERR, "};\n");
  //   fclose(ERR);
  //   printf("Error = %12.5e\n", sqrt(Error));
  //   fprintf(infoFile, "erreur interpolee = %12.5e\n", sqrt(Error));
  //   fclose(infoFile);
  // }

#else
  printf("In feIsotropicAdapt : Error : Gmsh is required.\n");
#endif
}

feIsotropicAdapt::~feIsotropicAdapt() {
}