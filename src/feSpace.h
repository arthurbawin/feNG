#ifndef _FESPACE_
#define _FESPACE_

#include <vector>
#include <string>
#include <iostream>

#include "feMessage.h"
#include "feQuadrature.h"
#include "feFunction.h"
#include "feCncGeo.h"
#include "feNG.h"

#include "../contrib/Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EigenMat;

// Supported geometries and choices of degrees of freedom to define finite element spaces
typedef enum { POINT, LINE, TRI } geometryType;
typedef enum { LAGRANGE, LEGENDRE, NONCONFORMAL } shapeType;

class feMesh;
class feNumber;
class feSolution;

// Create a finite element space and perform safety checks.
// This is the recommended way of creating an FE space.
// Call within feCheck( ... ) to safely exit if a problem arises.
//
//                   space: pointer to the FE space, initially undefined and assigned
//                          during the call.
//                    mesh: valid pointer to the mesh
//               dimension: dimension of the connectivity on which the FE space is defined
//                geometry: geometric support of the finite element (see enum above)
//                   shape: the choice of degrees of freedom (see enum above)
//                  degree: degree (order) of the polynomial interpolation
//               fieldName: name of the field for which the space is defined
//        connectivityName: name of the geometric connectivity on which the FE space is defined
//        degreeQuadrature: maximum degree of the polynomial to integrate perfectly with
//                          a quadrature rule
//                     fct: field used to initialize the degrees of freedom
// useGlobalShapeFunctions: choose between local (reference space) or global (physical space)
//                          shape functions (global still experimental, do not use yet).
feStatus createFiniteElementSpace(feSpace *&space,
  feMesh *mesh,
  const int dimension,
  const geometryType geometry,
  const shapeType shape,
  const int degree,
  const std::string fieldName,
  const std::string connectivityName,
  const int degreeQuadrature,
  feFunction *fct,
  const bool useGlobalShapeFunctions = false);

// Abstract class defining a finite element space, handling the interpolation
// functions defined over the elements. A finite element space can be used 
// to interpolate either the geometry or the fields of interest.
//
// An FE space is defined on a single geometric connectivity (which 
// corresponds to a Physical Entity in the mesh).
class feSpace
{
protected:
  // Pointer to the mesh on which the space is defined
  feMesh *_mesh; 
  // Name of the associated field, e.g "U", "V", "P"
  std::string _fieldID;
  // Number (tag) of the associated field
  int _fieldTag;
  // Name of the geometric connectivity on which the space is defined
  std::string _cncGeoID;
  // Number (tag) of the geometric connectivity
  int _cncGeoTag;

  // Quadrature rule on the reference element:
  // number of quadrature nodes, weights and coordinates
  int _nQuad;
  std::vector<double> _wQuad;
  std::vector<double> _xQuad;
  std::vector<double> _yQuad;
  std::vector<double> _zQuad;

  // Number of interpolation functions (i.e. number of DOFs)
  int _nFunctions;
  // Coordinates of the collocation nodes on the reference element
  std::vector<double> _Lcoor;
  // Evaluation of the (local) shape functions and their derivatives at
  // the quadrature nodes.
  std::vector<double> _L;
  std::vector<double> _dLdr;
  std::vector<double> _dLds;
  std::vector<double> _dLdt;

  // Use global (physical) interpolation functions (experimental, do not use yet)
  bool _useGlobalShapeFunctions = false;
  // Global shape functions and derivatives evaluated at quadrature points on each element
  std::vector<std::vector<double> > _Lglob;
  std::vector<std::vector<double> > _dLdxglob;
  std::vector<std::vector<double> > _dLdyglob;

  // Scalar field used to initialize the DOFs
  feFunction *_fct;

public:
  // Do not use the abstract class constructor directly, call the derived classes.
  feSpace(feMesh *mesh = nullptr, const std::string &fieldID = "", const std::string &cncGeoID = "",
    feFunction *fct = nullptr, const bool useGlobalShapeFunctions = false);
  virtual ~feSpace(){};

  const std::string &getFieldID() { return _fieldID; }
  int getFieldTag() { return _fieldTag; }
  const std::string &getCncGeoID() { return _cncGeoID; }
  int getCncGeoTag() { return _cncGeoTag; }
  feCncGeo *getCncGeo();

  // Assign the tag of a geometric space after the mesh has been created
  void setCncGeoTag(int tag){ _cncGeoTag = tag; }
  // Assign the mesh pointer of a geometric space after the mesh has been created
  void setMeshPtr(feMesh *mesh){ _mesh = mesh; }

  // Return the attributes of the geometric connectivity on which the space is defined
  int getDim();
  int getNbElm();
  int getNbNodePerElem();

  bool isFctDefined() { return !(_fct == nullptr); }
  bool useGlobalFunctions() { return _useGlobalShapeFunctions; }
  void useGlobalFunctions(bool flag) { _useGlobalShapeFunctions = flag; }

  // Return the number of shape functions (DOF) on an element
  virtual int getNbFunctions() = 0;
  // Return highest degree of the polynomial basis
  virtual int getPolynomialDegree() = 0;

  EigenMat innerProductBasisFunctions(int iElm);
  double innerProductBasisFunctions(int iElm, int ex, int ey);

  const std::vector<double> &getLcoor() { return _Lcoor; }

  // Evaluate the local shape functions and derivatives at quadrature nodes
  virtual std::vector<double> L(double *r) = 0;
  virtual void L(double *r, double *L) = 0;
  virtual std::vector<double> dLdr(double *r) = 0;
  virtual std::vector<double> dLds(double *r) = 0;
  virtual std::vector<double> dLdt(double *r) = 0;

  // Evaluate the global shape functions at quadrature nodes
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy) = 0;

  // Return value of shape function (or derivatives) at quadrature node
  double getFunctionAtQuadNode(int iFun, int iQuadNode)
  {
    return _L[_nFunctions * iQuadNode + iFun];
  }
  double getdFunctiondrAtQuadNode(int iFun, int iQuadNode)
  {
    return _dLdr[_nFunctions * iQuadNode + iFun];
  }
  double getdFunctiondsAtQuadNode(int iFun, int iQuadNode)
  {
    return _dLds[_nFunctions * iQuadNode + iFun];
  }
  double getGlobalFunctionAtQuadNode(int iElm, int iFun, int iQuadNode)
  {
    return _Lglob[iElm][_nFunctions * iQuadNode + iFun];
  }
  double getdGlobalFunctiondxAtQuadNode(int iElm, int iFun, int iQuadNode)
  {
    return _dLdxglob[iElm][_nFunctions * iQuadNode + iFun];
  }
  double getdGlobalFunctiondyAtQuadNode(int iElm, int iFun, int iQuadNode)
  {
    return _dLdyglob[iElm][_nFunctions * iQuadNode + iFun];
  }

  // Set and get quadrature rule
  feStatus setQuadratureRule(feQuadrature *quadratureRule);
  int getNbQuadPoints() { return _nQuad; }
  std::vector<double> &getQuadratureWeights() { return _wQuad; }
  std::vector<double> &getRQuadraturePoints() { return _xQuad; }
  std::vector<double> &getSQuadraturePoints() { return _yQuad; }
  std::vector<double> &getTQuadraturePoints() { return _zQuad; }

  virtual void initializeNumberingUnknowns(feNumber *number) = 0;
  virtual void initializeNumberingEssential(feNumber *number) = 0;
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr) = 0;

  // Evaluate prescribed scalar field
  double evalFun(const double t, const std::vector<double> &x) { return _fct->eval(t, x); }

  //
  // Interpolation of fields and derivatives at reference node r = [r,s,t] or physical node x:
  //

  // Interpolate scalar field or derivatives at reference node r = [r,s,t]
  // using local shape functions (default)
  double interpolateField(std::vector<double> &field, double *r);
  double interpolateField_rDerivative(std::vector<double> &field, double *r);
  double interpolateField_sDerivative(std::vector<double> &field, double *r);

  // Interpolate scalar field or derivatives at physical node x on element iElm
  // using global shape functions
  double interpolateField(std::vector<double> &field, int iElm, std::vector<double> &x);
  double interpolateField_xDerivative(std::vector<double> &field, int iElm, std::vector<double> &x);
  double interpolateField_yDerivative(std::vector<double> &field, int iElm, std::vector<double> &x);

  // Interpolate field and also get the shape functions at r.
  // size(field) = size(shape) = fieldSize.
  // shape will be filled with values of the shape functions.
  // The result of the interpolation is stored in res.
  void interpolateField(double *field, int fieldSize, double *r, double *shape, double &res);

  double interpolateField(feNumber *number, feSolution *sol, std::vector<double> &x);
  void interpolateField_gradrs(feNumber *number, feSolution *sol, std::vector<double> &x,
                               std::vector<double> &grad);

  // Interpolate vector field or derivatives at reference node r = [r,s,t]
  // using local shape functions (default). Result is stored in res.
  void interpolateVectorField(std::vector<double> &field, double *r, std::vector<double> &res);
  void interpolateVectorField_rDerivative(std::vector<double> &field, double *r, std::vector<double> &res);
  void interpolateVectorField_sDerivative(std::vector<double> &field, double *r, std::vector<double> &res);

  //
  // Interpolation of fields and derivatives at quadrature node:
  //

  // Interpolate scalar field or derivatives at iNode-th quadrature node using local shape functions (default)
  double interpolateFieldAtQuadNode(std::vector<double> &field, int iNode);
  double interpolateFieldAtQuadNode_rDerivative(std::vector<double> &field, int iNode);
  double interpolateFieldAtQuadNode_sDerivative(std::vector<double> &field, int iNode);
  // Interpolate scalar field or derivatives at iNode-th quadrature node using global shape functions
  double interpolateFieldAtQuadNode(std::vector<double> &field, int iElm, int iNode);
  double interpolateFieldAtQuadNode_xDerivative(std::vector<double> &field, int iElm, int iNode);
  double interpolateFieldAtQuadNode_yDerivative(std::vector<double> &field, int iElm, int iNode);

  // Interpolate vector field or derivatives at iNode-th quadrature node using local shape functions (default)
  void interpolateVectorFieldAtQuadNode(std::vector<double> &field, int iNode,
                                        std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> &field,
    int iNode, std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> &field,
    int iNode, std::vector<double> &res);
};

//
// Lagrange interpolation functions of degree 0 on 1D reference element [-1,1]
//
class feSpace1DP0 : public feSpace
{
public:
  feSpace1DP0(const std::string &cncGeoID);
  feSpace1DP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP0(){};

  int getNbFunctions() { return 1; }
  int getPolynomialDegree() { return 0; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
    std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 1 on 1D reference element [-1,1]
//
class feSpace1DP1 : public feSpace
{
protected:
public:
  feSpace1DP1(std::string cncGeoID);
  feSpace1DP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP1(){};

  int getNbFunctions() { return 2; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Non-conformal interpolation functions of degree 1
//
class feSpace1DP1_nonConsistant : public feSpace
{
protected:
public:
  feSpace1DP1_nonConsistant(std::string cncGeoID);
  feSpace1DP1_nonConsistant(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP1_nonConsistant() {}

  int getNbFunctions() { return 1; }
  int getPolynomialDegree() { return 0; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 2 on 1D reference element [-1,1]
//
class feSpace1DP2 : public feSpace
{
protected:
public:
  feSpace1DP2(std::string cncGeoID);
  feSpace1DP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP2() {}

  int getNbFunctions() { return 3; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
    std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 3 on 1D reference element [-1,1]
//
class feSpace1DP3 : public feSpace
{
protected:
public:
  feSpace1DP3(std::string cncGeoID);
  feSpace1DP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP3() {}

  int getNbFunctions() { return 4; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 4 on 1D reference element [-1,1]
//
class feSpace1DP4 : public feSpace
{
protected:
public:
  feSpace1DP4(std::string cncGeoID);
  feSpace1DP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP4() {}

  int getNbFunctions() { return 5; }
  int getPolynomialDegree() { return 4; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Legendre interpolation functions of arbitrary degree on 1D reference element [-1,1]
//
class feSpace1D_Legendre : public feSpace
{
protected:
  int _degree;
public:
  feSpace1D_Legendre(const int degree, feMesh *mesh,
    const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1D_Legendre() {}

  int getNbFunctions() { return _nFunctions; }
  int getPolynomialDegree() { return _degree; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 1 on 2D reference triangle r = [0,1], s = [0,1-r]
//
class feSpaceTriP1 : public feSpace
{
protected:
public:
  feSpaceTriP1(std::string cncGeoID);
  feSpaceTriP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP1() {}

  int getNbFunctions() { return 3; }
  int getPolynomialDegree() { return 1; }
  std::vector<double> L(double *r);
  void L(double *r, double *L);
  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Non-conformal interpolation functions of degree 1 on 2D reference triangle r = [0,1], s = [0,1-r]
//
class feSpaceTriP1_nonConsistant : public feSpace
{
protected:
public:
  feSpaceTriP1_nonConsistant(std::string cncGeoID);
  feSpaceTriP1_nonConsistant(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                             feFunction *fct);
  ~feSpaceTriP1_nonConsistant() {}

  int getNbFunctions() { return 3; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 2 on 2D reference triangle r = [0,1], s = [0,1-r]
//
class feSpaceTriP2 : public feSpace
{
protected:
public:
  feSpaceTriP2(std::string cncGeoID);
  feSpaceTriP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP2() {}

  int getNbFunctions() { return 6; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Non-conformal interpolation functions of degree 2 on 2D reference triangle r = [0,1], s = [0,1-r]
//
class feSpaceTriP2_nonConsistant : public feSpace
{
protected:
public:
  feSpaceTriP2_nonConsistant(std::string cncGeoID);
  feSpaceTriP2_nonConsistant(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                             feFunction *fct);
  ~feSpaceTriP2_nonConsistant() {}

  int getNbFunctions() { return 7; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 3 on 2D reference triangle r = [0,1], s = [0,1-r]
//
class feSpaceTriP3 : public feSpace
{
protected:
public:
  feSpaceTriP3(std::string cncGeoID);
  feSpaceTriP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP3() {}

  int getNbFunctions() { return 10; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

//
// Lagrange interpolation functions of degree 4 on 2D reference triangle r = [0,1], s = [0,1-r]
//
class feSpaceTriP4 : public feSpace
{
protected:
public:
  feSpaceTriP4(std::string cncGeoID);
  feSpaceTriP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP4() {}

  int getNbFunctions() { return 15; }
  int getPolynomialDegree() { return 4; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);

  feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };

  void initializeNumberingUnknowns(feNumber *number);
  void initializeNumberingEssential(feNumber *number);
  void initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr);
};

#endif
