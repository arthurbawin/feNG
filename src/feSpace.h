#ifndef _FESPACE_
#define _FESPACE_

#include <vector>
#include <string>
#include <iostream>

#include "feMessage.h"
#include "feQuadrature.h"
#include "feFunction.h"
#include "feCncGeo.h"

#include "../contrib/Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EigenMat;

/* Supported geometries to define finite element spaces */
typedef enum { POINT, LINE, LINE_CR, TRI, TRI_CR } elemType;

class feMesh;
class feNumber;
class feSolution;

feStatus createFiniteElementSpace(feSpace *&space, feMesh *mesh, int dim, elemType type, int deg,
                                  std::string fieldID, std::string cncGeoID, int dQuad,
                                  feFunction *fct, bool useGlobalShapeFunctions = false);

/* Un feSpace est utilisé
  - soit pour définir les interpolants géométriques sur une partie d'un maillage (sur une cncGeo)
  - soit pour définir les interpolants utilisés pour la résolution (sur un field)
*/
class feSpace
{
protected:
  feMesh *_mesh; // Pointer to the mesh on which the space is defined

  std::string _fieldID; // Name of the associated field, e.g "U", "V", "P"
  int _fieldTag; // Number (tag) of the associated field
  std::string _cncGeoID; // Name of the geometric connectivity on which the space is defined
  int _cncGeoTag; // Number (tag) of the geometric connectivity

  // Quadrature
  int _nQuad;
  std::vector<double> _wQuad;
  std::vector<double> _xQuad;
  std::vector<double> _yQuad;
  std::vector<double> _zQuad;

  // Interpolants
  int _nFunctions; // nielm
  std::vector<double> _Lcoor;
  std::vector<double> _L;
  std::vector<double> _dLdr;
  std::vector<double> _dLds;
  std::vector<double> _dLdt;

  // Global (physical) interpolation functions
  bool _useGlobalShapeFunctions = false;
  // The global shape functions and derivatives evaluated at quadrature points on each element
  std::vector<std::vector<double> > _Lglob;
  std::vector<std::vector<double> > _dLdxglob;
  std::vector<std::vector<double> > _dLdyglob;

  std::vector<int> _adr;
  std::vector<double> _sol;
  std::vector<double> _soldot;

  feFunction *_fct;

public:
  feSpace(feMesh *mesh = nullptr, std::string fieldID = "", std::string cncGeoID = "",
          feFunction *fct = nullptr, bool useGlobalShapeFunctions = false);
  virtual ~feSpace()
  {
    _Lcoor.clear();
    _Lcoor.resize(0);
  }

  int getDim();
  int getNbElm();
  int getNbNodePerElem();

  std::string getFieldID() { return _fieldID; }
  int getFieldTag() { return _fieldTag; }
  std::string getCncGeoID() { return _cncGeoID; }
  int getCncGeoTag() { return _cncGeoTag; }
  void setCncGeoTag(int tag)
  {
    _cncGeoTag = tag;
  } // Used to assign the tag of a geometric space after the mesh has been created

  void setMeshPtr(feMesh *mesh)
  {
    _mesh = mesh;
  } // Used to assign the mesh pointer of a geometric space after the mesh has been created

  feCncGeo *getCncGeo();

  bool isFctDefined() { return !(_fct == nullptr); }

  bool useGlobalFunctions() { return _useGlobalShapeFunctions; }
  void useGlobalFunctions(bool flag) { _useGlobalShapeFunctions = flag; }

  // The number of degrees of freedom on an element
  virtual int getNbFunctions() = 0;
  // The highest degree of the polynomial basis
  virtual int getPolynomialDegree() = 0;

  EigenMat innerProductBasisFunctions(int iElm);
  double innerProductBasisFunctions(int iElm, int ex, int ey);

  const std::vector<double> &getLcoor() { return _Lcoor; }

  virtual std::vector<double> L(double r[3]) = 0;
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy) = 0;
  virtual std::vector<double> dLdr(double r[3]) = 0;
  virtual std::vector<double> dLds(double r[3]) = 0;
  virtual std::vector<double> dLdt(double r[3]) = 0;

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

  feStatus setQuadratureRule(feQuadrature *quad);
  int getNbQuadPoints() { return _nQuad; }
  std::vector<double> &getQuadratureWeights() { return _wQuad; }
  std::vector<double> &getRQuadraturePoints() { return _xQuad; }
  std::vector<double> &getSQuadraturePoints() { return _yQuad; }
  std::vector<double> &getTQuadraturePoints() { return _zQuad; }

  void initializeSolution(feSolution *sol);
  void initializeSolution(std::vector<double> &sol);
  void initializeSolutionDot(feSolution *sol);

  virtual void initializeNumberingUnknowns(feNumber *number) = 0;
  virtual void initializeNumberingEssential(feNumber *number) = 0;
  virtual void initializeAddressingVector(feNumber *number, int numElem) = 0;
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr) = 0;

  int getAddressingVectorAt(int node) { return _adr[node]; }

  double evalFun(const double t, const std::vector<double> &x) { return _fct->eval(t, x); }

  double interpolateField(std::vector<double> field, double r[3]);
  double interpolateFieldAtQuadNode(std::vector<double> field, int iNode);
  double interpolateField_rDerivative(std::vector<double> field, double r[3]);
  double interpolateFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode);
  double interpolateFieldAtQuadNode_sDerivative(std::vector<double> field, int iNode);
  double interpolateFieldAtQuadNode_xDerivative(std::vector<double> field, int iElm, int iNode);
  double interpolateFieldAtQuadNode_yDerivative(std::vector<double> field, int iElm, int iNode);

  double interpolateSolution(double r[3]);
  double interpolateSolution_rDerivative(double r[3]);
  double interpolateSolution_sDerivative(double r[3]);
  double interpolateSolution(int iElm, std::vector<double> &x);
  double interpolateSolution(feNumber *number, feSolution *sol, std::vector<double> &x);
  void interpolateSolution_gradrs(feNumber *number, feSolution *sol, std::vector<double> &x,
                                  std::vector<double> &grad);
  double interpolateSolution_xDerivative(int iElm, std::vector<double> &x);
  double interpolateSolution_yDerivative(int iElm, std::vector<double> &x);

  double interpolateSolutionAtQuadNode(int iNode);
  double interpolateSolutionAtQuadNode(int iElm, int iNode);
  double interpolateSolutionAtQuadNodeWithPhysicalBasisFunctions(int iElm, std::vector<double> &x);
  double interpolateSolutionAtQuadNode_rDerivative(int iNode);
  double interpolateSolutionAtQuadNode_sDerivative(int iNode);
  double interpolateSolutionAtQuadNode_xDerivative(int iElm, int iNode);
  double interpolateSolutionAtQuadNode_yDerivative(int iElm, int iNode);

  double interpolateSolutionDotAtQuadNode(int iNode);

  void interpolateVectorField(std::vector<double> field, double r[3], std::vector<double> &res);
  void interpolateVectorField_rDerivative(std::vector<double> field, double r[3],
                                          std::vector<double> &res);
  void interpolateVectorField_sDerivative(std::vector<double> field, double r[3],
                                          std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode(std::vector<double> field, int iNode,
                                        std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode,
                                                    std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> field, int iNode,
                                                    std::vector<double> &res);

  std::vector<double> &getSolutionReference() { return _sol; };
  std::vector<double> &getSolutionReferenceDot() { return _soldot; };
  void printL();
  void printdLdr();
};

// FESpace pour interpolant de Lagrange 1D de degre 0
class feSpace1DP0 : public feSpace
{
protected:
public:
  // Pour la geometrie : fieldID "GEO" est reserve
  feSpace1DP0(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
  {
    _nFunctions = 1;
    _nQuad = 1;
    _Lcoor = {1., 0., 0.};
  };
  // Pour la resolution
  feSpace1DP0(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct)
  {
    _nFunctions = 1;
    _adr.resize(_nFunctions);
    _nQuad = 1;
    _Lcoor = {1., 0., 0.};
  };
  virtual ~feSpace1DP0() {}

  virtual int getNbFunctions() { return 1; }
  virtual int getPolynomialDegree() { return 0; }
  virtual std::vector<double> L(double r[3]) { return {1.}; };
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };
  virtual std::vector<double> dLdr(double r[3]) { return {0.}; };
  virtual std::vector<double> dLds(double r[3]) { return {0.}; };
  virtual std::vector<double> dLdt(double r[3]) { return {0.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr);
};

// FESpace pour interpolant de Lagrange 1D de degre 1
class feSpace1DP1 : public feSpace
{
protected:
public:
  feSpace1DP1(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
  {
    _nFunctions = 2;
    _Lcoor = {-1., 0., 0., 1., 0., 0.};
  };
  feSpace1DP1(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct)
  {
    _nFunctions = 2;
    _adr.resize(_nFunctions);
    _Lcoor = {-1., 0., 0., 1., 0., 0.};
  };
  virtual ~feSpace1DP1() {}

  virtual int getNbFunctions() { return 2; }
  virtual int getPolynomialDegree() { return 1; }
  virtual std::vector<double> L(double r[3]) { return {(1. - r[0]) / 2., (1. + r[0]) / 2.}; };
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };
  virtual std::vector<double> dLdr(double r[3]) { return {-1. / 2., 1. / 2.}; };
  virtual std::vector<double> dLds(double r[3]) { return {0., 0.}; };
  virtual std::vector<double> dLdt(double r[3]) { return {0., 0.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr);
};

class feSpace1DP1_nonConsistant : public feSpace
{
protected:
public:
  feSpace1DP1_nonConsistant(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
  {
    _nFunctions = 1;
    _Lcoor = {0., 0., 0.};
  };
  feSpace1DP1_nonConsistant(feMesh *mesh, std::string fieldID, std::string cncGeoID,
                            feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct)
  {
    _nFunctions = 1;
    _adr.resize(_nFunctions);
    _Lcoor = {0., 0., 0.};
  };
  virtual ~feSpace1DP1_nonConsistant() {}

  virtual int getNbFunctions() { return 1; }
  virtual int getPolynomialDegree() { return 0; }
  virtual std::vector<double> L(double r[3]) { return {1.}; };
  virtual std::vector<double> dLdr(double r[3]) { return {0.}; };
  virtual std::vector<double> dLds(double r[3]) { return {0.}; };
  virtual std::vector<double> dLdt(double r[3]) { return {0.}; };

  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr);
};

// FESpace pour interpolant de Lagrange 1D de degre 2
class feSpace1DP2 : public feSpace
{
protected:
public:
  feSpace1DP2(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
  {
    _nFunctions = 3;
    _Lcoor = {-1., 0., 0., 1., 0., 0., 0., 0., 0.};
  };
  feSpace1DP2(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct)
  {
    _nFunctions = 3;
    _adr.resize(_nFunctions);
    _Lcoor = {-1., 0., 0., 1., 0., 0., 0., 0., 0.};
  };
  virtual ~feSpace1DP2() {}

  virtual int getNbFunctions() { return 3; }
  virtual int getPolynomialDegree() { return 2; }
  virtual std::vector<double> L(double r[3])
  {
    return {-r[0] * (1. - r[0]) / 2., r[0] * (1. + r[0]) / 2., -(r[0] + 1.) * (r[0] - 1.)};
  };
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };
  virtual std::vector<double> dLdr(double r[3])
  {
    return {(2. * r[0] - 1.) / 2., (2. * r[0] + 1.) / 2., -2. * r[0]};
  };
  virtual std::vector<double> dLds(double r[3]) { return {0., 0., 0.}; };
  virtual std::vector<double> dLdt(double r[3]) { return {0., 0., 0.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr);
};

// FESpace pour interpolant de Lagrange 1D de degre 3
class feSpace1DP3 : public feSpace
{
protected:
public:
  feSpace1DP3(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
  {
    _nFunctions = 4;
    _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 3., 0., 0., 1. / 3., 0., 0.}; // TODO : écrire en long
                                                                           // ?
  };
  feSpace1DP3(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct)
  {
    _nFunctions = 4;
    _adr.resize(_nFunctions);
    _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 3., 0., 0., 1. / 3., 0., 0.};
  };
  virtual ~feSpace1DP3() {}

  virtual int getNbFunctions() { return 4; }
  virtual int getPolynomialDegree() { return 3; }
  virtual std::vector<double> L(double r[3])
  {
    return {-9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.),
            9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.),
            27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.),
            -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.)};
  };
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };
  virtual std::vector<double> dLdr(double r[3])
  {
    return {r[0] * (9. / 8.) - r[0] * r[0] * (27. / 16.) + 1. / 16.,
            r[0] * (9. / 8.) + r[0] * r[0] * (27. / 16.) - 1. / 16.,
            r[0] * (-9. / 8.) + r[0] * r[0] * (81. / 16.) - 27. / 16.,
            r[0] * (-9. / 8.) - r[0] * r[0] * (81. / 16.) + 27. / 16.};
  };
  virtual std::vector<double> dLds(double r[3]) { return {0., 0., 0., 0.}; };
  virtual std::vector<double> dLdt(double r[3]) { return {0., 0., 0., 0.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr);
};

// FESpace pour interpolant de Lagrange 1D de degre 4
class feSpace1DP4 : public feSpace
{
protected:
public:
  feSpace1DP4(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
  {
    _nFunctions = 5;
    _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 2., 0., 0., 0., 0., 0., 1. / 2., 0., 0.};
  };
  feSpace1DP4(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct)
  {
    _nFunctions = 5;
    _adr.resize(_nFunctions);
    _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 2., 0., 0., 0., 0., 0., 1. / 2., 0., 0.};
  };
  virtual ~feSpace1DP4() {}

  virtual int getNbFunctions() { return 5; }
  virtual int getPolynomialDegree() { return 4; }
  virtual std::vector<double> L(double r[3])
  {
    return {r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
            r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
            r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (-8.0 / 3.0),
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0,
            r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (-8.0 / 3.0)};
  };
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    printf("Not implemented\n");
    exit(-1);
  };
  virtual std::vector<double> dLdr(double r[3])
  {
    return {r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (2.0 / 3.0) +
              r[0] * (r[0] - 1.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
              r[0] * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
              (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
            r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (2.0 / 3.0) +
              r[0] * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
              r[0] * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
              (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
            r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (-8.0 / 3.0) -
              r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0) -
              r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0) -
              (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0),
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * 4.0 +
              (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * 4.0 +
              (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0 +
              (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0,
            r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (-8.0 / 3.0) -
              r[0] * (r[0] - 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0) -
              r[0] * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0) -
              (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0)};
  };
  virtual std::vector<double> dLds(double r[3]) { return {0., 0., 0., 0., 0.}; };
  virtual std::vector<double> dLdt(double r[3]) { return {0., 0., 0., 0., 0.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
  virtual void initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr);
};

#endif
