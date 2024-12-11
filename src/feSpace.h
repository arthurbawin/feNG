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

#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EigenMat;

// Supported choices of degrees of freedom to define finite element spaces
enum class elementType {
  LAGRANGE,
  VECTOR_LAGRANGE,
  DG_LAGRANGE,
  LEGENDRE,
  HERMITE,
  CROUZEIX_RAVIART
};

enum class dofLocation {
  VERTEX,
  EDGE,
  ELEMENT
};

// How the degrees of freedom should be initialized from the provided field f:
//
//  PREVIOUS_SOL: dofs are not initialized (use only when using solution from a
//                previous simulation)
//      NODEWISE: dofs are initialized such that uh(x_i) = f(x_i)
//
//                                               /           /
// LEAST_SQUARES: dofs are initialized such that | uh v dx = | f v dx for all v
//                                               /           /
enum class dofInitialization { PREVIOUS_SOL, NODEWISE, LEAST_SQUARES, EXTRAPOLATED_EULER_0D };

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
//                 element: the choice of finite element and degrees of freedom (see enum above)
//                  degree: degree (order) of the polynomial interpolation
//               fieldName: name of the field for which the space is defined
//        connectivityName: name of the geometric connectivity on which the FE space is defined
//        degreeQuadrature: maximum degree of the polynomial to integrate perfectly with
//                          a quadrature rule
//                     fct: field used to initialize the degrees of freedom.
//                          Default initialization method is NODEWISE, i.e. uh(x_i) = f(x_i)
// useGlobalShapeFunctions: choose between local (reference space) or global (physical space)
//                          shape functions (global still experimental, do not use yet).
feStatus createFiniteElementSpace(feSpace *&space, feMesh *mesh, const elementType element,
                                  const int degree, const std::string fieldName,
                                  const std::string connectivityName, const int degreeQuadrature,
                                  void *fct, const bool useGlobalShapeFunctions = false);

// Abstract class defining a finite element space, handling the interpolation
// functions defined over the elements. A finite element space can be used
// to interpolate either the geometry or the fields of interest.
//
// An FE space is defined on a single geometric connectivity (which
// corresponds to a Physical Entity in the mesh).
class feSpace
{
protected:
  // The type of finite element
  elementType _elementType;
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

  bool _isGeometricInterpolant = false;

  bool _ignoreNegativeJacobianWarning = false;

  int _dim;

  int _nComponents = 1;
  bool _essentialComponents[3] = {false, false, false};

  // Quadrature rule on the reference element:
  // number of quadrature nodes, weights and coordinates
  int _nQuad;
  std::vector<double> _wQuad;
  std::vector<double> _xQuad;
  std::vector<double> _yQuad;
  std::vector<double> _zQuad;

  // Number of interpolation functions (i.e. number of DOFs)
  int _nFunctions;
  // rst-Coordinates of the collocation nodes on the reference element
  // Layout: [r1 s1 t1 r2 s2 t2 ... rn sn tn]
  std::vector<double> _Lcoor;
  // Evaluation of the (local) shape functions and their derivatives at
  // the quadrature nodes.
  std::vector<double> _L;
  std::vector<double> _dLdr;
  std::vector<double> _dLds;
  std::vector<double> _dLdt;
  std::vector<double> _d2Ldr2;
  std::vector<double> _d2Ldrs;
  std::vector<double> _d2Lds2;
  std::vector<double> _d2Ldt2;

  // Barycentric (linear) coordinates at quadrature nodes
  std::vector<double> _barycentricCoordinates;

  // Use global (physical) interpolation functions (experimental, do not use yet)
  bool _useGlobalShapeFunctions = false;
  // Global shape functions and derivatives evaluated at quadrature points on each element
  std::vector<std::vector<double> > _Lglob;
  std::vector<std::vector<double> > _dLdxglob;
  std::vector<std::vector<double> > _dLdyglob;

  // Scalar or vector field used to initialize the DOFs
  feFunction *_scalarFct;
  feVectorFunction *_vectorFct;
  // How the DOFs should be initialized, see enum above
  dofInitialization _DOFinitialization = dofInitialization::NODEWISE;
  std::vector<dofLocation> _dofLocations;

  // Ptr to the numbering of the field associated to this space
  feNumber *_numbering;
  // Ptr to the connectivity on which the space is defined
  feCncGeo *_cnc;

public:
  friend class feBilinearForm;
  friend class feNorm;
  friend class feCncGeo;

  // Do not use the abstract class constructor directly, call the derived classes.
  feSpace(const int dimension, feMesh *mesh = nullptr, const std::string &fieldID = "",
          const std::string &cncGeoID = "", feFunction *scalarField = nullptr,
          feVectorFunction *vectorField = nullptr, const bool useGlobalShapeFunctions = false);
  virtual ~feSpace();

  const std::string &getFieldID() const { return _fieldID; }
  int getFieldTag() { return _fieldTag; }
  const std::string &getCncGeoID() const { return _cncGeoID; }
  int getCncGeoTag() { return _cncGeoTag; }

  // Assign the tag of a geometric space after the mesh has been created
  void setCncGeoTag(int tag) { _cncGeoTag = tag; }
  // Assign the pointer to the connectivity after it has been created
  void setCncPtr(feCncGeo *cnc) { _cnc = cnc; }
  // Assign the mesh pointer of a geometric space after the mesh has been created
  void setMeshPtr(feMesh *mesh) { _mesh = mesh; }
  feMesh *getMeshPtr() { return _mesh; }
  // Assign the pointer to the field numbering after the numbering has been created
  void setNumberingPtr(feNumber *numbering) { _numbering = numbering; }

  // Return the number of field components (1 for scalar finite element, 1/2/3 for vector element)
  int getNumComponents() const { return _nComponents; };

  void setEssentialComponent(int iComponent, bool flag)
  {
    _essentialComponents[iComponent] = flag;
  };
  bool isEssentialComponent(int iComponent) { return _essentialComponents[iComponent]; };

  // Return the attributes of the geometric connectivity on which the space is defined
  feCncGeo *getCncGeo() { return _cnc; }
  int getDim() { return _dim; };
  int getNumElements();
  int getNumVerticesPerElem();

  elementType getElementType() { return _elementType; };
  dofInitialization getDOFInitialization() { return _DOFinitialization; };
  void setDOFInitialization(dofInitialization init) { _DOFinitialization = init; };

  const std::vector<dofLocation> &getDOFLocations() const { return _dofLocations; }
  dofLocation getDOFLocation(const int iDOF) const { return _dofLocations[iDOF]; }

  bool useGlobalFunctions() { return _useGlobalShapeFunctions; }
  void useGlobalFunctions(bool flag) { _useGlobalShapeFunctions = flag; }

  // Return the number of shape functions (DOF) on an element
  virtual int getNumFunctions() const = 0;
  // Return highest degree of the polynomial basis
  virtual int getPolynomialDegree() = 0;

  EigenMat innerProductBasisFunctions(int iElm);
  double innerProductBasisFunctions(int iElm, int ex, int ey);

  const std::vector<double> &getLcoor() const { return _Lcoor; }

  const std::vector<double> &getBarycentricCoordinatesAtQuadNode() { return _barycentricCoordinates; }

  // Evaluate the local shape functions and derivatives at quadrature nodes
  virtual std::vector<double> L(double *r) = 0;
  virtual void L(double *r, double *L) = 0;
  virtual std::vector<double> dLdr(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };
  virtual std::vector<double> dLds(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };
  virtual std::vector<double> dLdt(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };
  virtual std::vector<double> d2Ldr2(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };
  virtual std::vector<double> d2Ldrs(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };
  virtual std::vector<double> d2Lds2(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };
  virtual std::vector<double> d2Ldt2(double *r) { UNUSED(r); return std::vector<double>(_nFunctions, 0.); };

  // Evaluate the global shape functions at quadrature nodes
  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
  {
    UNUSED(iElm, x, L, dLdx, dLdy); 
    printf("Not implemented\n");
    exit(-1);
  };

  // Return value of shape function (or derivatives) at quadrature node
  void getFunctionsAtQuadNode(const int iQuadNode, std::vector<double> &phi);
  double getFunctionAtQuadNode(const int iFun, const int iQuadNode)
  {
    return _L[_nFunctions * iQuadNode + iFun];
  }
  double getdFunctiondrAtQuadNode(const int iFun, const int iQuadNode)
  {
    return _dLdr[_nFunctions * iQuadNode + iFun];
  }
  double getdFunctiondsAtQuadNode(const int iFun, const int iQuadNode)
  {
    return _dLds[_nFunctions * iQuadNode + iFun];
  }
  double getdFunctiondtAtQuadNode(const int iFun, const int iQuadNode)
  {
    return _dLdt[_nFunctions * iQuadNode + iFun];
  }
  double getd2Functiondr2AtQuadNode(const int iFun, const int iQuadNode)
  {
    return _d2Ldr2[_nFunctions * iQuadNode + iFun];
  }
  double getd2FunctiondrsAtQuadNode(const int iFun, const int iQuadNode)
  {
    return _d2Ldrs[_nFunctions * iQuadNode + iFun];
  }
  double getd2Functionds2AtQuadNode(const int iFun, const int iQuadNode)
  {
    return _d2Lds2[_nFunctions * iQuadNode + iFun];
  }

  // Get the gradient in physical coordinates of the basis functions at quadrature node:
  // gradPhi = [dPhi1dx ... dPhindx, dPhi1dy ... dPhindy, dPhi1dz ... dPhindz]
  // Size of gradPhi = dim x nFunctions
  void getFunctionsPhysicalGradientAtQuadNode(const int iQuadNode, const ElementTransformation &T,
                                              double *gradPhi);

  void getFunctionsPhysicalHessianAtQuadNode(const int iQuadNode, const ElementTransformation &T,
                                              double *hessPhi);

  double getGlobalFunctionAtQuadNode(const int iElm, const int iFun, const int iQuadNode)
  {
    return _Lglob[iElm][_nFunctions * iQuadNode + iFun];
  }
  double getdGlobalFunctiondxAtQuadNode(const int iElm, const int iFun, const int iQuadNode)
  {
    return _dLdxglob[iElm][_nFunctions * iQuadNode + iFun];
  }
  double getdGlobalFunctiondyAtQuadNode(const int iElm, const int iFun, const int iQuadNode)
  {
    return _dLdyglob[iElm][_nFunctions * iQuadNode + iFun];
  }

  // Set and get quadrature rule
  feStatus setQuadratureRule(feQuadrature *quadratureRule);
  int getNumQuadPoints() { return _nQuad; }
  std::vector<double> &getQuadratureWeights() { return _wQuad; }
  std::vector<double> &getRQuadraturePoints() { return _xQuad; }
  std::vector<double> &getSQuadraturePoints() { return _yQuad; }
  std::vector<double> &getTQuadraturePoints() { return _zQuad; }

  virtual void initializeNumberingUnknowns() = 0;
  virtual void initializeNumberingEssential() = 0;
  virtual void initializeAddressingVector(int numElem, std::vector<feInt> &adr) = 0;

  // Evaluate prescribed scalar field
  double evalFun(const double t, const std::vector<double> &x) { return _scalarFct->eval(t, x); }
  void evalFun(const double t, const std::vector<double> &x, std::vector<double> &res)
  {
    _vectorFct->eval(t, x, res);
  }

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

  double interpolateField(feSolution *sol, std::vector<double> &x);
  void interpolateField_gradrs(feSolution *sol, std::vector<double> &x, std::vector<double> &grad);

  // Interpolate vector field or derivatives at reference node r = [r,s,t]
  // using local shape functions (default). Result is stored in res.
  void interpolateVectorField(std::vector<double> &field, double *r, std::vector<double> &res);
  void interpolateVectorField(std::vector<double> &field, int nComponents, double *r,
                              std::vector<double> &res);
  double interpolateVectorFieldComponent(std::vector<double> &field, int iComponent, double *r);
  void interpolateVectorField_rDerivative(std::vector<double> &field, double *r,
                                          std::vector<double> &res);
  void interpolateVectorField_sDerivative(std::vector<double> &field, double *r,
                                          std::vector<double> &res);

  //
  // Interpolation of fields and derivatives at quadrature node:
  //

  // Interpolate scalar field or derivatives at iNode-th quadrature node using local shape functions
  // (default)
  double interpolateFieldAtQuadNode(std::vector<double> &field, int iNode);
  double interpolateFieldAtQuadNode_rDerivative(std::vector<double> &field, int iNode);
  double interpolateFieldAtQuadNode_sDerivative(std::vector<double> &field, int iNode);
  double interpolateFieldAtQuadNode_rrDerivative(std::vector<double> &field, int iNode);
  double interpolateFieldAtQuadNode_ssDerivative(std::vector<double> &field, int iNode);

  void interpolateFieldAtQuadNode_physicalGradient(std::vector<double> &field, const int iQuadNode,
                                                   const ElementTransformation &T, double *grad);

  // Interpolate scalar field or derivatives at iNode-th quadrature node using global shape
  // functions
  double interpolateFieldAtQuadNode(std::vector<double> &field, int iElm, int iNode);
  double interpolateFieldAtQuadNode_xDerivative(std::vector<double> &field, int iElm, int iNode);
  double interpolateFieldAtQuadNode_yDerivative(std::vector<double> &field, int iElm, int iNode);

  // Interpolate vector field or derivatives at iNode-th quadrature node using local shape functions
  // (default)

  // Interpolate vector valued function using scalar valued FE space
  void interpolateVectorFieldAtQuadNode(std::vector<double> &field, int iNode,
                                        std::vector<double> &res);
  // Interpolate vector valued function using vector valued FE space
  void interpolateVectorFieldAtQuadNode(std::vector<double> &field, int iNode,
                                        std::vector<double> &res, int nComponents);
  // Interpolate scalar component of vector valued function using vector valued FE space
  double interpolateVectorFieldComponentAtQuadNode(std::vector<double> &field, int iNode, int iComponent);
  double interpolateVectorFieldComponentAtQuadNode_fullField(std::vector<double> &field, int iNode, int iComponent);

  void interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> &field, int iNode,
                                                    std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> &field, int iNode,
                                                    double res[3]);
  void interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> &field, int iNode,
                                                    std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> &field, int iNode,
                                                    double res[3]);
  void interpolateVectorFieldAtQuadNode_tDerivative(std::vector<double> &field, int iNode,
                                                    std::vector<double> &res);
  void interpolateVectorFieldAtQuadNode_tDerivative(std::vector<double> &field, int iNode,
                                                    double res[3]);

  // FIXME: Add comments
  void interpolateVectorFieldAtQuadNode_physicalGradient(std::vector<double> &field,
                                                         const int nComponents, const int iQuadNode,
                                                         const ElementTransformation &T,
                                                         double *grad);

  void interpolateVectorFieldAtQuadNode_physicalHessian(std::vector<double> &field,
                                                         const int nComponents, const int iQuadNode,
                                                         const ElementTransformation &T,
                                                         double *hessian);
};

class feScalarSpace : public feSpace
{
protected:
public:
  feScalarSpace(const int dimension, feMesh *mesh = nullptr, const std::string &fieldID = "",
                const std::string &cncGeoID = "", feFunction *scalarField = nullptr,
                const bool useGlobalShapeFunctions = false);
  ~feScalarSpace(){};
};

class feVectorSpace : public feSpace
{
protected:
public:
  feVectorSpace(const int dimension, feMesh *mesh = nullptr, const std::string &fieldID = "",
                const std::string &cncGeoID = "", feVectorFunction *vectorField = nullptr,
                const bool useGlobalShapeFunctions = false);
  ~feVectorSpace(){};
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 0 on 0D point
// -----------------------------------------------------------------------------
class feSpace0DP0 : public feScalarSpace
{
public:
  feSpace0DP0(const std::string &cncGeoID);
  feSpace0DP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace0DP0(){};

  int getNumFunctions() const { return 1; }
  int getPolynomialDegree() { return 0; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// 0D Hermite point to fix boundary conditions
// -----------------------------------------------------------------------------
class feSpace0D_Hermite : public feScalarSpace
{
public:
  feSpace0D_Hermite(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                    feFunction *fct);
  ~feSpace0D_Hermite(){};

  int getNumFunctions() const { return 2; }
  int getPolynomialDegree() { return 0; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 0 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1DP0 : public feScalarSpace
{
public:
  feSpace1DP0(const std::string &cncGeoID);
  feSpace1DP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP0(){};

  int getNumFunctions() const { return 1; }
  int getPolynomialDegree() { return 0; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 1 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1DP1 : public feScalarSpace
{
protected:
public:
  feSpace1DP1(std::string cncGeoID);
  feSpace1DP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP1(){};

  int getNumFunctions() const { return 2; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 1 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecP1 : public feVectorSpace
{
protected:
public:
  feSpaceVecP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
               feVectorFunction *fct);
  ~feSpaceVecP1(){};

  int getNumFunctions() const { return 2 * dim; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Discontinuous Lagrange element of degree 1 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1D_DG_P1 : public feScalarSpace
{
protected:
public:
  feSpace1D_DG_P1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                  feFunction *fct);
  ~feSpace1D_DG_P1(){};

  int getNumFunctions() const { return 2; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Crouzeix-Raviart element of degree 0 on line
// Trivial constant function L = 1
// -----------------------------------------------------------------------------
class feSpace1D_CR0 : public feScalarSpace
{
protected:
public:
  feSpace1D_CR0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                feFunction *fct);
  ~feSpace1D_CR0() {}

  int getNumFunctions() const { return 1; }
  int getPolynomialDegree() { return 0; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 2 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1DP2 : public feScalarSpace
{
protected:
public:
  feSpace1DP2(std::string cncGeoID);
  feSpace1DP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP2() {}

  int getNumFunctions() const { return 3; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 2 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecP2 : public feVectorSpace
{
protected:
public:
  feSpaceVecP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
               feVectorFunction *fct);
  ~feSpaceVecP2(){};

  int getNumFunctions() const { return 3 * dim; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 3 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1DP3 : public feScalarSpace
{
protected:
public:
  feSpace1DP3(std::string cncGeoID);
  feSpace1DP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP3() {}

  int getNumFunctions() const { return 4; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 3 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecP3 : public feVectorSpace
{
protected:
public:
  feSpaceVecP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
               feVectorFunction *fct);
  ~feSpaceVecP3(){};

  int getNumFunctions() const { return 4 * dim; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Hermite element of degree 3 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1D_H3 : public feScalarSpace
{
protected:
public:
  feSpace1D_H3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
               feFunction *fct);
  ~feSpace1D_H3() {}

  int getNumFunctions() const { return 4; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> d2Ldr2(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 4 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1DP4 : public feScalarSpace
{
protected:
public:
  feSpace1DP4(std::string cncGeoID);
  feSpace1DP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct);
  ~feSpace1DP4() {}

  int getNumFunctions() const { return 5; }
  int getPolynomialDegree() { return 4; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 4 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecP4 : public feVectorSpace
{
protected:
public:
  feSpaceVecP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
               feVectorFunction *fct);
  ~feSpaceVecP4(){};

  int getNumFunctions() const { return 5 * dim; }
  int getPolynomialDegree() { return 4; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Legendre interpolation functions of arbitrary degree on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
class feSpace1D_Legendre : public feScalarSpace
{
protected:
  int _degree;

public:
  feSpace1D_Legendre(const int degree, feMesh *mesh, const std::string fieldID,
                     const std::string cncGeoID, feFunction *fct);
  ~feSpace1D_Legendre() {}

  int getNumFunctions() const { return _nFunctions; }
  int getPolynomialDegree() { return _degree; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
class feSpaceTriP1 : public feScalarSpace
{
protected:
public:
  feSpaceTriP1(std::string cncGeoID);
  feSpaceTriP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP1() {}

  int getNumFunctions() const { return 3; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> d2Ldr2(double *r);
  std::vector<double> d2Ldrs(double *r);
  std::vector<double> d2Lds2(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecTriP1 : public feVectorSpace
{
protected:
public:
  feSpaceVecTriP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                  feVectorFunction *fct, const bool useGlobalShapeFunctions = false);
  ~feSpaceVecTriP1() {}

  int getNumFunctions() const { return 3 * dim; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> d2Ldr2(double *r);
  std::vector<double> d2Ldrs(double *r);
  std::vector<double> d2Lds2(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Non-conforming Crouzeix-Raviart element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
class feSpaceTri_CR1 : public feScalarSpace
{
protected:
public:
  feSpaceTri_CR1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                 feFunction *fct);
  ~feSpaceTri_CR1() {}

  int getNumFunctions() const { return 3; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> d2Ldr2(double *r);
  std::vector<double> d2Ldrs(double *r);
  std::vector<double> d2Lds2(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 2 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
class feSpaceTriP2 : public feScalarSpace
{
protected:
public:
  feSpaceTriP2(std::string cncGeoID);
  feSpaceTriP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP2() {}

  int getNumFunctions() const { return 6; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> d2Ldr2(double *r);
  std::vector<double> d2Ldrs(double *r);
  std::vector<double> d2Lds2(double *r);

  virtual feStatus Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 2 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecTriP2 : public feVectorSpace
{
public:
  feSpaceVecTriP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                  feVectorFunction *fct, const bool useGlobalShapeFunctions = false);
  ~feSpaceVecTriP2() {}

  int getNumFunctions() const { return 6 * dim; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> d2Ldr2(double *r);
  std::vector<double> d2Ldrs(double *r);
  std::vector<double> d2Lds2(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Non-conforming Crouzeix-Raviart element of degree 2 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
class feSpaceTri_CR2 : public feScalarSpace
{
protected:
public:
  feSpaceTri_CR2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                 feFunction *fct);
  ~feSpaceTri_CR2() {}

  int getNumFunctions() const { return 7; }
  int getPolynomialDegree() { return 2; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 3 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
class feSpaceTriP3 : public feScalarSpace
{
protected:
public:
  feSpaceTriP3(std::string cncGeoID);
  feSpaceTriP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP3() {}

  int getNumFunctions() const { return 10; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 3 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecTriP3 : public feVectorSpace
{
public:
  feSpaceVecTriP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                  feVectorFunction *fct, const bool useGlobalShapeFunctions = false);
  ~feSpaceVecTriP3() {}

  int getNumFunctions() const { return 10 * dim; }
  int getPolynomialDegree() { return 3; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 4 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
class feSpaceTriP4 : public feScalarSpace
{
protected:
public:
  feSpaceTriP4(std::string cncGeoID);
  feSpaceTriP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTriP4() {}

  int getNumFunctions() const { return 15; }
  int getPolynomialDegree() { return 4; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);
  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 4 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim> class feSpaceVecTriP4 : public feVectorSpace
{
public:
  feSpaceVecTriP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                  feVectorFunction *fct, const bool useGlobalShapeFunctions = false);
  ~feSpaceVecTriP4() {}

  int getNumFunctions() const { return 15 * dim; }
  int getPolynomialDegree() { return 4; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

// -----------------------------------------------------------------------------
// Lagrange element of degree 1 on reference tetrahedron
// -----------------------------------------------------------------------------
class feSpaceTetP1 : public feScalarSpace
{
protected:
public:
  feSpaceTetP1(std::string cncGeoID);
  feSpaceTetP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
               const bool useGlobalShapeFunctions = false);
  ~feSpaceTetP1() {}

  int getNumFunctions() const { return 4; }
  int getPolynomialDegree() { return 1; }

  std::vector<double> L(double *r);
  void L(double *r, double *L);

  std::vector<double> dLdr(double *r);
  std::vector<double> dLds(double *r);
  std::vector<double> dLdt(double *r);
  std::vector<double> d2Ldr2(double *r);
  std::vector<double> d2Ldrs(double *r);
  std::vector<double> d2Lds2(double *r);
  // std::vector<double> d2Ldt2(double *r);
  // std::vector<double> d2Ldrt(double *r);
  // std::vector<double> d2Ldst(double *r);

  void initializeNumberingUnknowns();
  void initializeNumberingEssential();
  void initializeAddressingVector(int numElem, std::vector<feInt> &adr);
};

#endif
