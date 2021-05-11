#ifndef _FESPACE_
#define _FESPACE_

#include <vector>
#include <string>
#include <iostream>

#include "feQuadrature.h"
#include "feFunction.h"

class feMesh;
class feNumber;
class feSolution;

/* Un feSpace est utilisé
  - soit pour définir les interpolants géométriques sur une partie d'un maillage (sur une cncGeo)
  - soit pour définir les interpolants utilisés pour la résolution (sur un field)
*/
class feSpace{

protected:
  feMesh *_mesh;         // Pointer to the mesh on which the space is defined

  std::string _fieldID;  // Name of the associated field, e.g "U", "V", "P"
  int _fieldTag;         // Number (tag) of the associated field
  std::string _cncGeoID; // Name of the geometric connectivity on which the space is defined
  int _cncGeoTag;        // Number (tag) of the geometric connectivity
	// Quadrature
	int _nQuad;
	std::vector<double> _wQuad;
	std::vector<double> _xQuad;
  // Interpolants
	int _nFunctions; // nielm
  std::vector<double> _Lcoor;
	std::vector<double> _L;
  std::vector<double> _dLdr;

  std::vector<int>    _adr;
  std::vector<double> _sol;
  std::vector<double> _soldot;

  feFunction *_fct;

public:
	feSpace(class feMesh *mesh = nullptr, std::string fieldID = "", std::string cncGeoID = "", feFunction *fct = nullptr);
	virtual ~feSpace() {}

  int getDim();
  int getNbElm();
  int getNbNodePerElem();

  std::string getFieldID(){ return _fieldID; }
	int getFieldTag(){ return _fieldTag; }
  std::string getCncGeoID(){ return _cncGeoID; }
  int getCncGeoTag(){ return _cncGeoTag; }
  void setCncGeoTag(int tag){ _cncGeoTag = tag; } // Used to assign the tag of a geometric space after the mesh has been created

  virtual int getNbFunctions(){ return 0; }

  const std::vector<double>& getLcoor(){ return _Lcoor; }

  virtual std::vector<double>    L(double r[3]){ return {0.}; };
  virtual std::vector<double> dLdr(double r[3]){ return {0.}; };

  double getFunctionAtQuadNode(int iFun, int iQuadNode){ return _L[_nFunctions*iQuadNode+iFun]; }  
  double getdFunctiondrAtQuadNode(int iFun, int iQuadNode){ return _dLdr[_nFunctions*iQuadNode+iFun]; }  

  void initializeSolution(feSolution *sol);
  void initializeSolutionDot(feSolution *sol);

  virtual void initializeNumberingUnknowns(feNumber *number) {};
  virtual void initializeNumberingEssential(feNumber *number) {};
  virtual void initializeAddressingVector(feNumber *number, int numElem) {};

  int getAddressingVectorAt(int node){ return _adr[node]; }

  void setQuadratureRule(feQuadrature *quad);
  int getNbQuadPoints(){ return _nQuad; }
  std::vector<double> &getQuadratureWeights(){ return _wQuad; }
  std::vector<double> &getQuadraturePoints(){ return _xQuad; }

  double evalFun(const double t, const std::vector<double> &x){ return _fct->eval(t,x); }

  double interpolateField(std::vector<double> field, double r[3]);
  double interpolateFieldAtQuadNode(std::vector<double> field, int iNode);
  double interpolateSolutionAtQuadNode(int iNode);
  double interpolateSolutionDotAtQuadNode(int iNode);
  double interpolateField_rDerivative(std::vector<double> field, double r[3]);
  double interpolateFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode);
  double interpolateSolutionAtQuadNode_rDerivative(int iNode);

  void printL();
  void printdLdr();
};

// FESpace pour interpolant de Lagrange 1D de degre 0
class feSpace1DP0 : public feSpace{
protected:
public:
  // Pour la geometrie : fieldID "GEO" est reserve ?
  feSpace1DP0(std::string cncGeoID) 
    : feSpace(nullptr, "GEO", cncGeoID, nullptr){
    _nFunctions = 1;
    _nQuad = 1;
    _Lcoor = {1.};
  };
  // Pour la resolution
	feSpace1DP0(class feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
    : feSpace(mesh, fieldID, cncGeoID, fct){
    _nFunctions = 1;
    _adr.resize(_nFunctions);
    _nQuad = 1;
    _Lcoor = {1.};
	};
	virtual ~feSpace1DP0() {}

  // virtual std::string getMeshID(){ return _mesh->getID(); }
  virtual int getNbFunctions(){ return 1; }
  virtual std::vector<double>    L(double r[3]){ return {1.}; };
  virtual std::vector<double> dLdr(double r[3]){ return {0.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

// FESpace pour interpolant de Lagrange 1D de degre 1
class feSpace1DP1 : public feSpace{
protected:
public:
  feSpace1DP1(std::string cncGeoID) 
    : feSpace(nullptr, "GEO", cncGeoID, nullptr){
    _nFunctions = 2;
    _Lcoor = {-1., 1.}; // Coord entre 0 et 1 pour un P1 il me semble (même si ça change rien)
  };
  feSpace1DP1(class feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct) 
    : feSpace(mesh, fieldID, cncGeoID, fct){
    _nFunctions = 2;
    _adr.resize(_nFunctions);
    _Lcoor = {-1., 1.};
  };
  ~feSpace1DP1() {}

  virtual int getNbFunctions(){ return 2; }
  virtual std::vector<double>    L(double r[3]){ return {(1.-r[0])/2., (1.+r[0])/2.}; };
  virtual std::vector<double> dLdr(double r[3]){ return {      -1./2.,        1./2.}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

// FESpace pour interpolant de Lagrange 1D de degre 2
class feSpace1DP2 : public feSpace{
protected:
public:
  feSpace1DP2(std::string cncGeoID)
    : feSpace(nullptr, "GEO", cncGeoID, nullptr){
    _nFunctions = 3;
    _Lcoor = {-1., 1., 0.};
  };
  feSpace1DP2(class feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct) 
    : feSpace(mesh, fieldID, cncGeoID, fct){
    _nFunctions = 3;
    _adr.resize(_nFunctions);
    _Lcoor = {-1., 1., 0.};
  };
  ~feSpace1DP2() {}

  virtual int getNbFunctions(){ return 3; }
  virtual std::vector<double>    L(double r[3]){ return {-r[0]*(1.-r[0])/2. , r[0]*(1.+r[0])/2. , -(r[0]+1.)*(r[0]-1.)}; };
  virtual std::vector<double> dLdr(double r[3]){ return {   (2.*r[0]-1.)/2. ,   (2.*r[0]+1.)/2. ,             -2.*r[0]}; };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

// FESpace pour interpolant de Lagrange 1D de degre 3
class feSpace1DP3 : public feSpace{
protected:
public:
  feSpace1DP3(std::string cncGeoID)
    : feSpace(nullptr, "GEO", cncGeoID, nullptr){
    _nFunctions = 4;
    _Lcoor = {-1., 1., -1./3., 1./3.}; //TODO : écrire en long ?
  };
  feSpace1DP3(class feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct) 
    : feSpace(mesh, fieldID, cncGeoID, fct){
    _nFunctions = 4;
    _adr.resize(_nFunctions);
    _Lcoor =  {-1., 1., -1./3., 1./3.};
  };
  ~feSpace1DP3() {}

  virtual int getNbFunctions(){ return 4; }
  virtual std::vector<double>    L(double r[3]){
    return { -9./16.*(r[0]+1./3.)*(r[0]-1./3.)*(r[0]-1.), 
              9./16.*(r[0]+1./3.)*(r[0]-1./3.)*(r[0]+1.), 
             27./16.*(r[0]+1.)*(r[0]-1./3.)*(r[0]-1.),
            -27./16.*(r[0]+1.)*(r[0]+1./3.)*(r[0]-1.)}; 
  };
  virtual std::vector<double> dLdr(double r[3]){
    return { r[0]*( 9./8.) - r[0]*r[0]*(27./16.) +  1./16.,
             r[0]*( 9./8.) + r[0]*r[0]*(27./16.) -  1./16.,
             r[0]*(-9./8.) + r[0]*r[0]*(81./16.) - 27./16.,
             r[0]*(-9./8.) - r[0]*r[0]*(81./16.) + 27./16. };
  };

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};


#endif