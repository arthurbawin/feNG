#ifndef _FESPACETRIANGLE_
#define _FESPACETRIANGLE_

#include "feSpace.h"

class feSpaceTriP1 : public feSpace {
protected:
public:
  feSpaceTriP1(std::string cncGeoID);
  feSpaceTriP1(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct);
  ~feSpaceTriP1() {}

  virtual int getNbFunctions() { return 3; }
  virtual int getPolynomialDegree() { return 1; }
  virtual std::vector<double> L(double r[3]);
  virtual std::vector<double> dLdr(double r[3]);
  virtual std::vector<double> dLds(double r[3]);
  virtual std::vector<double> dLdt(double r[3]);

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

class feSpaceTriP1_nonConsistant : public feSpace {
protected:
public:
  feSpaceTriP1_nonConsistant(std::string cncGeoID);
  feSpaceTriP1_nonConsistant(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct);
  ~feSpaceTriP1_nonConsistant() {}

  virtual int getNbFunctions() { return 3; }
  virtual int getPolynomialDegree() { return 1; }
  virtual std::vector<double> L(double r[3]);
  virtual std::vector<double> dLdr(double r[3]);
  virtual std::vector<double> dLds(double r[3]);
  virtual std::vector<double> dLdt(double r[3]);

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

class feSpaceTriP2 : public feSpace {
protected:
public:
  feSpaceTriP2(std::string cncGeoID);
  feSpaceTriP2(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct);
  ~feSpaceTriP2() {}

  virtual int getNbFunctions() { return 6; }
  virtual int getPolynomialDegree() { return 2; }
  virtual std::vector<double> L(double r[3]);
  virtual std::vector<double> dLdr(double r[3]);
  virtual std::vector<double> dLds(double r[3]);
  virtual std::vector<double> dLdt(double r[3]);

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

class feSpaceTriP3 : public feSpace {
protected:
public:
  feSpaceTriP3(std::string cncGeoID);
  feSpaceTriP3(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct);
  ~feSpaceTriP3() {}

  virtual int getNbFunctions() { return 10; }
  virtual int getPolynomialDegree() { return 3; }
  virtual std::vector<double> L(double r[3]);
  virtual std::vector<double> dLdr(double r[3]);
  virtual std::vector<double> dLds(double r[3]);
  virtual std::vector<double> dLdt(double r[3]);

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

class feSpaceTriP4 : public feSpace {
protected:
public:
  feSpaceTriP4(std::string cncGeoID);
  feSpaceTriP4(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct);
  ~feSpaceTriP4() {}

  virtual int getNbFunctions() { return 15; }
  virtual int getPolynomialDegree() { return 4; }
  virtual std::vector<double> L(double r[3]);
  virtual std::vector<double> dLdr(double r[3]);
  virtual std::vector<double> dLds(double r[3]);
  virtual std::vector<double> dLdt(double r[3]);

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};

#endif