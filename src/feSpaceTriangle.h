#ifndef _FESPACETRIANGLE_
#define _FESPACETRIANGLE_

#include "feSpace.h"

class feSpaceTriP1 : public feSpace{
protected:
public:
  feSpaceTriP1(std::string cncGeoID) 
    : feSpace(nullptr, "GEO", cncGeoID, nullptr){
    _nFunctions = 3;
    _Lcoor = {0., 0., 0.,
              1., 0., 0.,
              0., 1., 0.};
  };
  feSpaceTriP1(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct) 
    : feSpace(mesh, fieldID, cncGeoID, fct){
    _nFunctions = 3;
    _adr.resize(_nFunctions);
    _Lcoor = {0., 0., 0.,
              1., 0., 0.,
              0., 1., 0.};
  };
  ~feSpaceTriP1() {}

  virtual int getNbFunctions(){ return 3; }
  virtual std::vector<double>    L(double r[3]);
  virtual std::vector<double> dLdr(double r[3]);
  virtual std::vector<double> dLds(double r[3]);
  virtual std::vector<double> dLdt(double r[3]);

  virtual void initializeNumberingUnknowns(feNumber *number);
  virtual void initializeNumberingEssential(feNumber *number);
  virtual void initializeAddressingVector(feNumber *number, int numElem);
};


#endif