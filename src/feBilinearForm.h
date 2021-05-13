#ifndef _FEBILINEARFORM_
#define _FEBILINEARFORM_

#include "feNG.h"
#include "feMesh.h"
#include "feSpace.h"
#include "feSysElm.h"
#include "feQuadrature.h"
#include "feNumber.h"
#include "feSolution.h"

class feBilinearForm{

protected:
  feSysElm *_sysElm;
  std::vector<feSpace*> _intSpace;
  
  int _cncGeoTag;
  feSpace *_geoSpace;
  int _nCoord;
  int _nGeoNodes;
  int _nGeoElm;
  std::vector<double> _geoCoord;

  int _nQuad;
  int _degQuad;
  std::vector<double> _w;
  std::vector<double> _x; 

  std::vector<int> _iVar;
  std::vector<int> _jVar;
  feInt _niElm;
  feInt _njElm;
  std::vector<int> _adrI;
  std::vector<int> _adrJ;

  double** _Ae;
  double* _Be;

public:
	feBilinearForm(std::vector<feSpace*> &space, feMesh *mesh, int degQuadrature, feSysElm *sysElm);
	~feBilinearForm();

  int getCncGeoTag(){ return _cncGeoTag; }

  int getNiElm(){ return _niElm; }
  int getNjElm(){ return _njElm; }
  std::vector<int> &getAdrI(){ return _adrI; }
  std::vector<int> &getAdrJ(){ return _adrJ; }

  double** getAe(){ return _Ae; }
  double*  getBe(){ return _Be; }

  void initialize_vadij_only(feMetaNumber *metaNumber, int numElem);
  void initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem);

  void computeMatrix(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem);
  void computeResidual(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem);

  void printInfo(){
    printf("============== Bilinear form ==============\n");
    printf("_cncGeoTag = %d\n", _cncGeoTag);
    printf("_nCoord = %d\n", _nCoord);
    printf("_nGeoNodes = %d\n", _nGeoNodes);
    printf("_nGeoElm = %d\n", _nGeoElm);
    printf("_geoCoord :"); for(auto val : _geoCoord) std::cout<<val<<" "; std::cout<<std::endl;
    printf("_iVar :"); for(auto val : _iVar) std::cout<<val<<" "; std::cout<<std::endl;
    printf("_jVar :"); for(auto val : _jVar) std::cout<<val<<" "; std::cout<<std::endl;
    printf("_adrI :"); for(auto val : _adrI) std::cout<<val<<" "; std::cout<<std::endl;
    printf("_adrJ :"); for(auto val : _adrJ)  std::cout<<val<<" "; std::cout<<std::endl;
    // printf("_Ae :"); for(auto val : _Ae) std::cout<<val<<" "; std::cout<<std::endl;
    // printf("_Be :"); for(auto val : _Be) std::cout<<val<<" "; std::cout<<std::endl;
  }
};

#endif