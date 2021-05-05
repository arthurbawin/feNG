#ifndef _FEQUADRATURE_
#define _FEQUADRATURE_

#include <vector>

class feQuadrature{

protected:
  int _nQuad;
  std::vector<double> _w;
  std::vector<double> _x;

public:
	feQuadrature(int nQuad) : _nQuad(nQuad){
    switch(nQuad){
      case 2:
        _w = {1.,1.};
        _x = {5.773502691896257e-01,-5.773502691896257e-01};
        break;
      case 3:
        _w = {5.555555555555544e-01, 8.888888888888888e-01, 5.555555555555544e-01};
        _x = {7.745966692414834e-01, 0.0                  ,-7.745966692414834e-01};
        break;
      default: 
        printf("Règle d'intégration non définie.\n");
    }
	};
	~feQuadrature() {}

  int getNQuad(){ return _nQuad; }
  std::vector<double> getWeights(){ return _w; }
  std::vector<double> getPoints(){ return _x; }
};

#endif