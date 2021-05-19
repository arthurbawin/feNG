#ifndef _FEQUADRATURE_
#define _FEQUADRATURE_

#include <vector>

class feQuadrature{

protected:
  int _dim;
  int _nQuad;
  std::vector<double> _w;
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _z;

public:
	feQuadrature(int nQuad, int dim = 1) : _dim(dim), _nQuad(nQuad){
    if(dim == 1){
      switch(nQuad){
        case 2:
          _w = {1.,1.};
          _x = {5.773502691896257e-01,-5.773502691896257e-01};
          _y = {0., 0.};
          _z = {0., 0.};
          break;
        case 3:
          _w = {5.555555555555544e-01, 8.888888888888888e-01, 5.555555555555544e-01};
          _x = {7.745966692414834e-01, 0.0                  ,-7.745966692414834e-01};
          _y = {0., 0., 0.};
          _z = {0., 0., 0.};
          break;
        default: 
          printf("Règle d'intégration non définie.\n");
      }
    } else if(dim == 2){
      switch(nQuad){
        case 3:
          _w = {0.166666666666667,0.166666666666667,0.166666666666667};
          _x = {0.166666666666667,0.666666666666667,0.166666666666667};
          _y = {0.166666666666667,0.166666666666667,0.666666666666667};
          _z = {0., 0., 0.};
          break;
        case 6:
          _w = {0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667};
          _x = {0.659027622374092, 0.659027622374092, 0.231933368553031, 0.231933368553031, 0.109039009072877, 0.109039009072877};
          _y = {0.231933368553031, 0.109039009072877, 0.659027622374092, 0.109039009072877, 0.659027622374092, 0.231933368553031};
          _z = {0., 0., 0., 0., 0., 0.};
          break;
        default: 
          printf("Règle d'intégration non définie.\n");
      }
    }
	};
	~feQuadrature() {}

  int getDim(){ return _dim; };
  int getNQuad(){ return _nQuad; }
  std::vector<double> getWeights(){ return _w; }
  std::vector<double> getPoints(){ return _x; }
  std::vector<double> getYPoints(){ return _y; }
  std::vector<double> getZPoints(){ return _z; }
};

#endif