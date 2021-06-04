#ifndef _FEQUADRATURE_
#define _FEQUADRATURE_

#include <vector>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <cmath>



class feQuadrature{

protected:
  int _degQuad;
  std::vector<double> _w;
  std::vector<double> _x;

public:
<<<<<<< HEAD
  feQuadrature(int degQuad) : _degQuad(degQuad){
    
  };
  virtual ~feQuadrature() {}
  int _nQuad = (_degQuad +2)/2;
  virtual int getNQuad(){ return _nQuad; }
  virtual std::vector<double> getWeights(){ return _w; }
  virtual std::vector<double> getPoints(){ return _x; }
};

class feQuadrature2 : public feQuadrature{
public:
    feQuadrature2 ( int degQuad) : feQuadrature(degQuad) {
      int _nQuad = (_degQuad +2)/2;
      _w.resize(_nQuad);
      _x.resize(_nQuad);
      calculateWeightAndRoot();
      
      
      /*for(int i=0; i<nQuad;++i){
        _w[i]=_w[i+1];
        _x[i]=_x[i+1];
=======
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
        case 9:
          _x.resize(9);
          _y.resize(9);
          _w.resize(9);

          _x[0]=1.884094059521e-01;
          _y[0]=7.876594617608e-01;
          _w[0]=1.939638330596e-02;
          _x[1]=5.239790677201e-01;
          _y[1]=4.094668644407e-01;
          _w[1]=6.367808509989e-02;
          _x[2]=8.086943856777e-01;
          _y[2]=8.858795951270e-02;
          _w[2]=5.581442048304e-02;
          _x[3]=1.061702691196e-01;
          _y[3]=7.876594617608e-01;
          _w[3]=3.103421328954e-02;
          _x[4]=2.952665677796e-01;
          _y[4]=4.094668644407e-01;
          _w[4]=1.018849361598e-01;
          _x[5]=4.557060202436e-01;
          _y[5]=8.858795951270e-02;
          _w[5]=8.930307277287e-02;
          _x[6]=2.393113228708e-02;
          _y[6]=7.876594617608e-01;
          _w[6]=1.939638330596e-02;
          _x[7]=6.655406783916e-02;
          _y[7]=4.094668644407e-01;
          _w[7]=6.367808509989e-02;
          _x[8]=1.027176548096e-01;
          _y[8]=8.858795951270e-02;
          _w[8]=5.581442048304e-02;
          _z.resize(9, 0.);
          // _w = {0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667, 0.16666666666666666667};
          // _x = {0.659027622374092, 0.659027622374092, 0.231933368553031, 0.231933368553031, 0.109039009072877, 0.109039009072877};
          // _y = {0.231933368553031, 0.109039009072877, 0.659027622374092, 0.109039009072877, 0.659027622374092, 0.231933368553031};
          // _z = {0., 0., 0., 0., 0., 0.};
          break;
        case 16:
          _x.resize(16);
          _y.resize(16);
          _w.resize(16);
          _z.resize(16, 0.);

          _x[0]=1.300560792168e-01;
          _y[0]=8.602401356562e-01;
          _w[0]=5.423225910525e-03;
          _x[1]=3.874974834067e-01;
          _y[1]=5.835904323689e-01;
          _w[1]=2.258404928237e-02;
          _x[2]=6.729468631505e-01;
          _y[2]=2.768430136381e-01;
          _w[2]=3.538806789809e-02;
          _x[3]=8.774288093305e-01;
          _y[3]=5.710419611452e-02;
          _w[3]=2.356836819338e-02;
          _x[4]=9.363778443733e-02;
          _y[4]=8.602401356562e-01;
          _w[4]=1.016725956448e-02;
          _x[5]=2.789904634965e-01;
          _y[5]=5.835904323689e-01;
          _w[5]=4.233972452175e-02;
          _x[6]=4.845083266304e-01;
          _y[6]=2.768430136381e-01;
          _w[6]=6.634421610705e-02;
          _x[7]=6.317312516411e-01;
          _y[7]=5.710419611452e-02;
          _w[7]=4.418508852236e-02;
          _x[8]=4.612207990645e-02;
          _y[8]=8.602401356562e-01;
          _w[8]=1.016725956448e-02;
          _x[9]=1.374191041346e-01;
          _y[9]=5.835904323689e-01;
          _w[9]=4.233972452175e-02;
          _x[10]=2.386486597314e-01;
          _y[10]=2.768430136381e-01;
          _w[10]=6.634421610705e-02;
          _x[11]=3.111645522444e-01;
          _y[11]=5.710419611452e-02;
          _w[11]=4.418508852236e-02;
          _x[12]=9.703785126946e-03;
          _y[12]=8.602401356562e-01;
          _w[12]=5.423225910525e-03;
          _x[13]=2.891208422439e-02;
          _y[13]=5.835904323689e-01;
          _w[13]=2.258404928237e-02;
          _x[14]=5.021012321137e-02;
          _y[14]=2.768430136381e-01;
          _w[14]=3.538806789809e-02;
          _x[15]=6.546699455501e-02;
          _y[15]=5.710419611452e-02;
          _w[15]=2.356836819338e-02;
          break;
        default: 
          printf("Règle d'intégration non définie.\n");
>>>>>>> origin/master
      }
      _w.resize(nQuad);
      _x.resize(nQuad);*/
    }
    ~feQuadrature2() {}



    const double EPSILON =1e-15;

    struct Result {
        double value;
        double derivative;

        Result() : value(0), derivative(0) {}
        Result(double val, double deriv) : value(val), derivative(deriv) {}
    };

    void calculateWeightAndRoot(); 

    Result calculatePolynomialValueAndDerivative(double x); 


    
};

class feQuadratureTriangle : public feQuadrature{

  protected: 
    // std::vector<double> _W;
    // std::vector<double> _xr;
    // std::vector<double> _yr;
    // std::vector<double> _zr;
    double _W[100];
    double _xr[100];
    double _yr[100];
    double _zr[100]; //TODO : construire la table avec la bonne dimension


  public:
    feQuadratureTriangle ( int degQuad) : feQuadrature(degQuad) {

        int _nQuad = (_degQuad +2)/2;
        int _dimQuad = 3;
        

        // int _nQuad = (_degQuad +2)/2;
        // _W.resize(pow(_nQuad,2));
        // _xr.resize(pow(_nQuad,2));
        // _yr.resize(pow(_nQuad,2));
        // _zr.resize(pow(_nQuad,2));

// 
        // _W.resize(_nQuad * _nQuad);
        // _xr.resize(_nQuad * _nQuad);
        // _yr.resize(_nQuad * _nQuad);
        // _zr.resize(_nQuad * _nQuad);


        if (_dimQuad==2) {calculateWeightAndRootTri();}
        else {calculateWeightAndRootTetra();}
  }
  void calculateWeightAndRootTri();
  void calculateWeightAndRootTetra();
  double getXPoints(int i){ return _xr[i]; }
  double getYPoints(int i ){ return _yr[i]; }
  double getZPoints(int i ){ return _zr[i]; }
  double get3DWeights(int i ){ return _W[i]; }

  // std::vector<double> getXPoints(){ return _xr; }
  // std::vector<double> getYPoints(){ return _yr; }
  // std::vector<double> getZPoints(){ return _zr; }
  // std::vector<double> get3DWeights(){ return _W; }

};




/*const double LegendrePolynomial::EPSILON = 1e-15;*/




/*std::cout<<"Test regle quadrature2"<<"\n"<<"_____________"<<"\n";
    LegendrePolynomial *poly = new LegendrePolynomial(0.0,1.0,5);
    std::vector<double> Poids=poly->mWeight;
    std::vector<double> Racine=poly->mRoot;
    int Taille=Poids.size();
    std::cout<<"Poids "<<"\n";
    for(int i=1;i<Taille;++i){
      std::cout<< Poids[i] <<",";
    }
    std::cout<<"\n";
    std::cout<<"Racine "<<"\n";
    for(int i=1;i<Taille;++i){
      std::cout<< Racine[i] <<",";
    }
    std::cout<<"\n";*/


#endif