#ifndef _FEQUADRATURE_
#define _FEQUADRATURE_

#include <vector>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>



class feQuadrature{

protected:
  int _degQuad;
  std::vector<double> _w;
  std::vector<double> _x;

public:
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