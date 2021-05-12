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
      
    }
	};
	virtual ~feQuadrature() {}

  virtual int getNQuad(){ return _nQuad; }
  virtual std::vector<double> getWeights(){ return _w; }
  virtual std::vector<double> getPoints(){ return _x; }
};

class feQuadrature2 : public feQuadrature{
public:
    feQuadrature2 ( int nQuad) : feQuadrature(nQuad) {
        calculateWeightAndRoot();
    }
    ~feQuadrature2() {}

    
  
    std::vector<double> _w;
    std::vector<double> _x;
    // std::vector<double> mWeight;
    // std::vector<double> mRoot;
    int _nQuad;


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