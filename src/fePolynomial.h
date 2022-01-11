#ifndef _FEPOLYNOMIAL_
#define _FEPOLYNOMIAL_

#include "feNG.h"
#include "feSpace.h"
#include "feMesh.h"

class Polynomial {
protected:
  std::vector<double> _c; // coefficients
  int _deg; // degree
  int _nMonomials;

  std::vector<int> _exponentsX;
  std::vector<int> _exponentsY;

public:
  Polynomial(int deg, std::vector<double> &coeff);
  Polynomial(const Polynomial &p);
  ~Polynomial() {}

  inline double coeff(int i) const { return _c[i]; }
  inline int deg() const { return _deg; }
  inline int nMonomials() const { return _nMonomials; }

  Polynomial operator+(Polynomial &other){
    Polynomial p(*this);
    if(this->deg() == other.deg()){
      for(int i = 0; i < _nMonomials; i++) _c[i] += other._c[i];
    } else{
      printf("In fePolynomial : Error - Cannot add polynomial of different degrees for now. Returning current polynomial\n");
    }
    return p;
  }
  Polynomial operator-(Polynomial &other){
    Polynomial p(*this);
    if(this->deg() == other.deg()){
      for(int i = 0; i < _nMonomials; i++) _c[i] -= other._c[i];
    } else{
      printf("In fePolynomial : Error - Cannot add polynomial of different degrees for now. Returning current polynomial\n");
    }
    return p;
  }
  Polynomial operator*(double &other){
    Polynomial p(*this);
    for(int i = 0; i < _nMonomials; i++) _c[i] *= other;
    return p;
  }
  Polynomial &operator+=(Polynomial &other){
    if(this->deg() == other.deg()){
      for(int i = 0; i < _nMonomials; i++) _c[i] += other._c[i];
    } else{
      printf("In fePolynomial : Error - Cannot add polynomial of different degrees for now. Returning current polynomial\n");
    }
    return *this;
  }
  Polynomial &operator-=(Polynomial &other){
    if(this->deg() == other.deg()){
      for(int i = 0; i < _nMonomials; i++) _c[i] -= other._c[i];
    } else{
      printf("In fePolynomial : Error - Cannot add polynomial of different degrees for now. Returning current polynomial\n");
    }
    return *this;
  }
  Polynomial &operator*=(double &other){
    for(int i = 0; i < _nMonomials; i++) _c[i] *= other;
    return *this;
  }

  double eval(std::vector<double> &x);
  double integrate(feMesh *mesh, std::string cncGeoID, int iElm);
  double innerProduct(Polynomial &other, feMesh *mesh, std::string cncGeoID, int iElm);

  void print(){
    for(int i = 0; i < _nMonomials; ++i)
      printf("%+-12.16e ", _c[i]);
    printf("\n");
  }
};


#endif