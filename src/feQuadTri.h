#ifndef _FEQUADRATURETRI_
#define _FEQUADRATURETRI_

#include <vector>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>

class feQuadratureTri{


protected: 
	int _degQuad;
  	int _nQuad;
  	std::vector<double> _w;
  	std::vector<double> _x;

public:
	feQuadrature(int nQuad) : _nQuad(nQuad){
    
	};
	virtual ~feQuadrature() {}

  virtual int getNQuad(){ return _nQuad; }
  virtual std::vector<double> getWeights(){ return _w; }
  virtual std::vector<double> getPoints(){ return _x; }
};
















}