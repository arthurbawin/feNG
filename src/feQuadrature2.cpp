#include "feQuadrature.h"

#include <cmath>

feQuadrature::feQuadrature(int degQuad, int dimQuad, std::string form  ) : _degQuad(degQuad) , _dimQuad(dimQuad)   {

    _nQuad1D = (_degQuad +2)/2;
    _w.resize (pow(_nQuad1D,_dimQuad));
    _xr.resize(pow(_nQuad1D,_dimQuad));
    _yr.resize(pow(_nQuad1D,_dimQuad));
    _zr.resize(pow(_nQuad1D,_dimQuad));
    _w1D.resize(_nQuad1D);
    _x1D.resize(_nQuad1D);
        
    if(_dimQuad == 1){
      calculateWeightAndRoot(_nQuad1D);
      _w=_w1D;
      _xr=_x1D;
    }

    else if (_dimQuad == 2 && form =="TriP1" ){
      // _w.resize (_degQuad*(_degQuad+1)/2 -1); //On redimensionne car version allegé de la quadrature
      // _xr.resize(_degQuad*(_degQuad+1)/2 -1);
      // _yr.resize(_degQuad*(_degQuad+1)/2 -1);
      // _zr.resize(_degQuad*(_degQuad+1)/2 -1);
      // _w.resize (_nQuad1D*(_nQuad1D+1)/2 -1); //On redimensionne car version allegé de la quadrature
      // _xr.resize(_nQuad1D*(_nQuad1D+1)/2 -1);
      // _yr.resize(_nQuad1D*(_nQuad1D+1)/2 -1);
      // _zr.resize(_nQuad1D*(_nQuad1D+1)/2 -1);
      calculateWeightAndRootTri();
    }
    else if (_dimQuad == 2 && form =="SquareP1" ){
      calculateWeightAndRootSquare();
    }

    else if (_dimQuad == 3 && form =="SquareP1" ){
      calculateWeightAndRootCube();
    }

    else {
      calculateWeightAndRootTetra();
    }
}

void feQuadrature::calculateWeightAndRoot(int nQuadLocal) {
    // std::cout<<"Computing 1D points for "<<nQuadLocal<<" points"<<std::endl;
    _w1D.resize(nQuadLocal);
    _x1D.resize(nQuadLocal);

    if(nQuadLocal == 1){
        _x1D[0] = 0.0;
        _w1D[0] = 2.0;
    } else{
        _nQuad1D = nQuadLocal ; //On met a jour le nombre de point que l'on veut pour la règle 

        for(int i = 0; i < _nQuad1D; i++) {
            double root = cos(M_PI * ((i+1)-0.25)/(_nQuad1D+0.5));

            calculatePolynomialValueAndDerivative(root, result);

            double newtonRaphsonRatio;
            do {
                newtonRaphsonRatio = result[0]/result[1];
                root -= newtonRaphsonRatio;
                calculatePolynomialValueAndDerivative(root, result);
            } while (fabs(newtonRaphsonRatio) > EPSILON);

            _x1D[_nQuad1D-1-i] = root;
            _w1D[_nQuad1D-1-i] = 2.0/((1-root*root)*result[1]*result[1]);
        }
        // std::cout<<"Done computing 1D :"<<std::endl;
        // for(auto val : _x1D) std::cout<<val<<" ";
        // std::cout<<std::endl;
    }
}

void feQuadrature::calculatePolynomialValueAndDerivative(double x, double *result) {
    result[0] =x; 
    result[1] =0.0;
    double value_minus_1 = 1;
    const double f = 1/(x*x-1);

    for(int i = 2; i <= _nQuad1D; i++) {
        const double value = ((2*i-1)*x*result[0]-(i-1)*value_minus_1)/i;
        result[1] = i*f*(x*value - result[0]);

        value_minus_1 = result[0];
        result[0] = value;
    }
}

// void feQuadrature::calculateWeightAndRootTri() {

//     calculateWeightAndRoot(_degQuad-1); 
//     std::vector<double> _xi = _x1D;
//     std::vector<double> _wi = _w1D;

//     int l=0;
//     for(int i = 0; i < _degQuad-1; ++i){

//         calculateWeightAndRoot(_degQuad-i);
//         std::vector<double> _eta = _x1D;
//         std::vector<double> _wj  = _w1D; 

//         for(int j = 0; j < _degQuad-i; ++j){
//             std::cout<<i<<" - "<<j<<std::endl;
//             _yr[l] = (1+_xi[i])/2.0;
//             _xr[l] = (1-_xi[i]) * (1+_eta[j]) / 4.0 ; 
//             _w[l] = _wj[j]*_wi[i]*(1-_xi[i])/8.0;
//             ++l;
//         }
//     }
// }

void feQuadrature::calculateWeightAndRootTri() {

    _nQuad1D = (_degQuad +2)/2;
    calculateWeightAndRoot(_nQuad1D); 
    int l=0;
    for(int i = 0; i < _nQuad1D; ++i){
        for(int j = 0; j < _nQuad1D; ++j){
            std::cout<<i<<" - "<<j<<std::endl;
            _yr[l] = (1+_x1D[i])/2.0;
            _xr[l] = (1-_x1D[i]) * (1+_x1D[j]) / 4.0 ; 
            _w[l] = _w1D[j]*_w1D[i]*(1-_x1D[i])/8.0;
            ++l;
        }
    }
}


void feQuadrature::calculateWeightAndRootTetra() {
        

    // feQuadrature2 Quad = feQuadrature2(_degQuad);
    // std::vector<double>_wi=Quad.getWeights();
    // std::vector<double>_xi=Quad.getPoints();

    // int l=0; 
    // for(int i = 0; i< _nQuad; ++i){

    // // feQuadrature2 Quad = feQuadrature2((_nQuad+1-i)*2);
    // // std::vector<double>_wj=Quad.getWeights();
    // // std::vector<double>_eta=Quad.getPoints();

    //     for(int j = 0; j < _nQuad; ++j){
    //         _yr[l] = (1-_xi[i])/2;
    //         _xr[l] = (1-_xi[i]) * (1-_xi[j]) / 4 ; 
    //         _W[l] = _wi[j]*_wi[i]*(1+_xi[i])/8;
    //         l=l+1;
    //     }

    // feQuadrature2 Quad = feQuadrature2(_nQuad);
    // std::vector<double>_w=Quad.getWeights();
    // std::vector<double>_x=Quad.getPoints();

    calculateWeightAndRoot(_nQuad1D);



    int l=0; 
    for(int i = 0; i< _nQuad1D; ++i){
        for(int j = 0; j < _nQuad1D; ++j){
            for(int k = 0; k < _nQuad1D; ++k){
                _xr[l] = (1-_x1D[i])/2;
                _yr[l] = (1+_x1D[i])*(1-_x1D[j]) /4 ; 
                _zr[l] = (1+_x1D[i])*(1+_x1D[j])*(1-_x1D[k])/8;
                _w[l] = _w1D[j]*_w1D[i]*_w1D[k]*(1+_x1D[i])*(1+_x1D[i])*(1+_x1D[j])/64;
                l=l+1;

            }
            
        }
    } 
        
        
        
        
}


void feQuadrature::calculateWeightAndRootSquare() {
        

    calculateWeightAndRoot(_nQuad1D);


    int l=0; 
    for(int i = 0; i< _nQuad1D; ++i){
        for(int j = 0; j < _nQuad1D; ++j){
            _xr[l] = _x1D[i];
            _yr[l] = _x1D[j]; 
            _w[l] = _w1D[j]*_w1D[i];
            l=l+1;
        }
            
    }
} 
        
        
void feQuadrature::calculateWeightAndRootCube() {
        


    calculateWeightAndRoot(_nQuad1D);



    int l=0; 
    for(int i = 0; i< _nQuad1D; ++i){
        for(int j = 0; j < _nQuad1D; ++j){
            for(int k = 0; k < _nQuad1D; ++k){
                _xr[l] = _x1D[i];
                _yr[l] = _x1D[j]; 
                _zr[l] = _x1D[k];
                _w[l] = _w1D[j]*_w1D[i]*_w1D[k];
                l=l+1;

            }
            
        }
    } 
        
        
        
        
}
        

