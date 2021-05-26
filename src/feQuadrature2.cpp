#include "feQuadrature.h"
#include "feMesh.h"
#include <cmath>


void feQuadrature2::calculateWeightAndRoot() {
        for(int step = 1; step <= _nQuad; step++) {
            double root = cos(M_PI * (step-0.25)/(_nQuad+0.5));
            Result result = calculatePolynomialValueAndDerivative(root);

            double newtonRaphsonRatio;
            do {
                newtonRaphsonRatio = result.value/result.derivative;
                root -= newtonRaphsonRatio;
                result = calculatePolynomialValueAndDerivative(root);
            } while (fabs(newtonRaphsonRatio) > EPSILON);

            _x[step-1] = root;
            _w[step-1] = 2.0/((1-root*root)*result.derivative*result.derivative);
        }
    }


feQuadrature2::Result feQuadrature2::calculatePolynomialValueAndDerivative(double x) {
        Result result(x, 0);

        double value_minus_1 = 1;
        const double f = 1/(x*x-1);
        for(int step = 2; step <= _nQuad; step++) {
            const double value = ((2*step-1)*x*result.value-(step-1)*value_minus_1)/step;
            result.derivative = step*f*(x*value - result.value);

            value_minus_1 = result.value;
            result.value = value;
        }

        return result;
    }



void feQuadratureTriangle::calculateWeightAndRootTri() {
        

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

        feQuadrature2 Quad = feQuadrature2((_nQuad-1)*2-1);
        std::vector<double>_wi=Quad.getWeights();
        std::vector<double>_xi=Quad.getPoints();

        int l=0; 
        for(int i = 0; i< _nQuad-1; ++i){

        feQuadrature2 Quad = feQuadrature2((_nQuad-i)*2-1);
        std::vector<double>_wj=Quad.getWeights();
        std::vector<double>_eta=Quad.getPoints();


            for(int j = 0; j < _nQuad-i; ++j){
                _yr[l] = (1-_xi[i])/2;
                _xr[l] = (1-_xi[i]) * (1-_eta[j]) / 4 ; 
                _W[l] = _wj[j]*_wi[i]*(1+_xi[i])/8;
                l=l+1;
            }
        } 
        
        
        
        
}


void feQuadratureTriangle::calculateWeightAndRootTetra() {
        

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

    feQuadrature2 Quad = feQuadrature2(_nQuad);
    std::vector<double>_w=Quad.getWeights();
    std::vector<double>_x=Quad.getPoints();

    int l=0; 
    for(int i = 0; i< _nQuad; ++i){
        for(int j = 0; j < _nQuad; ++j){
            for(int k = 0; k < _nQuad; ++k){
                _xr[l] = (1-_x[i])/2;
                _yr[l] = (1+_x[i])*(1-_x[j]) /4 ; 
                _zr[l] = (1+_x[i])*(1+_x[j])*(1-_x[k])/8;
                _W[l] = _w[j]*_w[i]*_w[k]*(1+_x[i])*(1+_x[i])*(1+_x[j])/64;
                l=l+1;

            }
            
        }
    } 
        
        
        
        
}

