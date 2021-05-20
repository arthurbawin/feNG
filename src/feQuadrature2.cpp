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
        

        feQuadrature2 Quad = feQuadrature2(_degQuad);
        _x=Quad.getWeights();
        _w=Quad.getPoints();

        if (1)
        {
           for(int i = 1; i<= _nQuad; i++){

                for(int j = 1; j <= _nQuad; j++){

                     _xr[(i-1)*_nQuad + j] = (1+_x[i]) * (1-_x[j]) / 4 ; 
                     _yr[(i-1)*_nQuad + j] = (1-_x[j])/2;

                    _W[(i-1)*_nQuad + j] = _w[i]*_w[j]*(1-_x[j])/8;
                }
            } 
        }
        else{   
            // on remplit la coordonnés z aussi 


        }
        
        
    }


