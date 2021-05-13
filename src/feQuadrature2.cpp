#include "feQuadrature.h"
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



