#include "feSpace.h"

void feVectorSpace::dotProductShapeShape(const int iNode, std::vector<double> &res)
{
  int offset = iNode * _nFunctions * _nComponents;
  double dotProd;
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      dotProd = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        dotProd += _L[offset + i * _nComponents + m] * _L[offset + j * _nComponents + m];
      }
      res[i * _nFunctions + j] = dotProd;
    }
  }
}

void feVectorSpace::dotProductShapeOther(const int iNode, const std::vector<double> &other,
                                         std::vector<double> &res)
{
  int offset = iNode * _nFunctions * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    res[i] = 0.;
    for(int m = 0; m < _nComponents; ++m) {
      res[i] += other[m] * _L[offset + i * _nComponents + m];
    }
  }
}

void feVectorSpace::dotProductShapeShapeOtherSpace(const int iNode, const feSpace *other,
                                                   std::vector<double> &res)
{
  int nFunctionsOther = other->getNumFunctions();
  int nComponentsOther = other->getNumComponents();
  int offset = iNode * _nFunctions * _nComponents;
  // !!! Assumes that this space and other have the same quadrature rule (using iNode to offset)
  int offsetOther = iNode * nFunctionsOther * nComponentsOther;
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < nFunctionsOther; ++j) {
      res[i * nFunctionsOther + j] = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        res[i * nFunctionsOther + j] += _L[offset + i * _nComponents + m] *
                                        other->shape()[offsetOther + j * nComponentsOther + m];
      }
    }
  }
}

void feVectorSpace::vectorDotGradShapeDotShape(const int iNode, const std::vector<double> &gradPhi,
                                               const std::vector<double> &other,
                                               std::vector<double> &res)
{
  int offset = iNode * _nFunctions * _nComponents;
  int offsetGrad = _dim * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i * _nFunctions + j] = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        for(int n = 0; n < _nComponents; ++n) {
          res[i * _nFunctions + j] +=
            other[m] * gradPhi[j * offsetGrad + m * _dim + n] * _L[offset + i * _nComponents + n];
        }
      }
    }
  }
}

void feVectorSpace::shapeDotTensorDotShape(const int iNode, const std::vector<double> &other,
                                           std::vector<double> &res)
{
  int offset = iNode * _nFunctions * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i * _nFunctions + j] = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        for(int n = 0; n < _nComponents; ++n) {
          res[i * _nFunctions + j] +=
            _L[offset + j * _nComponents + m] * other[m * _dim + n] * _L[offset + i * _nComponents + n];
        }
      }
    }
  } 
}

void feVectorSpace::doubleContractionGradShapeGradShape(const std::vector<double> &gradPhi,
                                                        std::vector<double> &res)
{
  int offset = _dim * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i * _nFunctions + j] = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        for(int n = 0; n < _nComponents; ++n) {
          res[i * _nFunctions + j] +=
            gradPhi[i * offset + m * _dim + n] * gradPhi[j * offset + m * _dim + n];
        }
      }
    }
  }
}

void feVectorSpace::doubleContractionGradShapeGradShapeTransposed(
  const std::vector<double> &gradPhi, std::vector<double> &res)
{
  int offset = _dim * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i * _nFunctions + j] = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        for(int n = 0; n < _nComponents; ++n) {
          res[i * _nFunctions + j] +=
            gradPhi[i * offset + m * _dim + n] * gradPhi[j * offset + n * _dim + m];
        }
      }
    }
  }
}

void feVectorSpace::doubleContractionGradShapeOther(const std::vector<double> &gradPhi,
                                                    const std::vector<double> &other,
                                                    std::vector<double> &res)
{
  int offset = _dim * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    res[i] = 0.;
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        res[i] += gradPhi[i * offset + m * _dim + n] * other[m * _dim + n];
      }
    }
  }
}

void feVectorSpace::divergence(const std::vector<double> &gradPhi, std::vector<double> &res)
{
  int offset = _dim * _nComponents;
  for(int i = 0; i < _nFunctions; ++i) {
    if(_dim == 2) {
      res[i] = gradPhi[i * offset + 0] + gradPhi[i * offset + 3];
    } else {
      res[i] = gradPhi[i * offset + 0] + gradPhi[i * offset + 4] + gradPhi[i * offset + 8];
    }
  }
}