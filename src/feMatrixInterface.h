#ifndef _FEMETRICINTERFACE_
#define _FEMETRICINTERFACE_

#include "feNG.h"
#include <memory>

// A small interface for fixed sizes/small matrices to avoid
// including Eigen everywhere and slowing down compilation
//
// IMPORTANT: Classes use the PIMPL idiom, so instantiated objects are *not* copyable
//
// "Reference" for copy/move destructors/assignment operators:
// https://stackoverflow.com/questions/71104545/constructor-
// and-destructor-in-c-when-using-the-pimpl-idiom

class SquareMatrix;

class MetricTensor
{
public:
  friend class SquareMatrix;

  // "Real" constructor, for us mere mortals
  MetricTensor(const double val = 1.0);
  MetricTensor(const double m11, const double m12, const double m22);
  MetricTensor(const double eigenvalues[2], const double eigenvector1[2], const double eigenvector2[2]);
  MetricTensor(const SquareMatrix &other);

  // Default, copy/move constructors/assignments
  // To make insertion in e.g. maps safe
  ~MetricTensor();
  MetricTensor(MetricTensor const &) = delete;
  MetricTensor &operator=(MetricTensor const &) = delete;
  MetricTensor(MetricTensor &&) noexcept;
  MetricTensor &operator=(MetricTensor &&) noexcept;

  // Return a new MetricTensor whose matrix is a copy of this
  MetricTensor copy() const;
  void assignMatrixFrom(const MetricTensor &other);
  void assignMatrixFrom(const double other[2][2]);

  double &operator()(int i, int j);
  double operator()(int i, int j) const;
  MetricTensor &operator*=(const double &val);
  MetricTensor operator*(const double &val) const;
  MetricTensor operator/(const double &val) const;
  MetricTensor operator+(const MetricTensor &other) const;
  MetricTensor operator-(const MetricTensor &other) const;
  MetricTensor operator*(const MetricTensor &other) const;
  SquareMatrix operator*(const SquareMatrix &other) const;

  double trace() const;
  double determinant() const;
  MetricTensor inverse() const;
  MetricTensor transpose() const;
  MetricTensor log() const;
  MetricTensor exp() const;
  MetricTensor pow(const double &power) const;
  MetricTensor sqrt() const;

  double frobeniusNorm() const;

  // Returns the dot product weighted by the metric x1^T * M * x2
  double dotProduct(const double x1[2], const double x2[2]) const;

  double maxCoeff() const;

  double getMinEigenvalue() const;
  double getMaxEigenvalue() const;

  void getEigenvectorsAndEigenvalues(double v1[2], double v2[2], double &lambda1, double &lambda2) const;

  // Set matrix to (m + m^T)/2
  void symmetrize();

  // Compute diagonalization M = R * |D| * Rt
  // with |D| = diag( |lambda1|, |lambda2| )
  MetricTensor absoluteValueEigen() const;
  // Compute diagonalization M = R * |D| * Rt
  // with |D| = diag( |lambda1|, |lambda2| )
  // and bound the eigenvalues between lMin and lMax
  MetricTensor boundEigenvalues(const double lMin, const double lMax) const;
  MetricTensor boundEigenvaluesOfAbs(const double lMin, const double lMax) const;
  MetricTensor boundEigenvaluesOfAbsIsotropic(const double lMin, const double lMax) const;
  MetricTensor setEigenvectorsAndBoundEigenvalues(const double ev1[2], const double ev2[2], const double lMin, const double lMax) const;
  MetricTensor limitAnisotropy(const double alpha) const;

  MetricTensor spanMetricInMetricSpace(const double gradation, const double pq[2]) const;
  MetricTensor spanMetricInPhysicalSpace(const double gradation, const double pq[2]) const;
  MetricTensor spanMetricInMixedSpace(const double gradation, const double pq[2], const double t) const;

  void print(const int numDigits = 8) const;

private:
  class MetricTensorImpl;
  std::unique_ptr<MetricTensorImpl> _impl;
};

class Vector
{
public:
  friend class Matrix;
  friend class SquareMatrix;

  // "Real" constructor, for us mere mortals
  Vector(const int size);

  // Default, copy/move constructors/assignments
  Vector();
  ~Vector();
  Vector(Vector const &) = delete;
  Vector &operator=(Vector const &) = delete;
  Vector(Vector &&) noexcept;
  Vector &operator=(Vector &&) noexcept;

  int getSize() const;

  double &operator()(int i);

  void print() const;

private:
  class VectorImpl;
  std::unique_ptr<VectorImpl> _impl;
};

class SquareMatrix
{
public:
  friend class MetricTensor;
  friend class Matrix;

  // "Real" constructor, for us mere mortals
  SquareMatrix(const int size);

  // Default, copy/move constructors/assignments
  SquareMatrix();
  ~SquareMatrix();
  SquareMatrix(SquareMatrix const &) = delete;
  SquareMatrix &operator=(SquareMatrix const &) = delete;
  SquareMatrix(SquareMatrix &&) noexcept;
  SquareMatrix &operator=(SquareMatrix &&) noexcept;

  int getSize() const;
  int rank() const;

  void getEigenvectorsAndEigenvalues(double v1[2], double v2[2], double &lambda1, double &lambda2) const;

  double &operator()(int i, int j);
  double operator()(int i, int j) const;
  Vector operator*(const Vector &v);
  SquareMatrix operator*(const MetricTensor &other) const;
  SquareMatrix operator*(const SquareMatrix &other) const;

  void print() const;

  double determinant() const;
  SquareMatrix inverse() const;
  void inverse(SquareMatrix &) const;
  SquareMatrix transpose() const;

private:
  class SquareMatrixImpl;
  std::unique_ptr<SquareMatrixImpl> _impl;
};

class Matrix
{
public:
  // "Real" constructor, for us mere mortals
  Matrix(const int sizeM, const int sizeN);

  // Default, copy/move constructors/assignments
  Matrix();
  ~Matrix();
  Matrix(Matrix const &) = delete;
  Matrix &operator=(Matrix const &) = delete;
  Matrix(Matrix &&) noexcept;
  Matrix &operator=(Matrix &&) noexcept;

  int getSizeM() const;
  int getSizeN() const;
  int rank() const;

  double &operator()(int i, int j);
  Vector operator*(const Vector &v) const;

  void print() const;

  // Compute A^T*A (dim = NxN if A is MxN)
  SquareMatrix multiplyByTranspose() const;
  // Compute (A^T*A)^-1 * A^T (dim = NxM if A is MxN)
  Matrix getLeastSquaresMatrix() const;

private:
  class MatrixImpl;
  std::unique_ptr<MatrixImpl> _impl;
};

#endif