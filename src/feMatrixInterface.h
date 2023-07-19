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

class MetricTensor
{
public:
  // "Real" constructor, for us mere mortals
  MetricTensor(const double val = 1.0);
  MetricTensor(const double eigenvalues[2], const double eigenvector1[2], const double eigenvector2[2]);

  // Default, copy/move constructors/assignments
  // To make insertion in e.g. maps safe
  ~MetricTensor();
  MetricTensor(MetricTensor const &) = delete;
  MetricTensor &operator=(MetricTensor const &) = delete;
  MetricTensor(MetricTensor &&) noexcept;
  MetricTensor &operator=(MetricTensor &&) noexcept;

  double &operator()(int i, int j);
  double operator()(int i, int j) const;
  MetricTensor &operator*=(const double &val);
  MetricTensor operator*(const double &val) const;
  MetricTensor operator/(const double &val) const;
  MetricTensor operator+(const MetricTensor &other) const;
  MetricTensor operator-(const MetricTensor &other) const;
  MetricTensor operator*(const MetricTensor &other) const;

  double determinant() const;
  MetricTensor inverse() const;
  MetricTensor transpose() const;
  MetricTensor log() const;
  MetricTensor exp() const;

  double maxCoeff() const;

  // Compute diagonalization M = R * |D| * Rt
  // with |D| = diag( |lambda1|, |lambda2| )
  MetricTensor absoluteValueEigen() const;
  // Compute diagonalization M = R * |D| * Rt
  // with |D| = diag( |lambda1|, |lambda2| )
  // and bound the eigenvalues between lMin and lMax
  MetricTensor boundEigenvaluesOfAbs(const double lMin, const double lMax) const;
  MetricTensor boundEigenvaluesOfAbsIsotropic(const double lMin, const double lMax) const;
  MetricTensor setEigenvectorsAndBoundEigenvalues(const double ev1[2], const double ev2[2], const double lMin, const double lMax) const;

  void print() const;

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

  double &operator()(int i, int j);
  Vector operator*(const Vector &v);

  void print() const;

  SquareMatrix inverse() const;

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