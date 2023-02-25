#ifndef _FEMETRICINTERFACE_
#define _FEMETRICINTERFACE_

#include <memory>

// A small interface for fixed sizes/small matrices to avoid
// including Eigen everywhere and slowing down compilation

// IMPORTANT: Classes use the PIMPL idiom and are not copyable

class MetricTensor
{
public:
  // "Real" constructor, for us mere mortals
  MetricTensor(const double val = 1.0);

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

  MetricTensor inverse() const;
  double determinant() const;

  // Compute diagonalization M = R * |D| * Rt
  // with |D| = diag( |lambda1|, |lambda2| )
  // and bound the eigenvalues between lMin and lMax
  MetricTensor boundEigenvaluesOfAbs(const double lMin, const double lMax) const;

  void print() const;

private:
  class MetricTensorImpl;
  std::unique_ptr<MetricTensorImpl> _impl;
};

class Vector
{
public:
  friend class SquareMatrix;

  // "Real" constructor, for us mere mortals
  Vector(const int size);

  // Default, copy/move constructors/assignments
  // To make insertion in e.g. maps safe
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
  // "Real" constructor, for us mere mortals
  SquareMatrix(const int size);

  // Default, copy/move constructors/assignments
  // To make insertion in e.g. maps safe
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

#endif