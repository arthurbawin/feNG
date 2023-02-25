#include "feMatrixInterface.h"

#include "../contrib/Eigen/Dense"

// ============================================================
// Implementation of MetricTensor
// ============================================================
class MetricTensor::MetricTensorImpl
{
public:
  Eigen::Matrix2d _m;

  // Default, copy/move constructor/assignment
  MetricTensorImpl() = default;
  ~MetricTensorImpl() = default;
  MetricTensorImpl(MetricTensorImpl const &) = default;
  MetricTensorImpl &operator=(MetricTensorImpl const &) = default;
  MetricTensorImpl(MetricTensorImpl &&) noexcept = default;
  MetricTensorImpl &operator=(MetricTensorImpl &&) noexcept = default;

  // "Real" constructor
  MetricTensorImpl(const double val)
  {
    _m = Eigen::Matrix2d::Identity();
    _m(0, 0) = val;
    _m(1, 1) = val;
  }

  void setMatrix(const Eigen::Matrix2d &other) { _m = other; };
};

// Default, copy/move constructor/assignment
MetricTensor::~MetricTensor() = default;
MetricTensor::MetricTensor(MetricTensor &&) noexcept = default;
MetricTensor &MetricTensor::operator=(MetricTensor &&) noexcept = default;

// "Real" constructor
MetricTensor::MetricTensor(const double val) : _impl(std::make_unique<MetricTensorImpl>(val)) {}

// =============== MetricTensor member functions ===============
double &MetricTensor::operator()(int i, int j) { return _impl->_m(i, j); }
double MetricTensor::operator()(int i, int j) const { return _impl->_m(i, j); }
MetricTensor &MetricTensor::operator*=(const double &val)
{
  _impl->_m *= val;
  return *this;
}

void MetricTensor::print() const { std::cout << _impl->_m << std::endl; }

MetricTensor MetricTensor::inverse() const
{
  MetricTensor inv(1.0);
  inv._impl->setMatrix(this->_impl->_m.inverse());
  return inv;
}

double MetricTensor::determinant() const { return _impl->_m.determinant(); }

MetricTensor MetricTensor::boundEigenvaluesOfAbs(const double lMin, const double lMax) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  Eigen::Matrix2d D = Eigen::Matrix2d::Identity();
  D(0, 0) = fmin(lMax, fmax(lMin, fabs(ev(0).real())));
  D(1, 1) = fmin(lMax, fmax(lMin, fabs(ev(1).real())));
  MetricTensor res(1.0);
  res._impl->setMatrix(es.eigenvectors().real() * D * es.eigenvectors().transpose().real());
  return res;
}

// ============================================================
// Implementation of Vector
// ============================================================
class Vector::VectorImpl
{
public:
  int _size;
  Eigen::VectorXd _v;

  // Default, copy/move constructor/assignment
  VectorImpl() = default;
  ~VectorImpl() = default;
  VectorImpl(VectorImpl const &) = default;
  VectorImpl &operator=(VectorImpl const &) = default;
  VectorImpl(VectorImpl &&) noexcept = default;
  VectorImpl &operator=(VectorImpl &&) noexcept = default;

  // "Real" constructor
  VectorImpl(const int size) : _size(size) { _v = Eigen::VectorXd::Zero(size); }

  void setVector(const Eigen::VectorXd &other) { _v = other; };
};

// Default, copy/move constructor/assignment
Vector::Vector() : _impl(new VectorImpl{}) {}
Vector::~Vector() = default;
Vector::Vector(Vector &&) noexcept = default;
Vector &Vector::operator=(Vector &&) noexcept = default;

// "Real" constructor
Vector::Vector(const int size) : _impl(std::make_unique<VectorImpl>(size)) {}

// =============== Vector member functions ===============
int Vector::getSize() const { return _impl->_size; }
double &Vector::operator()(int i) { return _impl->_v(i); }
void Vector::print() const { std::cout << _impl->_v << std::endl; }

// ============================================================
// Implementation of SquareMatrix
// ============================================================
class SquareMatrix::SquareMatrixImpl
{
public:
  int _size;
  Eigen::MatrixXd _m;

  // Default, copy/move constructor/assignment
  SquareMatrixImpl() = default;
  ~SquareMatrixImpl() = default;
  SquareMatrixImpl(SquareMatrixImpl const &) = default;
  SquareMatrixImpl &operator=(SquareMatrixImpl const &) = default;
  SquareMatrixImpl(SquareMatrixImpl &&) noexcept = default;
  SquareMatrixImpl &operator=(SquareMatrixImpl &&) noexcept = default;

  // "Real" constructor
  SquareMatrixImpl(const int size) : _size(size) { _m = Eigen::MatrixXd::Zero(size, size); }

  void setMatrix(const Eigen::MatrixXd &other) { _m = other; };
};

// Default, copy/move constructor/assignment
SquareMatrix::SquareMatrix() : _impl(new SquareMatrixImpl{}) {}
SquareMatrix::~SquareMatrix() = default;
SquareMatrix::SquareMatrix(SquareMatrix &&) noexcept = default;
SquareMatrix &SquareMatrix::operator=(SquareMatrix &&) noexcept = default;

// "Real" constructor
SquareMatrix::SquareMatrix(const int size) : _impl(std::make_unique<SquareMatrixImpl>(size)) {}

// =============== SquareMatrix member functions ===============
int SquareMatrix::getSize() const { return _impl->_size; }

double &SquareMatrix::operator()(int i, int j) { return _impl->_m(i, j); }

// Matrix-vector product
Vector SquareMatrix::operator*(const Vector &v)
{
#if defined(FENG_DEBUG)
  assert(_impl->_size == v.getSize());
#endif
  Vector res(v._impl->_size);
  res._impl->setVector(_impl->_m * v._impl->_v);
  return res;
}

void SquareMatrix::print() const { std::cout << _impl->_m << std::endl; }

SquareMatrix SquareMatrix::inverse() const
{
  SquareMatrix inv(_impl->_size);
  inv._impl->setMatrix(this->_impl->_m.inverse());
  return inv;
}
