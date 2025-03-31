#include "feMatrixInterface.h"

// #include "../contrib/Eigen/Dense"
// #include "../contrib/unsupported/Eigen/MatrixFunctions"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

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

  // "Real" constructors
  MetricTensorImpl(const double val)
  {
    _m = Eigen::Matrix2d::Identity();
    _m(0, 0) = val;
    _m(1, 1) = val;
  }
  MetricTensorImpl(const double m11, const double m12, const double m22)
  {
    _m = Eigen::Matrix2d();
    _m(0, 0) = m11;
    _m(0, 1) = m12;
    _m(1, 0) = m12;
    _m(1, 1) = m22;
  }
  MetricTensorImpl(const double eigenvalues[2], const double eigenvector1[2], const double eigenvector2[2])
  {
    // Matrix whose columns are the eigenvectors
    Eigen::Matrix2d V, D;
    V(0,0) = eigenvector1[0];
    V(1,0) = eigenvector1[1];
    V(0,1) = eigenvector2[0];
    V(1,1) = eigenvector2[1];
    D << eigenvalues[0], 0., 0., eigenvalues[1];
    _m = V * D * V.transpose();
  }
  MetricTensorImpl(const SquareMatrix &other)
  {
    if(other.getSize() != 2) {
      printf("Cannot create MetricTensor from SquareMatrix because SquareMatrix has dimension %d\n", other.getSize());
      exit(-1);
    }
    // if(fabs(other(0,1) - other(1,0))/fabs(other(0,1)) > 1e-10) {
    //   printf("Cannot create MetricTensor from SquareMatrix because SquareMatrix is not symmetric: M12 = %+-1.12e - M21 = %+-1.12e\n",
    //     other(0,1), other(1,0));
    //   exit(-1);
    // }

    _m = Eigen::Matrix2d::Identity();
    _m(0, 0) = other(0,0);
    _m(0, 1) = other(0,1);
    _m(1, 0) = other(0,1);
    _m(1, 1) = other(1,1);
  }

  void setMatrix(const Eigen::Matrix2d &other) { _m = other; };
};

// Default, copy/move constructor/assignment
MetricTensor::~MetricTensor() = default;
MetricTensor::MetricTensor(MetricTensor &&) noexcept = default;
MetricTensor &MetricTensor::operator=(MetricTensor &&) noexcept = default;

// "Real" constructors
MetricTensor::MetricTensor(const double val) : _impl(std::make_unique<MetricTensorImpl>(val)) {}
MetricTensor::MetricTensor(const double m11, const double m12, const double m22)
 : _impl(std::make_unique<MetricTensorImpl>(m11, m12, m22)) {}
MetricTensor::MetricTensor(const double eigenvalues[2], const double eigenvector1[2], const double eigenvector2[2])
  : _impl(std::make_unique<MetricTensorImpl>(eigenvalues, eigenvector1, eigenvector2)) {}
MetricTensor::MetricTensor(const SquareMatrix &other)
: _impl(std::make_unique<MetricTensorImpl>(other)) {}

// =============== MetricTensor member functions ===============
MetricTensor MetricTensor::copy() const
{
  MetricTensor tmp(1.0);
  tmp._impl->setMatrix(_impl->_m);
  return tmp;
}

void MetricTensor::assignMatrixFrom(const MetricTensor &other)
{
  this->_impl->_m(0,0) = other._impl->_m(0,0);
  this->_impl->_m(0,1) = other._impl->_m(0,1);
  this->_impl->_m(1,0) = other._impl->_m(1,0);
  this->_impl->_m(1,1) = other._impl->_m(1,1);
}

void MetricTensor::assignMatrixFrom(const double other[2][2])
{
  this->_impl->_m(0,0) = other[0][0];
  this->_impl->_m(0,1) = other[0][1];
  this->_impl->_m(1,0) = other[1][0];
  this->_impl->_m(1,1) = other[1][1];
}

void MetricTensor::setToZero()
{
  this->_impl->_m(0,0) = 0.;
  this->_impl->_m(0,1) = 0.;
  this->_impl->_m(1,0) = 0.;
  this->_impl->_m(1,1) = 0.;
}

void MetricTensor::setToIdentity()
{
  this->_impl->_m(0,0) = 1.;
  this->_impl->_m(0,1) = 0.;
  this->_impl->_m(1,0) = 0.;
  this->_impl->_m(1,1) = 1.;
}

void MetricTensor::incrementMatrix(const double alpha, const MetricTensor &other)
{
  this->_impl->_m(0,0) += alpha * other._impl->_m(0,0);
  this->_impl->_m(0,1) += alpha * other._impl->_m(0,1);
  this->_impl->_m(1,0) += alpha * other._impl->_m(1,0);
  this->_impl->_m(1,1) += alpha * other._impl->_m(1,1);
}

double &MetricTensor::operator()(int i, int j) { return _impl->_m(i, j); }
double MetricTensor::operator()(int i, int j) const { return _impl->_m(i, j); }
MetricTensor &MetricTensor::operator*=(const double &val)
{
  _impl->_m *= val;
  return *this;
}

MetricTensor MetricTensor::operator*(const double &val) const
{
  MetricTensor tmp(1.0);
  tmp._impl->setMatrix(_impl->_m * val);
  return tmp;
}

MetricTensor MetricTensor::operator/(const double &val) const
{
  MetricTensor tmp(1.0);
  tmp._impl->setMatrix(_impl->_m / val);
  return tmp;
}

MetricTensor MetricTensor::operator+(const MetricTensor &other) const
{
  MetricTensor tmp(1.0);
  tmp._impl->setMatrix(_impl->_m + other._impl->_m);
  return tmp;
}

MetricTensor MetricTensor::operator-(const MetricTensor &other) const
{
  MetricTensor tmp(1.0);
  tmp._impl->setMatrix(_impl->_m - other._impl->_m);
  return tmp;
}

MetricTensor MetricTensor::operator*(const MetricTensor &other) const
{
  MetricTensor tmp(1.0);
  tmp._impl->setMatrix(_impl->_m * other._impl->_m);
  return tmp;
}

void MetricTensor::print(const int numDigits) const { std::cout << std::setprecision(numDigits) <<_impl->_m << std::endl; }

double MetricTensor::trace() const { return _impl->_m.trace(); }

double MetricTensor::determinant() const { return _impl->_m.determinant(); }

MetricTensor MetricTensor::inverse() const
{
  MetricTensor inv(1.0);
  inv._impl->setMatrix(this->_impl->_m.inverse());
  return inv;
}

MetricTensor MetricTensor::transpose() const
{
  MetricTensor tr(1.0);
  tr._impl->setMatrix(this->_impl->_m.transpose());
  return tr;
}

// #define INLINED_LOG

MetricTensor MetricTensor::log() const
{
  MetricTensor logm(1.0);
  #if defined(INLINED_LOG)
  double a = this->_impl->_m(0,0);
  double b = this->_impl->_m(0,1);
  double c = this->_impl->_m(1,1);
  double r4 = sqrt(a*a - 2.*a*c + 4.*b*b + c*c);
  double r3 = std::log( (a + c + r4) / 2. );
  double r2 = std::log( (a + c - r4) / 2. );
  double r1 = -b * (r2 - r3) / r4;
  logm(0,0) = (r2 * (c - a + r4) + r3 * (a - c + r4)) / (2.*r4);
  logm(0,1) = r1;
  logm(1,0) = r1;
  logm(1,1) = (r2 * (a - c + r4) + r3 * (c - a + r4)) / (2.*r4);
  #else
  logm._impl->setMatrix(this->_impl->_m.log());
  #endif
  return logm;
}

MetricTensor MetricTensor::exp() const
{
  MetricTensor expm(1.0);
  expm._impl->setMatrix(this->_impl->_m.exp());
  return expm;
}

MetricTensor MetricTensor::pow(const double &power) const
{
  MetricTensor powm(1.0);
  powm._impl->setMatrix(this->_impl->_m.pow(power));
  return powm;
}

MetricTensor MetricTensor::sqrt() const
{
  MetricTensor sqrtm(1.0);
  sqrtm._impl->setMatrix(this->_impl->_m.sqrt());
  return sqrtm;
}

double MetricTensor::frobeniusNorm() const
{
  return this->_impl->_m.norm();
}

double MetricTensor::dotProduct(const double x1[2], const double x2[2]) const
{
  return x2[0] * (this->_impl->_m(0,0) * x1[0] + this->_impl->_m(0,1) * x1[1])
       + x2[1] * (this->_impl->_m(0,1) * x1[0] + this->_impl->_m(1,1) * x1[1]);
}

double MetricTensor::maxCoeff() const
{
  return this->_impl->_m.maxCoeff();
}

double MetricTensor::getMinEigenvalue() const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  return fmin(ev(0).real(), ev(1).real());
}

double MetricTensor::getMaxEigenvalue() const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  return fmax(ev(0).real(), ev(1).real());
}

void MetricTensor::getEigenvectorsAndEigenvalues(double v1[2], double v2[2], double &lambda1, double &lambda2) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  lambda1 = ev(0).real();
  lambda2 = ev(1).real();
  v1[0] = es.eigenvectors().real().col(0)[0];
  v1[1] = es.eigenvectors().real().col(0)[1];
  v2[0] = es.eigenvectors().real().col(1)[0];
  v2[1] = es.eigenvectors().real().col(1)[1];
}

void MetricTensor::symmetrize()
{
  this->_impl->setMatrix( 0.5 * (this->_impl->_m + this->_impl->_m.transpose()) );
}

MetricTensor MetricTensor::absoluteValueEigen() const
{
  if(_impl->_m.isIdentity(1e-12)) {
    MetricTensor res(1.0);
    return res;
  } else {
    Eigen::EigenSolver<Eigen::Matrix2d> es;
    Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
    es.compute(_impl->_m, true);
    ev = es.eigenvalues();
    Eigen::Matrix2d D = Eigen::Matrix2d::Identity();
    D(0, 0) = fabs(ev(0).real());
    D(1, 1) = fabs(ev(1).real());
    Eigen::MatrixXcd sol = es.eigenvectors() * D * es.eigenvectors().inverse();
    if(sol.imag().norm() > 1e-14) {
      printf("In MetricTensor::absoluteValueEigen: Matrix has nontrivial imaginary part!");
      std::cout << sol << std::endl;
      exit(-1);
    }
    MetricTensor res(1.0);
    res._impl->setMatrix(sol.real());
    return res;
  }
}

MetricTensor MetricTensor::boundEigenvalues(const double lMin, const double lMax) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  D(0, 0) = fmin(lMax, fmax(lMin, ev(0).real()));
  D(1, 1) = fmin(lMax, fmax(lMin, ev(1).real()));
  Eigen::MatrixXcd sol = es.eigenvectors() * D * es.eigenvectors().inverse();
  if(sol.imag().norm() > 1e-14) {
    printf("In MetricTensor::boundEigenvalues: Matrix has nontrivial imaginary part after bounding eigenvalues!");
    std::cout << sol << std::endl;
    exit(-1);
  }
  MetricTensor res(1.0);
  res._impl->setMatrix(sol.real());
  return res;
}

MetricTensor MetricTensor::boundEigenvaluesOfAbs(const double lMin, const double lMax) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  D(0, 0) = fmin(lMax, fmax(lMin, fabs(ev(0).real())));
  D(1, 1) = fmin(lMax, fmax(lMin, fabs(ev(1).real())));
  Eigen::MatrixXcd sol = es.eigenvectors() * D * es.eigenvectors().inverse();
  if(sol.imag().norm() > 1e-14) {
    printf("In MetricTensor::boundEigenvaluesOfAbs: Matrix has nontrivial imaginary part after bounding eigenvalues!");
    std::cout << sol << std::endl;
    exit(-1);
  }
  MetricTensor res(1.0);
  res._impl->setMatrix(sol.real());
  return res;
}

MetricTensor MetricTensor::boundEigenvaluesOfAbsIsotropic(const double lMin, const double lMax) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  // Take most critical (largest) eigenvalue
  double lambda = fmax(ev(0).real(), ev(1).real());
  D(0, 0) = fmin(lMax, fmax(lMin, fabs( lambda )));
  D(1, 1) = fmin(lMax, fmax(lMin, fabs( lambda )));
  if(D.imag().norm() > 1e-14) {
    printf("In MetricTensor::boundEigenvaluesOfAbsIsotropic: Matrix has nontrivial imaginary part after bounding eigenvalues!");
    std::cout << D << std::endl;
    exit(-1);
  }
  MetricTensor res(1.0);
  res._impl->setMatrix(D.real());
  return res;
}

MetricTensor MetricTensor::setEigenvectorsAndBoundEigenvalues(const double ev1[2],
  const double ev2[2],const double lMin, const double lMax) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  D(0, 0) = fmin(lMax, fmax(lMin, fabs(ev(0).real())));
  D(1, 1) = fmin(lMax, fmax(lMin, fabs(ev(1).real())));
  Eigen::Matrix2d newEigenvectors;
  newEigenvectors(0,0) = ev1[0];
  newEigenvectors(1,0) = ev1[1];
  newEigenvectors(0,1) = ev2[0];
  newEigenvectors(1,1) = ev2[1];
  MetricTensor res(1.0);
  res._impl->setMatrix(newEigenvectors * D.real() * newEigenvectors.transpose());
  return res;
}

MetricTensor MetricTensor::limitAnisotropy(const double alpha) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  double lambda = ev(0).real();
  double mu     = ev(1).real();
  if(lambda >= mu) {
    D(0, 0) = lambda;
    D(1, 1) = fmax(mu, lambda / (alpha*alpha));
  } else {
    D(0, 0) = fmax(lambda, mu / (alpha*alpha));
    D(1, 1) = mu;
  }
  Eigen::MatrixXcd sol = es.eigenvectors() * D * es.eigenvectors().inverse();
  if(sol.imag().norm() > 1e-14) {
    printf("In MetricTensor::limitAnisotropy: Matrix has nontrivial imaginary part after bounding eigenvalues!");
    std::cout << sol << std::endl;
    exit(-1);
  }
  MetricTensor res(1.0);
  res._impl->setMatrix(sol.real());
  return res;
}

MetricTensor MetricTensor::spanMetricInMetricSpace(const double gradation, const double pq[2]) const
{
  double eta = (1. + std::sqrt(this->dotProduct(pq, pq)) * std::log(gradation));
  eta = 1. / (eta * eta);
  MetricTensor res(1.0);
  res._impl->setMatrix(this->_impl->_m * eta);
  return res;
}

MetricTensor MetricTensor::spanMetricInPhysicalSpace(const double gradation, const double pq[2]) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(this->_impl->_m, true);
  ev = es.eigenvalues();
  double l1 = ev(0).real();
  double l2 = ev(1).real();
  double normPQ = std::sqrt(pq[0]*pq[0] + pq[1]*pq[1]);
  double eta1 = 1. + std::sqrt(l1) * normPQ * std::log(gradation); eta1 = 1. / (eta1*eta1);
  double eta2 = 1. + std::sqrt(l2) * normPQ * std::log(gradation); eta2 = 1. / (eta2*eta2);
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  D(0, 0) = eta1 * l1;
  D(1, 1) = eta2 * l2;
  Eigen::MatrixXcd sol = es.eigenvectors() * D * es.eigenvectors().inverse();
  if(sol.imag().norm() > 1e-14) {
    printf("In MetricTensor::spanMetricInPhysicalSpace: Matrix has nontrivial imaginary part after bounding eigenvalues!");
    std::cout << sol << std::endl;
    exit(-1);
  }
  MetricTensor res(1.0);
  res._impl->setMatrix(sol.real());
  return res;
}

MetricTensor MetricTensor::spanMetricInMixedSpace(const double gradation, const double pq[2], const double t) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(this->_impl->_m, true);
  ev = es.eigenvalues();
  double l1 = ev(0).real();
  double l2 = ev(1).real();
  double normPQ = std::sqrt(pq[0]*pq[0] + pq[1]*pq[1]);
  double etaMetric = 1. + std::sqrt(this->dotProduct(pq, pq)) * std::log(gradation);
  double etaPhysical1 = 1. + std::sqrt(l1) * normPQ * std::log(gradation);
  double etaPhysical2 = 1. + std::sqrt(l2) * normPQ * std::log(gradation);
  double eta1 = std::pow(etaMetric, 1.-t) * std::pow(etaPhysical1, t); eta1 = 1. / (eta1*eta1);
  double eta2 = std::pow(etaMetric, 1.-t) * std::pow(etaPhysical2, t); eta2 = 1. / (eta2*eta2);
  Eigen::Matrix2d D = Eigen::Matrix2d::Identity();  
  D(0, 0) = eta1 * l1;
  D(1, 1) = eta2 * l2;
  Eigen::MatrixXcd sol = es.eigenvectors() * D * es.eigenvectors().inverse();
  if(sol.imag().norm() > 1e-14) {
    printf("In MetricTensor::spanMetricInMixedSpace: Matrix has nontrivial imaginary part after bounding eigenvalues!");
    std::cout << sol << std::endl;
    exit(-1);
  }
  MetricTensor res(1.0);
  res._impl->setMatrix(sol.real());
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
double SquareMatrix::operator()(int i, int j) const { return _impl->_m(i, j); }

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

// Product of this (2x2) SquareMatrix with a MetricTensor (symmetric 2x2 matrix)
SquareMatrix SquareMatrix::operator*(const MetricTensor &other) const
{
  if(_impl->_size != 2) {
    printf("Cannot multiply this SquareMatrix by MetricTensor because SquareMatrix has dimension %d\n", _impl->_size);
    exit(-1);
  }
  SquareMatrix tmp(2);
  tmp._impl->setMatrix(_impl->_m * other._impl->_m);
  return tmp;
}

// Product of this SquareMatrix with another SquareMatrix of same dimensions
SquareMatrix SquareMatrix::operator*(const SquareMatrix &other) const
{
  if(_impl->_size != other._impl->_size) {
    printf("Cannot multiply two SquareMatrix because their dimensions are %d and %d\n",_impl-> _size, other._impl->_size);
    exit(-1);
  }
  SquareMatrix tmp(_impl->_size);
  tmp._impl->setMatrix(_impl->_m * other._impl->_m);
  return tmp;
}

SquareMatrix MetricTensor::operator*(const SquareMatrix &other) const
{
  if(other.getSize() != 2) {
    printf("Cannot multiply this MetricTensor by SquareMatrix because SquareMatrix has dimension %d\n", other.getSize());
    exit(-1);
  }
  SquareMatrix tmp(2);
  tmp._impl->setMatrix(_impl->_m * other._impl->_m);
  return tmp;
}

void SquareMatrix::print() const { std::cout << _impl->_m << std::endl; }

double SquareMatrix::determinant() const { return _impl->_m.determinant(); }

SquareMatrix SquareMatrix::inverse() const
{
  SquareMatrix inv(_impl->_size);
  inv._impl->setMatrix(this->_impl->_m.inverse());
  return inv;
}

void SquareMatrix::inverse(SquareMatrix &res) const
{
  res._impl->setMatrix(this->_impl->_m.inverse());
}

SquareMatrix SquareMatrix::transpose() const
{
  SquareMatrix tr(_impl->_size);
  tr._impl->setMatrix(this->_impl->_m.transpose());
  return tr;
}

int SquareMatrix::rank() const
{
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(this->_impl->_m);
  return lu_decomp.rank();
}
void SquareMatrix::getEigenvectorsAndEigenvalues(double v1[2], double v2[2], double &lambda1, double &lambda2) const
{
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  es.compute(_impl->_m, true);
  ev = es.eigenvalues();
  lambda1 = ev(0).real();
  lambda2 = ev(1).real();
  v1[0] = es.eigenvectors().real().col(0)[0];
  v1[1] = es.eigenvectors().real().col(0)[1];
  v2[0] = es.eigenvectors().real().col(1)[0];
  v2[1] = es.eigenvectors().real().col(1)[1];
}

// ============================================================
// Implementation of Matrix
// ============================================================
class Matrix::MatrixImpl
{
public:
  int _sizeM;
  int _sizeN;
  Eigen::MatrixXd _m;

  // Default, copy/move constructor/assignment
  MatrixImpl() = default;
  ~MatrixImpl() = default;
  MatrixImpl(MatrixImpl const &) = default;
  MatrixImpl &operator=(MatrixImpl const &) = default;
  MatrixImpl(MatrixImpl &&) noexcept = default;
  MatrixImpl &operator=(MatrixImpl &&) noexcept = default;

  // "Real" constructor
  MatrixImpl(const int sizeM, const int sizeN) : _sizeM(sizeM), _sizeN(sizeN)
    { _m = Eigen::MatrixXd::Zero(sizeM, sizeN); }

  void setMatrix(const Eigen::MatrixXd &other) { _m = other; };
};

// Default, copy/move constructor/assignment
Matrix::Matrix() : _impl(new MatrixImpl{}) {}
Matrix::~Matrix() = default;
Matrix::Matrix(Matrix &&) noexcept = default;
Matrix &Matrix::operator=(Matrix &&) noexcept = default;

// "Real" constructor
Matrix::Matrix(const int sizeM, const int sizeN) : _impl(std::make_unique<MatrixImpl>(sizeM, sizeN)) {}

// =============== Matrix member functions ===============
int Matrix::getSizeM() const { return _impl->_sizeM; }
int Matrix::getSizeN() const { return _impl->_sizeN; }

double &Matrix::operator()(int i, int j) { return _impl->_m(i, j); }

// Matrix-vector product
Vector Matrix::operator*(const Vector &v) const
{
#if defined(FENG_DEBUG)
  assert(_impl->_sizeN == v.getSize());
#endif
  Vector res(this->_impl->_sizeM);
  res._impl->setVector(_impl->_m * v._impl->_v);
  return res;
}

void Matrix::print() const { std::cout << _impl->_m << std::endl; }

SquareMatrix Matrix::multiplyByTranspose() const
{
  SquareMatrix res(_impl->_sizeN);
  res._impl->setMatrix(this->_impl->_m.transpose() * this->_impl->_m);
  return res;
}

Matrix Matrix::getLeastSquaresMatrix() const
{
  Matrix res(_impl->_sizeN, _impl->_sizeM);
  SquareMatrix inv(_impl->_sizeN);
  inv._impl->setMatrix(this->_impl->_m.transpose() * this->_impl->_m);
  inv._impl->setMatrix(inv._impl->_m.inverse());
  res._impl->setMatrix(inv._impl->_m * this->_impl->_m.transpose());
  return res;
}

int Matrix::rank() const
{
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(this->_impl->_m);
  return lu_decomp.rank();
}