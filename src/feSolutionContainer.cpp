#include "feSolutionContainer.h"
#include "feLinearSystem.h"

feSolutionContainer::feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber)
  : _nDofs(metaNumber->getNbDOFs()), _nSol(nSol)
{
  _t.resize(_nSol);
  _t[0] = tn;
  _sol.resize(_nSol);
  _fResidual.resize(_nSol);
  _cn.resize(3);
  _d.resize(_nDofs);
}

void feSolutionContainer::initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber)
{
  sol->setCurrentTime(_t[0]);
  sol->initializeUnknowns(mesh);
  sol->initializeEssentialBC(mesh);
  _sol[0] = sol->getSolutionCopy();
  _fResidual[0].resize(_nDofs);
  // if(_nSol>1){
  //   for(int i = 1; i < _nSol - 1; ++i) {
  //     _sol[i].resize(_nDofs, 0.0);
  //   }
  // }
}

void feSolutionContainer::rotate(double dt)
{
  for(int i = _nSol - 1; i > 0; --i) {
    _t[i] = _t[i - 1];
    _sol[i] = _sol[i - 1];
    _fResidual[i] = _fResidual[i - 1];
  }
  _t[0] += dt;
}

void feStationarySolution::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, 0.0);
}

void feSolutionBDF1::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i]);
}

void feSolutionBDF2::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i] + _cn[2] * _sol[2][i]);
}

void feSolutionDCF::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i] + _cn[2] * _sol[2][i]);
  linearSystem->applyCorrectionToResidual(-1.0, _d);
}

void feSolutionDC2F::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  // std::cout << "_sol[0][0]   " << _sol[0][0] << "    _sol[1][0]   " << _sol[1][0] << std::endl;
  for(int i = 0; i < _nDofs; ++i) {
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i]);
  }
  linearSystem->applyCorrectionToResidual(-1.0, _d);
}

inline void tableDD(std::vector<double> &t, std::vector<double> &v, std::vector<double> &table,
                    std::vector<double> &delta)
{
  size_t n = t.size();
  size_t nDelta = n - 1;
  for(size_t i = 0; i < n; ++i) {
    table[i] = t[i];
    table[n + i] = v[i];
  }
  size_t id;
  for(size_t i = 0; i < nDelta; ++i) {
    id = 2 + i;
    for(size_t j = 0; j < nDelta - i; ++j) {
      // printf("i = %d j = %d - %f\n", i, j, (t[1+i+j] - t[j]));
      table[n * id + j] =
        (table[n * (id - 1) + j + 1] - table[n * (id - 1) + j]) / (t[1 + i + j] - t[j]);
    }
  }
  for(size_t i = 0; i < n; ++i)
    for(size_t j = 0; j < n - 1; ++j) delta[n * j + i] = table[n * (j + 2) + i];
}

inline void tableToCoeffBDF(std::vector<double> &t, std::vector<double> &coeff)
{
  size_t n = t.size();
  std::vector<double> facteur(n - 1, 0.);
  for(size_t i = 0; i < n - 1; ++i) {
    facteur[i] = 1;
    if(i > 0) {
      for(size_t j = 0; j < i; ++j) {
        facteur[i] *= t[0] - t[1 + j];
      }
    }
  }
  for(size_t i = 0; i < n; ++i) {
    coeff[i] = 0.0;
    std::vector<double> v(n, 0.0);
    v[i] = 1.0;
    std::vector<double> table(n * (n + 1), 0.0);
    std::vector<double> delta(n * (n - 1), 0.0);
    tableDD(t, v, table, delta);
    for(size_t j = 0; j < n - 1; ++j) coeff[i] += table[n * (2 + j)] * facteur[j];
  }
}

void initializeBDF2(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2)
{
  std::vector<double> t = solBDF2->getTime();
  int nLvl = 1;
  int orderBDF = 2, orderP1 = orderBDF + 1, nTMin = orderP1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeBDF2 : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) {
      cn[nLvl * k + i] = coeff[i];
    }
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solBDF2->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solBDF2->_cn[0] = cn[0];
  solBDF2->_cn[1] = cn[1];
  solBDF2->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solBDF2);
  for(int i = 0; i < nDOF; ++i) {
  }
}

void initializeDC2F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1, feSolutionDC2F *solDC2F)
{
  std::vector<double> t = solBDF1->getTime();
  int nLvl = 1;
  int orderBDF = 1, orderP1 = orderBDF + 1, nTMin = orderP1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeDC2F : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  double k1 = tt[0] - tt[1];
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) {
      cn[nLvl * k + i] = coeff[i];
    }
  }
  std::vector<double> sub(tt.begin(), tt.begin() + 3);
  std::vector<double> f(nLvl, 0.0);
  size_t n = sub.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d2u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solBDF1->_fResidual[0][k];
    f[1] = solBDF1->_fResidual[1][k];
    tableDD(sub, f, table, delta);
    d2u = delta[0];
    solDC2F->_d[k] = d2u * k1 / 2.0;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solBDF1->_sol[1][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solDC2F->_cn[0] = cn[0];
  solDC2F->_cn[1] = cn[1];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solDC2F);
}

void initializeBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1)
{
  std::vector<double> t = solBDF1->getTime();
  int nLvl = 1;
  int orderBDF = 1, orderP1 = orderBDF + 1, nTMin = orderP1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeBDF1 : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) {
      cn[nLvl * k + i] = coeff[i];
    }
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solBDF1->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solBDF1->_cn[0] = cn[0];
  solBDF1->_cn[1] = cn[1];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solBDF1);
}

void initializeDC3(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionBDF2 *solBDF2,
                   feSolutionDCF *solDC3)
{
  std::vector<double> t = solBDF2->getTime();
  int nLvl = 3;
  int orderBDF = 2, orderP1 = orderBDF + 1, nTMin = orderP1 + 2;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeDC3F : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  double k1 = tt[0] - tt[1];
  double k2 = tt[1] - tt[2];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) {
      cn[nLvl * k + i] = coeff[i];
    }
  }
  std::vector<double> sub(tt.begin(), tt.begin() + 3);
  std::vector<double> f(nLvl, 0.0);
  size_t n = sub.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solBDF2->_fResidual[0][k];
    f[1] = solBDF2->_fResidual[1][k];
    f[2] = solBDF2->_fResidual[2][k];
    tableDD(sub, f, table, delta);
    d3u = 2.0 * delta[n];
    solDC3->_d[k] = d3u / 6.0 * k1 * (k1 + k2);
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC3->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solDC3->_cn[0] = cn[0];
  solDC3->_cn[1] = cn[1];
  solDC3->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solDC3);
}

void initializeDC3F_centered(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3)
{
  std::vector<double> t = solDC2F->getTime();
  int nLvl = 1;
  int orderBDF = 1, orderP1 = orderBDF + 1, nTMin = orderP1 + 1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeDC3F : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) {
    tt[i] = t[i];
  }
  double tn = tt[1];
  double k1 = tt[0] - tt[1];
  double k2 = tt[1] - tt[2];
  // std::cout << "k1 vaut "<< k1 <<std::endl;
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) {
      cn[nLvl * k + i] = coeff[i];
    }
  }
  std::vector<double> sub(tt.begin(), tt.begin() + 3);
  std::vector<double> f(nLvl, 0.0);
  size_t n = sub.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  double d2u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC2F->_fResidual[0][k];
    f[1] = solDC2F->_fResidual[1][k];
    f[2] = solDC2F->_fResidual[2][k];
    tableDD(sub, f, table, delta);
    ;
    d3u = 2.0 * delta[n];
    d2u = delta[0] + delta[n] * (-k2); // l'intervalle vaut k2 pas k1
    solDC3->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC3->_sol[1][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solDC3->_cn[0] = cn[0];
  solDC3->_cn[1] = cn[1];
  // solDC3->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solDC3);
}

void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3)
{
  std::vector<double> t = solDC2F->getTime();
  int nLvl = 1;
  int orderBDF = 1, orderP1 = orderBDF + 1, nTMin = orderP1 + 1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeDC3F : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  double k1 = tt[0] - tt[1];
  double k2 = tt[1] - tt[2];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) {
      cn[nLvl * k + i] = coeff[i];
    }
  }
  std::vector<double> sub(tt.begin(), tt.begin() + 3);
  std::vector<double> f(nLvl, 0.0);
  size_t n = sub.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  double d2u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC2F->_fResidual[0][k];
    f[1] = solDC2F->_fResidual[1][k];
    f[2] = solDC2F->_fResidual[2][k];
    solDC3->_fResidual[0][k] = solDC2F->_fResidual[0][k];
    solDC3->_fResidual[1][k] = solDC2F->_fResidual[1][k];
    solDC3->_fResidual[2][k] = solDC2F->_fResidual[2][k];
    tableDD(sub, f, table, delta);
    d3u = 2.0 * delta[n];
    d2u = delta[0] + delta[n] * (k1 - k2);
    d2u = delta[0] + delta[n] * k1;
    solDC3->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC3->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solDC3->_cn[0] = cn[0];
  solDC3->_cn[1] = cn[1];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solDC3);
}

void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3, feSolutionDCF *solDC4)
{
}
