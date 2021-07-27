#include "feSolutionContainer.h"
#include "feLinearSystem.h"

feSolutionContainer::feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber)
  : _nDofs(metaNumber->getNbDOFs()), _nSol(nSol) {
  _t.resize(_nSol);
  _t[0] = tn;
  _sol.resize(_nSol);
  _fResidual.resize(_nSol);
  _cn.resize(3);
  _d.resize(_nDofs);
}

void feSolutionContainer::initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber) {
  sol->setCurrentTime(_t[0]);
  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);
  _sol[0] = sol->getSolutionCopy();
  _fResidual[0].resize(_nDofs);
}

void feSolutionContainer::rotate(double dt) {
  for(int i = _nSol - 1; i > 0; --i) {
    _t[i] = _t[i - 1];
    _sol[i] = _sol[i - 1];
    _fResidual[i] = _fResidual[i - 1];
  }
  _t[0] += dt;
}

void feStationarySolution::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem) {
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, 0.0);
}

void feSolutionBDF2::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem) {
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i] + _cn[2] * _sol[2][i]);
  // std::cout<<"coeffs"<<_cn[0]<<","<<_cn[1]<<","<<_cn[2]<<"timestep"<<_t[0]<<_t[1]<<std::endl;
}

void feSolutionBDF1::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem) {
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i]);

  // std::cout<<"coeffs"<<_cn[0]<<","<<_cn[1]<<","<<_cn[2]<<"timestep"<<_t[0]<<_t[1]<<std::endl;
}

void feSolutionDCF::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem) {
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i] + _cn[2] * _sol[2][i]);
  linearSystem->applyCorrectionToResidual(-1.0, _d);
}

inline void tableDD(std::vector<double> &t, std::vector<double> &v, std::vector<double> &table,
                    std::vector<double> &delta) {
  size_t n = t.size();
  size_t nDelta = n - 1;
  for(size_t i = 0; i < n; ++i) {
    table[i] = t[i];
    table[n + i] = v[i];
  }
  size_t id;
  for(size_t i = 0; i < nDelta; ++i) {
    for(size_t j = 0; j < nDelta - i; ++j) {
      id = 2 + i;
      // printf("i = %d j = %d - %f\n", i, j, (t[1+i+j] - t[j]));
      table[n * id + j] =
        (table[n * (id - 1) + j + 1] - table[n * (id - 1) + j]) / (t[1 + i + j] - t[j]);
    }
  }
  for(size_t i = 0; i < n; ++i)
    for(size_t j = 0; j < n - 1; ++j) delta[n * j + i] = table[n * (j + 2) + i];
}

inline void tableToCoeffBDF(std::vector<double> &t, std::vector<double> &coeff) {
  size_t n = t.size();
  std::vector<double> facteur(n - 1, 0.);
  for(size_t i = 0; i < n - 1; ++i) {
    facteur[i] = 1;
    if(i > 0) {
      for(size_t j = 0; j < i; ++j) { facteur[i] *= t[0] - t[1 + j]; }
    }
  }
  // for(auto val : facteur) std::cout<<val<<" "; std::cout<<std::endl;
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
                    feSolutionBDF2 *solBDF2) {
  std::vector<double> t = solBDF2->getTime();
  int nLvl = 1;
  int orderBDF = 2, orderP1 = orderBDF + 1, nTMin = orderP1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeBDF2 : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  // double k1 = tt[0] - tt[1];
  // double k2 = tt[1] - tt[2];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) { cn[nLvl * k + i] = coeff[i]; }
  }
  // for(auto val : cn) std::cout<<val<<" "; std::cout<<std::endl;
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
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solBDF2);
  for(int i = 0; i < nDOF; ++i) {}
}

void initializeBDF2withBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                            feSolutionBDF1 *solBDF1, feSolutionBDF2 *solBDF2) {
  std::vector<double> t = solBDF2->getTime();
  int nLvl = 1;
  int orderBDF = 2, orderP1 = orderBDF + 1, nTMin = orderP1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeBDF2 : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  // double k1 = tt[0] - tt[1];
  // double k2 = tt[1] - tt[2];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) { cn[nLvl * k + i] = coeff[i]; }
  }
  // for(auto val : cn) std::cout<<val<<" "; std::cout<<std::endl;
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solBDF1->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solBDF2->_cn[0] = cn[0];
  solBDF2->_cn[1] = cn[1];
  solBDF2->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solBDF1);
  for(int i = 0; i < nDOF; ++i) {
    std::cout << "sol dans l'ini" << sol->getSolAtDOF(i) << std::endl;
  }
}

void initializeBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1) {
  std::vector<double> t = solBDF1->getTime();
  int nLvl = 1;
  int orderBDF = 1, orderP1 = orderBDF + 1, nTMin = orderP1;
  std::vector<double> cn(nLvl * 3);
  if((int)t.size() < nTMin) printf("In feSolutionContainer::initializeBDF1 : t.size < nTMin\n");
  std::vector<double> tt(nTMin, 0.);
  for(int i = 0; i < nTMin; ++i) tt[i] = t[i];
  double tn = tt[0];
  // double k1 = tt[0] - tt[1];
  // double k2 = tt[1] - tt[2];
  // Coeffs BDF
  for(int k = 0; k < nLvl; ++k) {
    std::vector<double> sub(tt.begin() + (1 + k - 1), tt.begin() + (orderP1 + k - 1) + 1);
    std::vector<double> coeff(sub.size(), 0.);
    tableToCoeffBDF(sub, coeff);
    for(int i = 0; i < orderP1; ++i) { cn[nLvl * k + i] = coeff[i]; }
  }
  // for(auto val : cn) std::cout<<val<<" "; std::cout<<std::endl;
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solBDF1->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  // solBDF2->_cn[0] = cn[0];
  // solBDF2->_cn[1] = cn[1];
  // solBDF2->_cn[2] = cn[2];
  solBDF1->_cn[0] = cn[0];
  solBDF1->_cn[1] = cn[1];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solBDF1);
  // for(int i = 0; i < nDOF; ++i){
  //   std::cout<<sol->getSolAtDOF(i)<<std::endl;
  // }
}

void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3) {
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
    for(int i = 0; i < orderP1; ++i) { cn[nLvl * k + i] = coeff[i]; }
  }
  std::vector<double> sub(tt.begin(), tt.begin() + 3);
  // std::cout<<"sub.size() = "<<sub.size()<<std::endl;
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
    d3u =
      2.0 * delta[n]; // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
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
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC3);
}

void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3, feSolutionDCF *solDC4) {}
