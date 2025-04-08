#include "feSolutionContainer.h"
#include "feLinearSystem.h"

void feSolutionContainer::NaNify()
{
  double myNan = nan("");

  _nDofs = myNan;
  _nSol = myNan;
  for(auto &val : _t)
    val = myNan;
  for(auto &val : _deltaT)
    val = myNan;
  for(auto &val : _cn)
    val = myNan;
  for(auto &val : _d)
    val = myNan;
  for(auto &vec : _sol) {
    for(auto &val : vec)
      val = myNan;
  }
  for(auto &vec : _solDot) {
    for(auto &val : vec)
      val = myNan;
  }
  for(auto &vec : _fResidual) {
    for(auto &val : vec)
      val = myNan;
  }
}

feSolutionContainer::feSolutionContainer(int nSol, double tn, int nDOF)
  : _nDofs(nDOF), _nSol(nSol)
{
  _t.resize(_nSol);
  _deltaT.resize(_nSol);
  _t[0] = tn;
  _sol.resize(_nSol);
  _solDot.resize(_nSol);
  _fResidual.resize(_nSol);

  for(int i = 0; i < _nSol; ++i) {
    _sol[i].resize(_nDofs, 0.);
    _solDot[i].resize(_nDofs, 0.);
    _fResidual[i].resize(_nDofs, 0.);
  }

  _cn.resize(3, 0.);
  _d.resize(_nDofs);
}

feSolutionContainer::feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber)
  : _nDofs(metaNumber->getNbDOFs()), _nSol(nSol)
{
  _t.resize(_nSol);
  _deltaT.resize(_nSol);
  _t[0] = tn;
  _sol.resize(_nSol);
  _solDot.resize(_nSol);
  _fResidual.resize(_nSol);

  for(int i = 0; i < _nSol; ++i) {
    _sol[i].resize(_nDofs, 0.);
    _solDot[i].resize(_nDofs, 0.);
    _fResidual[i].resize(_nDofs, 0.);
  }

  _cn.resize(3, 0.);
  _d.resize(_nDofs);
}

void feSolutionContainer::setCurrentSolution(const feSolution *other)
{
  _sol[0] = other->getSolutionCopy();
}

void feSolutionContainer::setCurrentSolutionDot(const feSolution *other)
{
  _solDot[0] = other->getSolutionDotCopy();
}

void feSolutionContainer::computeBDFCoefficients(const int order, const std::vector<double> &deltaT)
{
  if(order == 0)
  {
    _cn[0] = 0.0;
  }
  else if(order == 1)
  {
    _cn[0] =  1.0/deltaT[0];
    _cn[1] = -1.0/deltaT[0];
  }
  else if(order == 2)
  {
    _cn[0] =  1.0/deltaT[0] + 1.0/(deltaT[0]+deltaT[1]);
    _cn[1] = -1.0/deltaT[0] - 1.0/(deltaT[1]);
    _cn[2] =  deltaT[0]/deltaT[1]* 1.0/(deltaT[0]+deltaT[1]);
  }
}

void feSolutionContainer::setBDFCoefficients(const std::vector<double> &coeff)
{
  size_t size = fmin(_cn.size(), coeff.size());
  for(size_t i = 0; i < size; ++i) {
    _cn[i] = coeff[i];
  } 
}

void feSolutionContainer::computeCurrentSolDot(feLinearSystem* /* linearSystem */)
{
  // Estimate time derivative with BDF coefficients
  for(int i = 0; i < _nDofs; ++i) {
    _solDot[0][i] = 0.;
    for(size_t j = 0; j <  _cn.size(); j++) {
      _solDot[0][i] += _cn[j]*_sol[j][i];
    }
  }

  // Uncomment this when DC methods are reimplemented
  // if(_correctionType == "SOLDOT") {
  //   for(int i = 0; i < _nDofs; ++i) {
  //     _solDot[0][i] += _d[i];
  //   }
  // } else if (_correctionType == "RESIDUAL") {
  //   linearSystem->applyCorrectionToResidual(-1.0, _d);
  // } 
}

//////////////////////////////////////////////////////////////////////////
// Everything below should be reworked

void feSolutionContainer::initialize(feSolution *sol, feMesh *mesh)
{
  sol->setCurrentTime(_t[0]);
  sol->initializeUnknowns(mesh);
  sol->initializeEssentialBC(mesh);
  _sol[0] = sol->getSolutionCopy();
  _fResidual[0].resize(_nDofs);
}

void feSolutionContainer::rotate(double dt)
{
  for(int i = _nSol - 1; i > 0; --i) {
    _t[i]         = _t[i - 1];
    _deltaT[i]    = _deltaT[i - 1];
    _sol[i]       = _sol[i - 1];
    for(int j = 0; j < _nDofs; ++j) {
      _sol[i][j]       = _sol[i - 1][j];
    }
    _solDot[i]    = _solDot[i - 1];
    _fResidual[i] = _fResidual[i - 1];
  }
  _t[0] += dt;
  _deltaT[0] = _t[0] - _t[1];
}

void feSolutionContainer::rotateWithoutTime()
{
  for(int i = _nSol - 1; i > 0; --i) {
    _t[i]         = _t[i - 1];
    _sol[i]       = _sol[i-1];
    _solDot[i]    = _solDot[i-1];
    _fResidual[i] = _fResidual[i-1];
  }
}

void feSolutionContainer::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  UNUSED(linearSystem);
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, 0.0);
}

void BDFContainer::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  UNUSED(linearSystem);
  for(int i = 0; i < _nDofs; ++i) {
    _solDot[0][i] = 0.;
    for(size_t j = 0; j <  _sol.size(); j++) {
      _solDot[0][i] += _cn[j] * _sol[j][i];
    }
    sol->setSolDotAtDOF(i, _solDot[0][i]);
  }
}

void feStationarySolution::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  UNUSED(linearSystem);
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, 0.0);
}

void feSolutionBDF1::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  UNUSED(linearSystem);
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i]);
}

void feSolutionBDF2::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  UNUSED(linearSystem);
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

  // std::vector<double> t = solBDF2->getTime();
  // int orderBDF = 2, NbCoeffBDF = orderBDF + 1;
  // int orderDC = 3;
  // t.resize(orderDC);
  // std::vector<double> cn(NbCoeffBDF, 0.);
  // std::vector<double> tt(NbCoeffBDF, 0.);
  // for(int i = 0; i < NbCoeffBDF; ++i) {
  //   tt[i] = t[i];
  // }
  // // std::cout<< " t "<< tt[0] << tt[1] << tt[2] <<std::endl;
  // double tn = t[0];
  // // Coeffs BDF
  // tableToCoeffBDF(tt, cn);

  // // Init FESOL
  // int nDOF = metaNumber->getNbDOFs();
  // for(int i = 0; i < nDOF; ++i) {
  //   sol->setSolAtDOF(i, solBDF2->_sol[0][i]);
  //   sol->setSolDotAtDOF(i, 0.);
  // }
  // sol->setC0(cn[0]);
  // sol->setCurrentTime(tn);
  // solBDF2->_cn[0] = cn[0];
  // solBDF2->_cn[1] = cn[1];
  // solBDF2->_cn[2] = cn[2];
  // // std::cout<<" cn 0 " << cn[0]<<std::endl;
  // sol->initializeEssentialBC(mesh, metaNumber);
  // sol->initializeEssentialBC(mesh, metaNumber, solBDF2);
}

void initializeDC2F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1, feSolutionDC2F *solDC2F)
{
  std::vector<double> t = solBDF1->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 2;
  t.resize(orderDC);
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = 0; i < NbCoeffBDF; ++i) {
    tt[i] = t[i];
  }
  double tn = t[0];
  double k1 = t[0] - t[1];
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);
  // std::vector<double> sub(tt.begin(), tt.begin() + 3);
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d2u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solBDF1->_fResidual[0][k];
    f[1] = solBDF1->_fResidual[1][k];
    // std::cout<<"les residus valent "<<solBDF1->_fResidual[0][k]<< " et "<<
    // solBDF1->_fResidual[1][k] <<std::endl;
    tableDD(t, f, table, delta);
    d2u = delta[0]; // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
    // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
    solDC2F->_d[k] = d2u * k1 / 2.0; 
    // std::cout<< "la correction du DC2F vaut " << d2u * k1 / 2.0 << std::endl;
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
  // std::cout<<"c0 = "<<cn[0]<<"c1 = "<<cn[1]<<std::endl;
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solDC2F);
}

void initializeBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1)
{
  std::vector<double> t = solBDF1->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 2;
  t.resize(orderDC);
  std::vector<double> cn(NbCoeffBDF, 0.);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = 0; i < NbCoeffBDF; ++i) {
    tt[i] = t[i];
  }
  // std::cout<< " t "<< tt[0] << tt[1] << tt[2] <<std::endl;
  double tn = t[0];
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);

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
  int orderBDF = 2, NbCoeffBDF = orderBDF + 1;
  int orderDC = 3;
  t.resize(orderDC);
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = 0; i < NbCoeffBDF; ++i) {
    tt[i] = t[i];
  }
  double tn = t[0];
  double k1 = t[0] - t[1];
  double k2 = t[1] - t[2];
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);
  
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  int nInc = metaNumber->getNbUnknowns();
  // std::cout<<"nInc = " <<nInc <<std::endl;
  for(int k = 0; k < nInc; ++k) {
    f[0] = solBDF2->_fResidual[0][k];
    f[1] = solBDF2->_fResidual[1][k];
    f[2] = solBDF2->_fResidual[2][k];
    // std::cout<<"les residus valent "<<f[0]<< " ,  "<< f[1]<< " et " << f[2] <<std::endl;
    tableDD(t, f, table, delta);
    // std::cout << "delta11 = " << delta[0] << " delta21 = " << delta[1]
    // << " et delta12 = " << delta[n] << std::endl;
    d3u =
      2.0 * delta[n]; // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
    // std::cout << "le derivee troisieme vaut " << d3u << std::endl;
    solDC3->_d[k] = d3u / 6.0 * k1 * (k1 + k2);
    // std::cout<<"les residus valent "<<f[0]<< " ,  "<< f[1]<< " et " << f[2] <<std::endl;
    // std::cout << "la correction du DC3 vaut " << d3u / 6.0 * k1 * (k1 + k2) << std::endl;
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
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 3;
  int shift_T= 1;
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = shift_T; i < NbCoeffBDF +shift_T; ++i) {
    tt[i-shift_T] = t[i];
    // std::cout<<"tt"<<tt[i]<<std::endl;
  }
  double tn = t[1];
  // double k1 = t[0] - t[1];
  // double k2 = t[1] - t[2];
  // double k1 = t[1] - t[2];
  double k2 = t[0] - t[1];


  // Coeffs BDF

  tableToCoeffBDF(tt, cn);
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  double d2u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC2F->_fResidual[0][k];
    f[1] = solDC2F->_fResidual[1][k];
    f[2] = solDC2F->_fResidual[2][k];
    // std::cout<<"les residus valent "<<f[0]<< "et "<< f[1]<< "et" << f[2] <<std::endl;
    tableDD(t, f, table, delta);
    // std::cout << "delta11 vaut "<< delta[0] << "et delta12 vaut quant a lui "<< delta[n] <<
    // std::endl;
    d3u =
      2.0 * delta[n]; // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
    d2u = delta[0] + delta[n] * (-k2); // l'intervalle vaut k2 pas k1
    // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
    // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
    solDC3->_d[k] = d2u * k2 / 2.0 - d3u / 6.0 * k2 * k2;
    // std::cout<< "la correction du DC3F vaut " << d2u * k2 / 2.0 - d3u / 6.0 * k2*k2 << std::endl;
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
  // std::cout<<"c_0 = "<<cn[0]<<" c_1 = "<< cn[1]<<std::endl;
  // solDC3->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh);
  sol->initializeEssentialBC(mesh, solDC3);
}

void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3)
{
  UNUSED(sol, metaNumber, mesh, solDC2F, solDC3);
  feErrorMsg(FE_STATUS_ERROR, "Réécrire les méthodes DC à partir de la version de Baptiste...");
  // std::vector<double> t = solDC2F->getTime();
  // int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  // int orderDC = 3;
  // std::vector<double> cn(NbCoeffBDF);
  // std::vector<double> tt(NbCoeffBDF, 0.);
  // for(int i = 0; i < orderDC; ++i) {
  //   tt[i] = t[i];
  // }
  // double tn = tt[0];
  // double k1 = tt[0] - tt[1];
  // double k2 = tt[1] - tt[2];
 
  // // Coeffs BDF
  // tableToCoeffBDF(tt, cn);
  // std::vector<double> f(orderDC, 0.0);
  // size_t n = t.size();
  // std::vector<double> table(n * (n + 1), 0.0);
  // std::vector<double> delta(n * (n - 1), 0.0);
  // double d3u;
  // double d2u;
  // int nInc = metaNumber->getNbUnknowns();
  // for(int k = 0; k < nInc; ++k) {
  //   f[0] = solDC2F->_fResidual[0][k];
  //   f[1] = solDC2F->_fResidual[1][k];
  //   f[2] = solDC2F->_fResidual[2][k];
  //   // solDC3->_fResidual[0][k] = solDC2F->_fResidual[0][k];
  //   // solDC3->_fResidual[1][k] = solDC2F->_fResidual[1][k];
  //   // solDC3->_fResidual[2][k] = solDC2F->_fResidual[2][k];
  //   // printf(" %12s \t %24.18e \n", "f[0]", f[0]);
  //   // printf(" %12s \t %24.18e \n", "f[1]", f[1]);
  //   // printf(" %12s \t %24.18e \n", "f[2]", f[2]);
  //   // std::cout<<"les residus valent "<<f[0]<< "  et  "<< f[1]<< "  et  " << f[2] <<std::endl;
  //   tableDD(t, f, table, delta);
  //   // std::cout << "delta11 vaut "<< delta[0] << "et delta12 vaut quant a lui "<< delta[n] <<
  //   // std::endl;
  //   d3u =
  //     2.0 * delta[n]; // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
  //   d2u = delta[0] + delta[n] * (k1 - k2);
  //   d2u = delta[0] + delta[n] * k1;
  //   // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
  //   // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
  //   solDC3->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1;
  //   // std::cout<< "la correction du DC3F vaut " << d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 <<
  //   // std::endl;
  // }
  // // Init FESOL
  // int nDOF = metaNumber->getNbDOFs();
  // for(int i = 0; i < nDOF; ++i) {
  //   sol->setSolAtDOF(i, solDC3->_sol[0][i]);
  //   sol->setSolDotAtDOF(i, 0.);
  // }
  // sol->setC0(cn[0]);
  // sol->setCurrentTime(tn);
  // solDC3->_cn[0] = cn[0];
  // solDC3->_cn[1] = cn[1];
  // // std::cout << "les coeffs sont "<< cn[0] <<" et "<< cn[1]<<std::endl;
  // // solDC3->_cn[2] = cn[2];
  // sol->initializeEssentialBC(mesh);
  // sol->initializeEssentialBC(mesh, solDC3);
}

void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3, feSolutionDCF *solDC4)
{
  UNUSED(sol, metaNumber, mesh, solBDF2, solDC3, solDC4);
}
