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
  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);
  _sol[0] = sol->getSolutionCopy();
  for(int i = 0; i < _nSol - 1; ++i) {
    _sol[i].resize(_nDofs, 0.0);
  }
  _fResidual[0].resize(_nDofs);
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

void feSolutionBDF2::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  printf(" %12s \t %24.18e \n", "sol[0][0]", _sol[0][0]);
  printf(" %12s \t %24.18e \n", "sol[0][1]", _sol[1][0]);
  printf(" %12s \t %24.18e \n", "sol[0][2]", _sol[2][0]);
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i] + _cn[2] * _sol[2][i]);
}

void feSolutionBDF1::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  // printf(" %12s \t %24.18e \n", "sol[0][0]", _sol[0][0]);
  // printf(" %12s \t %24.18e \n", "sol[0][1]", _sol[1][0]);
  // std::cout<<"nb de sol" << _nDofs << std::endl;
  for(int i = 0; i < _nDofs; ++i) sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i]);
}

void feSolutionDCF::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{
  printf(" %12s \t %24.18e \n", "sol[0][0]", _sol[0][0]);
  printf(" %12s \t %24.18e \n", "sol[0][1]", _sol[1][0]);
  printf(" %12s \t %24.18e \n", "sol[0][2]", _sol[2][0]);
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i] + _cn[2] * _sol[2][i]);
  linearSystem->applyCorrectionToResidual(-1.0, _d);
}

void feSolutionDC2F::computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem)
{ 
  // printf(" %12s \t %24.18e \n", "sol[0][0]", _sol[0][0]);
  // printf(" %12s \t %24.18e \n", "sol[0][1]", _sol[1][0]);
  // printf(" %12s \t %24.18e \n", "sol[0][2]", _sol[2][0]);
  // printf(" %12s \t %24.18e \n", "sol[0][3]", _sol[3][0]);
  // printf(" %12s \t %24.18e \n", "sol[0][4]", _sol[4][0]);
  for(int i = 0; i < _nDofs; ++i) {
    sol->setSolDotAtDOF(i, _cn[0] * _sol[0][i] + _cn[1] * _sol[1][i]);
  }
  linearSystem->applyCorrectionToResidual(-1.0, _d);
}

// used for Newton’s Divided Difference Interpolation Formula
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

//get the coeff of BDF using tableDD
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

//Set the coeff for BDF2
void initializeBDF2(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2)
{
  std::vector<double> t = solBDF2->getTime();
  int orderBDF = 2, NbCoeffBDF = orderBDF + 1;
  int orderDC = 3;
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
    sol->setSolAtDOF(i, solBDF2->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solBDF2->_cn[0] = cn[0];
  solBDF2->_cn[1] = cn[1];
  solBDF2->_cn[2] = cn[2];
  // std::cout<<" cn 0 " << cn[0]<<std::endl;
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solBDF2);
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
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC2F);
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
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solBDF1);
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
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC3);
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
  double k1 = t[1] - t[2];
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
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC3);
}

// void initializeDC4F_t1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
//                              feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4)
// {
//   std::vector<double> t = solDC3F->getTime();
//   int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
//   int orderDC = 4;
//   int shift_T= 2;
//   std::vector<double> cn(NbCoeffBDF);
//   std::vector<double> tt(NbCoeffBDF, 0.);
//     for(int i = shift_T; i < NbCoeffBDF+shift_T; ++i) {
//     tt[i-shift_T] = t[i];
//   }
//   double tn = t[2];
//   double k1 = t[0] - t[1];
//   double k2 = t[1] - t[2];
//   double k3 = t[2] - t[3];
//   // std::cout<<" k1 = "<<k1 << " k2 = "<< k2 << "k3 = "<<k3<<std::endl;
 
//   // Coeffs BDF
//   tableToCoeffBDF(tt, cn);
//   std::vector<double> f(orderDC, 0.0);
//   size_t n = t.size();
//   // std::cout<<"=============size of t  = "<<n<<std::endl;
//   std::vector<double> table(n * (n + 1), 0.0);
//   std::vector<double> delta(n * (n - 1), 0.0);
//   double d2u;
//   double d3u;
//   double d4u;
//   int nInc = metaNumber->getNbUnknowns();
//   for(int k = 0; k < nInc; ++k) {
//     f[0] = solDC3F->_fResidual[0][k];
//     f[1] = solDC3F->_fResidual[1][k];
//     f[2] = solDC3F->_fResidual[2][k];
//     f[3] = solDC3F->_fResidual[3][k];
//     // std::cout<<"les residus valent "<<f[0]<< " "<< f[1]<< "," << f[2] <<"et" << f[3]<<std::endl;
//     tableDD(t, f, table, delta);
//     // std::cout << "delta11 = "<< delta[0] << " delta12 = "<< delta[n]<< " et delta13 ="<< delta[2*n] <<
//     // std::endl;
//     d4u=6*delta[2*n];
//     d3u =
//       2.0 * delta[n] + 2*delta[2*n]*(-k1 - 2*k2); // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
//     d2u = delta[0] + delta[n]*(-k1 - 2*k2) + delta[2*n]*k2*(k1+k2);
//     // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
//     // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
//     // std::cout<< "le derivee quatrieme vaut " << d4u << std::endl;
//     solDC4->_d[k] = d2u * k3 / 2.0 - d3u / 6.0 * k3 * k3 + d4u / 24. * k3*k3*k3 ;
//     // std::cout<< "la correction du DC4F vaut " << d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 + d4u / 24. * k1*k1*k1  << std::endl;
//   }
//   // Init FESOL
//   int nDOF = metaNumber->getNbDOFs();
//   for(int i = 0; i < nDOF; ++i) {
//     sol->setSolAtDOF(i, solDC4->_sol[1][i]);
//     sol->setSolDotAtDOF(i, 0.);
//   }
//   sol->setC0(cn[0]);
//   sol->setCurrentTime(tn);
//   solDC4->_cn[0] = cn[0];
//   solDC4->_cn[1] = cn[1];
//   // std::cout<<"c_0 = "<<cn[0]<<" c_1 = "<< cn[1]<<std::endl;
//   // solDC3->_cn[2] = cn[2];
//   sol->initializeEssentialBC(mesh, metaNumber);
//   sol->initializeEssentialBC(mesh, metaNumber, solDC4);
// }

// void initializeDC4F_t2(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
//                              feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4)
// {
//   std::vector<double> t = solDC3F->getTime();
//   int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
//   int orderDC = 4;
//   int shift_T = 1;
//   std::vector<double> cn(NbCoeffBDF);
//   std::vector<double> tt(NbCoeffBDF, 0.);
//     for(int i = shift_T; i < NbCoeffBDF+shift_T; ++i) {
//       // std::cout<<"coucou"<<t[i]<<std::endl;
//       tt[i-shift_T] = t[i];
//   }
//   double tn = t[1];
//   double k1 = t[0] - t[1];
//   double k2 = t[1] - t[2];
//   double k3 = t[2] - t[3];
//   // std::cout<<"=============k2 = "<<k2<<std::endl;
 
//   // Coeffs BDF
//   tableToCoeffBDF(tt, cn);
//   std::vector<double> f(orderDC, 0.0);
//   size_t n = t.size();
//   std::vector<double> table(n * (n + 1), 0.0);
//   std::vector<double> delta(n * (n - 1), 0.0);
//   double d3u;
//   double d2u;
//   double d4u;
//   int nInc = metaNumber->getNbUnknowns();
//   for(int k = 0; k < nInc; ++k) {
//     f[0] = solDC3F->_fResidual[0][k];
//     f[1] = solDC3F->_fResidual[1][k];
//     f[2] = solDC3F->_fResidual[2][k];
//     f[3] = solDC3F->_fResidual[3][k];
//     // std::cout<<"les residus valent "<<f[0]<< "et "<< f[1]<< "et" << f[2] <<std::endl;
//     tableDD(t, f, table, delta);
//     // std::cout << "delta11 = "<< delta[0] << " delta12 = "<< delta[n]<< " et delta13 ="<< delta[2*n] <<
//     // std::endl;
//     d4u=6*delta[2*n];
//     d3u =
//       2.0 * delta[n] + 2*delta[2*n]*(-k1 + k2); // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
//     d2u = delta[0] + delta[n]*(-k1) + delta[2*n]*k2*(-k1); 
//     // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
//     // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
//     // std::cout<< "le derivee quatrieme vaut " << d4u << std::endl;
//     solDC4->_d[k] = d2u * k2 / 2.0 - d3u / 6.0 * k2* k2 + d4u / 24. * k2*k2*k2;
//     // std::cout<< "la correction du DC4F vaut " << d2u * k2 / 2.0 - d3u / 6.0 * k2* k2 + d4u / 24. * k2*k2*k2  << std::endl;
//   }
//   // Init FESOL
//   int nDOF = metaNumber->getNbDOFs();
//   for(int i = 0; i < nDOF; ++i) {
//     sol->setSolAtDOF(i, solDC4->_sol[1][i]);
//     sol->setSolDotAtDOF(i, 0.);
//   }
//   sol->setC0(cn[0]);
//   sol->setCurrentTime(tn);
//   solDC4->_cn[0] = cn[0];
//   solDC4->_cn[1] = cn[1];
//   // std::cout<<"c_0 = "<<cn[0]<<" c_1 = "<< cn[1]<<std::endl;
//   // solDC3->_cn[2] = cn[2];
//   sol->initializeEssentialBC(mesh, metaNumber);
//   sol->initializeEssentialBC(mesh, metaNumber, solDC4);
// }

void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3)
{ 
  std::vector<double> t = solDC2F->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 3;
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = 0; i < orderDC; ++i) {
    tt[i] = t[i];
  }
  double tn = tt[0];
  double k1 = tt[0] - tt[1];
  double k2 = tt[1] - tt[2];
 
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
    // solDC3->_fResidual[0][k] = solDC2F->_fResidual[0][k];
    // solDC3->_fResidual[1][k] = solDC2F->_fResidual[1][k];
    // solDC3->_fResidual[2][k] = solDC2F->_fResidual[2][k];
    // printf(" %12s \t %24.18e \n", "f[0]", f[0]);
    // printf(" %12s \t %24.18e \n", "f[1]", f[1]);
    // printf(" %12s \t %24.18e \n", "f[2]", f[2]);
    // std::cout<<"les residus valent "<<f[0]<< "  et  "<< f[1]<< "  et  " << f[2] <<std::endl;
    tableDD(t, f, table, delta);
    // std::cout << "delta11 vaut "<< delta[0] << "et delta12 vaut quant a lui "<< delta[n] <<
    // std::endl;
    d3u =
      2.0 * delta[n]; // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
    d2u = delta[0] + delta[n] * (k1 - k2);
    d2u = delta[0] + delta[n] * k1;
    // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
    // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
    solDC3->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1;
    // std::cout<< "la correction du DC3F vaut " << d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 <<
    // std::endl;
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
  // std::cout << "les coeffs sont "<< cn[0] <<" et "<< cn[1]<<std::endl;
  // solDC3->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC3);
}

void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4)
{ 
  std::vector<double> t = solDC3F->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 4;
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
    for(int i = 0; i < NbCoeffBDF; ++i) {
    tt[i] = t[i];
  }
  double tn = t[0];
  double k1 = t[0] - t[1];
  double k2 = t[1] - t[2];
  // std::cout<<"k1 = "<<k1<<"k2 = "<<k2<< std::endl;

 
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  double d2u;
  double d4u;
  int nInc = metaNumber->getNbUnknowns();
 for(int k = 0; k < nInc; ++k) {
    f[0] = solDC3F->_fResidual[0][k];
    f[1] = solDC3F->_fResidual[1][k];
    f[2] = solDC3F->_fResidual[2][k];
    f[3] = solDC3F->_fResidual[3][k];
    // printf(" %12s \t %24.18e \n", "f[0]", f[0]);
    // printf(" %12s \t %24.18e \n", "f[1]", f[1]);
    // printf(" %12s \t %24.18e \n", "f[2]", f[2]);
    // printf(" %12s \t %24.18e \n", "f[3]", f[3]);
    tableDD(t, f, table, delta);
    d4u=6*delta[2*n];
    d3u =
      2.0 * delta[n] + 2*delta[2*n]*(2*k1 + k2); // Indexé par delta[i][j] = delta[n*j+i] : delta(1,2) = delta[0][1] = delta[n]
    d2u = delta[0] + delta[n]*(k1) + delta[2*n]*(k1*(k1+k2)); 
    solDC4->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 + d4u / 24. * k1*k1*k1;
    // std::cout<< "la correction du DC4F vaut " << d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 + d4u / 24. * k1*k1*k1 << std::endl;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC4->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solDC4->_cn[0] = cn[0];
  solDC4->_cn[1] = cn[1];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC4);
}

void initializeDC4(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionBDF2 *solBDF2,
                   feSolutionDCF *solDC3, feSolutionDCF *solDC4)
{
  std::vector<double> t = solBDF2->getTime();
  int orderBDF = 2, NbCoeffBDF = orderBDF + 1;
  int orderDC = 4;
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
  double d3u, d4u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC3->_fResidual[0][k];
    f[1] = solDC3->_fResidual[1][k];
    f[2] = solDC3->_fResidual[2][k];
    f[3] = solDC3->_fResidual[3][k];
    tableDD(t, f, table, delta);
    // std::cout << "delta12 = " << delta[n] << " delta13 = " << delta[2*n] << std::endl;
    d3u = 2.0 * delta[n] + 2 * delta[2 * n] * (2 * k1 + k2);
    d4u = 6.0 * delta[2 * n];
    // std::cout << "le derivee trosieme vaut " << d3u << std::endl;
    // std::cout << "le derivee qutrieme vaut " << d4u << std::endl;
    solDC4->_d[k] = d3u / 6.0 * k1 * (k1 + k2) - d4u / 24.0 * k1 * (2 * k1 * k1 + 3 * k1 * k2 + k2 * k2);
    // std::cout << "la correction du DC4 vaut " << d3u / 6.0 * k1 * (k1 + k2) - d4u / 24.0 * k1 * (2 * k1 * k1 + 3 * k1 * k2 + k2 * k2) << std::endl;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC4->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);

  solDC4->_cn[0] = cn[0];
  solDC4->_cn[1] = cn[1];
  solDC4->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC4);
}

void initializeDC5(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionBDF2 *solBDF2,
                   feSolutionDCF *solDC3, feSolutionDCF *solDC4, feSolutionDCF *solDC5)
{
  std::vector<double> t = solBDF2->getTime();
  int orderBDF = 2, NbCoeffBDF = orderBDF + 1; //we need 3 points for the BDF2 
  int orderDC = 5; 
  t.resize(orderDC);
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = 0; i < NbCoeffBDF; ++i) {
    tt[i] = t[i];
  }
  double tn = t[0];
  double k1 = t[0] - t[1];
  double k2 = t[1] - t[2];
  double k3 = t[2] - t[3]; //in our situation, k1=k3
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);

  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  // std::cout<<"===n==="<<n<<std::endl;
  double d3u, d4u, d5u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC4->_fResidual[0][k];
    f[1] = solDC4->_fResidual[1][k];
    f[2] = solDC4->_fResidual[2][k];
    f[3] = solDC4->_fResidual[3][k];
    f[4] = solDC4->_fResidual[4][k];
    // printf(" %12s \t %24.18e \n", "f[0]", f[0]);
    // printf(" %12s \t %24.18e \n", "f[1]", f[1]);
    // printf(" %12s \t %24.18e \n", "f[2]", f[2]);
    // printf(" %12s \t %24.18e \n", "f[3]", f[3]);
    // printf(" %12s \t %24.18e \n", "f[4]", f[4]);
    // std::cout<<"les residus valent "<<f[0]<< " ,  "<< f[1]<<" ,  "<< f[2]<<" ,  "<< f[3]<< " et "
    // << f[4] <<std::endl;
    tableDD(t, f, table, delta);
    // for(int j=0 ; j<delta.size(); j++) {
    //   std::cout<<"delta"<<j;
    //   printf(" %24.18e \n", delta[j]);
    // }
    // std::cout << "delta12 = " << delta[n] << " delta13 = " << delta[2*n] << " delta14 = " <<
    // delta[3*n]<< std::endl;
    d3u = 2.0 * delta[n] + 2.0 * delta[2 * n] * (2.0 * k1 + k2) +
          2.0 * delta[3 * n] * (3.0 * k1 * k1 + 4.0 * k1 * k2 + 2.0 * k3 * k1 + k2 * k2 + k3 * k2);
    d4u = 6.0 * delta[2 * n] + 6.0 * delta[3 * n] * (3.0 * k1 + 2.0 * k2 + k3);
    d5u = 24.0 * delta[3 * n];
    // std::cout << "le derivee trosieme vaut " << d3u << std::endl;
    // std::cout << "le derivee qutrieme vaut " << d4u << std::endl;
    // std::cout << "le derivee cinqueme vaut " << d5u << std::endl;
    // std::cout<<" pas de temps k1 = "<<k1<<" k2 = "<<k1<< " k3 = "<<k3 << std::endl;
    solDC5->_d[k] = d3u / 6.0 * k1 * (k1 + k2) -
                    d4u / 24.0 * k1 * (2 * k1 * k1 + 3 * k1 * k2 + k2 * k2) +
                    d5u / 120. * (-pow(k1, 4.0) - pow(k1, 5.0) / k2 + (k1 / k2) * pow(k1 + k2, 4.));
    // std::cout << "la correction du DC5 vaut " << d3u / 6.0 * k1 * (k1 + k2) - d4u/ 24.0 * k1 *
    // (2*k1*k1 + 3*k1*k2 + k2*k2) + d5u/120. * (-pow(k1,4.0)-pow(k1,5.0)/k2 +
    // (k1/k2)*pow(k1+k2,4.)) << std::endl;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC5->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);

  solDC5->_cn[0] = cn[0];
  solDC5->_cn[1] = cn[1];
  solDC5->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC5);
}
///////////////////////////////////////////////////////////////////////////
// Set the correction for the first value of DC4F u^1, u^2, u^3 of order 4 .
void initializeDC4F_begining(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4 , std::string time )
{
  std::vector<double> t = solDC3F->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 4;
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  double t1 = t[2];
  double t2 = t[1];
  double k1 = t[0] - t[1];
  double k2 = t[1] - t[2];
  double k3 = t[2] - t[3];
  double k4 = t[3] - t[4];

  // Coeffs BDF
  tableToCoeffBDF(tt, cn);
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d3u;
  double d2u;
  double d4u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC3F->_fResidual[0][k];
    f[1] = solDC3F->_fResidual[1][k];
    f[2] = solDC3F->_fResidual[2][k];
    f[3] = solDC3F->_fResidual[3][k];
    // std::cout<<"les residus valent "<<f[0]<< "et "<< f[1]<< "et" << f[2] <<std::endl;
    tableDD(t, f, table, delta);
    // std::cout << "delta11 = "<< delta[0] << " delta12 = "<< delta[n]<< " et delta13 ="<< delta[2*n] <<
    // std::endl;
    if(time=="T1"){
      int shift_T = 2;
      for(int i = shift_T; i < NbCoeffBDF +shift_T; ++i) {
        tt[i-shift_T] = t[i];
        // std::cout<<" tt pour t1"<<tt[i-shift_T]<<std::endl;
      }
      // Coeffs BDF
      tableToCoeffBDF(tt, cn);
      d4u=6*delta[2*n];
      d3u =
      2.0 * delta[n] + 2*delta[2*n]*(-k1 - 2*k2); 
      d2u = delta[0] + delta[n]*(-k1 - 2*k2) + delta[2*n]*k2*(k1+k2);
      // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
      // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
      // std::cout<< "le derivee quatrieme vaut " << d4u << std::endl;
      solDC4->_d[k] = d2u * k3 / 2.0 - d3u / 6.0 * k3 * k3 + d4u / 24. * k3*k3*k3;
      // std::cout<< "la correction du DC4F vaut " << d2u * k3 / 2.0 - d3u / 6.0 * k3 * k3 + d4u / 24. * k3*k3*k3 << std::endl;
      sol->setCurrentTime(t1); 
    }
    if(time=="T2"){
      int shift_T = 1;
      for(int i = shift_T; i < NbCoeffBDF +shift_T; ++i) {
        tt[i-shift_T] = t[i];
        // std::cout<<" tt pour t2 "<<tt[i-shift_T]<<std::endl;
      }
      // Coeffs BDF
      tableToCoeffBDF(tt, cn);
      d4u=6*delta[2*n];
      d3u =
        2.0 * delta[n] + 2*delta[2*n]*(-k1 + k2); 
      d2u = delta[0] + delta[n]*(-k1) + delta[2*n]*k2*(-k1); 
      // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
      // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
      // std::cout<< "le derivee quatrieme vaut " << d4u << std::endl;
      solDC4->_d[k] = d2u * k2 / 2.0 - d3u / 6.0 * k2 * k2 + d4u / 24. * k2*k2*k2;
      // std::cout<< "la correction du DC4F vaut " << d2u * k2 / 2.0 - d3u / 6.0 * k2 * k2 + d4u / 24. * k2*k2*k2 << std::endl;
      sol->setCurrentTime(t2); 
    }
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC4->_sol[1][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  solDC4->_cn[0] = cn[0];
  solDC4->_cn[1] = cn[1];
  // std::cout<<"c_0 = "<<cn[0]<<" c_1 = "<< cn[1]<<std::endl;
  // solDC3->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC4);
}
/////////////////////////////////////////////////////////////////////////////////////
// Set the correction for the first value of DC5F u^1, u^2, u^3, u^4 of order 5 .
void initializeDC5F_begining(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionDC2F *solDC3F, 
                             feSolutionDC2F *solDC4F, feSolutionDC2F *solDC5, std::string time )
{ 
  std::vector<double> t = solDC3F->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 5;
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
    for(int i = 0; i < NbCoeffBDF; ++i) {
    tt[i] = t[i];
  }
  double t1 = t[3];
  double t2 = t[2];
  double t3 = t[1];
  double k1 = t[0] - t[1];
  double k2 = t[1] - t[2];
  double k3 = t[2] - t[3];
  double k4 = t[3] - t[4];
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d2u, d3u, d4u, d5u;
  int nInc = metaNumber->getNbUnknowns();
 for(int k = 0; k < nInc; ++k) {
    f[0] = solDC4F->_fResidual[0][k];
    f[1] = solDC4F->_fResidual[1][k];
    f[2] = solDC4F->_fResidual[2][k];
    f[3] = solDC4F->_fResidual[3][k];
    f[4] = solDC4F->_fResidual[4][k];
    // printf(" %12s \t %24.18e \n", "f[0]", f[0]);
    // printf(" %12s \t %24.18e \n", "f[1]", f[1]);
    // printf(" %12s \t %24.18e \n", "f[2]", f[2]);
    // printf(" %12s \t %24.18e \n", "f[3]", f[3]);
    // printf(" %12s \t %24.18e \n", "f[4]", f[4]);
    tableDD(t, f, table, delta);
    if(time=="T1"){
      int shift_T = 3 ;
      for(int i = shift_T; i < NbCoeffBDF +shift_T; ++i) {
        tt[i-shift_T] = t[i];
      }
      // Coeffs BDF
      tableToCoeffBDF(tt, cn);
      d5u =24*delta[3*n];
      d4u=6*delta[2*n] + 6*delta[3*n]*(-(k1+k2+k3)-(k2+k3)-k3);
      d3u =2.0 * delta[n] + 2*delta[2*n]*(-(k1+k2+k3) -(k2+k3)-k3) + 2*delta[3*n]*((k1+k2+k3)*(k2+k3) + (k1+k2+k3)*k3 + (k2+k3)*k3);
      d2u = delta[0] + delta[n]*(-(k1+k2+k3)-(k2+k3)) + delta[2*n]*(-(k1+k2+k3)*-(k2+k3) + (k1+k2+k3)*k3 + (k2+k3)*k3) - delta[3*n] *(k1+k2+k3)*(k2+k3)*(k3); 
      solDC5->_d[k] = d2u * k4 / 2.0 - d3u / 6.0 * k4 * k4 + d4u / 24. * k4*k4*k4  - d5u/120 *k4*k4*k4*k4;
      sol->setCurrentTime(t1); 
    }
    if(time=="T2"){
      int shift_T = 2;
      for(int i = shift_T; i < NbCoeffBDF +shift_T; ++i) {
        tt[i-shift_T] = t[i];
      }
      // Coeffs BDF
      tableToCoeffBDF(tt, cn);
      d5u =24*delta[3*n];
      d4u=6*delta[2*n] + 6*delta[3*n]*(-(k1+k2) -k2 +k3);
      d3u =2.0 * delta[n] + 2*delta[2*n]*(-(k1+k2) -k2) + 2*delta[3*n]*((k1+k2)*k2 - (k1+k2)*k3 - k2*k3);
      d2u = delta[0] + delta[n]*(-(k1+k2)-k2) + delta[2*n]*(-(k1+k2)*(-k2)) + delta[3*n] *(-(k1+k2))*(-k2)*(k3); 
      solDC5->_d[k] = d2u * k3 / 2.0 - d3u / 6.0 * k3 * k3 + d4u / 24. * k3*k3*k3  - d5u/120 *k3*k3*k3*k3;
      sol->setCurrentTime(t2); 
    }
    if(time=="T3"){
      int shift_T = 1;
      for(int i = shift_T; i < NbCoeffBDF + shift_T; ++i) {
        tt[i-shift_T] = t[i];
      }
      // Coeffs BDF
      tableToCoeffBDF(tt, cn);  
      d5u =24*delta[3*n];
      d4u=6*delta[2*n] + 6*delta[3*n]*(-k1 +2*k2 +k3 );
      d3u =2.0 * delta[n] + 2*delta[2*n]*(-k1+k2) + 2*delta[3*n]*(-k1*k2 -k1*(k2+k3) +k2*(k2+k3));
      d2u = delta[0] + delta[n]*(-k1) + delta[2*n]*(-k1*k2) + delta[3*n] *(-k1)*k2*(k2+k3); 
      solDC5->_d[k] = d2u * k2 / 2.0 - d3u / 6.0 * k2 * k2 + d4u / 24. * k2*k2*k2  - d5u/120 *k2*k2*k2*k2;
      sol->setCurrentTime(t3);
    }
    
    // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
    // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
    // std::cout<< "le derivee quatrieme vaut " << d4u << std::endl;
    // std::cout<< "le derivee cinquieme vaut " << d5u << std::endl;
    // solDC5->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 + d4u / 24. * k1*k1*k1  - d5u/120 *k1*k1*k1*k1;
    // std::cout<< "la correction du DC5F vaut " << solDC5->_d[k]<< std::endl;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC5->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  
  solDC5->_cn[0] = cn[0];
  solDC5->_cn[1] = cn[1];
  // std::cout << "les coeffs sont "<< cn[0] <<" et "<< cn[1]<<std::endl;
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC5);
}

void initializeDC5F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC3F,feSolutionDC2F *solDC4F, feSolutionDC2F *solDC5)
{ 
  std::vector<double> t = solDC3F->getTime();
  int orderBDF = 1, NbCoeffBDF = orderBDF + 1;
  int orderDC = 5;
  std::vector<double> cn(NbCoeffBDF);
  std::vector<double> tt(NbCoeffBDF, 0.);
  for(int i = 0; i < NbCoeffBDF+1; ++i) {
    tt[i] = t[i];
  }
  double tn = t[0];
  double k1 = t[0] - t[1];
  double k2 = t[1] - t[2];
  double k3 = t[2] - t[3];

 
  // Coeffs BDF
  tableToCoeffBDF(tt, cn);
  std::vector<double> f(orderDC, 0.0);
  size_t n = t.size();
  std::vector<double> table(n * (n + 1), 0.0);
  std::vector<double> delta(n * (n - 1), 0.0);
  double d2u, d3u , d4u, d5u;
  int nInc = metaNumber->getNbUnknowns();
  for(int k = 0; k < nInc; ++k) {
    f[0] = solDC4F->_fResidual[0][k];
    f[1] = solDC4F->_fResidual[1][k];
    f[2] = solDC4F->_fResidual[2][k];
    f[3] = solDC4F->_fResidual[3][k];
    f[4] = solDC4F->_fResidual[4][k];
    // printf(" %12s \t %24.18e \n", "f[0]", f[0]);
    // printf(" %12s \t %24.18e \n", "f[1]", f[1]);
    // printf(" %12s \t %24.18e \n", "f[2]", f[2]);
    // printf(" %12s \t %24.18e \n", "f[3]", f[3]);
    // printf(" %12s \t %24.18e \n", "f[4]", f[4]);
    tableDD(t, f, table, delta);
     // std::cout << "delta11 = "<< delta[0] << " delta12 = "<< delta[n]<< " et delta13 ="<< delta[2*n] <<"et delta14 ="<< delta[3*n] <<
    // std::endl;
    d5u =24*delta[3*n];
    d4u=6*delta[2*n] + 6*delta[3*n]*(k1+k1+k2+k1+k2+k3 );
    d3u =2.0 * delta[n] + 2*delta[2*n]*(k1+k1+k2) + 2*delta[3*n]*(k1*(k1+k2) + k1*(k1+k2+k3) + (k1+k2)*(k1+k2+k3));
    d2u = delta[0] + delta[n]*(k1) + delta[2*n]*(k1*(k1+k2)) + delta[3*n] *k1*(k1+k2)*(k1+k2+k3); 
    // std::cout<< "le derivee deuxieme vaut " << d2u << std::endl;
    // std::cout<< "le derivee troisieme vaut " << d3u << std::endl;
    // std::cout<< "le derivee quatrieme vaut " << d4u << std::endl;
    // std::cout<< "le derivee cinquieme vaut " << d5u << std::endl;
    solDC5->_d[k] = d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 + d4u / 24. * k1*k1*k1  - d5u/120 *k1*k1*k1*k1;
    // std::cout<< "la correction du DC5F vaut " << d2u * k1 / 2.0 - d3u / 6.0 * k1 * k1 + d4u / 24. * k1*k1*k1  - d5u/120 *k1*k1*k1*k1<< std::endl;
  }
  // Init FESOL
  int nDOF = metaNumber->getNbDOFs();
  for(int i = 0; i < nDOF; ++i) {
    sol->setSolAtDOF(i, solDC5->_sol[0][i]);
    sol->setSolDotAtDOF(i, 0.);
  }
  sol->setC0(cn[0]);
  sol->setCurrentTime(tn);
  solDC5->_cn[0] = cn[0];
  solDC5->_cn[1] = cn[1];
  // std::cout << "les coeffs sont "<< cn[0] <<" et "<< cn[1]<<std::endl;
  // solDC5->_cn[2] = cn[2];
  sol->initializeEssentialBC(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber, solDC5);
}