#include "feSysElm.h"
#include "feBilinearForm.h"

void feSysElm_0D_StiffSpring::createElementarySystem(std::vector<feSpace *> &space)
{
  _idX = 0;
  _idV = 1;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_StiffSpring::computeAe(feBilinearForm *form)
{
  double m = _par[0];
  double c = _par[1];
  double k = _par[2];
  form->_Ae[0][0] = form->_c0;
  form->_Ae[0][1] = -1.;
  form->_Ae[1][0] = k;
  form->_Ae[1][1] = form->_c0 * m + c;
}

void feSysElm_0D_StiffSpring::computeBe(feBilinearForm *form)
{
  double m = _par[0];
  double c = _par[1];
  double k = _par[2];
  double xDot = form->_intSpace[_idX]->interpolateFieldAtQuadNode(form->_solDot[_idX], 0);
  double vDot = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], 0);
  double x = form->_intSpace[_idX]->interpolateFieldAtQuadNode(form->_sol[_idX], 0);
  double v = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], 0);
  form->_Be[0] -= xDot - v;
  form->_Be[1] -= m * vDot + x * k + c * v;
}

void feSysElm_0D_Stiff2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idX = 0;
  _idY = 1;
  _idZ = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Stiff2::computeAe(feBilinearForm *form)
{
  double omega = _par;
  form->_Ae[0][0] = form->_c0 + 1.;
  form->_Ae[0][1] = 0.;
  form->_Ae[0][2] = 0.;

  form->_Ae[1][0] = 0.;
  form->_Ae[1][1] = form->_c0;
  form->_Ae[1][2] = -omega;

  form->_Ae[2][0] = 0.;
  form->_Ae[2][1] = +omega;
  form->_Ae[2][2] = form->_c0;
}

void feSysElm_0D_Stiff2::computeBe(feBilinearForm *form)
{
  double omega = _par;
  
  double xDot = form->_intSpace[_idX]->interpolateFieldAtQuadNode(form->_solDot[_idX], 0);
  double yDot = form->_intSpace[_idY]->interpolateFieldAtQuadNode(form->_solDot[_idY], 0);
  double zDot = form->_intSpace[_idZ]->interpolateFieldAtQuadNode(form->_solDot[_idZ], 0);
  double x = form->_intSpace[_idX]->interpolateFieldAtQuadNode(form->_sol[_idX], 0);
  double y = form->_intSpace[_idY]->interpolateFieldAtQuadNode(form->_sol[_idY], 0);
  double z = form->_intSpace[_idZ]->interpolateFieldAtQuadNode(form->_sol[_idZ], 0);
  form->_Be[0] -= xDot + x;
  form->_Be[1] -= yDot - omega * z;
  form->_Be[2] -= zDot + omega * y;
}

void feSysElm_0D_Stiff3::createElementarySystem(std::vector<feSpace *> &space)
{
  _idX = 0;
  _idY = 1;
  _idZ = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Stiff3::computeAe(feBilinearForm *form)
{
  double omega = _par;
  form->_Ae[0][0] = form->_c0 + 1.;
  form->_Ae[0][1] = -1.;
  form->_Ae[0][2] = -omega;

  form->_Ae[1][0] = 0.;
  form->_Ae[1][1] = form->_c0;
  form->_Ae[1][2] = -omega;

  form->_Ae[2][0] = 0.;
  form->_Ae[2][1] = +omega;
  form->_Ae[2][2] = form->_c0;
}

void feSysElm_0D_Stiff3::computeBe(feBilinearForm *form)
{
  double omega = _par;
  
  double xDot = form->_intSpace[_idX]->interpolateFieldAtQuadNode(form->_solDot[_idX], 0);
  double yDot = form->_intSpace[_idY]->interpolateFieldAtQuadNode(form->_solDot[_idY], 0);
  double zDot = form->_intSpace[_idZ]->interpolateFieldAtQuadNode(form->_solDot[_idZ], 0);
  double x = form->_intSpace[_idX]->interpolateFieldAtQuadNode(form->_sol[_idX], 0);
  double y = form->_intSpace[_idY]->interpolateFieldAtQuadNode(form->_sol[_idY], 0);
  double z = form->_intSpace[_idZ]->interpolateFieldAtQuadNode(form->_sol[_idZ], 0);
  form->_Be[0] -= xDot + x - y - omega * z;
  form->_Be[1] -= yDot - omega * z;
  form->_Be[2] -= zDot + omega * y;
}

void feSysElm_0D_weakBC::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC::computeAe(feBilinearForm *form)
{
  form->_Ae[0][1] = 1.;
  form->_Ae[1][0] = 1.;
}

void feSysElm_0D_weakBC::computeBe(feBilinearForm *form)
{
  std::vector<double> x(3);
  form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, 0, x);
  double gamma = _fct->eval(form->_tn, x);
  double lambda = form->_intSpace[_idL]->interpolateFieldAtQuadNode(form->_sol[_idL], 0);
  double u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], 0);
  form->_Be[0] -= lambda;
  form->_Be[1] -= (u - gamma);
}

void feSysElm_0D_weakBC_edo1::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC_edo1::computeAe(feBilinearForm *form)
{
  form->_Ae[0][1] = 1.;
  form->_Ae[1][0] = 1.;
  form->_Ae[1][2] = -1.;
  form->_Ae[2][2] = form->_c0;
}

void feSysElm_0D_weakBC_edo1::computeBe(feBilinearForm *form)
{
  std::vector<double> x(3);
  form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, 0, x);
  double gammaDot = _fct->eval(form->_tn, x);
  double lambda = form->_intSpace[_idL]->interpolateFieldAtQuadNode(form->_sol[_idL], 0);
  double u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], 0);
  double v = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], 0);
  double vDot = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], 0);
  // std::cout<< "v vaut "<< v << std::endl;
  // std::cout<<"gammaDot= "<<gammaDot<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
  // printf("%16.6e \n" , v);
  form->_Be[0] -= lambda;
  form->_Be[1] -= (u - v);
  form->_Be[2] -= (vDot - gammaDot);
}

void feSysElm_0D_weakBC_edo1_V2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC_edo1_V2::computeAe(feBilinearForm *form)
{
  form->_Ae[0][1] = 1.;
  form->_Ae[1][0] = 1.;
  form->_Ae[1][2] = -1.;
  form->_Ae[2][2] = form->_c0;
}

void feSysElm_0D_weakBC_edo1_V2::computeBe(feBilinearForm *form)
{
  std::vector<double> x(3);
  form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, 0, x);
  double gamma = _fct->eval(form->_tn, x);
  double lambda = form->_intSpace[_idL]->interpolateFieldAtQuadNode(form->_sol[_idL], 0);
  double u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], 0);
  double v = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], 0);
  double vDot = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], 0);
  // std::cout<<"gamma= "<<gamma<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
  printf("%16.6e \n", v);
  form->_Be[0] -= lambda;
  form->_Be[1] -= (u - v - gamma);
  form->_Be[2] -= (vDot);
}

void feSysElm_0D_weakBC_edo2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _idW = 3;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC_edo2::computeAe(feBilinearForm *form)
{
  form->_Ae[0][1] = 1.;
  form->_Ae[1][0] = 1.;
  form->_Ae[1][2] = -1.;
  form->_Ae[2][3] = -1.;
  form->_Ae[2][2] = form->_c0;
  form->_Ae[3][3] = form->_c0;
}

void feSysElm_0D_weakBC_edo2::computeBe(feBilinearForm *form)
{
  std::vector<double> x(3);
  form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, 0, x);
  double gammaDotDot = _fct->eval(form->_tn, x);
  double lambda = form->_intSpace[_idL]->interpolateFieldAtQuadNode(form->_sol[_idL], 0);
  double u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], 0);
  double v = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], 0);
  double vDot = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], 0);
  double w = form->_intSpace[_idW]->interpolateFieldAtQuadNode(form->_sol[_idW], 0);
  double wDot = form->_intSpace[_idW]->interpolateFieldAtQuadNode(form->_solDot[_idW], 0);
  form->_Be[0] -= lambda;
  form->_Be[1] -= (u - v);
  form->_Be[2] -= (vDot - w);
  form->_Be[3] -= (wDot - gammaDotDot);
}

void feSysElm_0D_Masse::createElementarySystem(std::vector<feSpace *> &space)
{
  // _idU = 0;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Masse::computeAe(feBilinearForm *form)
{
  // int                nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  // int        nFunctions = form->_intSpace[_idU]->getNbFunctions();

  for(int i = 0; i < form->_intSpace.size(); ++i) {
    form->_Ae[i][i] = rho * form->_c0;
  }
}

void feSysElm_0D_Masse::computeBe(feBilinearForm *form)
{
  // int                nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  // int        nFunctions = form->_intSpace[_idU]->getNbFunctions();

  double uDot;

  for(int i = 0; i < form->_intSpace.size(); ++i) {
    uDot = form->_intSpace[i]->interpolateFieldAtQuadNode(form->_solDot[i], 0);
    form->_Be[i] -= rho * uDot;
  }
}

void feSysElm_0D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _feU.resize(space[_idU]->getNbFunctions());
}

void feSysElm_0D_Source::computeBe(feBilinearForm *form)
{
  std::vector<double> x(3, 0.0);
  double S = _fct->eval(form->_tn, x);
  form->_Be[0] -= S;
}

void feSysElm_0D_Source_crossed::createElementarySystem(std::vector<feSpace *> &space)
{
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Source_crossed::computeBe(feBilinearForm *form)
{
  // int                nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  // int        nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> x(3, 0.0);
  std::vector<double> f(form->_intSpace.size(), 0);
  if(_fct != nullptr) _fct->eval(form->_tn, x, f);
  for(int i = 0; i < form->_intSpace.size(); i++) {
    form->_Be[i] = -f[i];
  }
}