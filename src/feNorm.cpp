#include "feNorm.h"
#include <cmath>
#include "feQuadrature.h"
feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int degQuad, feFunction *solRef,
               feVectorFunction *VecSolRef, std::string TypeNorm)
  : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()),
    geoSpace(mesh->getGeometricSpace(cncGeoTag)), nElmGeo(mesh->getNbElm(cncGeoTag)),
    dim(intSpace->getDim()), nNodePerElem(intSpace->getNbNodePerElem()), _solRef(solRef),
    _degQuad(degQuad), _TypeNorm(TypeNorm)
{
  feQuadrature *rule =
    new feQuadrature(_degQuad, intSpace->getDim(), intSpace->getCncGeo()->getForme());
  w = rule->getWeights();
  x = rule->getXPoints();
  _nQuad = rule->getNQuad();
  y = rule->getYPoints();
  z = rule->getZPoints();
  _VecfeSpace = {_intSpace};
  NbFields = _VecfeSpace.size();
  if(VecSolRef != nullptr) {
    _VecSolRef = VecSolRef;
  }
  // Attention : ça changle la règle de quadrature pour les interpolants avant la résolution !
  // Il faudrait choisir idépendamment les deux sans conséquences
  // geoSpace->setQuadratureRule(rule);
  // intSpace->setQuadratureRule(rule);
  delete rule;
}

feNorm::feNorm(std::vector<feSpace *> &VecfeSpace, feMesh *mesh, int degQuad, feFunction *solRef,
               feVectorFunction *VecSolRef, std::string TypeNorm)
  : _VecfeSpace(VecfeSpace), cncGeoTag(VecfeSpace[0]->getCncGeoTag()),
    geoSpace(mesh->getGeometricSpace(cncGeoTag)), nElmGeo(mesh->getNbElm(cncGeoTag)),
    dim(VecfeSpace[0]->getDim()), nNodePerElem(VecfeSpace[0]->getNbNodePerElem()), _solRef(solRef),
    _degQuad(degQuad), _TypeNorm(TypeNorm)
{
  feQuadrature *rule =
    new feQuadrature(_degQuad, VecfeSpace[0]->getDim(), VecfeSpace[0]->getCncGeo()->getForme());
  w = rule->getWeights();
  x = rule->getXPoints();
  _nQuad = rule->getNQuad();
  y = rule->getYPoints();
  z = rule->getZPoints();
  NbFields = _VecfeSpace.size();
  if(VecSolRef != nullptr) {
    _VecSolRef = VecSolRef;
  }
  // Attention : ça changle la règle de quadrature pour les interpolants avant la résolution !
  // Il faudrait choisir idépendamment les deux sans conséquences
  VecfeSpace[0]->setQuadratureRule(rule);
  geoSpace->setQuadratureRule(rule);
  delete rule;
}

void feNorm::computeL2Norm0D(feSolution *sol)
{
  double normL2 = 0.0, solInt, solRef, t = sol->getCurrentTime();

  // solRef = _intSpace[0]->evalFun(t, x);
  solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
  // std::cout << "la solution de ref est " << solRef << " au temps " << t << std::endl;
  solInt = sol->getSolAtDOF(0);
  // std::cout << "la solution de numérique est " << solInt << " au temps " << t << std::endl;
  // printf("Solution at (%f,%f,%f) = %10.10f - ref = %10.10f\n", x[0], x[1], x[2], solInt,
  // solRef);
  normL2 += (abs(solInt) - abs(solRef)) * (abs(solInt) - abs(solRef));

  norm = sqrt(normL2);
}

void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double normL2 = 0.0, solInt, solRef, J, t = sol->getCurrentTime();
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);

    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);

    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      if(dim == 0) {
        J = 1.;
      } else if(dim == 1) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
      } else if(dim == 2) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      }

      solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
      // std::cout<< "========"<<std::endl;
      // std::cout<<"    SolInt   "<<solInt<<std::endl;
      // std::cout<<"    SolRef   "<<solRef<<std::endl;
      normL2 += (solInt - solRef) * (solInt - solRef) * J * w[k];
    }
  }
  norm = sqrt(normL2);
}

void feNorm::computeL2NormVec(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double normL2 = 0.0, solIntU, solIntV, J, t = sol->getCurrentTime();
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<int> adr1(_VecfeSpace[1]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> sol1(adr1.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    _VecfeSpace[1]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm, adr1);
    for(size_t i = 0; i < adr1.size(); ++i) {
      sol1[i] = solVec[adr1[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solIntU = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _VecfeSpace[1]->interpolateFieldAtQuadNode(sol1, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      if(dim == 0) {
        J = 1.;
      } else if(dim == 1) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
      } else if(dim == 2) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      } else {
      }

      // solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
      std::vector<double> solRef(6, 0);
      if(_VecSolRef != nullptr) {
        _VecSolRef->eval(t, x, solRef);
      }
      // std::cout<< "========"<<std::endl;
      // std::cout<<"    SolU   "<<solIntU<<std::endl;
      // std::cout<<"    SolRefU   "<<solRef[0]<<std::endl;
      normL2 += ((solIntU - solRef[0]) * (solIntU - solRef[0]) +
                 (solIntV - solRef[1]) * (solIntV - solRef[1])) *
                J * w[k];
    }
  }
  norm = sqrt(normL2);
}

void feNorm::computeH1NormVec(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double normH1 = 0.0, solIntU, solIntV, J, t = sol->getCurrentTime(), E11, E22, E12;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<int> adr1(_VecfeSpace[1]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> sol1(adr1.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    _VecfeSpace[1]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm, adr1);
    for(size_t i = 0; i < adr1.size(); ++i) {
      sol1[i] = solVec[adr1[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solIntU = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _VecfeSpace[1]->interpolateFieldAtQuadNode(sol1, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      double drdx = dxds[1] / J;
      double drdy = -dxds[0] / J;
      double dsdx = -dxdr[1] / J;
      double dsdy = dxdr[0] / J;

      double dudx, dudy, dvdx, dvdy;
      dudx = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdx +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdx;
      dudy = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdy +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdy;
      dvdx = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdx +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdx;
      dvdy = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdy +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdy;

      std::vector<double> solRef(6, 0);
      if(_VecSolRef != nullptr) {
        _VecSolRef->eval(t, x, solRef);
      }
      E11 = 2 * dudx - 2 * solRef[2]; // solRef[2] = dUref/dx
      E22 = 2 * dvdy - 2 * solRef[5]; // solRef[2] = dVref/dy
      E12 = (dudy + dvdx) - (solRef[3] + solRef[4]);

      normH1 +=
        ((solIntU - solRef[0]) * (solIntU - solRef[0]) +
         (solIntV - solRef[1]) * (solIntV - solRef[1]) + (E11 * E11 + 2 * E12 * E12 + E22 * E22)) *
        J * w[k];
    }
  }
  norm = sqrt(normH1);
}

void feNorm::computeSemiH1NormVec(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double SeminormH1 = 0.0, solIntU, solIntV, J, t = sol->getCurrentTime(), E11, E22, E12;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<int> adr1(_VecfeSpace[1]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> sol1(adr1.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    _VecfeSpace[1]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm, adr1);
    for(size_t i = 0; i < adr1.size(); ++i) {
      sol1[i] = solVec[adr1[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solIntU = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _VecfeSpace[1]->interpolateFieldAtQuadNode(sol1, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      double drdx = dxds[1] / J;
      double drdy = -dxds[0] / J;
      double dsdx = -dxdr[1] / J;
      double dsdy = dxdr[0] / J;

      double dudx, dudy, dvdx, dvdy;
      dudx = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdx +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdx;
      dudy = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdy +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdy;
      dvdx = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdx +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdx;
      dvdy = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdy +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdy;

      std::vector<double> solRef(6, 0);
      if(_VecSolRef != nullptr) {
        _VecSolRef->eval(t, x, solRef);
      }

      E11 = 2 * dudx - 2 * solRef[2]; // solRef[2] = dUref/dx
      E22 = 2 * dvdy - 2 * solRef[5]; // solRef[2] = dVref/dy
      E12 = (dudy + dvdx) - (solRef[3] + solRef[4]);

      SeminormH1 += (E11 * E11 + 2 * E12 * E12 + E22 * E22) * J * w[k];
    }
  }
  norm = sqrt(SeminormH1);
}

void feNorm::computeNormLambda(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double normL2 = 0.0, solInt, solRef, J, t = sol->getCurrentTime();
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      if(dim == 0) {
        J = 1.;
      } else if(dim == 1) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
      } else {
      }

      solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
      double Kd;
      std::vector<double> gradf(3, 0);
      if(_VecSolRef != nullptr) {
        _VecSolRef->eval(t, x, gradf);
        Kd = _VecSolRef->getParam()[0];
      }
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;

      if(_VecSolRef != nullptr) {
        if(_TypeNorm == "IntFlux") {
          normL2 += (solInt + Kd * (gradf[0] * Nx + gradf[1] * Ny)) * J * w[k];
        } else if(_TypeNorm == "NormL2") {
          normL2 += (solInt + Kd * (gradf[0] * Nx + gradf[1] * Ny)) *
                    (solInt + Kd * (gradf[0] * Nx + gradf[1] * Ny)) * J * w[k];
        }
      }
    }
  }
  if(_TypeNorm == "NormL2") {
    norm = sqrt(normL2);
  } else if(_TypeNorm == "IntFlux") {
    norm = normL2;
  }
}

void feNorm::computeNormLambdaNS(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double Integral = 0.0, solInt, p, J, t = sol->getCurrentTime(), T11, T12, T22;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      if(dim == 0) {
        J = 1.;
      } else if(dim == 1) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
      } else {
      }

      p = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0; // here is the pressure
      double mu;
      std::vector<double> gradf(3, 0);
      if(_VecSolRef != nullptr) {
        _VecSolRef->eval(t, x,
                         gradf); // gradf[0] = dudx gradf[1] = dudy gradf[2] = dvdx gradf[3] = dvdy
        mu = _VecSolRef->getParam()[0];
      }
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;

      T11 = 2 * gradf[0];
      T12 = gradf[1] + gradf[2];
      T22 = 2 * gradf[3];
      // std::cout<<"solInt   "<<solInt<<std::endl;
      // std::cout<<"J   "<<J<<std::endl;
      // Integral +=  (solInt + mu * (T11*Nx + T12*Ny)) * J * w[k];
      Integral += (-solInt - p * Nx + mu * (T11 * Nx + T12 * Ny)) * J * w[k];
    }
  }
  norm = Integral;
}

// void feNorm::computeIntFluxVeloc(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh) {
//   double Integral = 0.0, solIntU, solIntV, J;
//   int nElm = _VecfeSpace[0]->getNbElm();
//   std::cout<<"fespace ID" << _VecfeSpace[0]->getFieldID() << std::endl;

//   for(int iElm = 0; iElm < nElm; ++iElm) {
//     _VecfeSpace[0]->initializeAddressingVector(metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()),
//     iElm); _VecfeSpace[0]->initializeSolution(sol);
//     _VecfeSpace[1]->initializeAddressingVector(metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()),
//     iElm); _VecfeSpace[1]->initializeSolution(sol); geoCoord = mesh->getCoord(cncGeoTag, iElm);
//     for(int k = 0; k < _nQuad; ++k) {
//       std::cout<<"nb de nQuad "<<_nQuad<<"  pour un nb d'element  "<<nElm<<std::endl;
//       std::cout<<"nb fonction U "<<_VecfeSpace[0]->getNbFunctions()<<"  nb fonction V   "<<
//       _VecfeSpace[1]->getNbFunctions()<<std::endl; solIntU =
//       _VecfeSpace[0]->interpolateSolutionAtQuadNode(k); solIntV =
//       _VecfeSpace[1]->interpolateSolutionAtQuadNode(k);

//       std::vector<double> x(3, 0.0);
//       geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

//       std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//       geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//       J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);

//       geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//       double Nx = dxdr[1];
//       double Ny = -dxdr[0];
//       double N = sqrt(Nx * Nx + Ny * Ny);
//       Nx /= N;
//       Ny /= N;

//       Integral += (solIntU *Nx + solIntV*Ny ) * J * w[k];

//     }
//   }
//   norm = Integral;
// }
void feNorm::computeIntFluxNx(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double Integral = 0.0, solInt, J;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);

      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;

      Integral += solInt * Nx * J * w[k];
    }
  }
  norm = Integral;
}
void feNorm::computeIntFluxNy(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double Integral = 0.0, solInt, J;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);

      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;

      Integral += solInt * Ny * J * w[k];
    }
  }
  norm = Integral;
}

void feNorm::computeIntegralNum(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double Integral = 0.0, solInt, p, J, t = sol->getCurrentTime(), T11, T12, T22;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      if(dim == 0) {
        J = 1.;
      } else if(dim == 1) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        J = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
      } else {
      }
      Integral += solInt * J * w[k];
    }
  }
  norm = Integral;
}

void feNorm::computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double area = 0.0, solInt, J;
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      // std::vector<double> j(3, 0.0);
      // geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete
      // for dim > 1 J = j[0];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      area += J * w[k];
    }
  }
  norm = area;
}

void feNorm::computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh,
                             feFunction *fun)
{
  double area = 0.0, solInt, J, t = sol->getCurrentTime();
  int nElm = _VecfeSpace[0]->getNbElm();
  std::vector<int> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      // std::vector<double> j(3, 0.0);
      // geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete
      // for dim > 1 J = j[0];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      // area += fun->eval(t, x) * J * w[k];
      area += solInt * J * w[k];
    }
  }
  norm = area;
}

/* Estimates the L2 norm of the error taking an external solution as the reference solution. */
feStatus feNorm::computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                      feMesh *mesh, feMetaNumber *refMN,
                                                      feSolution *refSol, feMesh *refMesh,
                                                      const std::vector<feSpace *> refSpaces)
{
  double normL2 = 0.0, solInt, solRef;
  int nElm = _intSpace->getNbElm();
  std::vector<int> adr(_intSpace->getNbFunctions());
  std::vector<double> sol0(adr.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  std::vector<double> &solVecS = refSol->getSolutionReference();

  if(fabs(sol->getCurrentTime() - refSol->getCurrentTime()) > 1e-3) {
    feWarning("Solutions time differ by  more than 1e-3 : tRef "
              "= %f - tSol = %f.\n",
              refSol->getCurrentTime(), sol->getCurrentTime());
    feWarning("Maybe the solution file does not match the current solution.\n");
  }

  std::vector<double> x(3, 0.0);
  std::vector<double> r(3, 0.0);
  std::vector<double> J = mesh->getCncGeoByTag(cncGeoTag)->getJacobians();

  // Identify the feSpace from among the set of spaces from the reference solution (we should use a
  // map) Aussi : la numérotation doit être recréée pour assigner correctement les ddl aux espaces.
  // Donc soit elle est la même que pour le problème courant, ce qui suppose que les fespaces sont
  // les memes, ce qui a du sens, soit il faut redéfinir un nouveau set de feSpace à partir d'une
  // autre simulation, ce qui semble hasardeux.
  bool matchingSpace = false;
  for(feSpace *fS : refSpaces) {
    std::vector<int> adrS(fS->getNbFunctions());
    std::vector<double> solS(adrS.size());
    if(fS->getFieldID() == _intSpace->getFieldID() &&
       fS->getCncGeoID() == _intSpace->getCncGeoID()) {
      matchingSpace = true;

      feNumber *n = refMN->getNumbering(fS->getFieldID());

      for(int iElm = 0; iElm < nElm; ++iElm) {
        _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()),
                                              iElm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          sol0[i] = solVec[adr[i]];
        }

        geoCoord = mesh->getCoord(cncGeoTag, iElm);
        for(int k = 0; k < _nQuad; ++k) {
          geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
          // Current solution
          solInt = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
          // Reference solution evaluated at quad node
          int elm = -1;
          bool isFound = static_cast<feMesh2DP1 *>(refMesh)->locateVertex(x, elm, r);
          if(!isFound) {
            feWarning("Point (%f, %f, %f) was not found in the mesh.\n", x[0], x[1], x[2]);
            solRef = solInt; // Points outside the mesh do not contribute
          } else {
            fS->initializeAddressingVector(n, elm, adrS);
            for(size_t i = 0; i < adr.size(); ++i) {
              solS[i] = solVecS[adrS[i]];
            }
            solRef = fS->interpolateField(solS, r.data());
            normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * w[k];
          }

          // normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * w[k];
        }
      }

      norm = sqrt(normL2);

      break;
    }
  }

  if(!matchingSpace)
    return feErrorMsg(FE_STATUS_ERROR, "No finite element space in the reference set "
                                       "matches the target space on this connectivity.\n");

  return FE_STATUS_OK;
}
