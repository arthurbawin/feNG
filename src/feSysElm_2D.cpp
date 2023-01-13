#include "feSysElm.h"
#include "feBilinearForm.h"

void feSysElm_2D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _x.resize(3);
}

void feSysElm_2D_Source::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpaces[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();
  bool globalFunctions = form->_intSpaces[_idU]->useGlobalFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _x);

    double S = _fct->eval(form->_tn, _x);

    for(int i = 0; i < nFunctions; ++i) {
      if(globalFunctions) {
        _feU[i] = form->_intSpaces[_idU]->getGlobalFunctionAtQuadNode(form->_numElem, i, k);
      } else {
        _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      }
      form->_Be[i] -= _feU[i] * S * jac * w[k];
      // #pragma omp critical
      // {
      // printf("%d;%d;%d;",form->_numElem,i,k);std::cout<<form->_Be[i]<<std::endl;
      // }
    }
  }
}

void feSysElm_2D_Diffusion::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Diffusion::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = form->_intSpaces[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();
  // bool globalFunctions = form->_intSpaces[_idU]->useGlobalFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);

    // if(!globalFunctions) {
      // jac = _dxdr[0]*_dxds[1]-_dxdr[1]*_dxds[0];
      jac = J[nG * form->_numElem + k];
      double drdx = _dxds[1] / jac;
      double drdy = -_dxds[0] / jac;
      double dsdx = -_dxdr[1] / jac;
      double dsdy = _dxdr[0] / jac;

      for(int i = 0; i < nFunctions; ++i) {
        _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                    form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
        _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                    form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      }

      for(int i = 0; i < nFunctions; ++i) {
        for(int j = 0; j < nFunctions; ++j) {
          form->_Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * jac * w[k];
        }
      }

    // } else {
    //   // Using global interpolation functions
    //   for(int i = 0; i < nFunctions; ++i) {
    //     _feUdx[i] = form->_intSpaces[_idU]->getdGlobalFunctiondxAtQuadNode(form->_numElem, i, k);
    //     _feUdy[i] = form->_intSpaces[_idU]->getdGlobalFunctiondyAtQuadNode(form->_numElem, i, k);
    //   }

    //   for(int i = 0; i < nFunctions; ++i) {
    //     for(int j = 0; j < nFunctions; ++j) {
    //       form->_Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * jac * w[k];
    //     }
    //   }
    // }
  }

  // if(false) {
  //   // 1D quadrature rule to integrate over edges
  //   feQuadrature rule(20, 1, "");
  //   std::vector<double> x1D = rule.getXPoints();
  //   std::vector<double> w1D = rule.getWeights();
  //   std::vector<double> r1D(x1D.size());
  //   int n1D = rule.getNQuad();
  //   // 1D quad nodes from [-1,1] to [0,1]
  //   for(int i = 0; i < n1D; ++i) {
  //     r1D[i] = (x1D[i] + 1.0) / 2.0;
  //   }

  //   double drdt[3] = {1.0, -1.0, 0.0};
  //   double dsdt[3] = {0.0, 1.0, -1.0};

  //   // Il faut recalculer les fonctions globales aux noeuds de quadrature 1D
  //   // Pas tres efficace : elles ont ete precalculees mais seulement aux noeuds 2D
  //   // TODO : stocker les matrices de coeff
  //   std::vector<double> l(nFunctions, 0.0);
  //   std::vector<double> dldx(nFunctions, 0.0);
  //   std::vector<double> dldy(nFunctions, 0.0);

  //   std::string name = "elem" + std::to_string(form->_numElem) + ".pos";
  //   FILE *f = fopen(name.c_str(), "w");
  //   fprintf(f, "View \"elem\" {\n");

  //   // Integrale sur les bords (devra devenir une autre forme faible si ça fonctionne)
  //   for(int iEdge = 0; iEdge < 3; ++iEdge) {
  //     double longueur = 0.0;
  //     double longueur2 = 0.0;
  //     for(int k = 0; k < n1D; ++k) {
  //       std::vector<double> x(3, 0.0);
  //       std::vector<double> _dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  //       std::vector<double> _dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]

  //       double jac1D;

  //       switch(iEdge) {
  //         case 0: {
  //           double r[3] = {r1D[k], 0., 0.};
  //           // These interpolations use local functions because geometry is defined with Lagrange
  //           // local polynomials
  //           form->_geoSpace->interpolateVectorField(form->_geoCoord, r, x);
  //           form->_geoSpace->interpolateVectorField_rDerivative(form->_geoCoord, r, _dxdr);
  //           form->_geoSpace->interpolateVectorField_sDerivative(form->_geoCoord, r, _dxds);

  //           double tx = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //           double ty = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //           // printf("Tangential vector on edge 12 at phys point (%+-4.4f, %+-4.4f) - ref point
  //           // (%4.4f, %4.4f) =
  //           // (%+-4.4f, %+-4.4f)\n", x[0], x[1], r[0], r[1], tx, ty);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,tx,ty,0.0);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,ty,-tx,0.0);
  //           break;
  //         }
  //         case 1: {
  //           double r[3] = {r1D[k], 1.0 - r1D[k], 0.};
  //           form->_geoSpace->interpolateVectorField(form->_geoCoord, r, x);
  //           form->_geoSpace->interpolateVectorField_rDerivative(form->_geoCoord, r, _dxdr);
  //           form->_geoSpace->interpolateVectorField_sDerivative(form->_geoCoord, r, _dxds);
  //           double tx = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //           double ty = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //           // printf("Tangential vector on edge 23 at phys point (%+-4.4f, %+-4.4f) - ref point
  //           // (%4.4f, %4.4f) =
  //           // (%+-4.4f, %+-4.4f)\n",
  //           //   x[0], x[1], r[0], r[1], tx, ty);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,tx,ty,0.0);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,ty,-tx,0.0);
  //           break;
  //         }
  //         case 2: {
  //           double r[3] = {0.0, r1D[k], 0.};
  //           form->_geoSpace->interpolateVectorField(form->_geoCoord, r, x);
  //           form->_geoSpace->interpolateVectorField_rDerivative(form->_geoCoord, r, _dxdr);
  //           form->_geoSpace->interpolateVectorField_sDerivative(form->_geoCoord, r, _dxds);
  //           double tx = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //           double ty = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //           // printf("Tangential vector on edge 31 at phys point (%+-4.4f, %+-4.4f) - ref point
  //           // (%4.4f, %4.4f) =
  //           // (%+-4.4f, %+-4.4f)\n",
  //           //   x[0], x[1], r[0], r[1], tx, ty);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,tx,ty,0.0);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,ty,-tx,0.0);
  //           break;
  //         }
  //       }

  //       double dxdt = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //       double dydt = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //       // Le 1/2 vient de la transformation de [-1,1] à [0,1], qui est ensuite envoyé dans l'espace
  //       // physique
  //       jac1D = sqrt(dxdt * dxdt + dydt * dydt) / 2.0;

  //       longueur += jac1D * w1D[k];
  //       longueur2 += J * w1D[k];

  //       double nx = dydt;
  //       double ny = -dxdt;
  //       double N = 2.0 * jac1D;
  //       nx /= N;
  //       ny /= N;
  //       fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., -ny, nx, 0.0);
  //       fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., nx, ny, 0.0);

  //       form->_intSpaces[_idU]->Lphys(form->_numElem, x, l, dldx, dldy);

  //       for(int i = 0; i < nFunctions; ++i) {
  //         _feU[i] = l[i];
  //         _feUdx[i] = dldx[i];
  //         _feUdy[i] = dldy[i];
  //       }

  //       for(int i = 0; i < nFunctions; ++i) {
  //         for(int j = 0; j < nFunctions; ++j) {
  //           form->_Ae[i][j] -= (_feU[i] * (_feUdx[j] * nx + _feUdy[j] * ny)) * kD * jac1D * w1D[k];
  //         }
  //       }
  //     }
  //     // printf("Longueur = %4.4f et %4.4f\n", longueur, longueur2);
  //   }

  //   fprintf(f, "};");
  //   fclose(f);
  // }
}

void feSysElm_2D_Diffusion::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = form->_intSpaces[_idU]->getNbFunctions();
  bool globalFunctions = form->_intSpaces[_idU]->useGlobalFunctions();

  double J, dudx, dudy;
  for(int k = 0; k < nG; ++k) {

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    J = _dxdr[0]*_dxds[1]-_dxdr[1]*_dxds[0];

    if(!globalFunctions) {
      double drdx = _dxds[1] / J;
      double drdy = -_dxds[0] / J;
      double dsdx = -_dxdr[1] / J;
      double dsdy = _dxdr[0] / J;

      dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
             form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
      dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
             form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;

      for(int i = 0; i < nFunctions; ++i) {
        _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                    form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
        _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                    form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
        form->_Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * J * w[k];
      }
    } else {

      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
    //   dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_xDerivative(sol, form->_numElem, k);
    //   dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_yDerivative(sol, form->_numElem, k);

    //   for(int i = 0; i < nFunctions; ++i) {
    //     _feUdx[i] = form->_intSpaces[_idU]->getdGlobalFunctiondxAtQuadNode(form->_numElem, i, k);
    //     _feUdy[i] = form->_intSpaces[_idU]->getdGlobalFunctiondyAtQuadNode(form->_numElem, i, k);
    //     form->_Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * J * w[k];
      // }
    }
  }

  // if(false) {
  //   // 1D quadrature rule to integrate over edges
  //   feQuadrature rule(30, 1, "");
  //   std::vector<double> x1D = rule.getXPoints();
  //   std::vector<double> w1D = rule.getWeights();
  //   std::vector<double> r1D(x1D.size());
  //   int n1D = rule.getNQuad();
  //   // 1D quad nodes from [-1,1] to [0,1]
  //   for(int i = 0; i < n1D; ++i) {
  //     r1D[i] = (x1D[i] + 1.0) / 2.0;
  //   }

  //   double drdt[3] = {1.0, -1.0, 0.0};
  //   double dsdt[3] = {0.0, 1.0, -1.0};

  //   // std::string name = "elem" + std::to_string(form->_numElem) + ".pos";
  //   // FILE *f = fopen(name.c_str(), "w");
  //   // fprintf(f, "View \"elem\" {\n");

  //   // Il faut recalculer les fonctions globales aux noeuds de quadrature 1D
  //   // Pas tres efficace : elles ont ete precalculees mais seulement aux noeuds 2D
  //   // TODO : stocker les matrices de coeff
  //   std::vector<double> l(nFunctions, 0.0);
  //   std::vector<double> dldx(nFunctions, 0.0);
  //   std::vector<double> dldy(nFunctions, 0.0);

  //   // std::vector<double> xc(3, 0.0);
  //   // double rc[3] = {1. / 3., 1. / 3., 1. / 3.};
  //   // form->_geoSpace->interpolateVectorField(form->_geoCoord, rc, xc);

  //   // Integrale sur les bords (devra devenir une autre forme faible si ça fonctionne)
  //   for(int iEdge = 0; iEdge < 3; ++iEdge) {
  //     double longueur = 0.0;
  //     double longueur2 = 0.0;

  //     double intEdge = 0.0;

  //     for(int k = 0; k < n1D; ++k) {
  //       std::vector<double> x(3, 0.0);
  //       // std::vector<double> _dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  //       // std::vector<double> _dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]

  //       double jac1D;

  //       switch(iEdge) {
  //         case 0: {
  //           double r[3] = {r1D[k], 0., 0.};
  //           // These interpolations use local functions because geometry is defined with Lagrange
  //           // local polynomials
  //           form->_geoSpace->interpolateVectorField(form->_geoCoord, r, x);
  //           form->_geoSpace->interpolateVectorField_rDerivative(form->_geoCoord, r, _dxdr);
  //           form->_geoSpace->interpolateVectorField_sDerivative(form->_geoCoord, r, _dxds);
  //           break;
  //         }
  //         case 1: {
  //           double r[3] = {1.0 - r1D[k], r1D[k], 0.};
  //           form->_geoSpace->interpolateVectorField(form->_geoCoord, r, x);
  //           form->_geoSpace->interpolateVectorField_rDerivative(form->_geoCoord, r, _dxdr);
  //           form->_geoSpace->interpolateVectorField_sDerivative(form->_geoCoord, r, _dxds);
  //           break;
  //         }
  //         case 2: {
  //           double r[3] = {0.0, 1.0 - r1D[k], 0.};
  //           form->_geoSpace->interpolateVectorField(form->_geoCoord, r, x);
  //           form->_geoSpace->interpolateVectorField_rDerivative(form->_geoCoord, r, _dxdr);
  //           form->_geoSpace->interpolateVectorField_sDerivative(form->_geoCoord, r, _dxds);
  //           break;
  //         }
  //       }

  //       double dxdt = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //       double dydt = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //       // Le 1/2 vient de la transformation de [-1,1] à [0,1], qui est ensuite envoyé dans l'espace
  //       // physique
  //       jac1D = sqrt(dxdt * dxdt + dydt * dydt);

  //       longueur += jac1D * w1D[k];
  //       // longueur2 += J * w1D[k];

  //       double nx = dydt;
  //       double ny = -dxdt;
  //       double N = sqrt(dxdt * dxdt + dydt * dydt);
  //       nx /= N;
  //       ny /= N;
  //       // fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., -ny, nx, 0.0);
  //       // fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., nx, ny, 0.0);

  //       // x[0] -= xc[0];
  //       // x[1] -= xc[1];
  //       // x[2] -= xc[2];

  //       form->_intSpaces[_idU]->Lphys(form->_numElem, x, l, dldx, dldy);

  //       dudx = form->_intSpaces[_idU]->interpolateField_xDerivative(sol, form->_numElem, x);
  //       dudy = form->_intSpaces[_idU]->interpolateField_yDerivative(sol, form->_numElem, x);

  //       // double sumPhi = 0.0;
  //       for(int i = 0; i < nFunctions; ++i) {
  //         _feU[i] = l[i];
  //         form->_Be[i] += _feU[i] * (dudx * nx + dudy * ny) * kD * jac1D * w1D[k];
  //         intEdge += (_feU[i] * (dudx * nx + dudy * ny)) * kD * jac1D * w1D[k];
  //         // sumPhi += _feU[i];
  //         // sumPhi += dldx[i];
  //         // sumPhi += dldy[i];
  //         // std::cout<<_feU[i]<<" - "<<dudx<<" - "<<dudy<<" - "<<nx<<" - "<<ny<<" - "<<jac1D<<" -
  //         // "<<w1D[k]<<std::endl;
  //       }
  //       // printf("Sumphi = %+-10.10e\n", sumPhi);
  //     }
  //     // printf("Integrale sur l'arete %d de l'elm %d = %+-10.10e\n", iEdge, form->_numElem, intEdge);
  //     // printf("Longueur = %10.16e\n", longueur);
  //   }
  //   // fprintf(f, "};");
  //   // fclose(f);
  // }
}

void feSysElm_2D_Masse::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Masse::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctionsU = form->_intSpaces[_idU]->getNbFunctions();

  double jac, u, dudt;
  for(int k = 0; k < nG; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < nFunctionsU; ++i) _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[i][j] += _feU[i] * rho * form->_c0 * _feU[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_2D_Masse::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctionsU = form->_intSpaces[_idU]->getNbFunctions();

  double jac, dudt;
  for(int k = 0; k < nG; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * rho * dudt * jac * w[k];
    }
  }
}

void feSysElm_2D_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Advection::computeAe(feBilinearForm *form)
{
  // Calculee par differences finies
}

void feSysElm_2D_Advection::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpaces[_idU]->getNbFunctions();

  double Jac, u, dudt, dudx, dudy;
  for(int k = 0; k < nG; ++k) {

    // Evaluate the exterior velocity field
    std::vector<double> v(2, 0.0);
    std::vector<double> x(3, 0.0);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    _fct->eval(form->_tn, x, v);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    Jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / Jac;
    double drdy = -_dxds[0] / Jac;
    double dsdx = -_dxdr[1] / Jac;
    double dsdy = _dxdr[0] / Jac;

    // Compute SUPG parameter
    // double tau1 = 0.0, tau2 = 0.0, tau3 = 0.0, c = 0.0, cT = 0.0, kT = 0.0;
    // for(int i = 0; i < nFunctions; ++i) {
    //   dudt = form->_intSpaces[_idU]->interpolateSolutionDotAtQuadNode(k);
    //   dudx = form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
    //          form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    //   dudy = form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
    //          form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    //   _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
    //   _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
    //               form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
    //   _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
    //               form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    //   c  += J * w[k] * _feU[i] * (v[0] * dudx + v[1] * dudy);
    //   cT += J * w[k] * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * dudt;
    //   kT += J * w[k] * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (v[0] * dudx + v[1] * dudy);
    // }
    double dx01 = fabs(form->_geoCoord[3]-form->_geoCoord[0]);
    double dx02 = fabs(form->_geoCoord[6]-form->_geoCoord[0]);
    double dx12 = fabs(form->_geoCoord[6]-form->_geoCoord[3]);
    double dy01 = fabs(form->_geoCoord[4]-form->_geoCoord[1]);
    double dy02 = fabs(form->_geoCoord[7]-form->_geoCoord[1]);
    double dy12 = fabs(form->_geoCoord[7]-form->_geoCoord[4]);
    double h01 = sqrt(dx01*dx01 + dy01*dy01);
    double h02 = sqrt(dx02*dx02 + dy02*dy02);
    double h12 = sqrt(dx12*dx12 + dy12*dy12);
    double h = fmax(h12,fmax(h01,h02));
    double normeU = sqrt(v[0]*v[0]+v[1]*v[1]);
    // double Peh = normeU * h / 2.0;

    // double tau = 1.0/( sqrt( 4.0/(form->_dt*form->_dt) + 4.0*normeU*normeU/h/h + 9.0*16.0/(h*h*h*h)) );
    double tau = 1.0/( sqrt( 4.0/(form->_dt*form->_dt) + 4.0*normeU*normeU/h/h ) );
    // tau = 0.01;
    // double tau = h/2.0/sqrt(v[0]*v[0]+v[1]*v[1]) * (cosh(Peh)/sinh(Peh) - 1.0/Peh);
    // double tau = (cosh(Peh)/sinh(Peh) - 1.0/Peh);
    // double tau = 0.01;
    // std::cout<<tau<<std::endl;
    // double deltaSUPG = c/kT;
    // double Re = (v[0]*v[0]+v[1]*v[1]) * c/kT;
    // printf("Re = %+-4.4e - c = %+-4.4e - cT = %+-4.4e - kT = %+-4.4e - delta = %+-4.4e\n", Re, c, cT, kT, deltaSUPG);

    for(int i = 0; i < nFunctions; ++i) {
      u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
      dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU],k);
      dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
             form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
      dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
             form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;

      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      // form->_Be[i] -= (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * u * J * w[k];
      // SUPG ?
      // double delta = form->_geoCoord[3]-form->_geoCoord[0];
      // double deltaSUPG = 0.01;
      // form->_Be[i] -= ((v[0] * dudx + v[1] * dudy) * _feU[i] + delta * (v[0] * dudx + v[1] * dudy) * (v[0] * _feUdx[i] + v[1] * _feUdy[i])) * J * w[k];
      // form->_Be[i] -= Jac * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + deltaSUPG * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
      form->_Be[i] -= Jac * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + tau       * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
    }
  }
}

void feSysElm_2D_Stokes::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idV = 1;
  _idP = 2;
  _fieldsLayoutI = {_idU, _idV, _idP};
  _fieldsLayoutJ = {_idU, _idV, _idP};
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _feV.resize(space[_idV]->getNbFunctions());
  _feVdx.resize(space[_idV]->getNbFunctions());
  _feVdy.resize(space[_idV]->getNbFunctions());
  _feP.resize(space[_idP]->getNbFunctions());
  _dxdr.resize(3); // [dx/dr, dy/dr, dz/dr]
  _dxds.resize(3); // [dx/ds, dy/ds, dz/ds]
}

void feSysElm_2D_Stokes::computeBe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNbFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNbFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNbFunctions();

  double J, u, v, p, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3);
  std::vector<double> f(2); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    J = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / J;
    double drdy = -_dxds[0] / J;
    double dsdx = -_dxdr[1] / J;
    double dsdy = _dxdr[0] / J;

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    Sxx = -p + 2. * mu * dudx;
    Sxy = mu * (dudy + dvdx);
    Syx = mu * (dudy + dvdx);
    Syy = -p + 2. * mu * dvdy;

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(
        i, k); // Fonction test de u : uniquement pour les forces volumiques
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (Sxx * _feUdx[i] + Sxy * _feUdy[i] - f[0] * _feU[i]) * J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(
        i, k); // Fonction test de v : uniquement pour les forces volumiques
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (Syx * _feVdx[i] + Syy * _feVdy[i] - f[1] * _feV[i]) * J * w[k];
    }
    // Residu pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
      form->_Be[cnt++] += _feP[i] * (dudx + dvdy) * J * w[k];
    }
  }
}

void feSysElm_2D_Stokes::computeAe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNbFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNbFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNbFunctions();

  double jac, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3);
  std::vector<double> f(2, 0.); // Body forces

  for(int k = 0; k < nG; ++k) {
    // Body forces
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    int I = 0, J;

    // Set up interpolation functions
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
    }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feVdx[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feVdy[i] * mu * 2. * _feVdy[j] + _feVdx[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] += 0.0;
      }
      ++I;
    }
  }
}

void feSysElm_2D_NavierStokes::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idV = 1;
  _idP = 2;
  _fieldsLayoutI = {_idU, _idV, _idP};
  _fieldsLayoutJ = {_idU, _idV, _idP};
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _feV.resize(space[_idV]->getNbFunctions());
  _feVdx.resize(space[_idV]->getNbFunctions());
  _feVdy.resize(space[_idV]->getNbFunctions());
  _feP.resize(space[_idP]->getNbFunctions());
  _dxdr.resize(3); // [dx/dr, dy/dr, dz/dr]
  _dxds.resize(3); // [dx/ds, dy/ds, dz/ds]
}

void feSysElm_2D_NavierStokes::computeBe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNbFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNbFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNbFunctions();

  double J, u, v, p, dudt, dvdt, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3, 0.);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    J = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / J;
    double drdy = -_dxds[0] / J;
    double dsdx = -_dxdr[1] / J;
    double dsdy = _dxdr[0] / J;

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);

    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);
    dvdt = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], k);

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    Sxx = -p + 2. * mu * dudx;
    Sxy = mu * (dudy + dvdx);
    Syx = mu * (dudy + dvdx);
    Syy = -p + 2. * mu * dvdy;

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (_feU[i] * rho * (dudt + u * dudx + v * dudy) + Sxx * _feUdx[i] +
                    Sxy * _feUdy[i] - f[0] * _feU[i]) *
                   J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (_feV[i] * rho * (dvdt + u * dvdx + v * dvdy) + Syx * _feVdx[i] +
                    Syy * _feVdy[i] - f[1] * _feV[i]) *
                   J * w[k];
    }
    // Residu pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
      form->_Be[cnt++] += _feP[i] * (dudx + dvdy) * J * w[k];
    }
  }
}

void feSysElm_2D_NavierStokes::computeAe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNbFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNbFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNbFunctions();

  double jac, u, v, dudt, dvdt, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);

    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);
    dvdt = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], k);

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    int I = 0, J = 0;

    // Set up interpolation functions
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
    }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] +=
          jac * w[k] *
          (rho * _feU[i] * (form->_c0 * _feU[j] + u * _feUdx[j] + v * _feUdy[j] + _feU[j] * dudx) +
           _feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (rho * _feU[i] * (_feV[j] * dudy) + _feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (rho * _feV[i] * (_feU[j] * dvdx) + _feVdx[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] +=
          jac * w[k] *
          (rho * _feV[i] * (form->_c0 * _feV[j] + u * _feVdx[j] + v * _feVdy[j] + _feV[j] * dvdy) +
           _feVdy[i] * mu * 2. * _feVdy[j] + _feVdx[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] += 0.0;
      }
      ++I;
    }
  }
}