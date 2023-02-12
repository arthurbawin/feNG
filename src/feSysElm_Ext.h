#include "feSysElm.h"

// 
// Contains undocumented 0D and 1D weak forms
//

class feSysElm_0D_Source_crossed : public feSysElm
{
  protected:
    feVectorFunction *_fct;
    double _par; // Parametre
    int _idU;
    std::vector<double> _feU;

  public:
    feSysElm_0D_Source_crossed(double par, feVectorFunction *fct)
      : feSysElm(0, 1, SOURCE_CROSSED_0D, false), _fct(fct), _par(par){};
    ~feSysElm_0D_Source_crossed() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_1D_weakBC_edo1 : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idU;
    int _idL;
    int _idV;
    std::vector<double> _feU;
    std::vector<double> _feV;
    std::vector<double> _feL;

  public:
    feSysElm_1D_weakBC_edo1(double par, feFunction *fct) : feSysElm(1, 3, WEAKBC_EDO1_1D, true), _fct(fct), _par(par){};
    ~feSysElm_1D_weakBC_edo1() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_StiffSpring : public feSysElm
{
  protected:
    feFunction *_fct;
    std::vector<double> _par; // Parametre
    int _idX; // pos
    int _idV; // vit
    std::vector<double> _feX;
    std::vector<double> _feV;

  public:
    feSysElm_0D_StiffSpring(std::vector<double> par, feFunction *fct)
      : feSysElm(0, 2, STIFFSPRING_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_StiffSpring() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_Stiff2 : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idX;
    int _idY;
    int _idZ;
    std::vector<double> _feX;
    std::vector<double> _feY;
    std::vector<double> _feZ;

  public:
    feSysElm_0D_Stiff2(double par, feFunction *fct) : feSysElm(0, 3, STIFF2_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_Stiff2() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_Stiff3 : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idX;
    int _idY;
    int _idZ;
    std::vector<double> _feX;
    std::vector<double> _feY;
    std::vector<double> _feZ;

  public:
    feSysElm_0D_Stiff3(double par, feFunction *fct) : feSysElm(0, 3, STIFF3_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_Stiff3() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_weakBC : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idU;
    int _idL;
    std::vector<double> _feU;
    std::vector<double> _feL;

  public:
    feSysElm_0D_weakBC(double par, feFunction *fct) : feSysElm(0, 2, WEAKBC_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_weakBC() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_weakBC_edo1 : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idU;
    int _idL;
    int _idV;
    std::vector<double> _feU;
    std::vector<double> _feL;

  public:
    feSysElm_0D_weakBC_edo1(double par, feFunction *fct) : feSysElm(0, 3, WEAKBC_EDO1_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_weakBC_edo1() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_weakBC_edo1_V2 : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idU;
    int _idL;
    int _idV;
    std::vector<double> _feU;
    std::vector<double> _feL;

  public:
    feSysElm_0D_weakBC_edo1_V2(double par, feFunction *fct) : feSysElm(0, 3, WEAKBC_EDO1_V2_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_weakBC_edo1_V2() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_weakBC_edo2 : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par; // Parametre
    int _idU;
    int _idL;
    int _idV;
    int _idW;
    std::vector<double> _feU;
    std::vector<double> _feL;

  public:
    feSysElm_0D_weakBC_edo2(double par, feFunction *fct) : feSysElm(0, 4, WEAKBC_EDO2_0D, true), _fct(fct), _par(par){};
    ~feSysElm_0D_weakBC_edo2() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

  // 2D Diffusion using global functions:
  // Matrix:

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

// Residual :
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