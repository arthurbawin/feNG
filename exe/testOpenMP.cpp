#include "feAPI.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par)  //solution exacte  
{
  return pow(x[0], 6);
}

double fSource(const double t, const std::vector<double> pos, const std::vector<double> par)  //terme source  //par: parametre
{
  double k = par[0];
  return k * 30. * (pow(pos[0], 4) + pow(pos[1], 4));
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.0;
}

int main(int argc, char **argv) {

  // Set the default parameters.
  const char *meshFile = "square.msh";
  int verbosity = 2;
  int order = 2;                   //par defaut 1 sinon lors de la compilation -o --order XX  ? 
  int degreeQuadrature = 10;       //par defaut 10 sinon lors de la compilation -dquad --degreeQuadrature XX  ?

  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
  options.addOption(&order, "-o", "--order", "Finite element space order");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  feCheck(options.parse());

  setVerbose(verbosity);

  feMesh2DP1 mesh(meshFile);

  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  int dim;
  feSpace *uBord, *uDomaine, *vBord, *vDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, order, "U", "Bord",
                                   degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Domaine",
                                   degreeQuadrature, funSol));

  std::vector<feSpace *> spaces = {uBord, uDomaine};
  std::vector<feSpace *> essentialSpaces = {uBord};

  feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);
  feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

  // feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));        
  // feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));

  #pragma omp parallel
  {
    feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
    #pragma omp for
    for(int i = 0; i < mesh.getNbInteriorElems(); ++i){
      diffU.initialize(&metaNumber, &mesh, &sol, i);
      feInfo("From thread %d/%d : Lecture elm %d : [%3d %3d %3d %3d %3d %3d]", 
        omp_get_thread_num(),
        omp_get_num_threads(), i,
        diffU._adr[0],
        diffU._adr[1],
        diffU._adr[2],
        diffU._adr[3],
        diffU._adr[4],
        diffU._adr[5]);
    }
  }

  // {
    // int ID = omp_get_thread_num();
    // printf("hello(%d/%d)\n", ID, omp_get_num_threads());
    // printf("world(%d/%d)\n", ID, omp_get_num_threads());

    
  // }

  // int cnt = 0;
  // #pragma omp parallel for private(cnt)
  // for(int i = 0; i < 100; ++i){
  //   // printf("Printing %3d from thread %d\n", i, omp_get_thread_num());
  //   printf("Thread %d has printed %2d times\n", omp_get_thread_num(), cnt++);
  // }

  return 0;
}