
#include "feAPI.h"

double fTest(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return pos[0] * pos[0];
}

int main(int argc, char **argv)
{
  feFunction *funTest = new feFunction(fTest, {});

  int dim, deg = 1, degreeQuadrature = 10;
  double xa = 0.;
  double xb = 1.;
  int nElm = 2;
  feMesh1DP1 mesh(xa, xb, nElm, "BXA", "BXB", "Domaine");

  feSpace *uDomaine;
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 1, LINE, LAGRANGE, deg, "U", "Domaine", degreeQuadrature, funTest));

  std::vector<feSpace*> spaces = {uDomaine};
  std::vector<feSpace*> essentialSpaces = {};

  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
  sol.initializeUnknowns(&mesh);

  int nInteriorPlotNodes = 40;
  feBasicViewer viewer("Testing the graphic window", mesh.getNbInteriorElems(), nInteriorPlotNodes);

  double xLim[2] = {0., 1.};
  double yLim[2] = {0., 1.};
  viewer.setAxesLimits(xLim, yLim);

  do{
    viewer.reshapeWindowBox(1.2, mesh, sol);
    viewer.draw1DCurve(mesh, numbering, sol, uDomaine, funTest);
    viewer.windowUpdate();
  } while(!viewer.windowShouldClose());

  return 0;
}