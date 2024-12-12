
#include "feAPI.h"

double fSol(const double /* t */, const std::vector<double> &pos, const std::vector<double> & /* par */)
{
  double x = pos[0];
  double y = pos[1];
  return x*y;
}

int main(int /* argc */, char ** /* argv */)
{
  feFunction fun(fSol);

  // Define a mesh from 0 to 1 with 10 elements of equal size
  // feMesh1DP1 mesh(0., 1., 10, "xa", "xb", "dom");
  feMesh2DP1 mesh("../data/square1.msh");

  // Define a 1D quadrature rule to integrate polynomials of order up to 2
  feQuadrature rule(2, geometryType::TRI);
  std::vector<double> wQ = rule.getWeights();
  std::vector<double> xsiQ = rule.getXPoints();
  std::vector<double> etaQ = rule.getYPoints();

  double I = 0.;
  std::vector<double> xPhys(3, 0.);

  // Loop over 2D connectivities (mesh entities)
  for(feCncGeo *cnc : mesh.getCncGeo())
  {
    if(cnc->getDim() == 2){
      std::vector<double> elemCoord(3 * cnc->getNumVerticesPerElem());

      // Set the quadrature rule of the FE space used to interpolate
      // the geometry (here: P1 Lagrange). This triggers the computation
      // of the jacobians evaluated at quadrature points, although 
      // here they are constant an equal to 2/h.
      feCheck(cnc->setQuadratureRule(&rule));

      // Get the jacobians
      const std::vector<double> &jacobians = cnc->getJacobians();

      // Loop over elements
      for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm)
      {
        // Get the element's physical coordinates
        mesh.getCoord(cnc, iElm, elemCoord);

        for(int j = 0; j < cnc->getNumVerticesPerElem(); ++j) {
          // feInfo("Local coordinates of element %d on cnc %s = %f - %f - %f",
          //   iElm, cnc->getID().data(),
          //   elemCoord[3 * j + 0],
          //   elemCoord[3 * j + 1],
          //   elemCoord[3 * j + 2]);
        }

        double x0 = elemCoord[3 * 0 + 0];
        double y0 = elemCoord[3 * 0 + 1];
        // double z0 = elemCoord[3 * 0 + 2];

        double x1 = elemCoord[3 * 1 + 0];
        double y1 = elemCoord[3 * 1 + 1];
        // double z1 = elemCoord[3 * 1 + 2];

        double x2 = elemCoord[3 * 2 + 0];
        double y2 = elemCoord[3 * 2 + 1];
        // double z2 = elemCoord[3 * 2 + 2];

        // Integrate over the element and sum
        double Iloc = 0.;
        double nQuad = rule.getNumQuadPoints();
        for(int i = 0; i < nQuad; ++i){

          // Map quadrature nodes to physical coordinates using reference to physical element transformation
          xPhys[0] = x0 * (1. - xsiQ[i] - etaQ[i]) + x1 * xsiQ[i] + x2 * etaQ[i];
          xPhys[1] = y0 * (1. - xsiQ[i] - etaQ[i]) + y1 * xsiQ[i] + y2 * etaQ[i];

          Iloc += fun(0., xPhys) * jacobians[nQuad * iElm + i] * wQ[i];
        }
        I += Iloc;
      }
    }
  }
  
  // Alternate way using a FE solution:
  // Define a P2 Lagrange FE space on the domain to represent fSol = xy exactly
  feSpace *u;
  int order = 2;
  int degreeQuadrature = 2;
  feCheck(createFiniteElementSpace(u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &fun));

  std::vector<feSpace*> allFEspaces = {u};
  std::vector<feSpace*> essentialSpaces = {u};

  // Define a numbering, solution and initialize essential DOFs
  feMetaNumber numbering(&mesh, allFEspaces, essentialSpaces);
  feSolution sol(numbering.getNbDOFs(), allFEspaces, essentialSpaces);
  sol.initializeEssentialBC(&mesh);

  // Define a solution-based norm
  feNorm norm(INTEGRAL, {u}, &sol);
  
  // Alternate way: use the feNorm made for user-defined fields
  feNorm norm2(mesh.getCncGeoByName("Domaine"), &fun);

  feInfo("Integral of xy on [0,1] x [0,1] by hand                         =  %+-1.10e", I);
  feInfo("Integral of xy on [0,1] x [0,1] using a FE discretization       =  %+-1.10e", norm.compute());
  feInfo("Integral of xy on [0,1] x [0,1] using integral of user function =  %+-1.10e", norm2.compute(INTEGRAL_F));

  return 0;
}