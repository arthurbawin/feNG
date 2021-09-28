#include <iostream>
#include <fstream>

#include "gmsh.h"
#include <set>

double mySizeField(int dim, int tag, double x, double y, double z) {
  return (x > 1.5 && x < 1.8) ? 0.02 : 2.0;
}

int main(int argc, char **argv) {
  gmsh::initialize();
  gmsh::model::add("myBeautifulModel");
  double lc = 0.2, xc = 3., yc = 2., r = 0.5;
  // Points
  gmsh::model::occ::addPoint(0, 0, 0, lc, 1);
  gmsh::model::occ::addPoint(10, 0, 0, lc, 2);
  gmsh::model::occ::addPoint(10, 4, 0, lc, 3);
  gmsh::model::occ::addPoint(0, 4, 0, lc, 4);
  // Courbes : un nouveau tag est renvoye si aucun n'est specifie
  int l1 = gmsh::model::occ::addLine(1, 2);
  int l2 = gmsh::model::occ::addLine(2, 3);
  int l3 = gmsh::model::occ::addLine(3, 4);
  int l4 = gmsh::model::occ::addLine(4, 1);
  int c1 = gmsh::model::occ::addCircle(xc, yc, 0, r);
  int boundary = gmsh::model::occ::addCurveLoop({l1, l2, l3, l4});
  int circle = gmsh::model::occ::addCurveLoop({c1});
  // Surface
  gmsh::model::occ::addPlaneSurface({boundary, circle}, 1);
  // Creation des structures Gmsh associees a la geometrie
  gmsh::model::occ::synchronize();
  // Entites physiques
  int entree = gmsh::model::addPhysicalGroup(1, {4});
  int bord = gmsh::model::addPhysicalGroup(1, {1, 3});
  gmsh::model::addPhysicalGroup(1, {2}, 3);
  gmsh::model::addPhysicalGroup(1, {5}, 4);
  int surface = gmsh::model::addPhysicalGroup(2, {1});
  gmsh::model::setPhysicalName(1, entree, "Entree");
  gmsh::model::setPhysicalName(1, bord, "Bord");
  gmsh::model::setPhysicalName(1, 3, "Sortie");
  gmsh::model::setPhysicalName(1, 4, "Cylindre");
  gmsh::model::setPhysicalName(2, surface, "Surface");
  // Champ Distance : specification des options via setNumber()/setNumbers()
  gmsh::model::mesh::field::add("Distance", 1); // tag > 0
  gmsh::model::mesh::field::setNumbers(1, "PointsList", {1});
  gmsh::model::mesh::field::setNumbers(1, "CurvesList", {3});
  gmsh::model::mesh::field::setNumber(1, "NumPointsPerCurve", 100);
  // Champ Threshold
  gmsh::model::mesh::field::add("Threshold", 2);
  gmsh::model::mesh::field::setNumber(2, "InField", 1);
  gmsh::model::mesh::field::setNumber(2, "SizeMin", lc / 30);
  gmsh::model::mesh::field::setNumber(2, "SizeMax", lc);
  gmsh::model::mesh::field::setNumber(2, "DistMin", 0.15);
  gmsh::model::mesh::field::setNumber(2, "DistMax", 1.5);
  // Champ MathEval
  gmsh::model::mesh::field::add("MathEval", 3);
  gmsh::model::mesh::field::setString(3, "F", "0.05 + Abs(Cos(2*pi*x/5))/5");
  // Champ Min
  gmsh::model::mesh::field::add("Min", 4);
  gmsh::model::mesh::field::setNumbers(4, "FieldsList", {2, 3});

  gmsh::model::mesh::field::setAsBackgroundMesh(4);

  // Carte de taille externe (callback)
  // double
  gmsh::model::mesh::setSizeCallback(mySizeField);

  // Ne pas tenir compte des frontieres pour evaluer la taille de maille
  gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

  gmsh::model::mesh::setAlgorithm(2, 1, 5);

  gmsh::model::mesh::generate(2);
  gmsh::write("foo.msh");

  // Launch the GUI to see the results:
  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
}