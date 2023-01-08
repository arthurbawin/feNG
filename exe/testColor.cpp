#include "feAPI.h"   //bibliotheque avec tous les references des fonctions

#if defined(HAVE_METIS)
  #include "metis.h"
#endif

//dans cet exemple, on connait déjà la forme de la solution (u(x,y)=x⁶), on calcule alors le terme source avec l'équation de la diffusion(f=k*u''=30*k*x⁴)

double fSol(const double t, const std::vector<double> x, const std::vector<double> par)  //solution exacte  
{
  return pow(x[0], 6);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par)  //terme source  //par: parametre
{
  double k = par[0];
  return k * 30. * pow(x[0], 4);
}

// int main(){
// #if defined(HAVE_METIS)

//     idx_t nVertices = 6;
//     idx_t nEdges    = 7;
//     idx_t nWeights  = 1;
//     idx_t nParts    = 2;

//     idx_t objval;

//     idx_t part[6];


//     // Indexes of starting points in adjacent array
//     idx_t xadj[7] = {0,2,5,7,9,12,14};

//     // Adjacent vertices in consecutive index order
//     idx_t adjncy[14] = {1,3,0,4,2,1,5,0,4,3,1,5,4,2};

//     // Weights of vertices
//     // if all weights are equal then can be set to NULL
//     idx_t vwgt[6];
    

//     // int ret = METIS_PartGraphRecursive(&nVertices,& nWeights, xadj, adjncy,
//     //               NULL, NULL, NULL, &nParts, NULL,
//     //               NULL, NULL, &objval, part);

//     int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj, adjncy,
//                NULL, NULL, NULL, &nParts, NULL,
//                NULL, NULL, &objval, part);

//     feInfo("Return value = %d", ret);
    
//     for(unsigned part_i = 0; part_i < nVertices; part_i++){
//   std::cout << part_i << " " << part[part_i] << std::endl;
//     }

// #endif
//     return 0;
// }

int main(int argc, char **argv)
{
  //petscInitialize(argc, argv);

  // Set the default parameters.
  const char *meshFile = "square.msh";
  int verbosity = 2;
  int order = 1;                   //par defaut 1 sinon lors de la compilation -o --order XX  ? 
  int degreeQuadrature = 10;       //par defaut 10 sinon lors de la compilation -dquad --degreeQuadrature XX  ?

  // Create an option parser and parse the command line arguments.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each "addOption" adds an optional field, although it can be made required by
  // setting the 5th argument of addOption to "true".
  feOptionsParser options(argc, argv);                                                        //options pour la compilation ou l'execution? 
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
  options.addOption(&order, "-o", "--order", "Finite element space order");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  feCheck(options.parse());

  // Set the global verbosity level :
  // - 0 : No information messages, only print warnings and errors
  // - 1 : Moderate information messages (default)
  // - 2 : All information messages
  setVerbose(verbosity);

  // Create a mesh structure from a Gmsh mesh file (version 2.2, 4.1+)
  feMesh2DP1 mesh(meshFile);   
  
  // std::vector <int> list=mesh.getList();  
  // for (int i=0;i<list.size();++i){
  //   feInfo("%d - %d",i,list[i]);
  // }
  // feInfo("%d",mesh.getNbInteriorElems());

  // std::vector <int> nbElmPerColor=mesh.getNbElmPerColor();  
  // std::vector<int> index=mesh.getIndexStartColorInList();

  // for (int i=0;i<nbElmPerColor.size();++i){
  //   feInfo("%d = %d - %d",i,nbElmPerColor[i],index[i]);
  // }
  // feInfo("%d",mesh.getNbColor());

  // int nbColor=mesh.getNbColor();
  // std::vector <int> colorElm=mesh.getColorElm();
  // for (int i=0;i<colorElm.size();++i){
  //   feInfo("%d = %d",i,colorElm[i]);
  // }
  // feInfo("%d",mesh.getNbInteriorElems());
  
  // int nbColor=mesh.getNbColor();
  // feInfo("%d",nbColor);
  

  // FILE *f = fopen("myBeautifulColors.pos", "w");
  // fprintf(f, "View \"coucou\" {\n");

  // feCncGeo *cnc = mesh.getCncGeoByTag(1);

  // int nElm = cnc->getNbElm();
  // for(int i = 0; i < nElm; ++i){
  //   fprintf(f, "ST(");
  //   for(int j = 0; j < 3; ++j){
  //     int iVert = cnc->getNodeConnectivity(i, j);
  //     Vertex *v = mesh.getVertex(iVert);
  //     fprintf(f, "%4.4f, %4.4f, %4.4f", v->x(), v->y(), v->z());
  //     if(j < 2)
  //       fprintf(f, ",");
  //   }
  //   fprintf(f, "){%d, %d, %d};\n", colorElm[i], colorElm[i], colorElm[i]);
  // }

  // fprintf(f, "};");
  // fclose(f);

  return 0;
}






















