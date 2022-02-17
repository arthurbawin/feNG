#include "feAPI.h"   //bibliotheque avec tous les references des fonctions

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
  feMesh2DP1 mesh(meshFile);    //P1 car order=1 ? //on peut inclure d'autres pararametre (exemple bool curved), si pas inclut prend valeurs par defaut definies dans .h                                                                                          //que contient mesh ?    
  
  std::vector<feCncGeo *> CncGeo=mesh.getCncGeo(); 
  
  feCncGeo *ConnecElmInt=CncGeo[1]; //CncGeoElmInt=connect des elements

  int nElm=ConnecElmInt->getNbElm();  //dejà dans feMesh (_nInteriorElm)
  int nNode=ConnecElmInt->getNbNodes();  //dejà dans feMesh (_nNod)
  int nNodePerElm =ConnecElmInt->getNbNodePerElem();
  
  std::vector<int> connecNodes =ConnecElmInt->getNodeConnectivityRef();     
  std::vector<int> nbOccurence (nNode);                                        

  
  for (int j=0;j<nNode;++j){     //checked
      nbOccurence[j]=std::count(connecNodes.begin(),connecNodes.end(),j); 
  }



  std::vector<std::vector<int>> ListeElmParNode(nNode, std::vector<int> (0,0));   //utiliser fonction patch dans feRecovery
  for (int i=0;i<nElm;++i){
    for(int j=0;j<nNodePerElm;j++){
      int nds=connecNodes[i*nNodePerElm+j];    
      ListeElmParNode[nds].push_back(i);
    }
  }
 
  std::vector<int> colorElm1(nElm);  

  for (int i=0;i<nElm;++i){  //checked
    colorElm1[i]=-1;
  }

  int nbColor=0;
  bool noColor=false;

  while (noColor==false){

    for (int k=0;k<nElm;++k){
      if (colorElm1[k]==-1){

        for(int i=0;i<nNodePerElm;++i){
          int s=connecNodes[k*nNodePerElm+i];
          for (int j=0;j<nbOccurence[s];++j){
            int NumElm=ListeElmParNode[s][j];
            if(NumElm !=k){
              if(colorElm1[NumElm]<0){
                colorElm1[NumElm]=-2;
              }
            }
          }
        }

      }
    }

    nbColor=nbColor+1;
    noColor=true;
    for (int k=0;k<nElm;++k){
      if(colorElm1[k]==-1){
        colorElm1[k]=nbColor;
      }
      if(colorElm1[k]==-2){
        colorElm1[k]=-1;
        noColor=false;
      }
    }
  }

  for (int k=0;k<nElm;++k){
    feInfo("%d-%d",k,colorElm1[k]);
  }

  std::vector<int> colorElm=mesh.meshColoring(1);  

  for (int k=0;k<colorElm.size();++k){
    feInfo("%d-%d",k,colorElm[k]);
  }

  // Free the used memory
  // delete solver;
  // delete exporter;
  // delete uBord;
  // delete uDomaine;
  // delete funSol;
  // delete funSource;

  // petscFinalize();


  return 0;
}



//recuperer le nombre d'element de la conncetivté
//recuperer le nombre de noeuds de la conncectivité

//parcourir les elements de la connectivité
//recuperer le nombre de noeuds par element
//recuperer le ni





















