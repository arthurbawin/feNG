#include "feNG.h"
#include "feMessage.h"
#if defined(HAVE_PETSC)
  #include "petscsys.h"
#endif
#if defined(HAVE_MPI)
  #include "mpi.h"
#endif
#if defined(HAVE_GMSH)
 #include "gmsh.h"
#endif

#include <stack>

using namespace std::chrono_literals;

// To allow nested calls to tic()/toc()
std::stack<std::chrono::steady_clock::time_point> timeStamps;

double tic(int mode)
{
  // static std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point begin;

  if(mode == 0){
    begin = std::chrono::steady_clock::now();
    timeStamps.push(begin);
    return 0.;
  }
  else {
    if(timeStamps.empty()) {
      feErrorMsg(FE_STATUS_ERROR, "toc() was called without calling tic() first");
      exit(-1);
    }
    auto end = std::chrono::steady_clock::now();
    auto lastBegin = timeStamps.top();
    timeStamps.pop();
    // std::chrono::duration<double, std::milli> d = end - begin;
    std::chrono::duration<double, std::milli> d = end - lastBegin;
    double timeInSeconds  = d / 1.0s;
    return timeInSeconds;
  }
}

double toc() { return tic(1); }

bool wasInitialized = false;
#if defined(HAVE_GMSH)
bool gmshWasInitialized = false;
extern int FE_VERBOSE;
#endif

void initialize(int argc, char **argv)
{
#if defined(HAVE_PETSC)
  if(!wasInitialized) {
    wasInitialized = true;
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, (char *)0, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    feErrorMsg(FE_STATUS_ERROR, "PETSc was already initialized\n");
    return;
  }
#elif defined(HAVE_MPI)
  if(!wasInitialized) {
    wasInitialized = true;
    MPI_Init(&argc, &argv);
  } else {
    feErrorMsg(FE_STATUS_ERROR, "MPI was already initialized\n");
    return;
  }
#else
  // Do nothing
  UNUSED(argc, argv);
#endif
}

void finalize()
{
#if defined(HAVE_PETSC)
  PetscFinalize();
#elif defined(HAVE_MPI)
  MPI_Finalize();
#else
  // Do nothing
#endif
}

#if defined(HAVE_GMSH)
void initializeGmsh()
{
  if(!gmshWasInitialized) {

    // Get max threads *before* initializing gmsh
    #if defined(HAVE_OMP)
    int maxNumThreads = omp_get_max_threads();
    #endif

    gmsh::initialize();

    #if defined(HAVE_OMP)
    gmsh::option::setNumber("General.NumThreads", maxNumThreads);
    #endif

    gmshWasInitialized = true;
  }

  if(FE_VERBOSE == VERBOSE_NONE) gmsh::option::setNumber("General.Verbosity", 2);
}

// Straight outta Gmsh's StringUtils.cpp
static std::vector<std::string> SplitFileName(const std::string &fileName)
{
  std::vector<std::string> s;
  s.resize(3);
  if(fileName.size()) {
    // returns [path, baseName, extension]
    int idot = (int)fileName.find_last_of('.');
    int islash = (int)fileName.find_last_of("/\\");
    if(idot == (int)std::string::npos) idot = -1;
    if(islash == (int)std::string::npos) islash = -1;
    if(idot > 0) s[2] = fileName.substr(idot);
    if(islash > 0) s[0] = fileName.substr(0, islash + 1);
    s[1] =
      fileName.substr(s[0].size(), fileName.size() - s[0].size() - s[2].size());
  }
  return s;
}

void openGmshModel(std::string &filename)
{
  // Extract the model name
  std::vector<std::string> splitName = SplitFileName(filename);
  std::string modelName = splitName[1];
  feInfo("Opening file name %s", filename.data());
  feInfo("Model        name %s", modelName.data());

  if(gmshWasInitialized) {
    std::vector<std::string> modelList;
    gmsh::model::list(modelList);

    feInfo("");
    feInfo("Checking in list of models:");
    for(const std::string &s : modelList) {
      feInfo("%s", s.data());
    }
    feInfo("");

    gmsh::fltk::run();

    for(const std::string &s : modelList) {
      if(s == modelName) {
        feInfo("Setting existing %s as current Gmsh model", modelName.data());
        gmsh::model::setCurrent(modelName);
        return;
      }
    }

    // Model was not already loaded, open it
    feInfo("Opening new model %s", filename.data());
    gmsh::open(filename);
  }
}

void finalizeGmsh()
{
  if(gmshWasInitialized) {
    gmsh::clear();
    gmsh::finalize();
    gmshWasInitialized = false;
  }
}
#endif