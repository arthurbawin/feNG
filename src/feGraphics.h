#ifndef _FEGRAPHICS_
#define _FEGRAPHICS_

/* Basic graphic features */

#include "feAPI.h"
#include "glfem.h"

class feBasicViewer
{
protected:
  GLFWwindow* _window;

  int _windowWidth;
  int _windowHeight;

  double _startingTime;
  double _endTime;
  double _appTime;
  double _simulationTime;

  double _framesPerIteration;

public:
  feBasicViewer(std::string windowTitle);

  void reshapeWindowBox(double xMin, double xMax, double yMin, double yMax);
  void reshapeWindowBox(double scaleFactor, feMesh &mesh, feSolution &solution);
  void windowUpdate();
  void getWindowBox();
  bool windowShouldClose();
};

// Draw the 1D curve discretized by continuous or discontinuous elements stored in solution on the active OpenGL window
void draw1DCurve(feMesh &mesh, feMetaNumber &numbering, feSolution &solution, feSpace *space, feFunction *analyticSolution = nullptr);
void draw1DCurveDG(feMesh &mesh, feMetaNumber &numbering, feSolution &solution, feSpace *space, feFunction *analyticSolution = nullptr);

#endif 