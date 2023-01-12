#ifndef _FEGRAPHICS_
#define _FEGRAPHICS_

/*! \brief Basic graphic features to visualize 1D solutions
 *         in an OpenGL window.
 *
 *  Detailed description starts here.
 */
#include "feAPI.h"
#include "glfem.h"

class feBasicViewer
{
protected:
#if defined(HAVE_GLFW)
  GLFWwindow* _window;
#endif

  int _windowWidth;
  int _windowHeight;

  int _nElm;
  int _nInteriorPlotNodes;
  int _nPlotNodes;

  bool _drawAxes;
  double _xLim[2];
  double _yLim[2];

  double _startingTime;
  double _endTime;
  double _appTime;
  double _simulationTime;

  double _framesPerIteration;

  double *_x;
  double *_u;
  double *_uh;

public:
  feBasicViewer(std::string windowTitle, int nElm, int nInteriorPlotNodes);
  ~feBasicViewer();

  void drawMessage(std::string &message);
  void setAxesLimits(double *xLim, double *yLim);
  void drawAxes();
  void draw1DCurve(feMesh &mesh,
                   feMetaNumber &numbering,
                   feSolution &solution,
                   feSpace *space,
                   feFunction *analyticSolution = nullptr);
  void reshapeWindowBox(double xMin, double xMax, double yMin, double yMax);
  void reshapeWindowBox(double scaleFactor, feMesh &mesh, feSolution &solution);
  void windowUpdate();
  void getWindowBox();
  bool windowShouldClose();
};

// Draw the 1D curve discretized by continuous or discontinuous elements stored in solution on the active OpenGL window
// void draw1DCurveDG(feMesh &mesh, feMetaNumber &numbering, feSolution &solution, feSpace *space, feFunction *analyticSolution = nullptr);

#endif 