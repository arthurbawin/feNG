#ifndef _FEGRAPHICS_
#define _FEGRAPHICS_

#include "feAPI.h"
#include "glfem.h"

class feBasicViewer
{
protected:
#if defined(HAVE_GLFW)
  GLFWwindow *_window;
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

  double *_rho;
  double *_mach;
  double *_p;

public:
  // Create an OpenGL window for 1D curves visualization. See testGraphic.cpp
  // for example use including resizing of the window, drawing and buffer update.
  //
  // Adapted from the material of MECA1120 Intro to FEM - UCLouvain by V. Legat
  //
  //               nElm: the number of line elements in the 1D mesh
  // nInteriorPlotNodes: the number of additional interior nodes used to plot
  //                     a smooth curve. 0 means the plot only connects the
  //                     extremities of each element.
  feBasicViewer(const std::string &windowTitle, int nElm, int nInteriorPlotNodes);
  ~feBasicViewer();

  void drawMessage(std::string &message);
  void setAxesLimits(double *xLim, double *yLim);
  void drawAxes();

  // Draw the 1D curve discretized by continuous or discontinuous
  // elements stored in solution on the active OpenGL window.
  void draw1DCurve(feMesh &mesh, feMetaNumber &numbering, feSolution &solution, feSpace *space,
                   feFunction *analyticSolution = nullptr, double yScaling = 1.0, int color = 0);

  void drawEulerNozzle(feMesh &mesh, feMetaNumber &numbering, feSolution &solution,
                       std::vector<feSpace *> &spaces, // rhoA, rhouA and eA
                       double yScaling[3], feFunction *analyticSolution = nullptr);

  void reshapeWindowBox(double xMin, double xMax, double yMin, double yMax);
  void reshapeWindowBox(double scaleFactorX, double scaleFactorY, feMesh &mesh,
                        feSolution &solution);
  void windowUpdate();
  void getWindowBox();
  bool windowShouldClose();
};

#endif