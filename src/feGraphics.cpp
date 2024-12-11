#include "feGraphics.h"

feBasicViewer::feBasicViewer(const std::string &windowTitle, int nElm, int nInteriorPlotNodes)
  : _nElm(nElm)
{
#if defined(HAVE_GLFW)
  _window = glfemInit(windowTitle.data());
  glfwMakeContextCurrent(_window);
  glfwGetFramebufferSize(_window, &_windowWidth, &_windowHeight);

  _drawAxes = false;

  _nInteriorPlotNodes = fmax(nInteriorPlotNodes, 0);
  // 2 boundary + nInterior per element
  _nPlotNodes = nElm * (_nInteriorPlotNodes + 2);

  // Allocate some vectors to plot
  _x = (double *)malloc(sizeof(double) * _nPlotNodes);
  _u = (double *)malloc(sizeof(double) * _nPlotNodes);
  _uh = (double *)malloc(sizeof(double) * _nPlotNodes);

  _rho = (double *)malloc(sizeof(double) * _nPlotNodes);
  _mach = (double *)malloc(sizeof(double) * _nPlotNodes);
  _p = (double *)malloc(sizeof(double) * _nPlotNodes);
#else
  UNUSED(windowTitle, nElm, nInteriorPlotNodes);
  feWarning("Cannot create a display window without GLFW support :/");
#endif
}

feBasicViewer::~feBasicViewer()
{
#if defined(HAVE_GLFW)
  free(_x);
  free(_u);
  free(_uh);
  free(_rho);
  free(_mach);
  free(_p);
#endif
}

void feBasicViewer::drawMessage(std::string &message) { glfemMessage(&message[0]); }

void feBasicViewer::setAxesLimits(double *xLim, double *yLim)
{
  _drawAxes = true;
  _xLim[0] = xLim[0];
  _xLim[1] = xLim[1];
  _yLim[0] = yLim[0];
  _yLim[1] = yLim[1];
}
void feBasicViewer::drawAxes()
{
#if defined(HAVE_GLFW)
  double zero[2] = {0., 0.};
  glColor3f(0.0, 0.0, 0.0);
  glfemDrawCurve(_xLim, zero, 2);
  glfemDrawCurve(zero, _yLim, 2);
#endif
}

void feBasicViewer::reshapeWindowBox(double xMin, double xMax, double yMin, double yMax)
{
#if defined(HAVE_GLFW)
  glfemReshapeWindowsBox(xMin, xMax, yMin, yMax, _windowWidth, _windowHeight);
#else
  UNUSED(xMin, xMax, yMin, yMax);
#endif
}

// Reshape the window according to the min/max positions in the 1D mesh
// and the min/max field values in the solution, with a scale factor.
void feBasicViewer::reshapeWindowBox(double scaleFactorX, double scaleFactorY, feMesh &mesh,
                                     feSolution &solution)
{
#if defined(HAVE_GLFW)
  getWindowBox();
  double xMin = DBL_MAX;
  double yMin = DBL_MAX;
  double xMax = -DBL_MAX;
  double yMax = -DBL_MAX;
  for(auto v : mesh.getVertices()) {
    xMin = fmin(xMin, v.x());
    xMax = fmax(xMax, v.x());
  }
  for(auto val : solution.getSolutionReference()) {
    yMin = fmin(yMin, val);
    yMax = fmax(yMax, val);
  }
  double xC = (xMin + xMax) / 2.;
  double yC = (yMin + yMax) / 2.;
  double dx = xMax - xMin;
  double dy = yMax - yMin;
  this->reshapeWindowBox(xC - dx / 2. * scaleFactorX, xC + dx / 2. * scaleFactorX,
                         yC - dy / 2. * scaleFactorY, yC + dy / 2. * scaleFactorY);
#else
  UNUSED(scaleFactorX, scaleFactorY, mesh, solution);
#endif
}

void feBasicViewer::windowUpdate()
{
#if defined(HAVE_GLFW)
  glfwSwapBuffers(_window);
  glfwPollEvents();
#endif
}

void feBasicViewer::getWindowBox()
{
#if defined(HAVE_GLFW)
  glfwGetFramebufferSize(_window, &_windowWidth, &_windowHeight);
#endif
}

bool feBasicViewer::windowShouldClose()
{
#if defined(HAVE_GLFW)
  return glfwGetKey(_window, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwWindowShouldClose(_window) == 1;
#else
  return true;
#endif
}

void feBasicViewer::draw1DCurve(feMesh &mesh, feMetaNumber &numbering, feSolution &solution,
                                feSpace *space, feFunction *analyticSolution, double yScaling,
                                int color)
{
#if defined(HAVE_GLFW)
  std::vector<double> &sol = solution.getSolutionReference();
  std::vector<double> elmSolution(space->getNumFunctions());
  std::vector<double> xloc;
  std::vector<feInt> adr(space->getNumFunctions());

  int elm, cnt = 0;
  double r[3], eps = 1e-4;
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    Vertex *v0 = mesh.getVertex(iElm);
    Vertex *v1 = mesh.getVertex(iElm + 1);
    double dx = fabs(v1->x() - v0->x()) / (_nInteriorPlotNodes + 1);

    for(int iPlot = 0; iPlot < _nInteriorPlotNodes + 2; ++iPlot) {
      _x[cnt] = v0->x() + iPlot * dx;

      // A little offset to stay inside the elements
      if(iPlot == 0) {
        _x[cnt] += eps;
      }
      if(iPlot == _nInteriorPlotNodes - 1) {
        _x[cnt] -= eps;
      }

      // FE solution uh at plot point x
      xloc = {_x[cnt], 0., 0.};
      if(mesh.locateVertex(xloc.data(), elm, r)) {
        space->initializeAddressingVector(elm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          elmSolution[i] = sol[adr[i]];
        }
        _uh[cnt] = yScaling * space->interpolateField(elmSolution, r);
      } else {
        _uh[cnt] = 0.;
        feWarning("Could not find plot node %1.4e in the mesh :/", _x[cnt]);
      }

      // Exact solution
      if(analyticSolution != nullptr) {
        _u[cnt] = yScaling * analyticSolution->eval(solution.getCurrentTime(), xloc);
      }

      cnt++;
    }

    // std::vector<double> xxx(3);
    // std::vector<double> localCoord(3*2);
    // mesh.getCoord(space->getCncGeoID(), iElm, localCoord);
    // double rr[3] = {-1., 0., 0.};
    // mesh.getGeometricSpace(space->getCncGeoID())->interpolateVectorField(localCoord, rr, xxx);

    // space->initializeAddressingVector(numbering.getNumbering(space->getFieldID()), iElm, adr);
    // for(size_t i = 0; i < adr.size(); ++i) {
    //   elmSolution[i] = sol[adr[i]];
    // }
    // double myVal = space->interpolateField(elmSolution, rr);
    // feInfo("Interpolation (graphics) en x = %f : uh = %f - u = %f - adr = %d %d %d", xxx[0],
    // myVal, space->evalFun(0., xxx), adr[0], adr[1], adr[2]); feInfo("Interpolation (graphics) en
    // x = %f : uh = %f - u = %f", xxx[0], 0., space->evalFun(0., xxx));
  }

  // Plot FE solution
  switch(color) {
    case 0:
      glColor3f(0.0, 0.0, 1.0);
      break;
    case 1:
      glColor3f(1.0, 0.0, 0.0);
      break;
    case 2:
      glColor3f(0.0, 1.0, 0.0);
      break;
    case 3:
      glColor3f(0.0, 0.0, 0.0);
      break;
  }
  glfemDrawCurveDG(_x, _uh, _nElm, _nInteriorPlotNodes + 2);
  // Plot exact solution
  if(analyticSolution != nullptr) {
    glColor3f(1.0, 0.0, 0.0);
    glfemDrawCurve(_x, _u, _nPlotNodes);
  }

  if(_drawAxes) {
    drawAxes();
  }

#else
  UNUSED(mesh, numbering, solution,space, analyticSolution, yScaling, color);
#endif
}

#if defined(HAVE_GLFW)
static double area(double x)
{
  return 1. + 2.2 * (x - 1.5) * (x - 1.5);
  // return 1.398 + 0.347 * tanh(0.8*x-4.);
}
#endif

void feBasicViewer::drawEulerNozzle(feMesh &mesh, feMetaNumber &numbering, feSolution &solution,
                                    std::vector<feSpace *> &spaces, double yScaling[3],
                                    feFunction *analyticSolution)
{
#if defined(HAVE_GLFW)
  std::vector<double> &sol = solution.getSolutionReference();
  std::vector<double> elmSolution(spaces[0]->getNumFunctions());
  std::vector<double> xloc;
  std::vector<feInt> adr(spaces[0]->getNumFunctions());

  int elm, cnt = 0;
  double r[3], eps = 1e-4;
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    Vertex *v0 = mesh.getVertex(iElm);
    Vertex *v1 = mesh.getVertex(iElm + 1);
    double dx = fabs(v1->x() - v0->x()) / (_nInteriorPlotNodes + 1);

    for(int iPlot = 0; iPlot < _nInteriorPlotNodes + 2; ++iPlot) {
      _x[cnt] = v0->x() + iPlot * dx;

      // A little offset to stay inside the elements
      if(iPlot == 0) {
        _x[cnt] += eps;
      }
      if(iPlot == _nInteriorPlotNodes - 1) {
        _x[cnt] -= eps;
      }

      // FE solution uh at plot point x
      xloc = {_x[cnt], 0., 0.};
      double A = area(_x[cnt]);
      double rhoA, rhouA, eA, rho, u, p, e, c;
      if(mesh.locateVertex(xloc.data(), elm, r)) {
        // rho
        spaces[0]->initializeAddressingVector(elm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          elmSolution[i] = sol[adr[i]];
        }
        rhoA = spaces[0]->interpolateField(elmSolution, r);

        // rhouA
        spaces[1]->initializeAddressingVector(elm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          elmSolution[i] = sol[adr[i]];
        }
        rhouA = spaces[1]->interpolateField(elmSolution, r);

        // eA
        spaces[2]->initializeAddressingVector(elm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          elmSolution[i] = sol[adr[i]];
        }
        eA = spaces[2]->interpolateField(elmSolution, r);

        rho = rhoA / A;
        u = rhouA / rhoA;
        e = eA / A;
        p = (e - rho * u * u / 2.) * 0.4; // gamma-1 = 0.4
        c = sqrt(1.4 * p / rho);

        _rho[cnt] = yScaling[0] * rho;
        _mach[cnt] = yScaling[1] * u / c;
        _p[cnt] = yScaling[2] * p;
        // feInfo("scaled rho = %f", _rho[cnt]);
        // feInfo("scaled mach = %f", _mach[cnt]);
        // feInfo("scaled p = %f", _p[cnt]);

      } else {
        _uh[cnt] = 0.;
        feWarning("Could not find plot node %1.4e in the mesh :/", _x[cnt]);
      }

      // Exact solution
      if(analyticSolution != nullptr) {
        _u[cnt] = yScaling[0] * analyticSolution->eval(solution.getCurrentTime(), xloc);
      }

      cnt++;
    }
  }

  // Plot FE solution
  glColor3f(0.0, 0.0, 1.0);
  glfemDrawCurveDG(_x, _rho, _nElm, _nInteriorPlotNodes + 2);
  glColor3f(1.0, 0.0, 0.0);
  glfemDrawCurveDG(_x, _mach, _nElm, _nInteriorPlotNodes + 2);
  glColor3f(0.0, 1.0, 0.0);
  glfemDrawCurveDG(_x, _p, _nElm, _nInteriorPlotNodes + 2);

  // Plot exact solution
  if(analyticSolution != nullptr) {
    glColor3f(1.0, 0.0, 0.0);
    glfemDrawCurve(_x, _u, _nPlotNodes);
  }

  if(_drawAxes) {
    drawAxes();
  }

#else
  UNUSED(mesh, numbering, solution, spaces, yScaling, analyticSolution);
#endif
}