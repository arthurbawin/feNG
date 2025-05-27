
#include "GenericSolver.h"

#include "feNorm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

feStatus SolverBase::createSpaces(feMesh                 *mesh,
                                  std::vector<feSpace *> &spaces,
                                  std::vector<feSpace *> &essentialSpaces) const
{
  spaces.clear();
  essentialSpaces.clear();

  //
  // Unknown fields
  //
  for (const auto &d : _fieldDescriptors)
  {
    void *initFct =
      d._isScalar ? (void *)d._scalarField : (void *)d._vectorField;

    feSpace *s;
    feCheckReturn(createFiniteElementSpace(s,
                                           mesh,
                                           d._feType,
                                           d._degree,
                                           d._fieldName,
                                           d._physicalEntityName,
                                           d._degreeQuadrature,
                                           initFct));
    spaces.push_back(s);
  }

  //
  // Boundary conditions
  //
  for (const auto &b : _boundaryConditions)
  {
    const FEDescriptor &d = b._descriptor;

    void *initFct =
      d._isScalar ? (void *)d._scalarField : (void *)d._vectorField;

    feSpace *s;
    feCheckReturn(createFiniteElementSpace(s,
                                           mesh,
                                           d._feType,
                                           d._degree,
                                           d._fieldName,
                                           d._physicalEntityName,
                                           d._degreeQuadrature,
                                           initFct));
    spaces.push_back(s);

    // Print a summary of the boundary condition
    feInfo("");

    //
    // Assign essential components
    //
    if (d._isEssential)
    {
      essentialSpaces.push_back(s);
    }
    else
    {
      for (int i = 0; i < 3; ++i)
      {
        if (d._essentialComponents[i])
        {
          s->setEssentialComponent(i, true);
        }
      }
    }
  }

  return FE_STATUS_OK;
}

feStatus
SolverBase::projectSolution(feMesh2DP1                  *currentMesh,
                            feMesh2DP1                  *nextMesh,
                            feSolutionContainer         *container,
                            feSolutionContainer         &bufferContainer,
                            feSolutionContainer        *&nextContainer,
                            feMetaNumber                *numbering,
                            const std::vector<feSpace *> spaces,
                            const std::vector<feSpace *> essentialSpaces) const
{
  // Copy current container into the buffer
  bufferContainer = *container;

  // For checks : last BDF solutions
  std::vector<feSolution> currentSolutions(bufferContainer.getNbSol());
  for (int i = 0; i < bufferContainer.getNbSol(); ++i)
  {
    currentSolutions[i] =
      feSolution(bufferContainer, i, spaces, essentialSpaces);
  }

  // Projection
  std::vector<feSpace *> nextSpaces, nextEssentialSpaces;
  feCheckReturn(this->createSpaces(nextMesh, nextSpaces, nextEssentialSpaces));
  feMetaNumber nextNumbering(nextMesh, nextSpaces, nextEssentialSpaces);
  feCheckReturn(currentMesh->transfer(nextMesh,
                                      numbering,
                                      &nextNumbering,
                                      &bufferContainer,
                                      spaces,
                                      essentialSpaces,
                                      nextSpaces));

  delete nextContainer;
  nextContainer = new feSolutionContainer(bufferContainer);

  // Check projection (visualization) and integral of each projected solution
  // in container
  for (int i = 0; i < bufferContainer.getNbSol(); ++i)
  {
    feInfo("Checking solution %d in bufferContainer :", i);
    feSolution projectedSolution(bufferContainer, i, spaces, essentialSpaces);
    // feExporter *nextExporter;
    // feCheck(createVisualizationExporter(nextExporter, VTK, &nextNumbering,
    // &projectedSolution, nextMesh, nextSpaces));
    // feCheck(nextExporter->writeStep("afterTransfer" + std::to_string(i) +
    // ".vtk")); delete nextExporter;

    // Integral of highest dimensionsal fields before and after transfer
    for (size_t j = 0; j < spaces.size(); ++j)
    {
      if (spaces[j]->getDim() == 2)
      {
        feNorm *integral, *integralNextMesh, *intDot, *intDotNextMesh;
        feCheck(
          createNorm(integral, INTEGRAL, {spaces[j]}, &currentSolutions[i]));
        feCheck(createNorm(
          integralNextMesh, INTEGRAL, {nextSpaces[j]}, &projectedSolution));
        feInfo("Integral of field %s on current mesh = %+-1.6e",
               spaces[j]->getFieldID().data(),
               integral->compute());
        feInfo("Integral of field %s on next    mesh = %+-1.6e",
               nextSpaces[j]->getFieldID().data(),
               integralNextMesh->compute());
        feCheck(
          createNorm(intDot, INTEGRAL_DT, {spaces[j]}, &currentSolutions[i]));
        feCheck(createNorm(
          intDotNextMesh, INTEGRAL_DT, {nextSpaces[j]}, &projectedSolution));
        feInfo("Integral of field %s dot on current mesh = %+-1.6e",
               spaces[j]->getFieldID().data(),
               intDot->compute());
        feInfo("Integral of field %s dot on next    mesh = %+-1.6e",
               nextSpaces[j]->getFieldID().data(),
               intDotNextMesh->compute());
        delete integral;
        delete integralNextMesh;
        delete intDot;
        delete intDotNextMesh;
      }
    }
  }

  for (feSpace *s : nextSpaces)
    delete s;

  return FE_STATUS_OK;
}

feStatus SolverBase::readReferenceTestCase(const std::string meshName,
                                           const std::string solutionFileName,
                                           feMesh          *&mesh,
                                           std::vector<feSpace *> &spaces,
                                           feMetaNumber          *&numbering,
                                           feSolution            *&sol) const
{
  mesh = new feMesh2DP1(meshName);

  // Essential spaces are not needed for error computation, and are not
  // returned
  std::vector<feSpace *> refEssentialSpaces;
  this->createSpaces(mesh, spaces, refEssentialSpaces);
  numbering = new feMetaNumber(mesh, spaces, refEssentialSpaces);

  std::filebuf fb;
  if (fb.open(solutionFileName, std::ios::in))
  {
    fb.close();

    // Read solution file if it exists
    sol = new feSolution(numbering->getNbDOFs(), solutionFileName);
  }
  else
  {
    return feErrorMsg(FE_STATUS_ERROR,
                      "Could not create reference test case because file "
                      "\"%s\" does not exist.",
                      solutionFileName.data());
  }
  return FE_STATUS_OK;
}