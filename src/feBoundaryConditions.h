#ifndef _FEBOUNDARYCONDITIONS_
#define _FEBOUNDARYCONDITIONS_

#include "feFunction.h"
#include "feSpace.h"

//
// The description required to create a finite element space
//
class FEDescriptor
{
public:
  std::string _fieldName;
  std::string _physicalEntityName;

  elementType _feType;

  int _degree;
  int _degreeQuadrature;

  bool                    _isScalar;
  const feFunction       *_scalarField;
  const feVectorFunction *_vectorField;

  bool              _isEssential;
  std::vector<bool> _essentialComponents;

public:
  FEDescriptor(const std::string  fieldName,
               const std::string  physicalEntityName,
               const elementType &feType,
               const int          polynomialOrder,
               const int          degreeQuadrature,
               const bool         isScalar,
               const feFunction        *scalarCallback,
               const feVectorFunction  *vectorCallback)
    : _fieldName(fieldName)
    , _physicalEntityName(physicalEntityName)
    , _feType(feType)
    , _degree(polynomialOrder)
    , _degreeQuadrature(degreeQuadrature)
    , _isScalar(isScalar)
    , _scalarField(scalarCallback)
    , _vectorField(vectorCallback)
    , _isEssential(false)
    , _essentialComponents({0, 0, 0})
  {}
};

class ScalarFEDescriptor : public FEDescriptor
{
public:
  ScalarFEDescriptor(const std::string  fieldName,
                     const std::string  physicalEntityName,
                     const elementType &feType,
                     const int          polynomialOrder,
                     const int          degreeQuadrature,
                     const feFunction        *initCallback)
    : FEDescriptor(fieldName,
                   physicalEntityName,
                   feType,
                   polynomialOrder,
                   degreeQuadrature,
                   true,
                   initCallback,
                   nullptr)
    {}
};

class VectorFEDescriptor : public FEDescriptor
{
public:
  VectorFEDescriptor(const std::string  fieldName,
                     const std::string  physicalEntityName,
                     const elementType &feType,
                     const int          polynomialOrder,
                     const int          degreeQuadrature,
                     const feVectorFunction        *initCallback)
    : FEDescriptor(fieldName,
                   physicalEntityName,
                   feType,
                   polynomialOrder,
                   degreeQuadrature,
                   false,
                   nullptr,
                   initCallback)
    {}
};

class BoundaryCondition
{
public:
  FEDescriptor _descriptor;

public:
  BoundaryCondition(const FEDescriptor &descriptor)
    : _descriptor(descriptor){};
};

class ScalarDirichletBC : public BoundaryCondition
{
public:
  ScalarDirichletBC(const ScalarFEDescriptor &descriptor)
  : BoundaryCondition(descriptor)
  {
    _descriptor._isEssential = true;
    _descriptor._essentialComponents = {1, 1, 1};
  }
};

class VectorDirichletBC : public BoundaryCondition
{
public:
  VectorDirichletBC(const VectorFEDescriptor &descriptor)
    : BoundaryCondition(descriptor)
  {
    _descriptor._isEssential = true;
    _descriptor._essentialComponents = {1, 1, 1};
  }
};

class DirichletSingleComponentBC : public BoundaryCondition
{
public:
  DirichletSingleComponentBC(const VectorFEDescriptor &descriptor,
                             const int               fixedComponent)
    : BoundaryCondition(descriptor)
  {
    _descriptor._isEssential = false;
    _descriptor._essentialComponents = {0, 0, 0};
    _descriptor._essentialComponents[fixedComponent] = 1;
  }
};

//
// Slip boundary condition.
// For vector-valued fields, typically the velocity.
//
class SlipX_BC : public DirichletSingleComponentBC
{
public:
  SlipX_BC(const VectorFEDescriptor &descriptor)
    : DirichletSingleComponentBC(descriptor, 1)
  {}
};

class SlipY_BC : public DirichletSingleComponentBC
{
public:
  SlipY_BC(const VectorFEDescriptor &descriptor)
    : DirichletSingleComponentBC(descriptor, 0)
  {}
};

#endif