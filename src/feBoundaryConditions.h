#ifndef _FEBOUNDARYCONDITIONS_
#define _FEBOUNDARYCONDITIONS_

#include "feFunction.h"
#include "feParameters.h"
#include "feSpace.h"

//
// The description required to create a finite element space
//
class FEDescriptor
{
public:
  std::string _fieldName;
  std::string _physicalEntityName;

  // The index in the vector of FE spaces once they have
  // been created.
  unsigned int _space_id;

  elementType _feType;

  int _degree;
  int _degreeQuadrature;

  bool                    _isScalar;
  const feFunction       *_scalarField;
  const feVectorFunction *_vectorField;

  bool             _isEssential;
  std::vector<int> _essentialComponents;

public:
  FEDescriptor(const std::string       fieldName,
               const std::string       physicalEntityName,
               const elementType      &feType,
               const int               polynomialOrder,
               const int               degreeQuadrature,
               const bool              isScalar,
               const feFunction       *scalarCallback,
               const feVectorFunction *vectorCallback)
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

  // A comparator which checks if two spaces:
  //
  // - are associated to the same field
  // - are defined on Physical groups with the same name
  //
  // Both spaces can be discretized by different finite elements.
  bool representsSameFieldAs(const FEDescriptor &other) const
  {
    bool same = true;
    same &= this->_fieldName == other._fieldName;
    same &= this->_physicalEntityName == other._physicalEntityName;
    return same;
  }
};

class ScalarFEDescriptor : public FEDescriptor
{
public:
  ScalarFEDescriptor(const std::string  fieldName,
                     const std::string  physicalEntityName,
                     const elementType &feType,
                     const int          polynomialOrder,
                     const int          degreeQuadrature,
                     const feFunction  *initCallback)
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
  VectorFEDescriptor(const std::string       fieldName,
                     const std::string       physicalEntityName,
                     const elementType      &feType,
                     const int               polynomialOrder,
                     const int               degreeQuadrature,
                     const feVectorFunction *initCallback)
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

namespace BoundaryConditions
{
  // FIXME: this is a mix of "fundamental" BC categories
  // (Dirichlet/Neumann/Robin), and of specific cases for
  // fluid flows BC (e.g. slip (mix Dirichlet-Neumann)
  // or no-slip (full Dirichlet).
  enum class Type
  {
    dirichlet,
    dirichlet_one_component,
    periodic,
    contactAngle
  };

  std::string typeToString(const Type &type);

  //
  // Base class for a boundary condition
  //
  class BoundaryCondition
  {
  public:
    // Description of the finite element space used for the trial space.
    // The trial space is the space affected by the boundary condition.
    FEDescriptor _descriptor_trialSpace;

    // Description of the finite element space used for the test space. The test
    // space is associated with the equation in which the (non-essential)
    // boundary condition appears, for instance a non-zero flux. Usually, the
    // test space will be identical to the trial space, unless the flux of a
    // variable appears in another equation, for instance atfer integrating by
    // parts the coupled diffusion terms in the Cahn-Hilliard system.
    FEDescriptor _descriptor_testSpace;

    Type _boundaryType;

  public:
    BoundaryCondition(const FEDescriptor &descriptor_trialSpace,
                      const FEDescriptor &descriptor_testSpace,
                      const Type         &boundaryType)
      : _descriptor_trialSpace(descriptor_trialSpace)
      , _descriptor_testSpace(descriptor_testSpace)
      , _boundaryType(boundaryType){};
    virtual ~BoundaryCondition() = default;
  };

  //
  // Scalar-valued Dirichlet boundary condition
  //
  class ScalarDirichlet : public BoundaryCondition
  {
  public:
    ScalarDirichlet(const ScalarFEDescriptor &descriptor)
      : BoundaryCondition(descriptor, descriptor, Type::dirichlet)
    {
      _descriptor_trialSpace._isEssential         = true;
      _descriptor_trialSpace._essentialComponents = {1, 1, 1};
    }
  };

  //
  // Vector-valued Dirichlet boundary condition
  //
  class VectorDirichlet : public BoundaryCondition
  {
  public:
    VectorDirichlet(const VectorFEDescriptor &descriptor)
      : BoundaryCondition(descriptor, descriptor, Type::dirichlet)
    {
      _descriptor_trialSpace._isEssential         = true;
      _descriptor_trialSpace._essentialComponents = {1, 1, 1};
    }
  };

  //
  // Vector-valued Dirichlet boundary condition
  // for a single component of the field.
  //
  class DirichletSingleComponent : public BoundaryCondition
  {
  public:
    DirichletSingleComponent(const VectorFEDescriptor &descriptor,
                             const int                 fixedComponent)
      : BoundaryCondition(descriptor, descriptor, Type::dirichlet_one_component)
    {
      _descriptor_trialSpace._isEssential                         = false;
      _descriptor_trialSpace._essentialComponents                 = {0, 0, 0};
      _descriptor_trialSpace._essentialComponents[fixedComponent] = 1;
    }
  };

  //
  // Slip boundary condition.
  // For vector-valued fields, typically the velocity.
  //
  class SlipX : public DirichletSingleComponent
  {
  public:
    SlipX(const VectorFEDescriptor &descriptor)
      : DirichletSingleComponent(descriptor, 1)
    {}
  };

  class SlipY : public DirichletSingleComponent
  {
  public:
    SlipY(const VectorFEDescriptor &descriptor)
      : DirichletSingleComponent(descriptor, 0)
    {}
  };

  //
  // Periodic boundary condition
  //
  // Matching periodic boundaries are detected for translation only.
  // The offset vector is the translation from a vertex on boundary 1
  // to its position on boundary 2, that is, v2 = v1 + offset.
  //
  class ScalarPeriodic : public BoundaryCondition
  {
  public:
    FEDescriptor _matching_descriptor;
    const std::vector<double> &_offset;
  public:
    ScalarPeriodic(const ScalarFEDescriptor &descriptor,
                   const ScalarFEDescriptor &matching_descriptor,
                   const std::vector<double> &offset)
      : BoundaryCondition(descriptor, descriptor, Type::periodic)
      , _matching_descriptor(matching_descriptor)
      , _offset(offset)
    {
      _descriptor_trialSpace._isEssential         = false;
      _descriptor_trialSpace._essentialComponents = {0, 0, 0};
    }
  };

  class VectorPeriodic : public BoundaryCondition
  {
  public:
    FEDescriptor _matching_descriptor;
    const std::vector<double> &_offset;
  public:
    VectorPeriodic(const VectorFEDescriptor &descriptor,
                   const VectorFEDescriptor &matching_descriptor,
                   const std::vector<double> &offset)
      : BoundaryCondition(descriptor, descriptor, Type::periodic)
      , _matching_descriptor(matching_descriptor)
      , _offset(offset)
    {
      _descriptor_trialSpace._isEssential         = false;
      _descriptor_trialSpace._essentialComponents = {0, 0, 0};
    }
  };

  //
  // Flux boundary conditions for scalar-valued fields,
  // where only the gradient angle is imposed, but not
  // the complete value of the flux:
  //
  //   ∇u . n = - |∇u| cos(contact_angle)
  //
  // The contactAngle parameter is given in DEGREES.
  //
  class FluxContactAngle : public BoundaryCondition
  {
  public:
    double _contactAngle = -1.;

  public:
    FluxContactAngle(const ScalarFEDescriptor &descriptor_trialSpace,
                     const ScalarFEDescriptor &descriptor_testSpace,
                     const double              contactAngle)
      : BoundaryCondition(descriptor_trialSpace,
                          descriptor_testSpace,
                          Type::contactAngle)
      , _contactAngle(contactAngle)
    {}

    double getContactAngle() const { return _contactAngle; }
  };
} // namespace BoundaryCondition

#endif