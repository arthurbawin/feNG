
#include "feBoundaryConditions.h"

std::string BoundaryConditions::typeToString(const Type &type)
{
  switch (type)
  {
    case Type::dirichlet:
      return "dirichlet";
    case Type::dirichlet_one_component:
      return "dirichlet_one_component";
    case Type::periodic:
      return "periodic";
    case Type::contactAngle:
      return "contactAngle";
    
    // Not adding a "default" clause allows to return an error
    // at compile time when a new type of BC was added to the
    // enum, but was not added here.  
  };

  return "Undefined BC";
}