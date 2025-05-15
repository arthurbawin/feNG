#include "feSpace.h"

template <int dim>
void feVectorSpace<dim>::dotProductShapeShape(const int            iNode,
                                              std::vector<double> &res) const
{
  int offset = iNode * _nFunctions * _nComponents;

  for (int i = 0; i < _nFunctions; ++i)
  {
    const double *phi_i = _Ldata + offset + i * _nComponents;
    for (int j = 0; j < _nFunctions; ++j)
    {
      const double *phi_j = _Ldata + offset + j * _nComponents;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        res[i * _nFunctions + j] = phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1];
      }

      if constexpr (dim == 3)
      {
        res[i * _nFunctions + j] =
          phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2];
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::dotProductShapeOther(const int                  iNode,
                                              const std::vector<double> &other,
                                              std::vector<double> &res) const
{
  int           offset  = iNode * _nFunctions * _nComponents;
  const double *other_p = other.data();
  const double *phi_i   = _Ldata + offset;
  double       *res_p   = res.data();
  const double  v0 = other_p[0], v1 = other_p[1];

  for (int i = 0; i < _nFunctions; ++i, phi_i += _nComponents)
  {
    static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

    if constexpr (dim == 2)
    {
      res_p[i] = phi_i[0] * v0 + phi_i[1] * v1;
    }

    if constexpr (dim == 3)
    {
      res_p[i] = phi_i[0] * v0 + phi_i[1] * v1 + phi_i[2] * other_p[2];
    }
  }
}

// Wrong, implemented with a scalar phi in mind. Modify for vector-valued phi.
// template <int dim>
// void feVectorSpace<dim>::dotProductGradShapeOther(const int                  iNode,
//                                                   const std::vector<double> &gradPhi,
//                                                   const std::vector<double> &other,
//                                                   std::vector<double> &res) const
// {
//   int           offsetGrad = _dim * _nComponents;
//   const double *gradPhi_p  = gradPhi.data();
//   const double *other_p    = other.data();
//   double       *res_p      = res.data();
//   const double  v0 = other_p[0], v1 = other_p[1];

//   const double *gphi_i = gradPhi_p;
//   for (int i = 0; i < _nFunctions; ++i, gphi_i += offsetGrad)
//   {
//     if constexpr (dim == 2)
//     {
//       res_p[i] = gphi_i[0] * v0 + gphi_i[1] * v1;
//     }

//     if constexpr (dim == 3)
//     {
//       res_p[i] = gphi_i[0] * v0 + gphi_i[1] * v1 + gphi_i[2] * other_p[2];
//     }
//   }
// }

template <int dim>
void feVectorSpace<dim>::dotProductShapeShapeOtherSpace(
  const int            iNode,
  const feSpace       *other,
  std::vector<double> &res) const
{
  int nFunctionsOther  = other->getNumFunctions();
  int nComponentsOther = other->getNumComponents();
  int offset           = iNode * _nFunctions * _nComponents;

  const double *otherShape = other->shape().data();

  // !!! Assumes that this space and other have the same quadrature rule (using
  // iNode to offset)
  int offsetOther = iNode * nFunctionsOther * nComponentsOther;

  for (int i = 0; i < _nFunctions; ++i)
  {
    const double *phi_i = _Ldata + offset + i * _nComponents;
    for (int j = 0; j < nFunctionsOther; ++j)
    {
      const double *phi_j_other =
        otherShape + offsetOther + j * nComponentsOther;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        res[i * nFunctionsOther + j] =
          phi_i[0] * phi_j_other[0] + phi_i[1] * phi_j_other[1];
      }

      if constexpr (dim == 3)
      {
        res[i * nFunctionsOther + j] = phi_i[0] * phi_j_other[0] +
                                       phi_i[1] * phi_j_other[1] +
                                       phi_i[2] * phi_j_other[2];
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::dotProductShapeGradShapeOtherSpace(
  const int                  iNode,
  const int                  nFunctionsOther,
  const std::vector<double> &gradOtherScalarShape,
  std::vector<double>       &res) const
{
  int           offset                 = iNode * _nFunctions * _nComponents;
  const double *gradOtherScalarShape_p = gradOtherScalarShape.data();

  for (int i = 0; i < _nFunctions; ++i)
  {
    const double *phi_i = _Ldata + offset + i * _nComponents;
    for (int j = 0; j < nFunctionsOther; ++j)
    {
      const double *gradOtherScalarShape_j =
        gradOtherScalarShape_p + j * _nComponents;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        res[i * nFunctionsOther + j] = phi_i[0] * gradOtherScalarShape_j[0] +
                                       phi_i[1] * gradOtherScalarShape_j[1];
      }

      if constexpr (dim == 3)
      {
        res[i * nFunctionsOther + j] = phi_i[0] * gradOtherScalarShape_j[0] +
                                       phi_i[1] * gradOtherScalarShape_j[1] +
                                       phi_i[2] * gradOtherScalarShape_j[2];
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::vectorDotGradShapeDotShape(
  const int                  iNode,
  const std::vector<double> &gradPhi,
  const std::vector<double> &other,
  std::vector<double>       &res) const
{
  int           offset     = iNode * _nFunctions * _nComponents;
  int           offsetGrad = _dim * _nComponents;
  const double *other_p    = other.data();
  const double *phi_i      = _Ldata + offset;
  const double *gradPhi_p  = gradPhi.data();
  double       *res_p      = res.data();

  for (int i = 0; i < _nFunctions;
       ++i, phi_i += _nComponents, res_p += _nFunctions)
  {
    const double *gphi_j = gradPhi_p;
    for (int j = 0; j < _nFunctions; ++j, gphi_j += offsetGrad)
    {
      if constexpr (dim == 2)
      {
        // res_p[i * _nFunctions + j] = other_p[0] * gphi_j[0] * phi_i[0]
        //                            + other_p[0] * gphi_j[1] * phi_i[1]
        //                            + other_p[1] * gphi_j[2] * phi_i[0]
        //                            + other_p[1] * gphi_j[3] * phi_i[1];
        res_p[j] = other_p[0] * (gphi_j[0] * phi_i[0] + gphi_j[1] * phi_i[1]) +
                   other_p[1] * (gphi_j[2] * phi_i[0] + gphi_j[3] * phi_i[1]);
      }

      if constexpr (dim == 3)
      {
        for (int m = 0; m < _nComponents; ++m)
        {
          for (int n = 0; n < _nComponents; ++n)
          {
            res_p[i * _nFunctions + j] +=
              other_p[m] * gphi_j[m * _dim + n] * phi_i[n];
          }
        }
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::shapeDotTensorDotShape(
  const int                  iNode,
  const std::vector<double> &other,
  std::vector<double>       &res) const
{
  int           offset  = iNode * _nFunctions * _nComponents;
  const double *other_p = other.data();
  const double *phi_i   = _Ldata + offset;
  double       *res_p   = res.data();

  for (int i = 0; i < _nFunctions;
       ++i, phi_i += _nComponents, res_p += _nFunctions)
  {
    // const double *phi_i = _Ldata + offset + i * _nComponents;
    const double *phi_j = _Ldata + offset;
    for (int j = 0; j < _nFunctions; ++j, phi_j += _nComponents)
    {
      // const double *phi_j = _Ldata + offset + j * _nComponents;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        // res_p[i * _nFunctions + j] = phi_j[0] * other_p[0] * phi_i[0]
        //                            + phi_j[0] * other_p[1] * phi_i[1]
        //                            + phi_j[1] * other_p[2] * phi_i[0]
        //                            + phi_j[1] * other_p[3] * phi_i[1];
        res_p[j] = phi_j[0] * (other_p[0] * phi_i[0] + other_p[1] * phi_i[1]) +
                   phi_j[1] * (other_p[2] * phi_i[0] + other_p[3] * phi_i[1]);
      }

      if constexpr (dim == 3)
      {
        double val = 0.;
        for (int m = 0; m < _nComponents; ++m)
        {
          for (int n = 0; n < _nComponents; ++n)
          {
            val += phi_j[m] * other_p[m * _dim + n] * phi_i[n];
          }
        }
        res_p[i * _nFunctions + j] = val;
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::gradOtherScalarShapeDotTensorDotShape(
  const int                  iNode,
  const int                  nFunctionsOther,
  const std::vector<double> &gradOtherScalarShape,
  const std::vector<double> &other,
  std::vector<double>       &res) const
{
  int           offset                 = iNode * _nFunctions * _nComponents;
  const double *other_p                = other.data();
  const double *gradOtherScalarShape_p = gradOtherScalarShape.data();

  for (int i = 0; i < _nFunctions; ++i)
  {
    const double *phi_i = _Ldata + offset + i * _nComponents;
    for (int j = 0; j < nFunctionsOther; ++j)
    {
      const double *gradOtherScalarShape_j =
        gradOtherScalarShape_p + j * _nComponents;

      // double val = 0.;
      // for(int m = 0; m < _nComponents; ++m) {
      //   for(int n = 0; n < _nComponents; ++n) {
      //     // res[i * nFunctionsOther + j] += gradOtherScalarShape[j *
      //     _nComponents + m] * other_p[m * _dim + n] * _L[offset + i *
      //     _nComponents + n]; val += gradOtherScalarShape_j[m] * other_p[m
      //     * _dim + n] * phi_i[n];
      //   }
      // }
      // res[i * nFunctionsOther + j] = val;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        res[i * nFunctionsOther + j] =
          gradOtherScalarShape_j[0] * other_p[0] * phi_i[0] +
          gradOtherScalarShape_j[0] * other_p[1] * phi_i[1] +
          gradOtherScalarShape_j[1] * other_p[2] * phi_i[0] +
          gradOtherScalarShape_j[1] * other_p[3] * phi_i[1];
      }

      if constexpr (dim == 3)
      {
        double val = 0.;
        for (int m = 0; m < _nComponents; ++m)
        {
          for (int n = 0; n < _nComponents; ++n)
          {
            val += gradOtherScalarShape_j[m] * other_p[m * _dim + n] * phi_i[n];
          }
        }
        res[i * nFunctionsOther + j] = val;
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::doubleContractionGradShapeGradShape(
  const std::vector<double> &gradPhi,
  std::vector<double>       &res) const
{
  int           offset    = _dim * _nComponents;
  const double *gradPhi_p = gradPhi.data();
  double       *res_p     = res.data();

  for (int i = 0; i < _nFunctions; ++i, res_p += _nFunctions)
  {
    const double *gradphi_i = gradPhi_p + i * offset;
    for (int j = 0; j < _nFunctions; ++j)
    {
      const double *gradphi_j = gradPhi_p + j * offset;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        res_p[j] = gradphi_i[0] * gradphi_j[0] + gradphi_i[1] * gradphi_j[1] +
                   gradphi_i[2] * gradphi_j[2] + gradphi_i[3] * gradphi_j[3];
      }

      if constexpr (dim == 3)
      {
        double val = 0.;
        for (int m = 0; m < _nComponents; ++m)
        {
          for (int n = 0; n < _nComponents; ++n)
          {
            // res[i * _nFunctions + j] += gradPhi[i * offset + m *
            // _dim + n] * gradPhi[j * offset + m * _dim + n];
            val += gradphi_i[m * _dim + n] * gradphi_j[m * _dim + n];
          }
        }
        res[i * _nFunctions + j] = val;
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::doubleContractionGradShapeGradShapeTransposed(
  const std::vector<double> &gradPhi,
  std::vector<double>       &res) const
{
  int           offset    = _dim * _nComponents;
  const double *gradphi_p = gradPhi.data();
  const double *gradphi_i = gradphi_p;
  double       *res_p     = res.data();

  for (int i = 0; i < _nFunctions;
       ++i, gradphi_i += offset, res_p += _nFunctions)
  {
    // const double *gradphi_i = gradPhi_p + i * offset;

    const double *gradphi_j = gradphi_p;
    for (int j = 0; j < _nFunctions; ++j, gradphi_j += offset)
    {
      // const double *gradphi_j = gradPhi_p + j * offset;

      // // res[i * _nFunctions + j] = 0.;
      // double val = 0.;
      // for(int m = 0; m < _nComponents; ++m) {
      //   for(int n = 0; n < _nComponents; ++n) {
      //     // res[i * _nFunctions + j] += gradPhi[i * offset + m * _dim +
      //     n] * gradPhi[j * offset + n * _dim + m]; val += gradphi_i[m *
      //     _dim + n] * gradphi_j[n * _dim + m];
      //   }
      // }
      // res[i * _nFunctions + j] = val;

      static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");

      if constexpr (dim == 2)
      {
        res_p[j] = gradphi_i[0] * gradphi_j[0] + gradphi_i[1] * gradphi_j[2] +
                   gradphi_i[2] * gradphi_j[1] + gradphi_i[3] * gradphi_j[3];
      }

      if constexpr (dim == 3)
      {
        double val = 0.;
        for (int m = 0; m < _nComponents; ++m)
        {
          for (int n = 0; n < _nComponents; ++n)
          {
            // res[i * _nFunctions + j] += gradPhi[i * offset + m *
            // _dim + n] * gradPhi[j * offset + n * _dim + m];
            val += gradphi_i[m * _dim + n] * gradphi_j[n * _dim + m];
          }
        }
        res[i * _nFunctions + j] = val;
      }
    }
  }
}

template <int dim>
void feVectorSpace<dim>::doubleContractionGradShapeOther(
  const std::vector<double> &gradPhi,
  const std::vector<double> &other,
  std::vector<double>       &res) const
{
  int           offset    = _dim * _nComponents;
  const double *other_p   = other.data();
  const double *gradphi_i = gradPhi.data();

  for (int i = 0; i < _nFunctions; ++i, gradphi_i += offset)
  {
    // const double *gradphi_i = gradPhi_p + i * offset;

    static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");
    if constexpr (dim == 2)
    {
      res[i] = gradphi_i[0] * other_p[0] + gradphi_i[1] * other_p[1] +
               gradphi_i[2] * other_p[2] + gradphi_i[3] * other_p[3];
    }
    if constexpr (dim == 3)
    {
      double val = 0.;
      for (int m = 0; m < _nComponents; ++m)
      {
        for (int n = 0; n < _nComponents; ++n)
        {
          // res[i] += gradPhi[i * offset + m * _dim + n] * other_p[m *
          // _dim + n];
          val += gradphi_i[m * _dim + n] * other_p[m * _dim + n];
        }
      }
      res[i] = val;
    }
  }
}

template <int dim>
void feVectorSpace<dim>::divergence(const std::vector<double> &gradPhi,
                                    std::vector<double>       &res) const
{
  int           offset    = _dim * _nComponents;
  const double *gradphi_i = gradPhi.data();

  for (int i = 0; i < _nFunctions; ++i, gradphi_i += offset)
  {
    static_assert(dim == 2 || dim == 3, "Dim should be 2 or 3");
    if constexpr (dim == 2)
    {
      res[i] = gradphi_i[0] + gradphi_i[3];
    }

    if constexpr (dim == 3)
    {
      res[i] = gradphi_i[0] + gradphi_i[4] + gradphi_i[8];
    }
  }
}

//
// Explicit instantiation
//
template void
feVectorSpace<2>::dotProductShapeShape(const int            iNode,
                                       std::vector<double> &res) const;

template void
              feVectorSpace<2>::dotProductShapeOther(const int                  iNode,
                                       const std::vector<double> &other,
                                       std::vector<double>       &res) const;
template void feVectorSpace<2>::dotProductShapeShapeOtherSpace(
  const int            iNode,
  const feSpace       *other,
  std::vector<double> &res) const;
template void feVectorSpace<2>::dotProductShapeGradShapeOtherSpace(
  const int                  iNode,
  const int                  nFunctionsOther,
  const std::vector<double> &gradOtherScalarShape,
  std::vector<double>       &res) const;
template void
feVectorSpace<2>::vectorDotGradShapeDotShape(const int                  iNode,
                                             const std::vector<double> &gradPhi,
                                             const std::vector<double> &other,
                                             std::vector<double> &res) const;
template void
              feVectorSpace<2>::shapeDotTensorDotShape(const int                  iNode,
                                         const std::vector<double> &other,
                                         std::vector<double>       &res) const;
template void feVectorSpace<2>::gradOtherScalarShapeDotTensorDotShape(
  const int                  iNode,
  const int                  nFunctionsOther,
  const std::vector<double> &gradOtherScalarShape,
  const std::vector<double> &other,
  std::vector<double>       &res) const;
template void feVectorSpace<2>::doubleContractionGradShapeGradShape(
  const std::vector<double> &gradPhi,
  std::vector<double>       &res) const;
template void feVectorSpace<2>::doubleContractionGradShapeGradShapeTransposed(
  const std::vector<double> &gradPhi,
  std::vector<double>       &res) const;
template void feVectorSpace<2>::doubleContractionGradShapeOther(
  const std::vector<double> &gradPhi,
  const std::vector<double> &other,
  std::vector<double>       &res) const;
template void feVectorSpace<2>::divergence(const std::vector<double> &gradPhi,
                                           std::vector<double> &res) const;