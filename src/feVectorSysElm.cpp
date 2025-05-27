#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feNumeric.h"

//
// SUPG element sizes and tau coefficient
//
static double hFortin(int dim, const double velocity[3], double normVelocity, int nFunctions,
                      const std::vector<double> &gradphi)
{
  // Compute h from (13.17) in Fortin & Garon
  double res = 0., dotprod;
  for(int i = 0; i < nFunctions; ++i) {
    dotprod = 0.;
    for(int iDim = 0; iDim < dim; ++iDim) {
      dotprod += velocity[iDim] / normVelocity * gradphi[i * dim + iDim];
    }
    res += dotprod * dotprod;
  }

  return sqrt(2.) / sqrt(res);
}

static double getInnerRadius(const std::vector<double> &triCoord)
{
  // radius of inscribed circle = 2 * Area / sum(Line_i)
  double dist[3], k = 0.;
  for(int i = 0; i < 3; i++) {
    double x0 = triCoord[3 * i + 0];
    double y0 = triCoord[3 * i + 1];
    double x1 = triCoord[3 * ((i + 1) % 3) + 0];
    double y1 = triCoord[3 * ((i + 1) % 3) + 1];
    dist[i] = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    k += 0.5 * dist[i];
  }
  double const area = std::sqrt(k * (k - dist[0]) * (k - dist[1]) * (k - dist[2]));
  return area / k;
}

static double getOuterRadius(const std::vector<double> &triCoord)
{
  // radius of circle circumscribing a triangle
  double dist[3], k = 0.0;
  for(int i = 0; i < 3; i++) {
    double x0 = triCoord[3 * i + 0];
    double y0 = triCoord[3 * i + 1];
    double x1 = triCoord[3 * ((i + 1) % 3) + 0];
    double y1 = triCoord[3 * ((i + 1) % 3) + 1];
    dist[i] = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    k += 0.5 * dist[i];
  }
  double const area = std::sqrt(k * (k - dist[0]) * (k - dist[1]) * (k - dist[2]));
  return dist[0] * dist[1] * dist[2] / (4 * area);
}

static double hCircumInscr(const std::vector<double> &triCoord)
{
  double hInscribed = 2. * getInnerRadius(triCoord);
  double hCircumscribed = 2. * getOuterRadius(triCoord);
  return (hInscribed + hCircumscribed) / 2.;
}

// Diameter of circle of equivalent area
inline double hEquivalentCircle(const double jac) { return 2. * sqrt(jac / (2. * M_PI)); }

template <int dim> 
static double compute_tau(const std::vector<double> &triCoord, double velocity[3],
                          double viscosity, int nFunctions, std::vector<double> &gradphi,
                          double jac)
{
  double normV = 0.;
  if constexpr(dim == 2) {
    normV = sqrt(velocity[0] * velocity[0]
               + velocity[1] * velocity[1]);
  }
  if constexpr(dim == 3) {
    normV = sqrt(velocity[0] * velocity[0]
               + velocity[1] * velocity[1]
               + velocity[2] * velocity[2]);
  }

  double h1 = hFortin(dim, velocity, normV, nFunctions, gradphi);
  if(normV < 1e-12) {
    return 0.;
  }
  double h2 = hCircumInscr(triCoord);
  double h3 = hEquivalentCircle(jac);

  UNUSED(h1, h2, h3);
  double h = h1;

  double r1 = 2. * normV / h;
  double r2 = 4. * viscosity / (h * h);

  return 1. / sqrt(r1 * r1 + 9. * r2 * r2);
}

// -----------------------------------------------------------------------------
// Linear form: vector-valued source term
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_VectorSource<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _SdotPhi.resize(_nFunctions);
}

template <int dim>
void feSysElm_VectorSource<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    (*_source)(form->_args, _S);

    uSpace->dotProductShapeOther(k, _S, _SdotPhi);

    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= _SdotPhi[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_VectorSource<0>;
// template class feSysElm_VectorSource<1>;
template class feSysElm_VectorSource<2>;
// template class feSysElm_VectorSource<3>;

// -----------------------------------------------------------------------------
// Linear form: weak gradient of scalar source term
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_GradSource<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradPhi.resize(_nComponents * _nComponents * _nFunctions);
  _divPhi.resize(_nFunctions);
}

template <int dim>
void feSysElm_GradSource<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, S;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    S = (*_source)(form->_args);

    // Get gradient and divergence of shape functions
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                         _gradPhi);
    uSpace->divergence(_gradPhi, _divPhi);

    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] += S * _divPhi[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_GradSource<0>;
// template class feSysElm_GradSource<1>;
template class feSysElm_GradSource<2>;
// template class feSysElm_GradSource<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mass matrix for vector-valued field
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_VectorMass<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _udotphi_i.resize(_nFunctions);
  _phi_idotphi_j.resize(_nFunctions*_nFunctions);
}

template <int dim>
void feSysElm_VectorMass<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    uSpace->dotProductShapeShape(k, _phi_idotphi_j);

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        form->_Ae[i][j] += coeff * _phi_idotphi_j[i*_nFunctions+j] * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_VectorMass<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    uSpace->dotProductShapeOther(k, _u, _udotphi_i);

    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= coeff * _udotphi_i[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_VectorMass<0>;
// template class feSysElm_VectorMass<1>;
template class feSysElm_VectorMass<2>;
// template class feSysElm_VectorMass<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed mass matrix for two vector-valued fields
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_MixedVectorMass<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nComponentsU = space[_idU]->getNumComponents();
  _u.resize(_nComponentsU);
  _udotphi_i.resize(_nFunctionsV);
  _phi_idotphi_j.resize(_nFunctionsV*_nFunctionsU);
}

template <int dim>
void feSysElm_MixedVectorMass<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *vSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idV]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    vSpace->dotProductShapeShapeOtherSpace(k, form->_intSpaces[_idU], _phi_idotphi_j);

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        form->_Ae[i][j] += coeff * _phi_idotphi_j[i*_nFunctionsU+j] * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_MixedVectorMass<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *vSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idV]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponentsU);
    vSpace->dotProductShapeOther(k, _u, _udotphi_i);

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * _udotphi_i[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedVectorMass<0>;
// template class feSysElm_MixedVectorMass<1>;
template class feSysElm_MixedVectorMass<2>;
// template class feSysElm_MixedVectorMass<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed scalar-vector mass matrix
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_MixedScalarVectorMass<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV};
  _fieldsLayoutJ = {_idU};
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponentsU = space[_idU]->getNumComponents();
  _phiV.resize(_nFunctionsV);
  _phiU.resize(_nFunctionsU);
  _u.resize(_nComponentsU);
}

template <int dim>
void feSysElm_MixedScalarVectorMass<dim>::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    for(int i = 0; i < _nFunctionsU; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
    for(int i = 0; i < _nFunctionsV; ++i)
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctionsV; ++i)
      for(int j = 0; j < _nFunctionsU; ++j)
        form->_Ae[i][j] += coeff * _phiV[i] * _phiU[j] * jac * _wQuad[k];
  }
}

template <int dim>
void feSysElm_MixedScalarVectorMass<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u,
                                                             _nComponentsU);

    for(int i = 0; i < _nFunctionsV; ++i) {
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * _phiV[i] * _u[i % _nComponentsU] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedScalarVectorMass<0>;
// template class feSysElm_MixedScalarVectorMass<1>;
template class feSysElm_MixedScalarVectorMass<2>;
// template class feSysElm_MixedScalarVectorMass<3>;

// -----------------------------------------------------------------------------
// Bilinear form: transient mass matrix for vector-valued field

template <int dim>// -----------------------------------------------------------------------------
void feSysElm_TransientVectorMass<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _dudt.resize(_nComponents);
  _dudtdotphi_i.resize(_nFunctions);
  _phi_idotphi_j.resize(_nFunctions*_nFunctions);
}

template <int dim>
void feSysElm_TransientVectorMass<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    uSpace->dotProductShapeShape(k, _phi_idotphi_j);
    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        form->_Ae[i][j] += coeff * form->_c0 * _phi_idotphi_j[i*_nFunctions+j] * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_TransientVectorMass<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt, _nComponents);
    uSpace->dotProductShapeOther(k, _dudt, _dudtdotphi_i);
    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= coeff * _dudtdotphi_i[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_TransientVectorMass<0>;
// template class feSysElm_TransientVectorMass<1>;
template class feSysElm_TransientVectorMass<2>;
// template class feSysElm_TransientVectorMass<3>;

// -----------------------------------------------------------------------------
// Bilinear form: diffusion of vector-valued field

template <int dim>// -----------------------------------------------------------------------------
void feSysElm_VectorDiffusion<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradu.resize(_nComponents * _nComponents);
  _gradPhi.resize(_nComponents * _nComponents * _nFunctions);
  _doubleContraction.resize(_nFunctions*_nFunctions);
  _doubleContraction_u.resize(_nFunctions);
}

template <int dim>
void feSysElm_VectorDiffusion<dim>::computeAe(feBilinearForm *form)
{ 
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff, diffusivity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar diffusivity
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    diffusivity = (*_diffusivity)(form->_args);
    coeff = (*_coeff)(form->_args);

    // Get grad phi
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,  _gradPhi);
    uSpace->doubleContractionGradShapeGradShape(_gradPhi, _doubleContraction);

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        form->_Ae[i][j] += coeff * diffusivity * _doubleContraction[i*_nFunctions+j] * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_VectorDiffusion<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff, diffusivity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar diffusivity k(t,x)
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    diffusivity = (*_diffusivity)(form->_args);
    coeff = (*_coeff)(form->_args);

    // Compute gradient of current solution
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());

    // Get gradient of shape functions
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhi);
    uSpace->doubleContractionGradShapeOther(_gradPhi, _gradu, _doubleContraction_u);

    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= coeff * diffusivity * _doubleContraction_u[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_VectorDiffusion<0>;
// template class feSysElm_VectorDiffusion<1>;
template class feSysElm_VectorDiffusion<2>;
// template class feSysElm_VectorDiffusion<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed weak gradient

template <int dim>// -----------------------------------------------------------------------------
void feSysElm_MixedGradient<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Vector-valued test functions
  _idU = 1; // Scalar unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nComponents = space[_idV]->getNumComponents();
  _phiU.resize(_nFunctionsU);
  _gradPhiV.resize(_nComponents * _nComponents * _nFunctionsV);
  _divPhiV.resize(_nFunctionsV);
}

template <int dim>
void feSysElm_MixedGradient<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *vSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idV]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    // Get phiU and gradPhiV
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    vSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV);
    vSpace->divergence(_gradPhiV, _divPhiV);

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        // div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
        form->_Ae[i][j] -= coeff * _phiU[j] * _divPhiV[i] * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_MixedGradient<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *vSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idV]);

  double jac, coeff, u;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    // Get u and gradPhiV
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    vSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV);
    vSpace->divergence(_gradPhiV, _divPhiV);

    for(int i = 0; i < _nFunctionsV; ++i) {
      // div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
      form->_Be[i] += coeff * u * _divPhiV[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedGradient<0>;
// template class feSysElm_MixedGradient<1>;
template class feSysElm_MixedGradient<2>;
// template class feSysElm_MixedGradient<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed weak gradient with coefficient depending on solution

template <int dim>// -----------------------------------------------------------------------------
void feSysElm_MixedGradientFieldDependentCoeff<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Vector-valued test functions
  _idP = 1; // Scalar unknown field
  _idW = 1; // Scalar unknown field on which the coefficient depends
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idV, _idP, _idW};
  _nComponents = space[_idV]->getNumComponents();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  // _phiU.resize(_nFunctionsU);
  _gradPhiV.resize(_nComponents * _nFunctionsV);
}

template <int dim>
void feSysElm_MixedGradientFieldDependentCoeff<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "feSysElm_MixedGradientFieldDependentCoeff : Solve with FD for now");
  exit(-1);
  // double jac, coeff, div_v;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Evaluate scalar coefficient
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
  //   coeff = (*_coeff)(form->_args);

  //   // Get phiU
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

  //   // Compute compacted grad(phiV) without the trivial zeros
  //   form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhiV.data());

  //   for(int i = 0; i < _nFunctionsV; ++i) {
  //     for(int j = 0; j < _nFunctionsU; ++j) {
  //       div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
  //       form->_Ae[i][j] -= coeff * _phiU[j] * div_v * jac * _wQuad[k];
  //     }
  //   }
  // }
}

template <int dim>
void feSysElm_MixedGradientFieldDependentCoeff<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, p, w, div_v;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Get current scalar solutions
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);
    w = form->_intSpaces[_idW]->interpolateFieldAtQuadNode(form->_sol[_idW], k);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    form->_args.u = w;
    coeff = (*_coeff)(form->_args);

    // Compute compacted grad(phiV) without the trivial zeros
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
      form->_Be[i] += coeff * p * div_v * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedGradientFieldDependentCoeff<0>;
// template class feSysElm_MixedGradientFieldDependentCoeff<1>;
template class feSysElm_MixedGradientFieldDependentCoeff<2>;
// template class feSysElm_MixedGradientFieldDependentCoeff<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed divergence

template <int dim>// -----------------------------------------------------------------------------
void feSysElm_MixedDivergence<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar test functions
  _idU = 1; // Vector-valued unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _gradu.resize(_nComponents * _nComponents);
  _phiV.resize(_nFunctionsV);
  _gradPhiU.resize(_nComponents * _nComponents * _nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_MixedDivergence<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get _phiV
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);
    // Get _gradPhiU
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);
    uSpace->divergence(_gradPhiU, _divPhiU);

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        // div_phiU = _gradPhiU[j * _nComponents + (j % _nComponents)];
        form->_Ae[i][j] += coeff * _phiV[i] * _divPhiU[j] * jac * _wQuad[k];
      }
    }
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

template <int dim>
void feSysElm_MixedDivergence<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, div_u;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Compute gradient of current solution
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU],
      _nComponents, k, form->_transformation, _gradu.data());

    div_u = 0.;
    for(int m = 0; m < _nComponents; ++m) {
      div_u += _gradu[(m + 1) * (m + 1) - 1];
    }

    // Get phiV
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * div_u * _phiV[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedDivergence<0>;
// template class feSysElm_MixedDivergence<1>;
template class feSysElm_MixedDivergence<2>;
// template class feSysElm_MixedDivergence<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed divergence with density function of phase marker (CHNS)
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_MixedDivergenceCHNS<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idP   = 0; // Pressure test functions
  _idU   = 1; // Vector-valued velocity
  _idPhi = 2; // Scalar phase marker
  _fieldsLayoutI = {_idP}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsP = space[_idP]->getNumFunctions();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradphi.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _phiU.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_MixedDivergenceCHNS<dim>::computeAe(feBilinearForm *form)
{
  double jac, coeff, phi, rho, drhodphi, div_phiU, div_rhophiU, gradPhiDotphiU;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    phi = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    form->_args.u = phi;
    rho = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);

    // _phiP
    form->_intSpaces[_idP]->getFunctionsAtQuadNode(k, _phiP);
    // _phiU
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    // _gradPhiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());
    // gradphi
    form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());

    for(int i = 0; i < _nFunctionsP; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        div_phiU = _gradPhiU[j * _nComponents + (j % _nComponents)];

        gradPhiDotphiU = _gradphi[j % _nComponents] * _phiU[j];

        div_rhophiU = drhodphi * gradPhiDotphiU + rho * div_phiU;

        form->_Ae[i][j] += coeff * _phiP[i] * div_rhophiU * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_MixedDivergenceCHNS<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, phi, drhodphi, rho, div_u, divrhou, gradphidotu;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    phi = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    form->_args.u = phi;
    rho = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);

    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // gradu = [grad(u_1) ... grad(u_nComponents)]
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // gradphi
    form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());
    // phiP
    form->_intSpaces[_idP]->getFunctionsAtQuadNode(k, _phiP);

    div_u = 0.;
    gradphidotu = 0.;
    for(int m = 0; m < _nComponents; ++m) {
      div_u += _gradu[(m + 1) * (m + 1) - 1];
      gradphidotu += _gradphi[m] * _u[m];
    }

    // Compute div(rho*u) = gradrho(phi) dot u + rho(phi) divu
    divrhou = drhodphi * gradphidotu + rho * div_u;

    for(int i = 0; i < _nFunctionsP; ++i) {
      form->_Be[i] -= coeff * divrhou * _phiP[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedDivergenceCHNS<0>;
// template class feSysElm_MixedDivergenceCHNS<1>;
template class feSysElm_MixedDivergenceCHNS<2>;
// template class feSysElm_MixedDivergenceCHNS<3>;

// -----------------------------------------------------------------------------
// SUPG stabilization for CHNS tracer equation
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_CHNS_Tracer_SUPG<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idPhi = 0;
  _idU   = 1;
  _idMu  = 2;
  _fieldsLayoutI = {_idPhi};
  _fieldsLayoutJ = {_idPhi, _idU, _idMu};
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsU   = space[_idU]->getNumFunctions();
  _nFunctionsMu  = space[_idMu]->getNumFunctions();
  // _nComponents = space[_idU]->getNumComponents();
  _u.resize(dim);
  _gradphi.resize(dim);
  _gradPhiphi.resize(dim * _nFunctionsPhi);
  if constexpr(dim == 2) _hessMu.resize(3); // xx xy yy
  if constexpr(dim == 3) _hessMu.resize(6); // xx xy yy xz yz zz
}

template <int dim>
void feSysElm_CHNS_Tracer_SUPG<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing matrix with finite differences");
  exit(-1);
}

template <int dim>
void feSysElm_CHNS_Tracer_SUPG<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, phi, M, R, uDotGradPhi, laplacian_mu, uDotGradPhiphi,
    dphidt, tau;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    phi = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    form->_args.u = phi;
    coeff = (*_coeff)(form->_args);
    M = (*_mobility)(form->_args);

    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, dim);
    // gradphi
    form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());
    // gradPhiphi
    form->_intSpaces[_idPhi]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiphi.data());
    // hessian of mu
    form->_intSpaces[_idMu]->interpolateFieldAtQuadNode_physicalHessian(form->_sol[_idMu], k, form->_transformation, _hessMu.data());

    if constexpr(dim == 2) {
      uDotGradPhi = _u[0] * _gradphi[0] + _u[1] * _gradphi[1];
      laplacian_mu = _hessMu[0] + _hessMu[2];
    }

    if constexpr(dim == 3) {
      uDotGradPhi = _u[0] * _gradphi[0] + _u[1] * _gradphi[1] + _u[2] * _gradphi[2];
      laplacian_mu = _hessMu[0] + _hessMu[2] + _hessMu[5];
    }

    // Residual of the tracer equation
    dphidt = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_solDot[_idPhi], k);
    R = dphidt + uDotGradPhi - M * laplacian_mu;

    // Which "viscosity" coefficient to take? Here the mobility M?
    tau = compute_tau<dim>(form->_geoCoord, _u.data(), M, _nFunctionsPhi, _gradPhiphi, jac);

    for(int i = 0; i < _nFunctionsPhi; ++i) {

      // Dot product u dot grad(shape_i)
      if constexpr(dim == 2) {
        uDotGradPhiphi = _u[0] * _gradPhiphi[i * dim]
                       + _u[1] * _gradPhiphi[i * dim + 1];
      }

      if constexpr(dim == 3) {
        uDotGradPhiphi = _u[0] * _gradPhiphi[i * dim]
                       + _u[1] * _gradPhiphi[i * dim + 1]
                       + _u[2] * _gradPhiphi[i * dim + 2];
      }

      form->_Be[i] -= coeff * R * tau * uDotGradPhiphi * jac * _wQuad[k];
    }
  }
}

template class feSysElm_CHNS_Tracer_SUPG<2>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed curl
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_MixedCurl<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar test functions
  _idU = 1; // Vector-valued unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = 2;
  _nFunctionsV = space[_idV]->getNumFunctions();
  _gradu.resize(_nComponents * _nComponents);
  _phiV.resize(_nFunctionsV);
}

template <int dim>
void feSysElm_MixedCurl<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

template <int dim>
void feSysElm_MixedCurl<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, curl_u;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Compute phi and gradu
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());

    // Only for 2D fields
    curl_u = _gradu[2] - _gradu[1];

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * curl_u * _phiV[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedCurl<0>;
// template class feSysElm_MixedCurl<1>;
template class feSysElm_MixedCurl<2>;
// template class feSysElm_MixedCurl<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed dot product
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_MixedDotProduct<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar test functions
  _idU = 1; // Vector-valued unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = 2;
  _u.resize(_nComponents);
  _fVal.resize(_nComponents);
  _phiV.resize(_nFunctionsV);
  _phiU.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_MixedDotProduct<dim>::computeAe(feBilinearForm *form)
{
  double jac, phiudotf;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate user-defined vector field
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    (*_f)(form->_args, _fVal);

    // _phiV
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);
    // _phiU
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {

        phiudotf = _phiU[j] * _fVal[j % _nComponents];

        form->_Ae[i][j] += _phiV[i] * phiudotf * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_MixedDotProduct<dim>::computeBe(feBilinearForm *form)
{
  double jac, udotf;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate user-defined vector field
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    (*_f)(form->_args, _fVal);

    // Get u and phiV
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    udotf = 0.;
    for(int iComp = 0; iComp < _nComponents; ++iComp) {
      udotf += _u[iComp] * _fVal[iComp];
    }

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= udotf * _phiV[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_MixedDotProduct<0>;
// template class feSysElm_MixedDotProduct<1>;
template class feSysElm_MixedDotProduct<2>;
// template class feSysElm_MixedDotProduct<3>;

// -----------------------------------------------------------------------------
// Bilinear form: weak divergence
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_ScalarVectorProduct<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar unknown, associated to test functions
  _idU = 1; // Vector unknown
  _fieldsLayoutI = {_idV};
  _fieldsLayoutJ = {_idV, _idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _u.resize(dim);
  _phiV.resize(_nFunctionsV);
  _gradPhiV.resize(_nFunctionsV * dim);
  _phiUdotGradPhiV.resize(_nFunctionsV * _nFunctionsU);
}

template <int dim>
void feSysElm_ScalarVectorProduct<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, jacw, coeff, v, uDotGradPhiV;
  for(int k = 0; k < _nQuad; ++k)
  {
    jac = form->_J[_nQuad * form->_numElem + k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate fields and coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get all relevant values and gradients
    // u
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, dim);
    // v
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    // phiV
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);
    // grad_phiV
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiV.data());

    uSpace->dotProductShapeGradShapeOtherSpaceTranspose(k, _nFunctionsV, _gradPhiV, _phiUdotGradPhiV);
    
    // Increment local matrix
    for(int i = 0; i < _nFunctionsV; ++i)
    {
      if constexpr(dim == 2) {
        uDotGradPhiV = _u[0] * _gradPhiV[i * dim + 0] + _u[1] *  _gradPhiV[i * dim + 1];
      }

      int J = 0;
      // phiV - V block : contribution of (delta v * u) cdot phiV_i
      for(int j = 0; j < _nFunctionsV; ++j, ++J)
      {
        form->_Ae[i][J] += coeff * (_phiV[j] * uDotGradPhiV) * jacw;
      }

      // phiV - U block : contribution of (v * delta u) cdot phiV_i
      for(int j = 0; j < _nFunctionsU; ++j, ++J)
      {
        form->_Ae[i][J] += coeff * (v * _phiUdotGradPhiV[i*_nFunctionsU+j]) * jacw;
      }
    }
  }
}

template <int dim>
void feSysElm_ScalarVectorProduct<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, v, uDotGradPhiV;
  
  for(int k = 0; k < _nQuad; ++k)
  {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate fields and coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Scalar v
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    // Vector u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, dim);

    // Gradient of shape functions of v
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i)
    {
      if constexpr(dim == 2) {
        uDotGradPhiV = _u[0] * _gradPhiV[i * dim + 0] + _u[1] *  _gradPhiV[i * dim + 1];
      }

      form->_Be[i] -= coeff * v * uDotGradPhiV * jac * _wQuad[k];
    }
  }
}

template class feSysElm_ScalarVectorProduct<2>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed weak gradient involving 2 scalar fields and the test
//                functions of a 3rd vector field
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_TripleMixedGradient<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idW = 0; // Vector-valued test functions
  _idU = 1; // Scalar unknown field
  _idV = 2; // Scalar unknown field
  _fieldsLayoutI = {_idW}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU, _idV};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nFunctionsW = space[_idW]->getNumFunctions();
  _nComponents = space[_idW]->getNumComponents();
  _gradV.resize(_nComponents); // Space dimension = nComponents of vector field
  _phiW.resize(_nFunctionsW);
}

template <int dim>
void feSysElm_TripleMixedGradient<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement local matrix for TripleMixedGradient");
  exit(-1);
}

template <int dim>
void feSysElm_TripleMixedGradient<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, u, gradVDotPhiW;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get u, gradV and phiW
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_intSpaces[_idV]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idV], k, form->_transformation, _gradV.data());
    form->_intSpaces[_idW]->getFunctionsAtQuadNode(k, _phiW);

    for(int i = 0; i < _nFunctionsW; ++i)
    {
      // _phiW contains the basis functions for each component of the vector field.
      // Take the i-th vector of _phiW = [phiW_x,1 phiW_y,1 phiW_z,1 ... phiW_x,N phiW_y,N phiW_z,N]
      gradVDotPhiW = dotProdN(_nComponents, _gradV.data(), _phiW.data() + i*_nComponents);
      form->_Be[i] -= coeff * u * gradVDotPhiW * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_TripleMixedGradient<0>;
// template class feSysElm_TripleMixedGradient<1>;
template class feSysElm_TripleMixedGradient<2>;
// template class feSysElm_TripleMixedGradient<3>;

// -----------------------------------------------------------------------------
// Bilinear form: CHNS Momentum equation

template <int dim>// -----------------------------------------------------------------------------
void feSysElm_CHNS_Momentum<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;
  _fieldsLayoutI = {_idU};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu = space[_idMu]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _f.resize(_nComponents);
  _u.resize(_nComponents);
  _dudt.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _uDotGradu.resize(_nComponents);
  _gradphi.resize(_nComponents);
  _gradmu.resize(_nComponents);
  _gradmuDotGradu.resize(_nComponents);
  _symmetricGradu.resize(_nComponents * _nComponents);

  // phiU is the test function of U, unrelated to the phase marker Phi
  _phiU.resize(_nFunctionsU);
  _gradPhiU.resize(_nComponents * _nComponents * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(_nComponents * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(_nComponents * _nFunctionsMu);

  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _gradPhiDotphiU.resize(_nFunctionsU);
  _gradMuDotgradUdotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);

  _phi_idotphi_j.resize(_nFunctionsU*_nFunctionsU);
  _u0DotGradPhiUDotPhiU.resize(_nFunctionsU*_nFunctionsU);
  _phiUDotGradu0DotPhiU.resize(_nFunctionsU*_nFunctionsU);
  _doubleContractionPhiPhi.resize(_nFunctionsU*_nFunctionsU);
  _doubleContractionPhiPhiT.resize(_nFunctionsU*_nFunctionsU);
  _gradMu0DotgradUdotphiU.resize(_nFunctionsU*_nFunctionsU);
  _gradMuDotgradU0DotphiU.resize(_nFunctionsU*_nFunctionsMu);
  _gradPhi0dotPhiU.resize(_nFunctionsU);
  _gradPhiPhiDotPhiU.resize(_nFunctionsU*_nFunctionsPhi);
  _symGraduDDotGradPhiU.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_CHNS_Momentum<dim>::computeAe(feBilinearForm *form)
{
  // feErrorMsg(FE_STATUS_ERROR, "Escaping");
  // exit(-1);

  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);
  const feSpace *pSpace = form->_intSpaces[_idP];
  const feSpace *fSpace = form->_intSpaces[_idPhi];
  const feSpace *mSpace = form->_intSpaces[_idMu];
  const double *Jptr = form->_J.data() + _nQuad * form->_numElem;

  double jac, jacw, rho0, mobility, drhodphi, coeffKorteweg, eta0, detadphi, phi0, mu;
  for(int k = 0; k < _nQuad; ++k)
  {
    // jac = form->_J[_nQuad * form->_numElem + k];
    jac = Jptr[k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    phi0 = fSpace->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    form->_args.u = phi0;
    rho0 = (*_density)(form->_args);
    mobility = (*_mobility)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta0 = (*_viscosity)(form->_args);
    detadphi = (*_dviscdphi)(form->_args);
    coeffKorteweg = (*_coeffKorteweg)(form->_args);
    (*_volumeForce)(form->_args, _f);

    // Get all relevant values and gradients
    // u
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // dudt
    uSpace->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt, _nComponents);
    // grad_u
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_phiU
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);

    // phiP
    pSpace->getFunctionsAtQuadNode(k, _phiP);

    // grad_phi
    fSpace->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());
    // phiPhi
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    // grad_phiPhi
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiPhi.data());

    // mu
    mu = mSpace->interpolateFieldAtQuadNode(form->_sol[_idMu], k);
    // grad_mu
    mSpace->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idMu], k, form->_transformation, _gradmu.data());
    // phiMu
    mSpace->getFunctionsAtQuadNode(k, _phiMu);
    // grad_phiMu
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiMu.data());

    // Compute u0 dot gradu0
    // std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    // for(int m = 0; m < _nComponents; ++m) {
    //   for(int n = 0; n < _nComponents; ++n) {
    //     // Contract _u and _gradmu with the lines of _gradu
    //     _uDotGradu[m] += _u[n] * _gradu[n * _nComponents + m];
    //     _symmetricGradu[m * _nComponents + n] =
    //       _gradu[m * _nComponents + n] + _gradu[n * _nComponents + m];
    //   }
    // }

    if constexpr(dim == 2) {
      _uDotGradu[0] = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1] = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];
    }

    if constexpr(dim == 3) {
      feErrorMsg(FE_STATUS_ERROR, "Add contractions in 3D");
      exit(-1);
    } 

    // Compute all the required dot products and contractions:

    // Acceleration, convective acceleration and source term
    uSpace->dotProductShapeShape(k, _phi_idotphi_j);
    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    // [(u0 dot gradu) + (u dot gradu0)] cdot phiU
    uSpace->vectorDotGradShapeDotShape(k, _gradPhiU, _u, _u0DotGradPhiUDotPhiU);
    uSpace->shapeDotTensorDotShape(k, _gradu, _phiUDotGradu0DotPhiU);
    // Diffusive flux
    // [(gradmu0 dot gradu) + (gradmu dot gradu0)] cdot phiU
    uSpace->vectorDotGradShapeDotShape(k, _gradPhiU, _gradmu, _gradMu0DotgradUdotphiU);
    uSpace->gradOtherScalarShapeDotTensorDotShape(k, _nFunctionsMu, _gradPhiMu, _gradu, _gradMuDotgradU0DotphiU);
    // Korteweg force
    uSpace->dotProductShapeOther(k, _gradphi, _gradPhi0dotPhiU);
    uSpace->dotProductShapeGradShapeOtherSpace(k, _nFunctionsPhi, _gradPhiPhi, _gradPhiPhiDotPhiU);
    // Viscous term
    uSpace->doubleContractionGradShapeGradShape(_gradPhiU, _doubleContractionPhiPhi);
    uSpace->doubleContractionGradShapeGradShapeTransposed(_gradPhiU, _doubleContractionPhiPhiT);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _symGraduDDotGradPhiU);
    // Pressure gradient
    uSpace->divergence(_gradPhiU, _divPhiU);

    UNUSED(mobility);
    UNUSED(drhodphi);
    
    // Increment local matrix
    for(int i = 0; i < _nFunctionsU; ++i)
    {
      int J = 0;
      // phiU - u block
      for(int j = 0; j < _nFunctionsU; ++j, ++J)
      {
        form->_Ae[i][J] += (
          // Acceleration
          rho0 * (form->_c0 * _phi_idotphi_j[i*_nFunctionsU+j]
                     + _u0DotGradPhiUDotPhiU[i*_nFunctionsU+j]
                     + _phiUDotGradu0DotPhiU[i*_nFunctionsU+j])
          // Diffusive flux
          + mobility * drhodphi * _gradMu0DotgradUdotphiU[i*_nFunctionsU+j]
          // Viscous term
          + eta0 * (_doubleContractionPhiPhi [i*_nFunctionsU+j]
                  + _doubleContractionPhiPhiT[i*_nFunctionsU+j])) * jacw;
      }

      // phiU - p block
      for(int j = 0; j < _nFunctionsP; ++j, ++J)
      {
        // Pressure gradient
        form->_Ae[i][J] += (- _phiP[j]) * _divPhiU[i] * jacw;
      }

      // phiU - phi block
      for(int j = 0; j < _nFunctionsPhi; ++j, ++J)
      {
        form->_Ae[i][J] += (
          // Acceleration and source term
          drhodphi * _phiPhi[j] * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])
          // Korteweg force
          + coeffKorteweg * mu * _gradPhiPhiDotPhiU[i*_nFunctionsPhi+j]
          // Viscous term
          + detadphi * _phiPhi[j] * _symGraduDDotGradPhiU[i]
          ) * jacw;
      }

      // phiU - mu block
      for(int j = 0; j < _nFunctionsMu; ++j, ++J)
      {
        form->_Ae[i][J] += (
          // Diffusive flux
          + mobility * drhodphi * _gradMuDotgradU0DotphiU[i*_nFunctionsMu+j]
          // Korteweg force
          + coeffKorteweg * _phiMu[j] * _gradPhi0dotPhiU[i]) * jacw;
      }
    }
  }
}

template <int dim>
void feSysElm_CHNS_Momentum<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, rho, mobility, drhodphi, coeffKorteweg, eta, p, phi, mu;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    phi = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    
    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    rho = (*_density)(form->_args);
    mobility = (*_mobility)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta = (*_viscosity)(form->_args);
    coeffKorteweg = (*_coeffKorteweg)(form->_args);
    (*_volumeForce)(form->_args, _f);
    
    // Get all relevant values and gradients
    // u
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // dudt
    uSpace->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt, _nComponents);
    // grad_u
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_phiU
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);

    // p
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);

    // grad_phi
    form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());

    // mu
    mu = form->_intSpaces[_idMu]->interpolateFieldAtQuadNode(form->_sol[_idMu], k);
    // grad_mu
    form->_intSpaces[_idMu]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idMu], k, form->_transformation, _gradmu.data());

    // Compute u dot gradu
    // Compute gradmu dot gradu
    // Compute symmetric gradient (2*d(u))
    // std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    // std::fill(_gradmuDotGradu.begin(), _gradmuDotGradu.end(), 0.);
    // for(int m = 0; m < _nComponents; ++m) {
    //   for(int n = 0; n < _nComponents; ++n) {
    //     // Contract _u and _gradmu with the lines of _gradu
    //     _uDotGradu[m] += _u[n] * _gradu[n * _nComponents + m];
    //     _gradmuDotGradu[m] += _gradmu[n] * _gradu[n * _nComponents + m];
    //     _symmetricGradu[m * _nComponents + n] =
    //       _gradu[m * _nComponents + n] + _gradu[n * _nComponents + m];
    //   }
    // }

    if constexpr(dim == 2) {
      _uDotGradu[0] = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1] = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _gradmuDotGradu[0] = _gradmu[0] * _gradu[0] + _gradmu[1] * _gradu[2];
      _gradmuDotGradu[1] = _gradmu[0] * _gradu[1] + _gradmu[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];
    }

    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _gradmuDotGradu, _gradMuDotgradUdotphiU);
    uSpace->dotProductShapeOther(k, _gradphi, _gradPhiDotphiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _doubleContraction);

    UNUSED(mobility);
    UNUSED(drhodphi);

    // Increment RHS
    for(int i = 0; i < _nFunctionsU; ++i)
    {
      form->_Be[i] -= jac * _wQuad[k] * (
        // Acceleration
        rho * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])
        // Diffusive flux
        + mobility * drhodphi * _gradMuDotgradUdotphiU[i]
        // Korteweg force
        + coeffKorteweg * mu * _gradPhiDotphiU[i]
        // Pressure gradient
        - p * _divPhiU[i]
        // Viscous term
        + eta * _doubleContraction[i]);
    }
  }
}

// template class feSysElm_CHNS_Momentum<0>;
// template class feSysElm_CHNS_Momentum<1>;
template class feSysElm_CHNS_Momentum<2>;
// template class feSysElm_CHNS_Momentum<3>;

// // -----------------------------------------------------------------------------
// // Bilinear form: Alternative CHNS Momentum equation
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_CHNS_Momentum_Alternative<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;
  _fieldsLayoutI = {_idU};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu = space[_idMu]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();

  _f.resize(_nComponents);
  _u.resize(_nComponents);
  _dudt.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _uDotGradu.resize(_nComponents);
  _gradphi.resize(_nComponents);
  _gradmu.resize(_nComponents);
  _gradmuDotGradu.resize(_nComponents);
  _symmetricGradu.resize(_nComponents * _nComponents);

  // phiU is the test function of U, unrelated to the phase marker Phi
  _phiU.resize(_nFunctionsU);
  _gradPhiU.resize(_nComponents * _nComponents * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(_nComponents * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(_nComponents * _nFunctionsMu);

  _uDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _gradPhiDotphiU.resize(_nFunctionsU);
  _gradMuDotphiU.resize(_nFunctionsU);
  _gradMuDotgradUdotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);

  _phi_idotphi_j.resize(_nFunctionsU*_nFunctionsU);
  _u0DotGradPhiUDotPhiU.resize(_nFunctionsU*_nFunctionsU);
  _phiUDotGradu0DotPhiU.resize(_nFunctionsU*_nFunctionsU);
  _doubleContractionPhiPhi.resize(_nFunctionsU*_nFunctionsU);
  _doubleContractionPhiPhiT.resize(_nFunctionsU*_nFunctionsU);
  _gradMu0DotgradUdotphiU.resize(_nFunctionsU*_nFunctionsU);
  _gradMuDotgradU0DotphiU.resize(_nFunctionsU*_nFunctionsMu);
  _gradPhi0dotPhiU.resize(_nFunctionsU);
  _gradPhiPhiDotPhiU.resize(_nFunctionsU*_nFunctionsPhi);
  _symGraduDDotGradPhiU.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_CHNS_Momentum_Alternative<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
  // const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);
  // const feSpace *pSpace = form->_intSpaces[_idP];
  // const feSpace *fSpace = form->_intSpaces[_idPhi];
  // const feSpace *mSpace = form->_intSpaces[_idMu];
  // const double *Jptr = form->_J.data() + _nQuad * form->_numElem;

  // double jac, jacw, rho0, mobility, drhodphi, coeffKorteweg, eta0, detadphi, phi0, mu;
  // for(int k = 0; k < _nQuad; ++k)
  // {
  //   // jac = form->_J[_nQuad * form->_numElem + k];
  //   jac = Jptr[k];
  //   jacw = jac * _wQuad[k];
  //   form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Evaluate scalar coefficients
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
  //   phi0 = fSpace->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
  //   form->_args.u = phi0;
  //   rho0 = (*_density)(form->_args);
  //   mobility = (*_mobility)(form->_args);
  //   drhodphi = (*_drhodphi)(form->_args);
  //   eta0 = (*_viscosity)(form->_args);
  //   detadphi = (*_dviscdphi)(form->_args);
  //   coeffKorteweg = (*_coeffKorteweg)(form->_args);
  //   (*_volumeForce)(form->_args, _f);

  //   // Get all relevant values and gradients
  //   // u
  //   uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
  //   // dudt
  //   uSpace->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt, _nComponents);
  //   // grad_u
  //   uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
  //   // grad_phiU
  //   uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);

  //   // phiP
  //   pSpace->getFunctionsAtQuadNode(k, _phiP);

  //   // grad_phi
  //   fSpace->interpolateFieldAtQuadNode_physicalGradient(
  //     form->_sol[_idPhi], k, form->_transformation, _gradphi.data());
  //   // phiPhi
  //   fSpace->getFunctionsAtQuadNode(k, _phiPhi);
  //   // grad_phiPhi
  //   fSpace->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhiPhi.data());

  //   // mu
  //   mu = mSpace->interpolateFieldAtQuadNode(form->_sol[_idMu], k);
  //   // grad_mu
  //   mSpace->interpolateFieldAtQuadNode_physicalGradient(
  //     form->_sol[_idMu], k, form->_transformation, _gradmu.data());
  //   // phiMu
  //   mSpace->getFunctionsAtQuadNode(k, _phiMu);
  //   // grad_phiMu
  //   mSpace->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhiMu.data());

  //   // Compute u0 dot gradu0
  //   // std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
  //   // for(int m = 0; m < _nComponents; ++m) {
  //   //   for(int n = 0; n < _nComponents; ++n) {
  //   //     // Contract _u and _gradmu with the lines of _gradu
  //   //     _uDotGradu[m] += _u[n] * _gradu[n * _nComponents + m];
  //   //     _symmetricGradu[m * _nComponents + n] =
  //   //       _gradu[m * _nComponents + n] + _gradu[n * _nComponents + m];
  //   //   }
  //   // }

  //   if constexpr(dim == 2) {
  //     _uDotGradu[0] = _u[0] * _gradu[0] + _u[1] * _gradu[2];
  //     _uDotGradu[1] = _u[0] * _gradu[1] + _u[1] * _gradu[3];
  //     _symmetricGradu[0] = _gradu[0] + _gradu[0];
  //     _symmetricGradu[1] = _gradu[1] + _gradu[2];
  //     _symmetricGradu[2] = _gradu[2] + _gradu[1];
  //     _symmetricGradu[3] = _gradu[3] + _gradu[3];
  //   }

  //   if constexpr(dim == 3) {
  //     feErrorMsg(FE_STATUS_ERROR, "Add contractions in 3D");
  //     exit(-1);
  //   } 

  //   // Compute all the required dot products and contractions:

  //   // Acceleration, convective acceleration and source term
  //   uSpace->dotProductShapeShape(k, _phi_idotphi_j);
  //   uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
  //   uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
  //   uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
  //   // [(u0 dot gradu) + (u dot gradu0)] cdot phiU
  //   uSpace->vectorDotGradShapeDotShape(k, _gradPhiU, _u, _u0DotGradPhiUDotPhiU);
  //   uSpace->shapeDotTensorDotShape(k, _gradu, _phiUDotGradu0DotPhiU);
  //   // Diffusive flux
  //   // [(gradmu0 dot gradu) + (gradmu dot gradu0)] cdot phiU
  //   uSpace->vectorDotGradShapeDotShape(k, _gradPhiU, _gradmu, _gradMu0DotgradUdotphiU);
  //   uSpace->gradOtherScalarShapeDotTensorDotShape(k, _nFunctionsMu, _gradPhiMu, _gradu, _gradMuDotgradU0DotphiU);
  //   // Korteweg force
  //   uSpace->dotProductShapeOther(k, _gradphi, _gradPhi0dotPhiU);
  //   uSpace->dotProductShapeGradShapeOtherSpace(k, _nFunctionsPhi, _gradPhiPhi, _gradPhiPhiDotPhiU);
  //   // Viscous term
  //   uSpace->doubleContractionGradShapeGradShape(_gradPhiU, _doubleContractionPhiPhi);
  //   uSpace->doubleContractionGradShapeGradShapeTransposed(_gradPhiU, _doubleContractionPhiPhiT);
  //   uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _symGraduDDotGradPhiU);
  //   // Pressure gradient
  //   uSpace->divergence(_gradPhiU, _divPhiU);

  //   UNUSED(mobility);
  //   UNUSED(drhodphi);
    
  //   // Increment local matrix
  //   for(int i = 0; i < _nFunctionsU; ++i)
  //   {
  //     int J = 0;
  //     // phiU - u block
  //     for(int j = 0; j < _nFunctionsU; ++j, ++J)
  //     {
  //       form->_Ae[i][J] += (
  //         // Acceleration
  //         rho0 * (form->_c0 * _phi_idotphi_j[i*_nFunctionsU+j]
  //                    + _u0DotGradPhiUDotPhiU[i*_nFunctionsU+j]
  //                    + _phiUDotGradu0DotPhiU[i*_nFunctionsU+j])
  //         // Diffusive flux
  //         // + mobility * drhodphi * _gradMu0DotgradUdotphiU[i*_nFunctionsU+j]
  //         // Viscous term
  //         + eta0 * (_doubleContractionPhiPhi [i*_nFunctionsU+j]
  //                 + _doubleContractionPhiPhiT[i*_nFunctionsU+j])) * jacw;
  //     }

  //     // phiU - p block
  //     for(int j = 0; j < _nFunctionsP; ++j, ++J)
  //     {
  //       // Pressure gradient
  //       form->_Ae[i][J] += (- _phiP[j]) * _divPhiU[i] * jacw;
  //     }

  //     // phiU - phi block
  //     for(int j = 0; j < _nFunctionsPhi; ++j, ++J)
  //     {
  //       form->_Ae[i][J] += (
  //         // Acceleration and source term
  //         drhodphi * _phiPhi[j] * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])
  //         // Korteweg force
  //         + coeffKorteweg * mu * _gradPhiPhiDotPhiU[i*_nFunctionsPhi+j]
  //         // Viscous term
  //         + detadphi * _phiPhi[j] * _symGraduDDotGradPhiU[i]
  //         ) * jacw;
  //     }

  //     // phiU - mu block
  //     for(int j = 0; j < _nFunctionsMu; ++j, ++J)
  //     {
  //       form->_Ae[i][J] += (
  //         // Diffusive flux
  //         // mobility * drhodphi * _gradMuDotgradU0DotphiU[i*_nFunctionsMu+j]
  //         // Korteweg force
  //         + coeffKorteweg * _phiMu[j] * _gradPhi0dotPhiU[i]) * jacw;
  //     }
  //   }
  // }
}

template <int dim>
void feSysElm_CHNS_Momentum_Alternative<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, rho, drhodphi, eta, p, phi, dphidt, divu, divRhoU, uDotGradRho;

  for(int k = 0; k < _nQuad; ++k)
  {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    phi = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    
    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    rho = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta = (*_viscosity)(form->_args);
    (*_volumeForce)(form->_args, _f);
    
    // Get all relevant values and gradients
    // u
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // dudt
    uSpace->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt, _nComponents);
    // grad_u
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_phiU
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);

    // p
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);

    // dphidt
    dphidt = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_solDot[_idPhi], k);
    // grad_phi
    form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());

    // grad_mu
    form->_intSpaces[_idMu]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idMu], k, form->_transformation, _gradmu.data());

    if constexpr(dim == 2) {
      _uDotGradu[0] = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1] = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];

      divu = _gradu[0] + _gradu[3];
      uDotGradRho = drhodphi * (_u[0] * _gradphi[0] + _u[1] * _gradphi[1]);
    }

    divRhoU = rho * divu + uDotGradRho;

    uSpace->dotProductShapeOther(k, _u, _uDotPhiU);
    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _gradmu, _gradMuDotphiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _doubleContraction);

    // Increment RHS
    for(int i = 0; i < _nFunctionsU; ++i)
    {
      form->_Be[i] -= jac * _wQuad[k] * (
        
        // Acceleration
        rho * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])
        
        // Mass conservation
        + 0.5 * _uDotPhiU[i] * (drhodphi * dphidt + divRhoU)

        // Korteweg force
        + phi * _gradMuDotphiU[i]
        
        // Pressure gradient
        - p * _divPhiU[i]
        
        // Viscous term : Symmetric gradient and compressible contribution
        + eta * (_doubleContraction[i]
         - (2. / (double) dim) * divu * _divPhiU[i]
         ));
    }
  }
}

template class feSysElm_CHNS_Momentum_Alternative<2>;

// -----------------------------------------------------------------------------
// Bilinear form: CHNS equations for volume-averaged velocity
// -----------------------------------------------------------------------------
template <int dim>
void CHNS_VolumeAveraged<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;

  _fieldsLayoutI = {_idU, _idP, _idPhi, _idMu};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};

  _nFunctionsU   = space[_idU]->getNumFunctions();
  _nFunctionsP   = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu  = space[_idMu]->getNumFunctions();

  // Volume force and source term
  _f.resize(dim);
  _Su.resize(dim);

  // Current fields and gradients
  _u.resize(dim);
  _dudt.resize(dim);
  _gradu.resize(dim * dim);
  _symmetricGradu.resize(dim * dim);
  _uDotGradu.resize(dim);
  _gradp.resize(dim);
  _gradphi.resize(dim);
  _gradmu.resize(dim);
  _gradmuDotGradu.resize(dim);

  // Test functions and gradients
  _gradPhiU.resize(dim * dim * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _gradPhiP.resize(dim * _nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(dim * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(dim * _nFunctionsMu);

  // Contractions
  _uDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _gradMuDotgradUdotphiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _SuDotPhiU.resize(_nFunctionsU);
  _gradPhiDotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);
}

template <int dim>
void CHNS_VolumeAveraged<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
}

template <int dim>
void CHNS_VolumeAveraged<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);
  const feSpace *pSpace = form->_intSpaces[_idP];
  const feSpace *fSpace = form->_intSpaces[_idPhi];
  const feSpace *mSpace = form->_intSpaces[_idMu];
  std::vector<double> &uSol    = form->_sol[_idU];
  std::vector<double> &pSol    = form->_sol[_idP];
  std::vector<double> &fSol    = form->_sol[_idPhi];
  std::vector<double> &mSol    = form->_sol[_idMu];
  std::vector<double> &uSolDot = form->_solDot[_idU];
  std::vector<double> &fSolDot = form->_solDot[_idPhi];
  ElementTransformation &T = form->_transformation;

  // const double hsize = hCircumInscr(form->_geoCoord);

  double jac, jacw, rho, eta, p, phi, mu, dphidt, div_u;
  double M, Sp, Sphi, Smu;
  double uDotGradPhi;

  for(int k = 0; k < _nQuad; ++k)
  {
    jac = form->_J[_nQuad * form->_numElem + k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, T);

    //
    // Get current fields and their derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSol, k, _u, dim);
    p   = pSpace->interpolateFieldAtQuadNode(pSol, k);
    phi = fSpace->interpolateFieldAtQuadNode(fSol, k);
    mu  = mSpace->interpolateFieldAtQuadNode(mSol, k);

    //
    // Current time derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSolDot, k, _dudt, dim);
    dphidt = fSpace->interpolateFieldAtQuadNode(fSolDot, k);

    // Gradients
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(uSol, dim, k, T, _gradu.data());
    pSpace->interpolateFieldAtQuadNode_physicalGradient(pSol, k, T, _gradp.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol, k, T, _gradphi.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol, k, T, _gradmu.data());

    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    //
    // Evaluate all scalar coefficients
    //
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    rho      = (*_density)(form->_args);
    eta      = (*_viscosity)(form->_args);
    M        = (*_mobility)(form->_args);
    (*_volumeForce)(form->_args, _f);

    (*_sourceU)(form->_args, _Su);
    Sp   = (*_sourceP)(form->_args);
    Sphi = (*_sourcePhi)(form->_args);
    Smu  = (*_sourceMu)(form->_args);
    
    // 
    // Test functions and gradients
    // Velocity test functions are accessed through the contraction functions
    //
    pSpace->getFunctionsAtQuadNode(k, _phiP);
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    mSpace->getFunctionsAtQuadNode(k, _phiMu);

    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiU);
    pSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiP.data());
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiPhi.data());
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiMu.data());

    if constexpr(dim == 2) {
      _uDotGradu[0] = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1] = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];
      div_u = _gradu[0] + _gradu[3];
      uDotGradPhi = _u[0] * _gradphi[0] + _u[1] * _gradphi[1];
    }

    uSpace->dotProductShapeOther(k, _u, _uDotPhiU);
    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _Su, _SuDotPhiU);
    uSpace->dotProductShapeOther(k, _gradphi, _gradPhiDotphiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _doubleContraction);

    //
    // Increment RHS
    //
    int I = 0;

    //
    // u test functions block : Momentum equation
    //
    for(int i = 0; i < _nFunctionsU; ++i, ++I)
    {
      form->_Be[I] -= jacw * (
        
        // Acceleration
        rho * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])

        // Diffusive flux
        // + drhodphi * M * _gradMuDotgradUdotphiU[i]
        
        // Korteweg force
        - mu * _gradPhiDotphiU[i]
        
        // Pressure gradient
        - p * _divPhiU[i]
        
        // Viscous term : Symmetric gradient
        + eta * _doubleContraction[i]

        // Source term if needed
        + _SuDotPhiU[i]

        );
    }

    //
    // p test functions block : continuity
    //
    for(int i = 0; i < _nFunctionsP; ++i, ++I)
    {
      form->_Be[I] -= jacw * (div_u * _phiP[i] + Sp * _phiP[i]);
    }

    //
    // phi test functions block : tracer equation
    //
    double  gradMu_dot_gradPhiPhi;
    for(int i = 0; i < _nFunctionsPhi; ++i, ++I)
    {
      if constexpr(dim == 2) {
        gradMu_dot_gradPhiPhi = _gradmu[0] * _gradPhiPhi[2 * i + 0]
                              + _gradmu[1] * _gradPhiPhi[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
        // Time derivative
        dphidt * _phiPhi[i]
        // Convective term
        + uDotGradPhi * _phiPhi[i]
        // Diffusive term
        + M * gradMu_dot_gradPhiPhi
        // Source term if needed
        + Sphi * _phiPhi[i]
        );
    }

    //
    // mu test functions block : potential equation
    //
    double gradPhi_dot_gradPhiMu;
    for(int i = 0; i < _nFunctionsMu; ++i, ++I)
    {
      if constexpr(dim == 2) {
        gradPhi_dot_gradPhiMu = _gradphi[0] * _gradPhiMu[2 * i + 0]
                              + _gradphi[1] * _gradPhiMu[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
        // Mu term
        mu * _phiMu[i]
        // Phi laplacian
        - _lambda  * gradPhi_dot_gradPhiMu
        // Double well potential
        - _lambda / (_epsilon * _epsilon) * phi * (phi * phi - 1.) * _phiMu[i]
        // Source term if needed
        + Smu * _phiMu[i]
        );
    }
  }
}

template class CHNS_VolumeAveraged<2>;

// // -----------------------------------------------------------------------------
// // Bilinear form: Alternative CHNS equations
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_CHNS_Alternative<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;

  _fieldsLayoutI = {_idU, _idP, _idPhi, _idMu};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};

  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu = space[_idMu]->getNumFunctions();

  _f.resize(dim);
  _Su.resize(dim);

  _u.resize(dim);
  _dudt.resize(dim);
  _gradu.resize(dim * dim);
  _symmetricGradu.resize(dim * dim);
  _uDotGradu.resize(dim);
  _gradp.resize(dim);
  _gradphi.resize(dim);
  _gradmu.resize(dim);
  // _gradmuDotGradu.resize(dim);

  // phiU is the test function of U, unrelated to the phase marker Phi
  // _phiU.resize(_nFunctionsU);
  _gradPhiU.resize(dim * dim * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _gradPhiP.resize(dim * _nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(dim * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(dim * _nFunctionsMu);

  _uDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _SuDotPhiU.resize(_nFunctionsU);
  _gradMuDotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);

  // _phi_idotphi_j.resize(_nFunctionsU*_nFunctionsU);
  // _u0DotGradPhiUDotPhiU.resize(_nFunctionsU*_nFunctionsU);
  // _phiUDotGradu0DotPhiU.resize(_nFunctionsU*_nFunctionsU);
  // _doubleContractionPhiPhi.resize(_nFunctionsU*_nFunctionsU);
  // _doubleContractionPhiPhiT.resize(_nFunctionsU*_nFunctionsU);
  // _gradMu0DotgradUdotphiU.resize(_nFunctionsU*_nFunctionsU);
  // _gradMuDotgradU0DotphiU.resize(_nFunctionsU*_nFunctionsMu);
  // _gradPhi0dotPhiU.resize(_nFunctionsU);
  // _gradPhiPhiDotPhiU.resize(_nFunctionsU*_nFunctionsPhi);
  // _symGraduDDotGradPhiU.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_CHNS_Alternative<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
}

template <int dim>
void feSysElm_CHNS_Alternative<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);
  const feSpace *pSpace = form->_intSpaces[_idP];
  const feSpace *fSpace = form->_intSpaces[_idPhi];
  const feSpace *mSpace = form->_intSpaces[_idMu];
  std::vector<double> &uSol    = form->_sol[_idU];
  std::vector<double> &pSol    = form->_sol[_idP];
  std::vector<double> &fSol    = form->_sol[_idPhi];
  std::vector<double> &mSol    = form->_sol[_idMu];
  std::vector<double> &uSolDot = form->_solDot[_idU];
  std::vector<double> &fSolDot = form->_solDot[_idPhi];
  ElementTransformation &T = form->_transformation;

  double jac, jacw, rho, drhodphi, eta, p, phi, mu, dphidt, div_u, divRhoU, uDotGradRho;
  double M, Sp, Sphi, Smu;

  for(int k = 0; k < _nQuad; ++k)
  {
    jac = form->_J[_nQuad * form->_numElem + k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, T);

    //
    // Get current fields and their derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSol, k, _u, dim);
    p   = pSpace->interpolateFieldAtQuadNode(pSol, k);
    phi = fSpace->interpolateFieldAtQuadNode(fSol, k);
    mu  = mSpace->interpolateFieldAtQuadNode(mSol, k);

    //
    // Current time derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSolDot, k, _dudt, dim);
    dphidt = fSpace->interpolateFieldAtQuadNode(fSolDot, k);

    // Gradients
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(uSol, dim, k, T, _gradu.data());
    pSpace->interpolateFieldAtQuadNode_physicalGradient(pSol, k, T, _gradp.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol, k, T, _gradphi.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol, k, T, _gradmu.data());

    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    //
    // Evaluate all scalar coefficients
    //
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    rho      = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta      = (*_viscosity)(form->_args);
    M        = (*_mobility)(form->_args);
    (*_volumeForce)(form->_args, _f);

    (*_sourceU)(form->_args, _Su);
    Sp   = (*_sourceP)(form->_args);
    Sphi = (*_sourcePhi)(form->_args);
    Smu  = (*_sourceMu)(form->_args);
    
    // 
    // Test functions and gradients
    // Velocity test functions are accessed through the contraction functions
    //
    pSpace->getFunctionsAtQuadNode(k, _phiP);
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    mSpace->getFunctionsAtQuadNode(k, _phiMu);

    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiU);
    pSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiP.data());
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiPhi.data());
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiMu.data());

    if constexpr(dim == 2) {
      _uDotGradu[0] = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1] = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];

      div_u = _gradu[0] + _gradu[3];
      uDotGradRho = drhodphi * (_u[0] * _gradphi[0] + _u[1] * _gradphi[1]);
    }

    divRhoU = rho * div_u + uDotGradRho;

    uSpace->dotProductShapeOther(k, _u, _uDotPhiU);
    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _Su, _SuDotPhiU);
    uSpace->dotProductShapeOther(k, _gradmu, _gradMuDotphiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _doubleContraction);

    //
    // Increment RHS
    //
    int I = 0;

    //
    // u test functions block : Momentum equation
    //
    for(int i = 0; i < _nFunctionsU; ++i, ++I)
    {
      form->_Be[I] -= jacw * (
        
        // Acceleration
        rho * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])
        
        // Mass conservation
        + 0.5 * _uDotPhiU[i] * (drhodphi * dphidt + divRhoU)

        // Korteweg force
        + phi * _gradMuDotphiU[i]
        
        // Pressure gradient
        - p * _divPhiU[i]
        
        // Viscous term : Symmetric gradient and compressible contribution
        + eta * (_doubleContraction[i] - (2. / (double) dim) * div_u * _divPhiU[i])

        // Source term if needed
        + _SuDotPhiU[i]

        );
    }

    //
    // p test functions block : continuity
    //
    double gradMu_dot_gradPhiP, gradP_dot_gradPhiP;
    for(int i = 0; i < _nFunctionsP; ++i, ++I)
    {
      if constexpr(dim == 2) {
        gradMu_dot_gradPhiP = _gradmu[0] * _gradPhiP[2 * i + 0]
                            + _gradmu[1] * _gradPhiP[2 * i + 1];
        gradP_dot_gradPhiP  =  _gradp[0] * _gradPhiP[2 * i + 0]
                            +  _gradp[1] * _gradPhiP[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
        // Div u
        div_u * _phiP[i]
        // Div (M (grad(mu + alpha * p))
        + _alpha * M * (gradMu_dot_gradPhiP + _alpha * gradP_dot_gradPhiP)
        // Source term if needed
        + Sp * _phiP[i]);
    }

    //
    // phi test functions block : tracer equation
    //
    double  gradMu_dot_gradPhiPhi, gradP_dot_gradPhiPhi, uDotGradPhiPhi;
    for(int i = 0; i < _nFunctionsPhi; ++i, ++I)
    {
      if constexpr(dim == 2) {
        uDotGradPhiPhi = _u[0] * _gradPhiPhi[2 * i + 0]
                       + _u[1] * _gradPhiPhi[2 * i + 1];
        gradMu_dot_gradPhiPhi = _gradmu[0] * _gradPhiPhi[2 * i + 0]
                              + _gradmu[1] * _gradPhiPhi[2 * i + 1];
        gradP_dot_gradPhiPhi  =  _gradp[0] * _gradPhiPhi[2 * i + 0]
                              +  _gradp[1] * _gradPhiPhi[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
        // Time derivative
        dphidt * _phiPhi[i]
        // Convective term
        - phi * uDotGradPhiPhi
        // Diffusive term
        + M * (gradMu_dot_gradPhiPhi + _alpha * gradP_dot_gradPhiPhi)
        // Source term if needed
        + Sphi * _phiPhi[i]
        );
    }

    //
    // mu test functions block : potential equation
    //
    double gradPhi_dot_gradPhiMu;
    for(int i = 0; i < _nFunctionsMu; ++i, ++I)
    {
      if constexpr(dim == 2) {
        gradPhi_dot_gradPhiMu = _gradphi[0] * _gradPhiMu[2 * i + 0]
                              + _gradphi[1] * _gradPhiMu[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
        // Mu term
        mu * _phiMu[i]
        // Phi laplacian
        - _tau * gradPhi_dot_gradPhiMu
        // Double well potential
        - _beta * (phi * phi * phi - phi) * _phiMu[i]
        // Source term if needed
        + Smu * _phiMu[i]
        );
    }
  }
}

template class feSysElm_CHNS_Alternative<2>;

// -----------------------------------------------------------------------------
// Bilinear form: convective acceleration of vector-valued field
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_VectorConvectiveAcceleration<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _uDotGradu.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  // _phiU.resize(_nFunctions);
  _gradPhiU.resize(_nComponents * _nComponents * _nFunctions);
  _uDotGraduDotPhiU.resize(_nFunctions);
  _u0DotGradPhiUDotPhiU.resize(_nFunctions*_nFunctions);
  _phiUDotGradu0DotPhiU.resize(_nFunctions*_nFunctions);
}

template <int dim>
void feSysElm_VectorConvectiveAcceleration<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get u, grad_u, phiU and gradPhiU
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // uSpace->getFunctionsAtQuadNode(k, _phiU);
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU);

    // Compute [(u0 dot gradu) + (u dot gradu0)] cdot phiU
    uSpace->vectorDotGradShapeDotShape(k, _gradPhiU, _u, _u0DotGradPhiUDotPhiU);
    uSpace->shapeDotTensorDotShape(k, _gradu, _phiUDotGradu0DotPhiU);

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        form->_Ae[i][j] += coeff * (_u0DotGradPhiUDotPhiU[i*_nFunctions+j]
                                  + _phiUDotGradu0DotPhiU[i*_nFunctions+j]) * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_VectorConvectiveAcceleration<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get u, grad_u and phiU
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // uSpace->getFunctionsAtQuadNode(k, _phiU);

    // Compute (u dot gradu)
    std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        // Contract _u with the lines of _gradu
        _uDotGradu[m] += _u[n] * _gradu[n * _nComponents + m];
      }
    }

    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);

    // Increment with (u dot gradu) dot v
    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= coeff * _uDotGraduDotPhiU[i] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_VectorConvectiveAcceleration<0>;
// template class feSysElm_VectorConvectiveAcceleration<1>;
template class feSysElm_VectorConvectiveAcceleration<2>;
// template class feSysElm_VectorConvectiveAcceleration<3>;

// -----------------------------------------------------------------------------------
// Bilinear form: convection of a scalar tracer C in the resolved velocity field u
// -----------------------------------------------------------------------------------
template <int dim>
void feSysElm_TracerConvection<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idC = 0; // Test functions of the tracer
  _idU = 1;
  _fieldsLayoutI = {_idC};
  _fieldsLayoutJ = {_idU, _idC};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsC = space[_idC]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _gradC.resize(_nComponents); // Dimension = number of velocity components
  _phiC.resize(_nFunctionsC);
  _gradPhiC.resize(_nComponents * _nFunctionsC);
  _phiUdotGradc0.resize(_nFunctionsU);
}

template <int dim>
void feSysElm_TracerConvection<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // u
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // gradC
    form->_intSpaces[_idC]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idC], k, form->_transformation, _gradC.data());
    // phiC
    form->_intSpaces[_idC]->getFunctionsAtQuadNode(k, _phiC);
    //gradPhiC
    form->_intSpaces[_idC]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiC.data());

    uSpace->dotProductShapeOther(k, _gradC, _phiUdotGradc0);

    for(int i = 0; i < _nFunctionsC; ++i)
    {
      int J = 0;
      for(int j = 0; j < _nFunctionsU; ++j, ++J) {
        // u dot gradC0 = phiU_j dot gradC0
        form->_Ae[i][J] += coeff * _phiC[i] * _phiUdotGradc0[j] * jac * _wQuad[k];
      }

      for(int j = 0; j < _nFunctionsC; ++j, ++J) {
        // Compute u0 dot gradC = u0 dot gradPhiC_j
        double u0DotGradC = 0.;
        for(int m = 0; m < _nComponents; ++m) {
          u0DotGradC += _u[m] * _gradPhiC[j * _nComponents + m];
        }
        form->_Ae[i][J] += coeff * _phiC[i] * u0DotGradC * jac * _wQuad[k];
      }

    }
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

template <int dim>
void feSysElm_TracerConvection<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get u, gradC and phiC
    uSpace->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idC]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idC], k, form->_transformation, _gradC.data());
    form->_intSpaces[_idC]->getFunctionsAtQuadNode(k, _phiC);

    // Compute u dot gradC
    double uDotGradC = 0.;
    for(int m = 0; m < _nComponents; ++m) {
      uDotGradC += _u[m] * _gradC[m];
    }

    // Increment elementwise residual with (u dot gradC) * phiC
    for(int i = 0; i < _nFunctionsC; ++i) {
      form->_Be[i] -= coeff * uDotGradC * _phiC[i] * jac * _wQuad[k];
    }
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local RHS on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     feInfo("B[%2d] = %+-1.10e", i, form->_Be[i]);
  // }
}

// template class feSysElm_TracerConvection<0>;
// template class feSysElm_TracerConvection<1>;
template class feSysElm_TracerConvection<2>;
// template class feSysElm_TracerConvection<3>;

// ------------------------------------------------------------------------------------------------------
// Bilinear form: convective acceleration of vector-valued field for the adjoint Navier-Stokes equations
// ------------------------------------------------------------------------------------------------------
template <int dim>
void feSysElm_VectorAdjointConvectiveAcceleration<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idUad = 1;
  _fieldsLayoutI = {_idUad}; // Local matrix is only the submatrix at (uAd, phi_uAD) ?
  _fieldsLayoutJ = {_idUad};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsUad = space[_idUad]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nFunctionsU);
  _uDotGraduAd.resize(_nComponents);
  _grad_uAd.resize(_nComponents * _nComponents);
  _phi_uAd.resize(_nFunctionsUad);
  _gradPhiU.resize(_nFunctionsU * _nFunctionsU);
}

template <int dim>
void feSysElm_VectorAdjointConvectiveAcceleration<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement FE matrix for VectorAdjointConvectiveAcceleration weak form.");
  exit(-1);
}

template <int dim>
void feSysElm_VectorAdjointConvectiveAcceleration<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get u, grad_uAd and phi_uAd
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idUad]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idUad], _nComponents, k, form->_transformation, _grad_uAd.data());
    form->_intSpaces[_idUad]->getFunctionsAtQuadNode(k, _phi_uAd);

    // Compute (u dot grad_uAd)
    std::fill(_uDotGraduAd.begin(), _uDotGraduAd.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGraduAd[m] += _u[n] * _grad_uAd[m * _nComponents + n];
        // ESSAI : grad transpose
        _uDotGraduAd[m] += _u[n] * _grad_uAd[n * _nComponents + m];
      }
    }

    // Increment with (u dot grad_uAd) dot v
    for(int i = 0; i < _nFunctionsUad; ++i) {
      form->_Be[i] -= coeff * _phi_uAd[i] * _uDotGraduAd[i % _nComponents] * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_VectorAdjointConvectiveAcceleration<0>;
// template class feSysElm_VectorAdjointConvectiveAcceleration<1>;
template class feSysElm_VectorAdjointConvectiveAcceleration<2>;
// template class feSysElm_VectorAdjointConvectiveAcceleration<3>;

// -----------------------------------------------------------------------------
// Bilinear form: weak divergence of newtonian stress tensor (mixed u and p)
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_DivergenceNewtonianStress<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU};
  _fieldsLayoutJ = {_idU, _idP};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradu.resize(_nComponents * _nComponents);
  _symmetricGradu.resize(_nComponents * _nComponents);
  _gradPhiU.resize(_nComponents * _nComponents * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);
  _doubleContractionPhiPhi.resize(_nFunctionsU*_nFunctionsU);
  _doubleContractionPhiPhiT.resize(_nFunctionsU*_nFunctionsU);
}

template <int dim>
void feSysElm_DivergenceNewtonianStress<dim>::computeAe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, mu, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    mu = (*_viscosity)(form->_args);
    coeff = (*_coeff)(form->_args);

    // Get relevant values and gradients
    // grad_phiU
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeGradShape(_gradPhiU, _doubleContractionPhiPhi);
    uSpace->doubleContractionGradShapeGradShapeTransposed(_gradPhiU, _doubleContractionPhiPhiT);
    // phiP
    form->_intSpaces[_idP]->getFunctionsAtQuadNode(k, _phiP);

    for(int i = 0; i < _nFunctionsU; ++i)
    {
      int J = 0;
      // phiU - u block (viscous part)
      for(int j = 0; j < _nFunctionsU; ++j)
      {
        form->_Ae[i][J++] += - coeff * mu * (_doubleContractionPhiPhi [i*_nFunctionsU+j]
                                           + _doubleContractionPhiPhiT[i*_nFunctionsU+j]) * jac * _wQuad[k];
      }
      // phiU - p block (pressure part)
      for(int j = 0; j < _nFunctionsP; ++j)
      {
        form->_Ae[i][J++] += coeff * _phiP[j] * _divPhiU[i] * jac * _wQuad[k];
      }
    }
  }
}

template <int dim>
void feSysElm_DivergenceNewtonianStress<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace = static_cast<feVectorSpace<dim>*>(form->_intSpaces[_idU]);

  double jac, coeff, mu, p;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    mu = (*_viscosity)(form->_args);
    coeff = (*_coeff)(form->_args);

    // Get p, grad_u and gradphiU
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU);

    // Symmetric velocity gradient
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _symmetricGradu[m * _nComponents + n] =
          _gradu[m * _nComponents + n] + _gradu[n * _nComponents + m];
      }
    }

    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _symmetricGradu, _doubleContraction);

    for(int i = 0; i < _nFunctionsU; ++i) {
      form->_Be[i] -= coeff * (p * _divPhiU[i] - mu * _doubleContraction[i]) * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_DivergenceNewtonianStress<0>;
// template class feSysElm_DivergenceNewtonianStress<1>;
template class feSysElm_DivergenceNewtonianStress<2>;
// template class feSysElm_DivergenceNewtonianStress<3>;

// -----------------------------------------------------------------------------
// Bilinear forms: Stabilization for the (Navier-)Stokes equations
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_Stokes_SUPG_PSPG<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

template <int dim>
void feSysElm_Stokes_SUPG_PSPG<dim>::computeAe(feBilinearForm *form)
{
  // Computed with finite differences
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Local matrix for SUPG/PSPG Stokes should be computed with finite differences");
  exit(-1);
}

template <int dim>
void feSysElm_Stokes_SUPG_PSPG<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, rho, mu, tau, dotProdU, dotProdP, residualQ[2];

  for(int k = 0, cnt; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    rho = (*_density)(form->_args);
    mu = (*_viscosity)(form->_args);
    (*_volumeForce)(form->_args, _f);

    // Get u, grad_u, grad_p, gradphiU, gradphiP and _hessu (2nd order derivatives)
    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // grad_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_p
    form->_intSpaces[_idP]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idP], k, form->_transformation, _gradp.data());
    // grad_phiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());
    // grad_phiP
    form->_intSpaces[_idP]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiP.data());
    // hess_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalHessian(
      form->_sol[_idU], _nComponents, k, form->_transformation, _hessu.data());
    // hessPhiU
    form->_intSpaces[_idU]->getFunctionsPhysicalHessianAtQuadNode(k, form->_transformation,
                                                                  _hessPhiU.data());

    // Stabilization parameter
    tau = compute_tau<dim>(form->_geoCoord, _u.data(), mu, _nFunctionsU, _gradPhiU, jac);

    // Residual for Navier-Stokes (laplacian form):
    double lapu[2] = {_hessu[0] + _hessu[2], _hessu[3] + _hessu[5]};
    _residual[0] = -_gradp[0] + mu * lapu[0] + rho * _f[0];
    _residual[1] = -_gradp[1] + mu * lapu[1] + rho * _f[1];

    cnt = 0;
    for(int iU = 0; iU < _nFunctionsU; ++iU) {
      double d2phidx2 = _hessPhiU[iU * 3 + 0];
      double d2phidy2 = _hessPhiU[iU * 3 + 2];
      double lap_phiU = d2phidx2 + d2phidy2;

      double residualV = mu * lap_phiU;

      // Dot product only has one non-zero term
      dotProdU = residualV * _residual[iU % _nComponents];

      form->_Be[cnt++] += coeff * tau * dotProdU * jac * _wQuad[k];
    }

    for(int iP = 0; iP < _nFunctionsP; ++iP) {
      residualQ[0] = -_gradPhiP[iP * _nComponents + 0];
      residualQ[1] = -_gradPhiP[iP * _nComponents + 1];

      dotProdP = _residual[0] * residualQ[0] + _residual[1] * residualQ[1];

      form->_Be[cnt++] -= coeff * tau * dotProdP * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_Stokes_SUPG_PSPG<0>;
// template class feSysElm_Stokes_SUPG_PSPG<1>;
template class feSysElm_Stokes_SUPG_PSPG<2>;
// template class feSysElm_Stokes_SUPG_PSPG<3>;

template <int dim>
void feSysElm_NS_SUPG_PSPG<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);
  _uDotGradu.resize(_nComponents);
  _uDotGradPhiu.resize(_nComponents);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

template <int dim>
void feSysElm_NS_SUPG_PSPG<dim>::computeAe(feBilinearForm *form)
{
  // Computed with finite differences
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Local matrix for SUPG/PSPG Navier-Stokes should be computed with finite differences");
  exit(-1);
}

template <int dim>
void feSysElm_NS_SUPG_PSPG<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, rho, mu, tau, dotProdU, dotProdP, residualQ[2];

  for(int k = 0, cnt; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);
    rho = (*_density)(form->_args);
    mu = (*_viscosity)(form->_args);
    (*_volumeForce)(form->_args, _f);

    // Get u, grad_u, grad_p, gradphiU, gradphiP and _hessu (2nd order derivatives)
    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // grad_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_p
    form->_intSpaces[_idP]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idP], k, form->_transformation, _gradp.data());
    // grad_phiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());
    // grad_phiP
    form->_intSpaces[_idP]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiP.data());
    // hess_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalHessian(
      form->_sol[_idU], _nComponents, k, form->_transformation, _hessu.data());
    // hessPhiU
    form->_intSpaces[_idU]->getFunctionsPhysicalHessianAtQuadNode(k, form->_transformation,
                                                                  _hessPhiU.data());

    // Stabilization parameter
    tau = compute_tau<dim>(form->_geoCoord, _u.data(), mu, _nFunctionsU, _gradPhiU, jac);

    // Compute (u dot gradu)
    std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGradu[m] += _u[n] * _gradu[m * _nComponents + n];
      }
    }

    // Residual for Navier-Stokes (laplacian form):
    double lapu[2] = {_hessu[0] + _hessu[2], _hessu[3] + _hessu[5]};
    _residual[0] = -_uDotGradu[0] - _gradp[0] + mu * lapu[0] + rho * _f[0];
    _residual[1] = -_uDotGradu[1] - _gradp[1] + mu * lapu[1] + rho * _f[1];

    cnt = 0;
    for(int iU = 0; iU < _nFunctionsU; ++iU) {
      // Compute (u dot grad_v)
      double uDotGradPhiU = 0.;
      for(int n = 0; n < _nComponents; ++n) {
        uDotGradPhiU += _u[n] * _gradPhiU[iU * _nComponents + n];
      }

      double d2phidx2 = _hessPhiU[iU * 3 + 0];
      double d2phidy2 = _hessPhiU[iU * 3 + 2];
      double lap_phiU = d2phidx2 + d2phidy2;

      double residualV = -uDotGradPhiU + mu * lap_phiU;

      // Dot product only has one non-zero term
      dotProdU = residualV * _residual[iU % _nComponents];

      form->_Be[cnt++] += coeff * tau * dotProdU * jac * _wQuad[k];
    }

    for(int iP = 0; iP < _nFunctionsP; ++iP) {
      residualQ[0] = -_gradPhiP[iP * _nComponents + 0];
      residualQ[1] = -_gradPhiP[iP * _nComponents + 1];

      dotProdP = _residual[0] * residualQ[0] + _residual[1] * residualQ[1];

      form->_Be[cnt++] -= coeff * tau * dotProdP * jac * _wQuad[k];
    }
  }
}

// template class feSysElm_NS_SUPG_PSPG<0>;
// template class feSysElm_NS_SUPG_PSPG<1>;
template class feSysElm_NS_SUPG_PSPG<2>;
// template class feSysElm_NS_SUPG_PSPG<3>;

template <int dim>
void feSysElm_GLS_Stokes<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

template <int dim>
void feSysElm_GLS_Stokes<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

template <int dim>
void feSysElm_GLS_Stokes<dim>::computeBe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

// template class feSysElm_GLS_Stokes<0>;
// template class feSysElm_GLS_Stokes<1>;
template class feSysElm_GLS_Stokes<2>;
// template class feSysElm_GLS_Stokes<3>;

template <int dim>
void feSysElm_GLS_NavierStokes<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);
  _uDotGradu.resize(_nComponents);
  _uDotGradPhiu.resize(_nComponents);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

template <int dim>
void feSysElm_GLS_NavierStokes<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

template <int dim>
void feSysElm_GLS_NavierStokes<dim>::computeBe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

// template class feSysElm_GLS_NavierStokes<0>;
// template class feSysElm_GLS_NavierStokes<1>;
template class feSysElm_GLS_NavierStokes<2>;
// template class feSysElm_GLS_NavierStokes<3>;