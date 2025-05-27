#ifndef _FECONSTANTS_
#define _FECONSTANTS_

//
// Useful scalar and vector constants
//
namespace scalarConstant {

	const feConstantFunction zero(0.);
  const feConstantFunction one(1.);
  const feConstantFunction minusOne(-1.);
}

namespace vectorConstant {
  const feConstantVectorFunction zero({0., 0., 0.});
  const feConstantVectorFunction one({1., 1., 1.});
  const feConstantVectorFunction one_x({1., 0., 0.});
  const feConstantVectorFunction one_y({0., 1., 0.});
  const feConstantVectorFunction one_z({0., 0., 1.});
  const feConstantVectorFunction minusOne_x({-1., 0., 0.});
  const feConstantVectorFunction minusOne_y({0., -1., 0.});
  const feConstantVectorFunction minusOne_z({0., 0., -1.});
}

#endif