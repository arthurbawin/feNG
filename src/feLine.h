#ifndef _FELINE_
#define _FELINE_

#include "feElement.h"

// Inspired from MFEM

class feLine : public feElement{

protected:
	int nodes[2];

public:
	feLine() : feElement(LINE) {}
	virtual ~feLine() {}

  virtual GeoType getGeoType() { return LINE; };
};

#endif