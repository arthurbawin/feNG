#ifndef _FEELEMENT_
#define _FEELEMENT_

// Inspired from MFEM

class feElement {
public:
  enum GeoType { POINT, LINE, TRIANGLE, QUADRILATERAL };

protected:
  int _physTag;
  GeoType _baseType;

public:
  feElement(GeoType baseType = POINT) : _physTag(-1), _baseType(baseType) {}
  virtual ~feElement() {}

  int getPhysicalTag() { return _physTag; };
  void setPhysicalTag(int tag) { _physTag = tag; };

  virtual GeoType getGeoType() = 0;
};

#endif