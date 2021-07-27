#ifndef _FEEDGE_
#define _FEEDGE_

#include <cmath>

#include "feVertex.h"

class Edge {
protected:
  Vertex *_v0;
  Vertex *_v1;
  int _tag;

public:
  Edge() {}
  Edge(Vertex *v0, Vertex *v1, int tag = -1) {
    _v0 = v0;
    _v1 = v1;
    _tag = tag;
  }
  ~Edge() {}

  inline Vertex *getVertex(int i) { return (i == 0) ? _v0 : _v1; }
  // Edge's own tag
  inline int getTag() const { return _tag; }
  // Tags of the vertices
  inline int getTag(int i) const { return (i == 0) ? _v0->getTag() : _v1->getTag(); }

  void print() { printf("Edge %4d : %4d - %4d\n", _tag, _v0->getTag(), _v1->getTag()); }
};

struct EdgeLessThan {
  bool operator()(const Edge e1, const Edge e2) const {
    int tag10 = e1.getTag(0), tag11 = e1.getTag(1);
    int tag20 = e2.getTag(0), tag21 = e2.getTag(1);
    int diffMin = fmin(tag10, tag11) - fmin(tag20, tag21);
    int diffMax = fmax(tag10, tag11) - fmax(tag20, tag21);
    if(diffMin < 0) return true;
    if(diffMin > 0) return false;
    if(diffMax < 0) return true;
    if(diffMax > 0) return false;
    return false;
  }
};

#endif