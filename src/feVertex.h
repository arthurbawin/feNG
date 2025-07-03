#ifndef _FEVERTEX_
#define _FEVERTEX_

// Inspired by MFEM

class Vertex
{
protected:
  double coord[3];
  int _tag; // The Gmsh numbering of the node, used to create edges

public:
  Vertex() {}
  Vertex(double x, int tag = -1)
  {
    _tag = tag;
    coord[0] = x;
  }
  Vertex(double x, double y, int tag = -1)
  {
    _tag = tag;
    coord[0] = x;
    coord[1] = y;
  }
  Vertex(double x, double y, double z, int tag = -1)
  {
    _tag = tag;
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
  }
  Vertex(const std::vector<double> &pos, int tag = -1)
  {
    _tag = tag;
    coord[0] = pos[0];
    coord[1] = pos[1];
    coord[2] = pos[2];
  }
  ~Vertex() {}

  double x() const { return coord[0]; }
  double y() const { return coord[1]; }
  double z() const { return coord[2]; }

  // Returns pointer to the coordinates of the vertex.
  inline double *operator()() const { return (double *)coord; }
  // Returns the i'th coordinate of the vertex.
  inline double &operator()(int i) { return coord[i]; }

  inline int getTag() const { return _tag; }

  // From Gmsh's MVertex.h
  double distance(const Vertex &v) const
  {
    double dx = coord[0] - v.coord[0];
    double dy = coord[1] - v.coord[1];
    double dz = coord[2] - v.coord[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  Vertex operator+(const Vertex &other) const
  {
    Vertex v(*this);
    v.coord[0] += other.coord[0];
    v.coord[1] += other.coord[1];
    v.coord[2] += other.coord[2];
    return v;
  }

  Vertex operator*(const double factor) const
  {
    Vertex v(*this);
    v.coord[0] *= factor;
    v.coord[1] *= factor;
    v.coord[2] *= factor;
    return v;
  }
};

#endif