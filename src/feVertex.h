#ifndef _FEVERTEX_
#define _FEVERTEX_

// Inspired by MFEM

class Vertex{

protected:
	double coord[3];
  int _tag; // The Gmsh numbering of the node, used to create edges
public:
	Vertex(){}
	Vertex(double x, int tag = -1){ _tag = tag; coord[0] = x; }
	Vertex(double x, double y, int tag = -1){ _tag = tag; coord[0] = x; coord[1] = y; }
	Vertex(double x, double y, double z, int tag = -1){ _tag = tag; coord[0] = x; coord[1] = y; coord[2] = z; }
	~Vertex() {}

  double x() const { return coord[0]; } 
  double y() const { return coord[1]; } 
  double z() const { return coord[2]; } 
  // Returns pointer to the coordinates of the vertex.
  inline double *operator() () const { return (double*) coord; }
  // Returns the i'th coordinate of the vertex.
  inline double &operator() (int i) { return coord[i]; }

  inline int getTag() const { return _tag; }
};

#endif