#ifndef _FECOLORING_
#define _FECOLORING_

//
// A coloring of the connectivity on a Physical Entity.
//
// Contributed by Baptiste Berlioux.
// Refactored by Arthur Bawin.
//
class feColoring
{
protected:

  // Choice of coloring algorithm:
  //  - 1 : Non homogeneous distribution of elements per color.
  //  - 2 : Homogeneous distribution of elements per color.
  //  - 3 : Homogeneous distribution of elements per color,
  //        does not require creating patches of elements.
  int _coloringAlgorithm = 1;

  int _numColors;
  std::vector<int> _elemToColor;
  std::vector<int> _numElmPerColor;
  std::vector<std::vector<int>> _elementsInColor;

public:
  //
  // Color a connectivity
  //
  feColoring(int coloringAlgorithm, int numVerticesPerElem, std::vector<int> &vertices, std::vector<int> &elements);
  ~feColoring() {}

  int getNbColor() const { return _numColors; };
  const std::vector<int> &getColorElm() const { return _elemToColor; };
  const std::vector<int> &getNbElmPerColor() const { return _numElmPerColor; };
  const std::vector<std::vector<int> > &getListElmPerColor() const { return _elementsInColor; };
  int getNbElmPerColorI(int i) const { return _numElmPerColor[i]; };
  const std::vector<int> &getListElmPerColorI(int i) const { return _elementsInColor[i]; };
  int getElmColored(int iColor, int iElmC) const { return _elementsInColor[iColor][iElmC]; };
};

#endif
