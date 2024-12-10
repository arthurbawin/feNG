#include "feNG.h"
#include "feColoring.h"
#include "feMessage.h"

extern int FE_VERBOSE;

static int colorChoice(std::vector<bool> availableColor, std::vector<int> &nbElmPerColor,
                       int &nbColor, std::vector<std::vector<int> > &listElmPerColor)
{
  int minColor = *max_element(nbElmPerColor.begin(), nbElmPerColor.end());
  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true) {
      if(minColor > nbElmPerColor[i]) {
        minColor = nbElmPerColor[i];
      }
    }
  }

  int activeColor = -1;
  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true && nbElmPerColor[i] == minColor) {
      activeColor = i;
    }
  }

  if(activeColor == -1) {
    activeColor = nbColor;
    nbColor = nbColor + 1;
    nbElmPerColor.push_back(0);
    listElmPerColor.push_back(std::vector<int>(0, 0));
  }
  return activeColor;
}

feColoring::feColoring(int coloringAlgorithm, int numVerticesPerElem,
  std::vector<int> &vertices, std::vector<int> &elements)
{
  feInfoCond(FE_VERBOSE > 0, "\t\tColoring connectivity...");
  
  tic();
  int nElements = elements.size();

  // Create the node patch :
  int size = *std::max_element(vertices.begin(), vertices.end()) + 1;
  std::vector<std::vector<int>> listElmPerNode(size);
  std::vector<int> numElemPerNode(size);

  for(int i = 0; i < nElements; ++i) {
    for(int j = 0; j < numVerticesPerElem; ++j) {
      int nds = vertices[i * numVerticesPerElem + j];
      numElemPerNode[nds] += 1;
      listElmPerNode[nds].push_back(i);
    }
  }

  switch(coloringAlgorithm) {
    case 1: // Returning a non homogeneous distribution of elements per color
    {
      _numColors = 0;
      _elemToColor.resize(nElements, -1);

      bool noColor = false;
      while(!noColor) {
        for(int k = 0; k < nElements; ++k) {
          if(_elemToColor[k] == -1) {
            for(int i = 0; i < numVerticesPerElem; ++i) {
              int s = vertices[k * numVerticesPerElem + i];
              for(int j = 0; j < numElemPerNode[s]; ++j) {
                int numElm = listElmPerNode[s][j];
                if(numElm != k) {
                  if(_elemToColor[numElm] < 0) {
                    _elemToColor[numElm] = -2;
                  }
                }
              }
            }
          }
        }

        _numElmPerColor.push_back(0);
        _elementsInColor.push_back(std::vector<int>(0, 0));

        noColor = true;
        for(int k = 0; k < nElements; ++k) {
          if(_elemToColor[k] == -1) {
            _elemToColor[k] = _numColors;
            _numElmPerColor[_numColors] = _numElmPerColor[_numColors] + 1;
            _elementsInColor[_numColors].push_back(k);
          } else if(_elemToColor[k] == -2) {
            _elemToColor[k] = -1;
            noColor = false;
          }
        }
        _numColors = _numColors + 1;

      } // while
      break;
    }

    case 2: // Returning a homogeneous distribution of elements per color
    {
      // Create the Elements patch :
      int sizeElm = *std::max_element(elements.begin(), elements.end()) + 1;
      std::vector<int> nbElmPerElm(sizeElm);
      std::vector<std::vector<int>> listElmPerElm(sizeElm);
      std::vector<std::set<int>> listElmPerElmSet(sizeElm);
      for(int i = 0; i < nElements; ++i) {
        int elm = elements[i];
        for(int j = 0; j < numVerticesPerElem; ++j) {
          int nds = vertices[i * numVerticesPerElem + j];
          for(int ngb : listElmPerNode[nds]) {
            if(ngb != i) {
              listElmPerElmSet[elm].insert(ngb);
            }
          }
        }
      }
      // Convert sets to vectors
      for(int i = 0; i < nElements; ++i) {
        int elm = elements[i];
        std::vector<int> v(listElmPerElmSet[elm].begin(), listElmPerElmSet[elm].end());
        nbElmPerElm[elm] = v.size();
        listElmPerElm[elm] = v;
      }

      // Finding the Elm with the most neighbours to asign nb Colors
      _numColors = *std::max_element(nbElmPerElm.begin(), nbElmPerElm.end()) + 1;

      _elemToColor.resize(nElements, -1);
      _numElmPerColor.resize(_numColors);
      _elementsInColor.resize(_numColors, std::vector<int>(0, 0));
      for(int iElm = 0; iElm < nElements; iElm++) {
        if(_elemToColor[iElm] == -1) {
          int elm = elements[iElm];
          // feInfo("%d",elm);
          std::vector<int> iPatchElm = listElmPerElm[elm];
          std::vector<bool> availableColor(_numColors, true);
          for(int i = 0; i < nbElmPerElm[elm]; i++) {
            if(_elemToColor[iPatchElm[i]] != -1 &&
               availableColor[_elemToColor[iPatchElm[i]]] == true) {
              availableColor[_elemToColor[iPatchElm[i]]] = false;
            }
          }

          int activeColor = colorChoice(availableColor, _numElmPerColor, _numColors, _elementsInColor);
          _elemToColor[iElm] = activeColor;
          _numElmPerColor[activeColor] = _numElmPerColor[activeColor] + 1;
          _elementsInColor[activeColor].push_back(iElm);
        }
      }
      break;
    }

    case 3: // Idem 2 but not using Patch Elm
    {
      // Finding the Node with the most neighbours to asign initial nb Colors
      _numColors = *std::max_element(numElemPerNode.begin(), numElemPerNode.end());
      // for(auto pair : numElemPerNode) {
      //   if(pair.second > _numColors) {
      //     _numColors = pair.second;
      //   }
      // }

      _elemToColor.resize(nElements, -1);
      _numElmPerColor.resize(_numColors);
      _elementsInColor.resize(_numColors, std::vector<int>(0, 0));

      bool noColor = false;
      while(noColor == false) {
        for(int k = 0; k < nElements; ++k) {
          if(_elemToColor[k] == -1) {
            std::vector<bool> availableColor(_numColors, true);
            for(int i = 0; i < numVerticesPerElem; ++i) {
              int s = vertices[k * numVerticesPerElem + i];
              for(int j = 0; j < numElemPerNode[s]; ++j) {
                int numElm = listElmPerNode[s][j];
                if(numElm != k) {
                  if(_elemToColor[numElm] < 0) {
                    _elemToColor[numElm] = -2;
                  } else if(availableColor[_elemToColor[numElm]] == true) {
                    availableColor[_elemToColor[numElm]] = false;
                  }
                }
              }
            }
            int activeColor =
              colorChoice(availableColor, _numElmPerColor, _numColors, _elementsInColor);
            _elemToColor[k] = activeColor;
            // feInfo("%d couleur pour elm %d",activeColor,k);
            _numElmPerColor[activeColor] = _numElmPerColor[activeColor] + 1;
            _elementsInColor[activeColor].push_back(k);
          }
        }

        noColor = true;
        for(int k = 0; k < nElements; ++k) {
          if(_elemToColor[k] == -2) {
            _elemToColor[k] = -1;
            noColor = false;
          }
        }
      } // while
      break;
    }
  }

  feInfoCond(FE_VERBOSE > 0, "\t\tDone in %f s", toc());

  feInfoCond(FE_VERBOSE > 0, "\t\tColoring statistics:");
  for(size_t i = 0; i < _numElmPerColor.size(); ++i) {
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of elements in color %d: %d", i, _numElmPerColor[i]);
  }
}