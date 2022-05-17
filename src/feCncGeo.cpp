#include "feNG.h"
#include "feCncGeo.h"
#include "feSpace.h"
#include "feMesh.h"

void feCncGeo::computeJacobians()
{
  int nQuad = _space->getNbQuadPoints();
  _J.resize(_nElm * nQuad);

  std::vector<double> geoCoord;

  if(_space->getDim() == 1) {
    std::vector<double> j(3, 0.0);
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      geoCoord = _mesh->getCoord(_tag, iElm);
      for(int k = 0; k < nQuad; ++k) {
        _space->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j);
        _J[nQuad * iElm + k] = sqrt(j[0] * j[0] + j[1] * j[1]);
      }
    }
  } else if(_space->getDim() == 2) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      geoCoord = _mesh->getCoord(_tag, iElm);
      for(int k = 0; k < nQuad; ++k) {
        _space->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        _space->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        _J[nQuad * iElm + k] = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
        if(_J[nQuad * iElm + k] < 0) {
          printf("In feCncGeo::computeJacobians : Error - Element jacobian = %+-12.12e\n",
                 _J[nQuad * iElm + k]);
          exit(-1);
        }
      }
    }
  } else if(_space->getDim() == 0) {
    _J[0] = 1.0;
  }

  else {
    printf("TODO : Implement jacobian for 3D elements.\n");
  }
}

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

  int activColor = -1;
  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true && nbElmPerColor[i] == minColor) {
      activColor = i;
    }
  }

  if(activColor == -1) {
    activColor = nbColor;
    nbColor = nbColor + 1;
    nbElmPerColor.push_back(0);
    listElmPerColor.push_back(std::vector<int>(0, 0));
  }
  return activColor;
}

void feCncGeo::colorElements(int coloringAlgorithm)
{
  feInfo("Coloring connectivity %s ...", _ID.c_str());
  tic();
  // feInfo("Coloring connectivity with %d elm", _nElm);
  // feInfo("Coloring connectivity with %d node per elm", _nNodPerElm);
  // feInfo("Coloring connectivity with %d sommets", _nNod);

  std::map<int, std::set<int> > listElmPerElm;

  // Create the node and element patches :
  // Nodes
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _nNodPerElm; ++j) {
      int nds = _connecNodes[i * _nNodPerElm + j];
      _listElmPerNode[nds].push_back(i);
    }
  }

  for(auto pair : _listElmPerNode) {
    _nbElmPerNode[pair.first] = pair.second.size();
  }

  // Elements
  // for(int i = 0; i < _nElm; ++i) {
  //   for(int j = 0; j < _nNodPerElm; ++j) {
  //     int nds = _connecNodes[i * _nNodPerElm + j];
  //     for(int elm : _listElmPerNode[nds]){
  //       if(elm != i)
  //         listElmPerElm[i].insert(elm);
  //     }
  //   }
  // }

  // for(auto pair : listElmPerElm) {
  //   _nbElmPerElm[pair.first] = pair.second.size();
  // }

  // // Convert sets to vectors
  // for(int i = 0; i < _nElm; ++i) {
  //   std::vector<int> v(listElmPerElm[i].begin(), listElmPerElm[i].end());
  //   _listElmPerElm[i] = v;
  // }

  switch(coloringAlgorithm) {
    case 1: // Returning a non homogeneous distribution of elements per color
    {
      _nbColor = 0;
      _elmToColor.resize(_nElm, -1);

      bool noColor = false;
      while(!noColor) {
        for(int k = 0; k < _nElm; ++k) {
          if(_elmToColor[k] == -1) {
            for(int i = 0; i < _nNodPerElm; ++i) {
              int s = _connecNodes[k * _nNodPerElm + i];
              for(int j = 0; j < _nbElmPerNode[s]; ++j) {
                int numElm = _listElmPerNode[s][j];
                if(numElm != k) {
                  if(_elmToColor[numElm] < 0) {
                    _elmToColor[numElm] = -2;
                  }
                }
              }
            }
          }
        }

        _nbElmPerColor.push_back(0);
        _listElmPerColor.push_back(std::vector<int>(0, 0));

        noColor = true;
        for(int k = 0; k < _nElm; ++k) {
          if(_elmToColor[k] == -1) {
            _elmToColor[k] = _nbColor;
            _nbElmPerColor[_nbColor] = _nbElmPerColor[_nbColor] + 1;
            _listElmPerColor[_nbColor].push_back(k);
          } else if(_elmToColor[k] == -2) {
            _elmToColor[k] = -1;
            noColor = false;
          }
        }
        _nbColor = _nbColor + 1;

      } // while
      feInfo("nbColor : %d", _nbColor);
      for(int i = 0; i < _nbColor; i++) {
        feInfo("%d Elm for color %d", _nbElmPerColor[i], i);
      }
      break;
    }

      // case 2:  //Returning a homogeneous distribution of elements per color
      // {
      //   //Finding the Elm with the most neighbours to asign nb Colors
      //   _nbColor=0;
      //   for(auto pair : _nbElmPerElm){
      //     if(pair.second>_nbColor){
      //       _nbColor=pair.second;
      //     }
      //   }

      //   _elmToColor.resize(_nElm, -1);
      //   _nbElmPerColor.resize(_nbColor);
      //   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

      //   for(int iElm = 0; iElm < _nElm; iElm++) {
      //     if(_elmToColor[iElm] == -1) {
      //       std::vector<int> iPatchElm = _listElmPerElm[iElm];
      //       std::vector<bool> availableColor(_nbColor, true);
      //       for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
      //         if(_elmToColor[iPatchElm[i]] != -1 && availableColor[_elmToColor[iPatchElm[i]]] ==
      //         true) {
      //           availableColor[_elmToColor[iPatchElm[i]]] = false;
      //         }
      //       }

      //       int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor,
      //       _listElmPerColor); _elmToColor[iElm] = activColor; _nbElmPerColor[activColor] =
      //       _nbElmPerColor[activColor] + 1; _listcdElmPerColor[activColor].push_back(iElm);
      //     }
      //   }
      //   feInfo("nbColor : %d", _nbColor);
      //   for (int i=0;i<_nbColor;i++){
      //     feInfo("%d Elm for color %d",_nbElmPerColor[i],i);
      //   }
      //   break;
      // }

    case 3: // Idem 2 but not using Patch Elm
    {
      // Finding the Node with the most neighbours to asign initial nb Colors
      _nbColor = 0;
      for(auto pair : _nbElmPerNode) {
        if(pair.second > _nbColor) {
          _nbColor = pair.second;
        }
      }

      _elmToColor.resize(_nElm, -1);
      _nbElmPerColor.resize(_nbColor);
      _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

      bool noColor = false;
      while(noColor == false) {
        for(int k = 0; k < _nElm; ++k) {
          if(_elmToColor[k] == -1) {
            std::vector<bool> availableColor(_nbColor, true);
            for(int i = 0; i < _nNodPerElm; ++i) {
              int s = _connecNodes[k * _nNodPerElm + i];
              for(int j = 0; j < _nbElmPerNode[s]; ++j) {
                int numElm = _listElmPerNode[s][j];
                if(numElm != k) {
                  if(_elmToColor[numElm] < 0) {
                    _elmToColor[numElm] = -2;
                  } else if(availableColor[_elmToColor[numElm]] == true) {
                    availableColor[_elmToColor[numElm]] = false;
                  }
                }
              }
            }
            int activColor =
              colorChoice(availableColor, _nbElmPerColor, _nbColor, _listElmPerColor);
            _elmToColor[k] = activColor;
            // feInfo("%d couleur pour elm %d",activColor,k);
            _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
            _listElmPerColor[activColor].push_back(k);
          }
        }

        noColor = true;
        for(int k = 0; k < _nElm; ++k) {
          if(_elmToColor[k] == -2) {
            _elmToColor[k] = -1;
            noColor = false;
          }
        }

      } // while

      feInfo("nbColor : %d", _nbColor);
      for(int i = 0; i < _nbColor; i++) {
        feInfo("%d Elm for color %d", _nbElmPerColor[i], i);
      }
      break;
    }
  }
  toc();
  feInfo("Done");
}
