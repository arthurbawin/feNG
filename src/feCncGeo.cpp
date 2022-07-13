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
        _J[nQuad * iElm + k] = sqrt(j[0]*j[0] + j[1]*j[1]);
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
        if(_J[nQuad * iElm + k] <= 0) {
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

static int colorChoice(std::vector<bool> availableColor, std::vector<int> nbElmPerColor,
                        int nbColor)
{
  int minColor = *max_element(nbElmPerColor.begin(), nbElmPerColor.end());
  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true) {
      if(minColor >= nbElmPerColor[i]) {
        minColor = nbElmPerColor[i];
      }
    }
  }

  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true && nbElmPerColor[i] == minColor) {
      return i;
    }
  }
}

void feCncGeo::colorElements(int coloringAlgorithm){

  // feInfo("Coloring connectivity %s", _ID.c_str());
  // feInfo("Coloring connectivity with %d elm", _nElm);
  // feInfo("Coloring connectivity with %d node per elm", _nNodPerElm);
  // feInfo("Coloring connectivity with %d sommets", _nNod);

  std::map<int, std::set<int>> listElmPerElm;

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
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _nNodPerElm; ++j) {
      int nds = _connecNodes[i * _nNodPerElm + j];
      for(int elm : _listElmPerNode[nds]){
        if(elm != i)
          listElmPerElm[i].insert(elm);
      }
    }
  }

  // _nbElmPerElm.resize(_nElm);
  // for(int i = 0; i < _nElm; ++i) {
  //   _nbElmPerElm[i] = _listElmPerElm[i].size();
  // }
  for(auto pair : listElmPerElm) {
    _nbElmPerElm[pair.first] = pair.second.size();
  }

  // Convert sets to vectors
  for(int i = 0; i < _nElm; ++i) {
    std::vector<int> v(listElmPerElm[i].begin(), listElmPerElm[i].end());
    _listElmPerElm[i] = v;
  }

  // for(int i = 0; i < _nElm; ++i) {
  //   for(int j = 0; j < _nNodPerElm; ++j) {
  //     int nds = _connecNodes[i * _nNodPerElm + j];
  //     for(int k = 0; k < _nbElmPerNode[nds]; ++k) {
  //       if(_listElmPerNode[nds][k] != i) {
  //         int check = 0;
  //         for(int iElm = 0; iElm < _nbElmPerElm[i]; ++iElm) {
  //           if(_listElmPerNode[nds][elm] == _listElmPerElm[i][iElm]) {
  //             check++;
  //             break;
  //           }
  //         }
  //         if(check == 0) {
  //           _nbElmPerElm[i]++;
  //           _listElmPerElm[i].push_back(_listElmPerNode[nds][elm]);
  //         }
  //       }
  //     }
  //   }
  // }

  // feInfo("LIST ELM PER NODE");
  // for(auto p : _listElmPerNode) {
  //   feInfo("Vector of size %d", p.second.size());
  //   feInfo("Il y a %d elements associes au noeud", _nbElmPerNode[p.first]);
  //   for(auto val : p.second){
  //     feInfo("val = %d", val);
  //   }
  // }
  // feInfo("LIST ELM PER ELM");
  // for(auto p : _listElmPerElm) {
  //   feInfo("Vector of size %d", p.second.size());
  //   feInfo("Il y a %d elements associes à l'elm", _nbElmPerElm[p.first]);
  //   for(auto val : p.second){
  //     feInfo("val = %d", val);
  //   }
  // }

  switch(coloringAlgorithm){
    case 1 :
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
      break;
    }
    // case 2 :
    // {
    //   _nbColor = *max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;
    //   _colorElm.resize(_nElm, -1);
    //   _nbElmPerColor.resize(_nbColor);
    //   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));
    //   std::vector<int> color(_nbColor);
    //   for(int i = 0; i < _nbColor; ++i) {
    //     color[i] = i;
    //   }

    //   for(int iElm = 0; iElm < _nElm; iElm++) {
    //     if(_colorElm[iElm] == -1) {
    //       std::vector<int> iElmPatch = _listElmPerElm[iElm];
    //       std::vector<bool> availableColor(_nbColor, true);
    //       for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
    //         if(_colorElm[iElmPatch[i]] != -1 && availableColor[_colorElm[iElmPatch[i]]] == true) {
    //           availableColor[_colorElm[iElmPatch[i]]] = false;
    //         }
    //       }

    //       for(int iC = 0; iC < _nbColor; iC++) {
    //         if(availableColor[iC] == true) {
    //           _colorElm[iElm] = iC;
    //           _nbElmPerColor[iC] = _nbElmPerColor[iC] + 1;
    //           _listElmPerColor[iC].push_back(iElm);
    //           break;
    //         }
    //       }
    //     }
    //   }
    //   break;
    // }
    // case 3:
    // {
    //   _nbColor = *max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;
    //   _colorElm.resize(_nElm, -1);
    //   _nbElmPerColor.resize(_nbColor);
    //   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

    //   for(int iElm = 0; iElm < _nElm; iElm++) {
    //     if(_colorElm[iElm] == -1) {
    //       std::vector<int> iPatchElm = _listElmPerElm[iElm];
    //       std::vector<bool> availableColor(_nbColor, true);
    //       for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
    //         if(_colorElm[iPatchElm[i]] != -1 && availableColor[_colorElm[iPatchElm[i]]] == true) {
    //           availableColor[_colorElm[iPatchElm[i]]] = false;
    //         }
    //       }

    //       int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor);
    //       _colorElm[iElm] = activColor;
    //       _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
    //       _listElmPerColor[activColor].push_back(iElm);
    //     }
    //   }
    //   break;
    // }
    // case 4:
    // {
    //   _nbColor = *max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;
    //   _colorElm.resize(_nElm, -1);
    //   _nbElmPerColor.resize(_nbColor);
    //   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

    //   bool noColor = false;
    //   while(noColor == false) {
    //     for(int iElm = 0; iElm < _nElm; iElm++) {
    //       if(_colorElm[iElm] == -1) {
    //         std::vector<int> iPatchElm = _listElmPerElm[iElm];
    //         std::vector<bool> availableColor(_nbColor, true);
    //         for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
    //           if(_colorElm[iPatchElm[i]] < 0) {
    //             _colorElm[iPatchElm[i]] = -2;
    //           } else if(availableColor[_colorElm[iPatchElm[i]]] == true) {
    //             availableColor[_colorElm[iPatchElm[i]]] = false;
    //           }
    //         }

    //         int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor);
    //         _colorElm[iElm] = activColor;
    //         _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
    //         _listElmPerColor[activColor].push_back(iElm);
    //       }
    //     }

    //     noColor = true;
    //     for(int iElm = 0; iElm < _nElm; iElm++) {
    //       if(_colorElm[iElm] == -2) {
    //         _colorElm[iElm] = -1;
    //         noColor = false;
    //       }
    //     }
    //     // feInfo("%d",_nbColor);
    //   }
    //   break;
    // }
    // case 5:
    // {
    //   _nbColor = 16;
    //   _colorElm.resize(_nElm, -1);
    //   _nbElmPerColor.resize(_nbColor);
    //   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

    //   bool noColor = false;
    //   while(noColor == false) {
    //     for(int k = 0; k < _nElm; ++k) {
    //       if(_colorElm[k] == -1) {
    //         std::vector<bool> availableColor(_nbColor, true);
    //         for(int i = 0; i < _nNodPerElm; ++i) {
    //           int s = _connecNodes[k * _nNodPerElm + i];
    //           for(int j = 0; j < _nbElmPerNode[s]; ++j) {
    //             int numElm = _listElmPerNode[s][j];
    //             if(numElm != k) {
    //               if(_colorElm[numElm] < 0) {
    //                 _colorElm[numElm] = -2;
    //               } else if(availableColor[_colorElm[numElm]] == true) {
    //                 availableColor[_colorElm[numElm]] = false;
    //               }
    //             }
    //           }
    //         }
    //         int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor);
    //         _colorElm[k] = activColor;
    //         _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
    //         _listElmPerColor[activColor].push_back(k);
    //       }
    //     }

    //     noColor = true;
    //     for(int k = 0; k < _nElm; ++k) {
    //       if(_colorElm[k] == -2) {
    //         _colorElm[k] = -1;
    //         noColor = false;
    //       }
    //     }

    //   } // while
    //   break;
    // }
  }
}

// void feCncGeo::color2(int cncGeoTag)
// {
//   // feCncGeo *cnc = _cncGeo[cncGeoTag];
//   // int nbElm = cnc->getNbElm();
//   // int _nNodPerElm = cnc->getNbNodePerElem();
//   // std::vector<int> &cncNodes = cnc->getNodeConnectivityRef();

//   // patchNode(nbElm, _nNodPerElm, cncNodes);
//   // patchElm(nbElm, _nNodPerElm, cncNodes);

//   // _nbColor = *max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;
//   // _colorElm.resize(nbElm, -1);
//   // _nbElmPerColor.resize(_nbColor);
//   // _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));
//   // std::vector<int> color(_nbColor);
//   // for(int i = 0; i < _nbColor; ++i) {
//   //   color[i] = i;
//   // }

//   // for(int iElm = 0; iElm < nbElm; iElm++) {
//   //   if(_colorElm[iElm] == -1) {
//   //     std::vector<int> iElmPatch = _listElmPerElm[iElm];
//   //     std::vector<bool> availableColor(_nbColor, true);
//   //     for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
//   //       if(_colorElm[iElmPatch[i]] != -1 && availableColor[_colorElm[iElmPatch[i]]] == true) {
//   //         availableColor[_colorElm[iElmPatch[i]]] = false;
//   //       }
//   //     }

//   //     for(int iC = 0; iC < _nbColor; iC++) {
//   //       if(availableColor[iC] == true) {
//   //         _colorElm[iElm] = iC;
//   //         _nbElmPerColor[iC] = _nbElmPerColor[iC] + 1;
//   //         _listElmPerColor[iC].push_back(iElm);
//   //         break;
//   //       }
//   //     }
//   //   }
//   // }

//   // for(int i=0;i<_nbColor;i++){
//   //   feInfo("%d Elm pour couleur %d",_nbElmPerColor[i],i);
//   // }

//   // feInfo("%d",_nbColor);
// }

// void feCncGeo::color3(int cncGeoTag)
// {
//   feCncGeo *cnc = _cncGeo[cncGeoTag];
//   int nbElm = cnc->getNbElm();
//   int _nNodPerElm = cnc->getNbNodePerElem();
//   std::vector<int> &cncNodes = cnc->getNodeConnectivityRef();

//   patchNode(nbElm, _nNodPerElm, cncNodes);
//   patchElm(nbElm, _nNodPerElm, cncNodes);

//   _nbColor = *max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;
//   _colorElm.resize(nbElm, -1);
//   _nbElmPerColor.resize(_nbColor);
//   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

//   for(int iElm = 0; iElm < nbElm; iElm++) {
//     if(_colorElm[iElm] == -1) {
//       std::vector<int> iPatchElm = _listElmPerElm[iElm];
//       std::vector<bool> availableColor(_nbColor, true);
//       for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
//         if(_colorElm[iPatchElm[i]] != -1 && availableColor[_colorElm[iPatchElm[i]]] == true) {
//           availableColor[_colorElm[iPatchElm[i]]] = false;
//         }
//       }

//       int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor);
//       _colorElm[iElm] = activColor;
//       _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
//       _listElmPerColor[activColor].push_back(iElm);
//     }
//   }

//   // for(int i=0;i<_nbColor;i++){
//   //   feInfo("%d Elm pour couleur %d",_nbElmPerColor[i],i);
//   // }

//   // feInfo("%d",_nbColor);
// }

// void feCncGeo::color4(int cncGeoTag)
// {
//   feCncGeo *cnc = _cncGeo[cncGeoTag];
//   int nbElm = cnc->getNbElm();
//   int _nNodPerElm = cnc->getNbNodePerElem();
//   std::vector<int> &cncNodes = cnc->getNodeConnectivityRef();

//   patchNode(nbElm, _nNodPerElm, cncNodes);
//   patchElm(nbElm, _nNodPerElm, cncNodes);

//   _nbColor = *max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;
//   _colorElm.resize(nbElm, -1);
//   _nbElmPerColor.resize(_nbColor);
//   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

//   bool noColor = false;
//   while(noColor == false) {
//     for(int iElm = 0; iElm < nbElm; iElm++) {
//       if(_colorElm[iElm] == -1) {
//         std::vector<int> iPatchElm = _listElmPerElm[iElm];
//         std::vector<bool> availableColor(_nbColor, true);
//         for(int i = 0; i < _nbElmPerElm[iElm]; i++) {
//           if(_colorElm[iPatchElm[i]] < 0) {
//             _colorElm[iPatchElm[i]] = -2;
//           } else if(availableColor[_colorElm[iPatchElm[i]]] == true) {
//             availableColor[_colorElm[iPatchElm[i]]] = false;
//           }
//         }

//         int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor);
//         _colorElm[iElm] = activColor;
//         _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
//         _listElmPerColor[activColor].push_back(iElm);
//       }
//     }

//     noColor = true;
//     for(int iElm = 0; iElm < nbElm; iElm++) {
//       if(_colorElm[iElm] == -2) {
//         _colorElm[iElm] = -1;
//         noColor = false;
//       }
//     }
//     // feInfo("%d",_nbColor);
//   }
//   // for(int i=0;i<_nbColor;i++){
//   //   feInfo("%d Elm pour couleur %d",_nbElmPerColor[i],i);
//   // }
// }

// void feCncGeo::color5(int cncGeoTag)
// {
//   feCncGeo *cnc = _cncGeo[cncGeoTag];
//   int nbElm = cnc->getNbElm();
//   int _nNodPerElm = cnc->getNbNodePerElem();
//   std::vector<int> &cncNodes = cnc->getNodeConnectivityRef();

//   patchNode(nbElm, _nNodPerElm, cncNodes);

//   _nbColor = 16;
//   _colorElm.resize(nbElm, -1);
//   _nbElmPerColor.resize(_nbColor);
//   _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

//   bool noColor = false;
//   while(noColor == false) {
//     for(int k = 0; k < nbElm; ++k) {
//       if(_colorElm[k] == -1) {
//         std::vector<bool> availableColor(_nbColor, true);
//         for(int i = 0; i < _nNodPerElm; ++i) {
//           int s = cncNodes[k * _nNodPerElm + i];
//           for(int j = 0; j < _nbElmPerNode[s]; ++j) {
//             int numElm = _listElmPerNode[s][j];
//             if(numElm != k) {
//               if(_colorElm[numElm] < 0) {
//                 _colorElm[numElm] = -2;
//               } else if(availableColor[_colorElm[numElm]] == true) {
//                 availableColor[_colorElm[numElm]] = false;
//               }
//             }
//           }
//         }
//         int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor);
//         _colorElm[k] = activColor;
//         _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
//         _listElmPerColor[activColor].push_back(k);
//       }
//     }

//     noColor = true;
//     for(int k = 0; k < nbElm; ++k) {
//       if(_colorElm[k] == -2) {
//         _colorElm[k] = -1;
//         noColor = false;
//       }
//     }

//   } // while
//   // for(int i=0;i<_nbColor;i++){
//   //   feInfo("%d Elm pour couleur %d",_nbElmPerColor[i],i);
//   // }
// }