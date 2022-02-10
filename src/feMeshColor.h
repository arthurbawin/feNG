#ifndef _FEMESHCOLOR_
#define _FEMESHCOLOR_


#include "feMesh.h"
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <set>
//....




class feMeshColor
{
protected:
	int _dim;
	int _nElm;
    int _nNodePerElem;
    int _nEdgePerElem;
    int _nbColor;

	

    std::vector<int> _connecElem;
    std::vector<int> _connecNodes;
    std::vector<int> _connecEdges;

    std::vector<int> _CouleurElm;
	std::vector<double> _coord;

public :
	feMeshColor(int nNod = 0, int dim = 0, int nCncGeo = 0, std::string ID = "")
    : _ID(ID), _dim(dim), _nNod(nNod), _nEdg(0), _nCncGeo(nCncGeo){};


}

