#include "feMesh.h"
#include "feElement.h"
#include "feTriangle.h"
#include "feLine.h"
#include "feSpaceTriangle.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "SBoundingBox3d.h"

extern int FE_VERBOSE;

// Number of nodes for each type of Gmsh elements, type is the index of the array + 1
int nodes_of_gmsh_element[] = {
  2, // 2-node line.
  3, // 3-node triangle.
  4, // 4-node quadrangle.
  4, // 4-node tetrahedron.
  8, // 8-node hexahedron.
  6, // 6-node prism.
  5, // 5-node pyramid.
  3, /* 3-node second order line (2 nodes associated with the
        vertices and 1 with the edge). */
  6, /* 6-node second order triangle (3 nodes associated with the
        vertices and 3 with the edges). */
  9, /* 9-node second order quadrangle (4 nodes associated with the
        vertices, 4 with the edges and 1 with the face). */
  10, /* 10-node second order tetrahedron (4 nodes associated with
         the vertices and 6 with the edges). */
  27, /* 27-node second order hexahedron (8 nodes associated with the
         vertices, 12 with the edges, 6 with the faces and 1 with the
         volume). */
  18, /* 18-node second order prism (6 nodes associated with the
         vertices, 9 with the edges and 3 with the quadrangular
         faces). */
  14, /* 14-node second order pyramid (5 nodes associated with the
         vertices, 8 with the edges and 1 with the quadrangular
         face). */
  1, // 1-node point.
  8, /* 8-node second order quadrangle (4 nodes associated with the
        vertices and 4 with the edges). */
  20, /* 20-node second order hexahedron (8 nodes associated with the
         vertices and 12 with the edges). */
  15, /* 15-node second order prism (6 nodes associated with the
         vertices and 9 with the edges). */
  13, /* 13-node second order pyramid (5 nodes associated with the
         vertices and 8 with the edges). */
  9, /* 9-node third order incomplete triangle (3 nodes associated
        with the vertices, 6 with the edges) */
  10, /* 10-node third order triangle (3 nodes associated with the
         vertices, 6 with the edges, 1 with the face) */
  12, /* 12-node fourth order incomplete triangle (3 nodes associated
         with the vertices, 9 with the edges) */
  15, /* 15-node fourth order triangle (3 nodes associated with the
         vertices, 9 with the edges, 3 with the face) */
  15, /* 15-node fifth order incomplete triangle (3 nodes associated
         with the vertices, 12 with the edges) */
  21, /* 21-node fifth order complete triangle (3 nodes associated
         with the vertices, 12 with the edges, 6 with the face) */
  4, /* 4-node third order edge (2 nodes associated with the
        vertices, 2 internal to the edge) */
  5, /* 5-node fourth order edge (2 nodes associated with the
        vertices, 3 internal to the edge) */
  6, /* 6-node fifth order edge (2 nodes associated with the
        vertices, 4 internal to the edge) */
  20, /* 20-node third order tetrahedron (4 nodes associated with the
         vertices, 12 with the edges, 4 with the faces) */
  35, /* 35-node fourth order tetrahedron (4 nodes associated with
         the vertices, 18 with the edges, 12 with the faces, and 1
         with the volume) */
  56, /* 56-node fifth order tetrahedron (4 nodes associated with the
         vertices, 24 with the edges, 24 with the faces, and 4 with
         the volume) */
  -1,   -1, /* unsupported tetrahedral types */
  -1,   -1, /* unsupported polygonal and polyhedral types */
  16, /* 16-node third order quadrilateral (4 nodes associated with
         the vertices, 8 with the edges, 4 wth the face) */
  25, /* 25-node fourth order quadrilateral (4 nodes associated with
         the vertices, 12 with the edges, 9 wth the face) */
  36, /* 36-node fifth order quadrilateral (4 nodes associated with
         the vertices, 16 with the edges, 16 wth the face) */
  -1,   -1, -1, /* unsupported quadrilateral types */
  28, /* 28-node sixth order complete triangle (3 nodes associated
         with the vertices, 15 with the edges, 10 with the face) */
  36, /* 36-node seventh order complete triangle (3 nodes associated
         with the vertices, 18 with the edges, 15 with the face) */
  45, /* 45-node eighth order complete triangle (3 nodes associated
         with the vertices, 21 with the edges, 21 with the face) */
  55, /* 55-node ninth order complete triangle (3 nodes associated
         with the vertices, 24 with the edges, 28 with the face) */
  66, /* 66-node tenth order complete triangle (3 nodes associated
         with the vertices, 27 with the edges, 36 with the face) */
  49, /* 49-node sixth order quadrilateral (4 nodes associated with
         the vertices, 20 with the edges, 25 wth the face) */
  64, /* 64-node seventh order quadrilateral (4 nodes associated with
         the vertices, 24 with the edges, 36 wth the face) */
  81, /* 81-node eighth order quadrilateral (4 nodes associated with
         the vertices, 28 with the edges, 49 wth the face) */
  100, /* 100-node ninth order quadrilateral (4 nodes associated with
          the vertices, 32 with the edges, 64 wth the face) */
  121, /* 121-node tenth order quadrilateral (4 nodes associated with
          the vertices, 36 with the edges, 81 wth the face) */
  -1,   -1, -1, -1, -1, /* unsupported triangular types */
  -1,   -1, -1, -1, -1, /* unsupported quadrilateral types */
  7, /* 7-node sixth order edge (2 nodes associated with the
        vertices, 5 internal to the edge) */
  8, /* 8-node seventh order edge (2 nodes associated with the
        vertices, 6 internal to the edge) */
  9, /* 9-node eighth order edge (2 nodes associated with the
        vertices, 7 internal to the edge) */
  10, /* 10-node ninth order edge (2 nodes associated with the
         vertices, 8 internal to the edge) */
  11, /* 11-node tenth order edge (2 nodes associated with the
         vertices, 9 internal to the edge) */
  -1, /* unsupported linear types */
  -1,   -1, -1, /* unsupported types */
  84, /* 84-node sixth order tetrahedron (4 nodes associated with the
         vertices, 30 with the edges, 40 with the faces, and 10 with
         the volume) */
  120, /* 120-node seventh order tetrahedron (4 nodes associated with
          the vertices, 36 with the edges, 60 with the faces, and 20
          with the volume) */
  165, /* 165-node eighth order tetrahedron (4 nodes associated with
          the vertices, 42 with the edges, 84 with the faces, and 35
          with the volume) */
  220, /* 220-node ninth order tetrahedron (4 nodes associated with
          the vertices, 48 with the edges, 112 with the faces, and 56
          with the volume) */
  286, /* 286-node tenth order tetrahedron (4 nodes associated with
          the vertices, 54 with the edges, 144 with the faces, and 84
          with the volume) */
  -1,   -1, -1, /* undefined types */
  -1,   -1, -1, -1, -1, /* unsupported tetrahedral types */
  -1,   -1, -1, -1, -1, -1, /* unsupported types */
  40, /* 40-node third order prism (6 nodes associated with the
         vertices, 18 with the edges, 14 with the faces, and 2 with
         the volume) */
  75, /* 75-node fourth order prism (6 nodes associated with the
         vertices, 27 with the edges, 33 with the faces, and 9 with
         the volume) */
  64, /* 64-node third order hexahedron (8 nodes associated with the
         vertices, 24 with the edges, 24 with the faces and 8 with
         the volume).*/
  125, /* 125-node fourth order hexahedron (8 nodes associated with
          the vertices, 36 with the edges, 54 with the faces and 27
          with the volume).*/
  216, /* 216-node fifth order hexahedron (8 nodes associated with the
          vertices, 48 with the edges, 96 with the faces and 64 with
          the volume).*/
  343, /* 343-node sixth order hexahedron (8 nodes associated with the
          vertices, 60 with the edges, 150 with the faces and 125 with
          the volume).*/
  512, /* 512-node seventh order hexahedron (8 nodes associated with
          the vertices, 72 with the edges, 216 with the faces and 216
          with the volume).*/
  729, /* 729-node eighth order hexahedron (8 nodes associated with
          the vertices, 84 with the edges, 294 with the faces and 343
          with the volume).*/
  1000, /* 1000-node ninth order hexahedron (8 nodes associated with
           the vertices, 96 with the edges, 384 with the faces and 512
           with the volume).*/
  -1,   -1, -1, -1, -1, -1, -1, /* unsupported hexahedron types */
  126, /* 126-node fifth order prism (6 nodes associated with the
          vertices, 36 with the edges, 60 with the faces, and 24 with
          the volume) */
  196, /* 196-node sixth order prism (6 nodes associated with the
          vertices, 45 with the edges, 95 with the faces, and 50 with
          the volume) */
  288, /* 288-node seventh order prism (6 nodes associated with the
          vertices, 54 with the edges, 138 with the faces, and 90 with
          the volume) */
  405, /* 405-node eighth order prism (6 nodes associated with the
          vertices, 63 with the edges, 189 with the faces, and 147
          with the volume) */
  550, /* 550-node ninth order prism (6 nodes associated with the
          vertices, 72 with the edges, 248 with the faces, and 224
          with the volume) */
  -1,   -1, -1, -1, -1, -1, -1, /* unsupported prism types */
  30, /* 30-node third order pyramid (5 nodes associated with the
         vertices, 16 with the edges and 8 with the faces, and 1 with
         the volume). */
  55, /* 55-node fourth order pyramid (5 nodes associated with the
         vertices, 24 with the edges and 21 with the faces, and 5
         with the volume). */
  91, /* 91-node fifth order pyramid (5 nodes associated with the
         vertices, 32 with the edges and 40 with the faces, and 14
         with the volume). */
  140, /* 140-node sixth order pyramid (5 nodes associated with the
          vertices, 40 with the edges and 65 with the faces, and 30
          with the volume). */
  204, /* 204-node seventh order pyramid (5 nodes associated with the
          vertices, 48 with the edges and 96 with the faces, and 55
          with the volume). */
  285, /* 285-node eighth order pyramid (5 nodes associated with the
          vertices, 56 with the edges and 133 with the faces, and 91
          with the volume). */
  385 /* 385-node ninth order pyramid (5 nodes associated with the
         vertices, 64 with the edges and 176 with the faces, and 140
         with the volume). */
};

// Dimension of the first elements of the list :
int dim_of_gmsh_element[] = {
  1, // 2-node line.
  2, // 3-node triangle.
  2, // 4-node quadrangle.
  3, // 4-node tetrahedron.
  3, // 8-node hexahedron.
  3, // 6-node prism.
  3, // 5-node pyramid.
  1, /* 3-node second order line (2 nodes associated with the
        vertices and 1 with the edge). */
  2, /* 6-node second order triangle (3 nodes associated with the
        vertices and 3 with the edges). */
  2, /* 9-node second order quadrangle (4 nodes associated with the
        vertices, 4 with the edges and 1 with the face). */
  3, /* 10-node second order tetrahedron (4 nodes associated with
         the vertices and 6 with the edges). */
  3, /* 27-node second order hexahedron (8 nodes associated with the
         vertices, 12 with the edges, 6 with the faces and 1 with the
         volume). */
  3, /* 18-node second order prism (6 nodes associated with the
         vertices, 9 with the edges and 3 with the quadrangular
         faces). */
  3, /* 14-node second order pyramid (5 nodes associated with the
         vertices, 8 with the edges and 1 with the quadrangular
         face). */
  0, // 1-node point.
  2, /* 8-node second order quadrangle (4 nodes associated with the
        vertices and 4 with the edges). */
  3, /* 20-node second order hexahedron (8 nodes associated with the
         vertices and 12 with the edges). */
  3, /* 15-node second order prism (6 nodes associated with the
         vertices and 9 with the edges). */
  3, /* 13-node second order pyramid (5 nodes associated with the
         vertices and 8 with the edges). */
  2, /* 9-node third order incomplete triangle (3 nodes associated
        with the vertices, 6 with the edges) */
  2, /* 10-node third order triangle (3 nodes associated with the
         vertices, 6 with the edges, 1 with the face) */
  2, /* 12-node fourth order incomplete triangle (3 nodes associated
         with the vertices, 9 with the edges) */
  2, /* 15-node fourth order triangle (3 nodes associated with the
         vertices, 9 with the edges, 3 with the face) */
  2, /* 15-node fifth order incomplete triangle (3 nodes associated
         with the vertices, 12 with the edges) */
  2, /* 21-node fifth order complete triangle (3 nodes associated
         with the vertices, 12 with the edges, 6 with the face) */
  1, /* 4-node third order edge (2 nodes associated with the
        vertices, 2 internal to the edge) */
  1, /* 5-node fourth order edge (2 nodes associated with the
        vertices, 3 internal to the edge) */
  1, /* 6-node fifth order edge (2 nodes associated with the
        vertices, 4 internal to the edge) */
  3, /* 20-node third order tetrahedron (4 nodes associated with the
         vertices, 12 with the edges, 4 with the faces) */
  3, /* 35-node fourth order tetrahedron (4 nodes associated with
         the vertices, 18 with the edges, 12 with the faces, and 1
         with the volume) */
  3, /* 56-node fifth order tetrahedron (4 nodes associated with the
         vertices, 24 with the edges, 24 with the faces, and 4 with
         the volume) */
  // TODO : Complete the table (-:
};

feStatus feMesh2DP1::readMsh2(std::istream &input, bool curved, mapType physicalEntitiesDescription)
{
  std::string buffer;
  int ph1; // Placeholder
  // std::map<int, int> _verticesMap; // Gmsh tag (which may include gaps) to sequential tag. Just
  // in case : not sure there are gaps in msh2...
  std::vector<int> numNodesInBlock;

  if(physicalEntitiesDescription.size() > 0) {
    for(auto pair : physicalEntitiesDescription) {
      // pair = { {dim,tag}, name }
      physicalEntity pE;
      pE.dim = pair.first.first;
      pE.tag = pair.first.second;
      pE.name = pair.second;
      pE.nElm = 0;
      pE.nNodePerElem = -1;
      pE.nEdgePerElem = -1;
      pE.cncID = "";
      _physicalEntities[pair.first] = pE;
    }
  }

  while(input >> buffer) {
    if(buffer == "$PhysicalNames") { // Read physical entities
      input >> _nPhysicalEntities;
      getline(input, buffer);
      for(int i = 0; i < _nPhysicalEntities; ++i) {
        physicalEntity pE;
        // dimension(ASCII int) - physicalTag(ASCII int) - "name"(127 characters max)
        input >> pE.dim >> pE.tag >> pE.name;
        if(pE.name == "") {
          return feErrorMsg(FE_STATUS_ERROR, "Empty name for physical entity.");
        }
        pE.name.pop_back();
        pE.name.erase(pE.name.begin());
        pE.nElm = 0;
        pE.nNodePerElem = -1;
        pE.nEdgePerElem = -1;
        pE.cncID = "";
        _physicalEntities[{pE.dim, pE.tag}] = pE;
      }

      // Save the description of the physical entities : this is used when looping and adapting the
      // mesh. Since mmg does not keep track of the name of the physical entities, they have to be
      // saved and loaded during the following adaptation cycles.
      // for(int i = 0; i < _nPhysicalEntities; ++i) {
      //   physicalEntity pE = _physicalEntities[i];
      //   std::pair<int, int> p = {pE.dim, pE.tag};
      //   _physicalEntitiesDescription[p] = pE.name;
      // }
      for(auto pair : _physicalEntities) {
        _physicalEntitiesDescription[pair.first] = _physicalEntities[pair.first].name;
      }
    }

    // There is no Entities block in msh2 : entities information is in the elements

    else if(buffer == "$Nodes") { // Read nodes
      int gmshNumber, countVertex = 0;
      double x, y, z;
      input >> _nNod;
      _vertices.resize(_nNod);
      getline(input, buffer);
      for(int iNode = 0; iNode < _nNod; ++iNode) {
        // node-number x-coord y-coord z-coord
        input >> gmshNumber >> x >> y >> z;
        _verticesMap[gmshNumber] = countVertex;
        _vertices[countVertex++] = Vertex(x, y, z, gmshNumber);
      }

      if((int)_verticesMap.size() != _nNod) {
        return feErrorMsg(FE_STATUS_ERROR, "Vertices indices are not unique.");
      }
    }

    else if(buffer == "$Elements") { // Read elements

      if(_physicalEntities.size() == 0) {
        return feErrorMsg(FE_STATUS_ERROR, "No physical entity defined in the mesh.");
      }

      int tagCount = 1, nEdges = 1; // To number the geometric entities
      int gmshElemTag, elemType, numEntities, physicalTag, geometricTag;
      getline(input, buffer);
      // number-of-elements
      input >> _nElm;
      _nTotalElm = _nElm;

      int serialNumberElm = 0; // To number the elements when skipping the Points
      for(int iElm = 0; iElm < _nElm; ++iElm) {
        getline(input, buffer);
        // elm-number elm-type number-of-tags < tag > ??? node-number-list
        input >> gmshElemTag >> elemType >> numEntities >> physicalTag;
        // printf("Read %d %d %d %d\n", gmshElemTag, elemType, numEntities, physicalTag);
        // Physical and geometric entities are not separated, and physical seem to be given first.
        // Since we only accept non-overlapping physical domains (i.e. each element belongs to a
        // single physical entity), let's just assume the first number is the physical domain and
        // the following are the geometric entities.

        int entityDim = dim_of_gmsh_element[elemType - 1];

        // Physical tag (number of the physical entity) should be strictly positive.
        // Meshes produced by mmg2d include points with physicalTag = 0, which we skip.
        if(physicalTag <= 0) continue;

        // numEntities should be 2 since an element belongs to max 1 physical entity
        // and should belong to a single geometric entity (?) : no need to loop
        // for(int i = 1; i < numEntities; ++i){
        input >> geometricTag;

        std::pair<int, int> p = {entityDim, geometricTag};

        // Meshes produced by mmg2d seem to include unwanted 0D physical entities (?).
        // If a description of the original physical entities is provided, we make sure
        // the read physical entities is indeed part of these entities.
        if(physicalEntitiesDescription.size() > 0) {
          bool isPartOfTheOriginalEntities = false;
          for(auto pair : physicalEntitiesDescription) {
            // pair = { {dim,tag}, name }
            // printf("Comparing entity dim = %d - tag = %d with stored entity %s - dim = %d - tag =
            // %d\n", entityDim, physicalTag, pair.second.c_str(), pair.first.first,
            // pair.first.second);
            if(entityDim == pair.first.first && physicalTag == pair.first.second)
              isPartOfTheOriginalEntities = true;
          }
          if(!isPartOfTheOriginalEntities) {
            continue;
          }
        }

        // If geometric entity doesn't exist, create it :
        if(_entities.find(p) == _entities.end()) {
          entity e;
          e.dim = entityDim;
          e.tag = tagCount++;
          e.tagGmsh = geometricTag;
          e.numPhysicalTags = 1;
          e.physicalTags.push_back(physicalTag);
          e.nElm = 0;
          // _physicalEntities[physicalTag - 1].listEntities.push_back(e.tag);
          std::pair<int, int> p = {e.dim, physicalTag};
          _physicalEntities[p].listEntities.push_back(e.tag);
          _entities[p] = e;
        }

        _entities[p].nElm++;
        // _entities[p].connecElem.push_back(gmshElemTag-1); // Gmsh elem tag starts at 1
        // The Gmsh elem tag takes the Points into account : trying with a counter
        _entities[p].connecElem.push_back(serialNumberElm++); // Gmsh elem tag starts at 1

        int nNodePerElem = nodes_of_gmsh_element[elemType - 1];
        std::vector<int> elemNodesGmsh(nNodePerElem);
        std::vector<int> elemNodes(nNodePerElem);

        std::map<int, int>::iterator it;
        for(int j = 0; j < nNodePerElem; ++j) {
          input >> ph1;
          // Verification
          it = _verticesMap.find(ph1);
          if(it == _verticesMap.end()) {
            return feErrorMsg(FE_STATUS_ERROR, "Element node does not exist /-:");
          }
          elemNodesGmsh[j] = ph1; // The Gmsh number stored in ph1 is used to create the edges
          elemNodes[j] =
            it->second; // The node number (0...nNode) is used to create the connectivity
          _entities[p].connecNodes.push_back(elemNodes[j]);
        }

        // Fill geometric entity data :
        switch(elemType) {
          case 1:
            [[gnu::fallthrough]]; //  2-node line
          case 8:
            if(curved) {
              _entities[p].cncID = "LineP2";
            }
            [[gnu::fallthrough]]; //  3-node line (2nd order)
          case 26:
            if(curved) {
              _entities[p].cncID = "LineP3";
            }
            [[gnu::fallthrough]]; //  4-node line (3rd order)
          case 27:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  5-node line (4th order)
          case 28:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  6-node line (5th order)
          case 62:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  7-node line (6th order)
          case 63:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  8-node line (7th order)
          case 64:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  9-node line (8th order)
          case 65:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; // 10-node line (9th order)
          case 66:
            if(curved) { /* Not supported */
            } // 11-node line (10th order)
            {
              // Default is linear interpolation for the geometry
              if(!curved) {
                _entities[p].cncID = "LineP1";
              }
              _entities[p].nNodePerElem = (curved) ? nNodePerElem : 2;
              _entities[p].nEdgePerElem = 0;
              break;
            }
          case 2:
            [[gnu::fallthrough]]; //  3-node triangle
          case 9:
            if(curved) {
              _entities[p].cncID = "TriP2";
            }
            [[gnu::fallthrough]]; //  6-node triangle (2nd order)
          case 21:
            if(curved) {
              _entities[p].cncID = "TriP3";
            }
            [[gnu::fallthrough]]; // 10-node triangle (3rd order)
          case 23:
            [[gnu::fallthrough]]; // 15-node triangle (4th order)
          case 25:
            [[gnu::fallthrough]]; // 21-node triangle (5th order)
          case 42:
            [[gnu::fallthrough]]; // 28-node triangle (6th order)
          case 43:
            [[gnu::fallthrough]]; // 36-node triangle (7th order)
          case 44:
            [[gnu::fallthrough]]; // 45-node triangle (8th order)
          case 45:
            [[gnu::fallthrough]]; // 55-node triangle (9th order)
          case 46: // 66-node triangle (10th order)
          {
            if(!curved) {
              _entities[p].cncID = "TriP1";
            }
            _entities[p].nNodePerElem = (curved) ? nNodePerElem : 3;
            _entities[p].nEdgePerElem = 3;
            // Construct the triangle edges :
            Vertex *v0, *v1;
            for(int k = 0; k < 3; ++k) {
              if(k == 2) {
                v0 = &_vertices[_verticesMap[elemNodesGmsh[2]]];
                v1 = &_vertices[_verticesMap[elemNodesGmsh[0]]];
              } else {
                v0 = &_vertices[_verticesMap[elemNodesGmsh[k]]];
                v1 = &_vertices[_verticesMap[elemNodesGmsh[k + 1]]];
              }
              std::pair<std::set<Edge, EdgeLessThan>::iterator, bool> ret;
              Edge e(v0, v1, nEdges, _entities[p].physicalTags[0]);
              ret = _edges.insert(e);
              if(ret.second) {
                // Edge was added to the set : nEdges is added to connecEdges
                _entities[p].connecEdges.push_back(nEdges++);
              } else {
                // Edge is already in the set : the negative is added to connecEdges
                // Assumes an edge is shared by only two triangles in 2D
                // More tests required in 3D, where an edge can be shared by N tets
                _entities[p].connecEdges.push_back(-ret.first->getTag());
              }
            }
            break;
          }
          case 3:
            [[gnu::fallthrough]]; //   4-node quadrangle
          case 10:
            [[gnu::fallthrough]]; //   9-node quadrangle (2nd order)
          case 36:
            [[gnu::fallthrough]]; //  16-node quadrangle (3rd order)
          case 37:
            [[gnu::fallthrough]]; //  25-node quadrangle (4th order)
          case 38:
            [[gnu::fallthrough]]; //  36-node quadrangle (5th order)
          case 47:
            [[gnu::fallthrough]]; //  49-node quadrangle (6th order)
          case 48:
            [[gnu::fallthrough]]; //  64-node quadrangle (7th order)
          case 49:
            [[gnu::fallthrough]]; //  81-node quadrangle (8th order)
          case 50:
            [[gnu::fallthrough]]; // 100-node quadrangle (9th order)
          case 51: // 121-node quadrangle (10th order)
          {
            return feErrorMsg(FE_STATUS_ERROR,
                              "Interpolant pas (encore) pris en charge pour la g??om??trie "
                              "de l'entit?? (dim = %d, tag = %d) (quad).\n",
                              entityDim, geometricTag);
          }
          case 4:
            [[gnu::fallthrough]]; //   4-node tetrahedron
          case 11:
            [[gnu::fallthrough]]; //  10-node tetrahedron (2nd order)
          case 29:
            [[gnu::fallthrough]]; //  20-node tetrahedron (3rd order)
          case 30:
            [[gnu::fallthrough]]; //  35-node tetrahedron (4th order)
          case 31:
            [[gnu::fallthrough]]; //  56-node tetrahedron (5th order)
          case 71:
            [[gnu::fallthrough]]; //  84-node tetrahedron (6th order)
          case 72:
            [[gnu::fallthrough]]; // 120-node tetrahedron (7th order)
          case 73:
            [[gnu::fallthrough]]; // 165-node tetrahedron (8th order)
          case 74:
            [[gnu::fallthrough]]; // 220-node tetrahedron (9th order)
          case 75: // 286-node tetrahedron (10th order)
          {
            return feErrorMsg(FE_STATUS_ERROR,
                              "Interpolant pas (encore) pris en charge pour la g??om??trie "
                              "de l'entit?? (dim = %d, tag = %d) (tet).\n",
                              entityDim, geometricTag);
          }
          case 5:
            [[gnu::fallthrough]]; //    8-node hexahedron
          case 12:
            [[gnu::fallthrough]]; //   27-node hexahedron (2nd order)
          case 92:
            [[gnu::fallthrough]]; //   64-node hexahedron (3rd order)
          case 93:
            [[gnu::fallthrough]]; //  125-node hexahedron (4th order)
          case 94:
            [[gnu::fallthrough]]; //  216-node hexahedron (5th order)
          case 95:
            [[gnu::fallthrough]]; //  343-node hexahedron (6th order)
          case 96:
            [[gnu::fallthrough]]; //  512-node hexahedron (7th order)
          case 97:
            [[gnu::fallthrough]]; //  729-node hexahedron (8th order)
          case 98: // 1000-node hexahedron (9th order)
          {
            return feErrorMsg(FE_STATUS_ERROR,
                              "Interpolant pas (encore) pris en charge pour la g??om??trie "
                              "de l'entit?? (dim = %d, tag = %d) (hex).\n",
                              entityDim, geometricTag);
          }
          case 15: // 1-node point
          {
            _entities[p].cncID = "Point0D";
            _entities[p].nNodePerElem = 1;
            _entities[p].nEdgePerElem = 0;
            break;
          }
          default:
            return feErrorMsg(FE_STATUS_ERROR, "Unsupported Gmsh element type.");
        } // switch(elemType)
      } // for iElm
    } // if buffer = "Elements"
  } // while input

  return FE_STATUS_OK;
}

feStatus feMesh2DP1::readMsh4(std::istream &input, bool curved, mapType physicalEntitiesDescription)
{
  std::string buffer;
  // Placeholders
  int ph1, ph2, ph3;
  double ph1D, ph2D, ph3D, ph4D, ph5D, ph6D;
  std::vector<int> numNodesInBlock;
  _nPhysicalEntities = 0;

  if(physicalEntitiesDescription.size() > 0) {
    for(auto pair : physicalEntitiesDescription) {
      // pair = { {dim,tag}, name }
      physicalEntity pE;
      pE.dim = pair.first.first;
      pE.tag = pair.first.second;
      pE.name = pair.second;
      pE.nElm = 0;
      pE.nNodePerElem = -1;
      pE.nEdgePerElem = -1;
      pE.cncID = "";
      _physicalEntities[pair.first] = pE;
      feInfoCond(FE_VERBOSE > 1, "Created physical entity \"%s\"", pE.name.c_str());
    }
  }

  while(input >> buffer) {
    if(buffer == "$PhysicalNames") { // Read physical entities
      input >> _nPhysicalEntities;
      getline(input, buffer);
      // The tag for physical entities does not necessarily start at 0
      int sequentialTagForPhysicals = 0;
      for(int i = 0; i < _nPhysicalEntities; ++i) {
        physicalEntity pE;
        // dimension(ASCII int) - physicalTag(ASCII int) - "name"(127 characters max)
        input >> pE.dim >> pE.tag >> pE.name;
        if(pE.name == "") {
          return feErrorMsg(FE_STATUS_ERROR, "Empty name for physical entity.");
        }
        pE.name.pop_back();
        pE.name.erase(pE.name.begin());
        pE.nElm = 0;
        pE.nNodePerElem = -1;
        pE.nEdgePerElem = -1;
        pE.cncID = "";

        _physicalEntities[{pE.dim, pE.tag}] = pE;
        feInfoCond(FE_VERBOSE > 1,
                   "Created Physical entity \"%s\" with dimension %d and gmsh tag %d",
                   pE.name.c_str(), pE.dim, pE.tag);
      }

      // Save the description of the physical entities : this is used when looping and adapting the
      // mesh. Since mmg does not keep track of the name of the physical entities, they have to be
      // saved and loaded during the following adaptation cycles.
      // for(int i = 0; i < _nPhysicalEntities; ++i) {
      //   physicalEntity pE = _physicalEntities[i];
      //   std::pair<int, int> p = {pE.dim, pE.tag};
      //   _physicalEntitiesDescription[p] = pE.name;
      // }
      for(auto pair : _physicalEntities) {
        _physicalEntitiesDescription[pair.first] = _physicalEntities[pair.first].name;
      }
    }

    else if(buffer == "$Entities") { // Read geometric entities
      // numPoints(size_t) numCurves(size_t) numSurfaces(size_t) numVolumes(size_t)
      input >> _nPoints >> _nCurves >> _nSurfaces >> _nVolumes;

      // Pour le moment : pas de cncGeo sur les Points
      numNodesInBlock.resize(_nCurves + _nSurfaces + _nVolumes);

      int tagCount = 1;
      int tagUnnumberedEntities = 0;
      // Points
      for(int i = 0; i < _nPoints; ++i) {
        entity e;
        e.dim = 0;
        e.tag = tagCount++;
        // We only allow a Point geometric to have a single node
        e.nElm = 0;
        //  pointTag(int) X(double) Y(double) Z(double)
        //  numPhysicalTags(size_t) physicalTag(int)
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> e.numPhysicalTags;
        if(e.numPhysicalTags > 1) {
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
        }
        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          if(_physicalEntities.find({e.dim, ph1}) != _physicalEntities.end()) {
            _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
            // if(FE_VERBOSE) printf("In readGmsh4 : Geometric entity with dim = %d - tag = %d is
            // part of Physical
            // \"%s\"\n", e.dim, e.tagGmsh, _physicalEntities[{e.dim, ph1}].name.c_str());
          } else if(_physicalEntities.find({e.dim, -ph1}) !=
                    _physicalEntities.end()) { // Physical entity numbered in reverse
            _physicalEntities[{e.dim, -ph1}].listEntities.push_back(e.tag);
            if(FE_VERBOSE > -1) {
              // printf("In readGmsh4 : Geometric entity with dim = %d - tag = %d is part of
              // Physical \"%s\"\n", e.dim, e.tagGmsh, _physicalEntities[{e.dim,
              // -ph1}].name.c_str());
              feWarning("Physical entity \"%s\" is negative on geometric entity (%d,%d) and is "
                        "thus numbered in reverse order.",
                        _physicalEntities[{e.dim, -ph1}].name.c_str(), e.dim, e.tagGmsh);
            }
          }
        }
        _entities.insert({{e.dim, e.tagGmsh}, e});
      }
      // Curves
      for(int i = 0; i < _nCurves; ++i) {
        entity e;
        e.dim = 1;
        e.tag = tagCount++;
        int numBoundingPoints;
        // curveTag(int) minX(double) minY(double) minZ(double)
        //  maxX(double) maxY(double) maxZ(double)
        //  numPhysicalTags(size_t) physicalTag(int) ...
        //  numBoundingPoints(size_t) pointTag(int) ... >>>>>>>>>>>>> Ignor??s pour le moment
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> ph4D >> ph5D >> ph6D >> e.numPhysicalTags;
        if(e.numPhysicalTags > 1) {
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
        }
        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          // It's possible physical groups were added to the model, but the geometric entities refer
          // to them with the physicalTag "0" (e.g. caviteP2WithGroups.msh). For now, the fix is to
          // sequentially number the physicals with tag 0, but it will conflict if some physicals
          // are numbered. An alternative would be to keep a distinct list of physical groups for
          // each dimension.
          if(ph1 != 0) {
            if(_physicalEntities.find({e.dim, ph1}) != _physicalEntities.end()) {
              _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
              // if(FE_VERBOSE) printf("In readGmsh4 (FE_VERBOSE) : Geometric entity with dim = %d -
              // tag = %d is part of Physical \"%s\"\n", e.dim, e.tagGmsh, _physicalEntities[{e.dim,
              // ph1}].name.c_str());
            } else if(_physicalEntities.find({e.dim, -ph1}) !=
                      _physicalEntities.end()) { // Physical entity numbered in reverse
              _physicalEntities[{e.dim, -ph1}].listEntities.push_back(e.tag);
              if(FE_VERBOSE > -1) {
                // printf("In readGmsh4 : Geometric entity with dim = %d - tag = %d is part of
                // Physical \"%s\"\n", e.dim, e.tagGmsh, _physicalEntities[{e.dim,
                // -ph1}].name.c_str());
                feWarning("Physical entity \"%s\" is negative on geometric entity (%d,%d) and is "
                          "thus numbered in "
                          "reverse order.",
                          _physicalEntities[{e.dim, -ph1}].name.c_str(), e.dim, e.tagGmsh);
              }
            }
          } else {
            _physicalEntities[{e.dim, tagUnnumberedEntities++}].listEntities.push_back(e.tag);
          }
        }
        input >> numBoundingPoints;
        for(int j = 0; j < numBoundingPoints; ++j) {
          input >> ph1;
          e.isBoundedBy.push_back(ph1);
        }
        _entities.insert({{e.dim, e.tagGmsh}, e});
        getline(input, buffer);
      }
      // Surfaces
      for(int i = 0; i < _nSurfaces; ++i) {
        entity e;
        e.dim = 2;
        e.tag = tagCount++;
        int numBoundingCurves;
        // surfaceTag(int) minX(double) minY(double) minZ(double)
        //   maxX(double) maxY(double) maxZ(double)
        //   numPhysicalTags(size_t) physicalTag(int) ...
        //   numBoundingCurves(size_t) curveTag(int) ... >>>>>>>>>>>>> Ignor??s pour le moment
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> ph4D >> ph5D >> ph6D >> e.numPhysicalTags;
        if(e.numPhysicalTags > 1) {
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
        }
        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          if(ph1 != 0) {
            if(_physicalEntities.find({e.dim, ph1}) != _physicalEntities.end()) {
              _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
              // if(FE_VERBOSE) printf("In readGmsh4 (FE_VERBOSE) : Geometric entity with dim = %d -
              // tag = %d is part of Physical \"%s\"\n", e.dim, ph1, _physicalEntities[{e.dim,
              // ph1}].name.c_str());
            } else if(_physicalEntities.find({e.dim, -ph1}) !=
                      _physicalEntities.end()) { // Physical entity numbered in reverse
              _physicalEntities[{e.dim, -ph1}].listEntities.push_back(e.tag);
              if(FE_VERBOSE > -1) {
                // printf("In readGmsh4 : Geometric entity with dim = %d - tag = %d is part of
                // Physical \"%s\"\n", e.dim, e.tagGmsh, _physicalEntities[{e.dim,
                // -ph1}].name.c_str());
                feWarning("Physical entity \"%s\" is negative on geometric entity (%d,%d) and is "
                          "thus numbered in "
                          "reverse order.\n",
                          _physicalEntities[{e.dim, -ph1}].name.c_str(), e.dim, e.tagGmsh);
              }
            }
          } else {
            _physicalEntities[{e.dim, tagUnnumberedEntities++}].listEntities.push_back(e.tag);
          }
        }
        input >> numBoundingCurves;
        for(int j = 0; j < numBoundingCurves; ++j) {
          input >> ph1;
          e.isBoundedBy.push_back(ph1);
        }
        _entities.insert({{e.dim, e.tagGmsh}, e});
        getline(input, buffer);
      }
      // Volumes
      for(int i = 0; i < _nVolumes; ++i) {
        entity e;
        e.dim = 3;
        e.tag = tagCount++;
        int numBoundingSurfaces;
        // volumeTag(int) minX(double) minY(double) minZ(double)
        //   maxX(double) maxY(double) maxZ(double)
        //   numPhysicalTags(size_t) physicalTag(int) ...
        //   numBoundingSurfaces(size_t) surfaceTag(int) ... >>>>>>>>>>>>> Ignor??s pour le moment
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> ph4D >> ph5D >> ph6D >> e.numPhysicalTags;
        if(e.numPhysicalTags > 1) {
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
        }
        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          if(ph1 != 0) {
            _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
          } else {
            _physicalEntities[{e.dim, tagUnnumberedEntities++}].listEntities.push_back(e.tag);
          }
        }
        input >> numBoundingSurfaces;
        for(int j = 0; j < numBoundingSurfaces; ++j) {
          input >> ph1;
          e.isBoundedBy.push_back(ph1);
        }
        _entities.insert({{e.dim, e.tagGmsh}, e});
        getline(input, buffer);
      }
    }

    else if(buffer == "$Nodes") { // Read nodes
      int numEntityBlocks, countVertex = 0;
      double x, y, z;
      // numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
      input >> numEntityBlocks >> _nNod; // != nVertices, mais alloue _vertices : ambigu (:
      getline(input, buffer);
      _vertices.resize(_nNod);
      _sequentialNodeToGmshNode.resize(_nNod);
      // A counter for the dim > 0 entities
      int entityCount = 0;
      for(int i = 0; i < numEntityBlocks; ++i) {
        int entityDim, numNodes;
        // entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
        input >> entityDim >> ph2 >> ph3 >> numNodes;
        if(entityDim > 0) numNodesInBlock[entityCount++] = numNodes;
        std::vector<int> gmshNodeTag(numNodes);
        for(int iNode = 0; iNode < numNodes; ++iNode) {
          input >> ph1;
          gmshNodeTag[iNode] = ph1;
          _sequentialNodeToGmshNode[countVertex] = gmshNodeTag[iNode];
          _verticesMap[gmshNodeTag[iNode]] = countVertex++;
        }
        countVertex -= numNodes;
        for(int iNode = 0; iNode < numNodes; ++iNode) {
          input >> x >> y >> z;
          _vertices[countVertex++] = Vertex(x, y, z, gmshNodeTag[iNode]);
        }
      }
      if((int)_verticesMap.size() != _nNod) {
        return feErrorMsg(FE_STATUS_ERROR, "Vertices indices are not unique.");
      }
    }

    else if(buffer == "$Elements") { // Read elements
      int numEntityBlocks, serialNumber, maxElementTag;
      // numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
      input >> numEntityBlocks >> _nElm >> ph1 >> maxElementTag;
      _nTotalElm = _nElm;

      int countElems = 0;
      _nNodesWithNoPhysical = 0;

      getline(input, buffer);
      int nEdges = 1; // Starting edge numbering at 0 to distinguish + or - in the connectivity

      for(int i = 0; i < numEntityBlocks; ++i) {
        int numElementsInBlock, elemType, entityTag, entityDim;
        // entityDim(int) entityTag(int) elementType(int) numElementsInBlock(size_t)
        input >> entityDim >> entityTag >> elemType >> numElementsInBlock;

        std::pair<int, int> p = {entityDim, entityTag};
        bool printNodeWarning = true;
        _entities[p].nElm = numElementsInBlock;
        // Element connectivity : there may be gaps in Gmsh's element numbering.
        // We use a serial numbering instead (countElems).
        _entities[p].connecElem.resize(numElementsInBlock);

        // Determine geometric feSpace and allocate nodes connectivity based on elemType
        switch(elemType) {
          case 1:
            [[gnu::fallthrough]]; //  2-node line
          case 8:
            if(curved && _entities[p].cncID == "") {
              _entities[p].cncID = "LineP2";
            }
            [[gnu::fallthrough]]; //  3-node line (2nd order)
          case 26:
            if(curved && _entities[p].cncID == "") {
              _entities[p].cncID = "LineP3";
            }
            [[gnu::fallthrough]]; //  4-node line (3rd order)
          case 27:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  5-node line (4th order)
          case 28:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  6-node line (5th order)
          case 62:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  7-node line (6th order)
          case 63:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  8-node line (7th order)
          case 64:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; //  9-node line (8th order)
          case 65:
            if(curved) { /* Not supported */
            }
            [[gnu::fallthrough]]; // 10-node line (9th order)
          case 66:
            if(curved) { /* Not supported */
            } // 11-node line (10th order)
            {
              // Default is linear interpolation for the geometry
              if(!curved) {
                _entities[p].cncID = "LineP1";
              }
              _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType - 1] : 2;
              _entities[p].nEdgePerElem = 0;
              _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
              break;
            }
          case 2:
            [[gnu::fallthrough]]; //  3-node triangle
          case 9:
            if(curved && _entities[p].cncID == "") {
              _entities[p].cncID = "TriP2";
            }
            [[gnu::fallthrough]]; //  6-node triangle (2nd order)
          case 21:
            [[gnu::fallthrough]]; // 10-node triangle (3rd order)
          case 23:
            [[gnu::fallthrough]]; // 15-node triangle (4th order)
          case 25:
            [[gnu::fallthrough]]; // 21-node triangle (5th order)
          case 42:
            [[gnu::fallthrough]]; // 28-node triangle (6th order)
          case 43:
            [[gnu::fallthrough]]; // 36-node triangle (7th order)
          case 44:
            [[gnu::fallthrough]]; // 45-node triangle (8th order)
          case 45:
            [[gnu::fallthrough]]; // 55-node triangle (9th order)
          case 46: // 66-node triangle (10th order)
          {
            if(!curved) {
              _entities[p].cncID = "TriP1";
            }
            _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType - 1] : 3;
            _entities[p].nEdgePerElem = 3;
            _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
            _entities[p].connecEdges.resize(numElementsInBlock * _entities[p].nEdgePerElem);
            break;
          }
          case 3:
            [[gnu::fallthrough]]; //   4-node quadrangle
          case 10:
            [[gnu::fallthrough]]; //   9-node quadrangle (2nd order)
          case 36:
            [[gnu::fallthrough]]; //  16-node quadrangle (3rd order)
          case 37:
            [[gnu::fallthrough]]; //  25-node quadrangle (4th order)
          case 38:
            [[gnu::fallthrough]]; //  36-node quadrangle (5th order)
          case 47:
            [[gnu::fallthrough]]; //  49-node quadrangle (6th order)
          case 48:
            [[gnu::fallthrough]]; //  64-node quadrangle (7th order)
          case 49:
            [[gnu::fallthrough]]; //  81-node quadrangle (8th order)
          case 50:
            [[gnu::fallthrough]]; // 100-node quadrangle (9th order)
          case 51: // 121-node quadrangle (10th order)
          {
            _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType - 1] : 4;
            _entities[p].nEdgePerElem = 4;
            _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
            _entities[p].connecEdges.resize(numElementsInBlock * _entities[p].nEdgePerElem);
            return feErrorMsg(FE_STATUS_ERROR,
                              "Interpolant pas pris en charge pour la g??om??trie de "
                              "l'entit?? %d (quad).\n",
                              entityTag);
            break;
          }
          case 4:
            [[gnu::fallthrough]]; //   4-node tetrahedron
          case 11:
            [[gnu::fallthrough]]; //  10-node tetrahedron (2nd order)
          case 29:
            [[gnu::fallthrough]]; //  20-node tetrahedron (3rd order)
          case 30:
            [[gnu::fallthrough]]; //  35-node tetrahedron (4th order)
          case 31:
            [[gnu::fallthrough]]; //  56-node tetrahedron (5th order)
          case 71:
            [[gnu::fallthrough]]; //  84-node tetrahedron (6th order)
          case 72:
            [[gnu::fallthrough]]; // 120-node tetrahedron (7th order)
          case 73:
            [[gnu::fallthrough]]; // 165-node tetrahedron (8th order)
          case 74:
            [[gnu::fallthrough]]; // 220-node tetrahedron (9th order)
          case 75: // 286-node tetrahedron (10th order)
          {
            return feErrorMsg(FE_STATUS_ERROR,
                              "Interpolant pas pris en charge pour la g??om??trie de "
                              "l'entit?? %d (tet).\n",
                              entityTag);
          }
          case 5:
            [[gnu::fallthrough]]; //    8-node hexahedron
          case 12:
            [[gnu::fallthrough]]; //   27-node hexahedron (2nd order)
          case 92:
            [[gnu::fallthrough]]; //   64-node hexahedron (3rd order)
          case 93:
            [[gnu::fallthrough]]; //  125-node hexahedron (4th order)
          case 94:
            [[gnu::fallthrough]]; //  216-node hexahedron (5th order)
          case 95:
            [[gnu::fallthrough]]; //  343-node hexahedron (6th order)
          case 96:
            [[gnu::fallthrough]]; //  512-node hexahedron (7th order)
          case 97:
            [[gnu::fallthrough]]; //  729-node hexahedron (8th order)
          case 98: // 1000-node hexahedron (9th order)
          {
            return feErrorMsg(FE_STATUS_ERROR,
                              "Interpolant pas pris en charge pour la g??om??trie de "
                              "l'entit?? %d (hex).\n",
                              entityTag);
          }
          case 15: // 1-node point
          {
            // We only allow a Point geometric to have a single node
            // if(numElementsInBlock > 1){
            //   printf("In readMesh : Warning - More than one node are present in geometric entity
            //   (dim = %d, gmshTag = %d)."
            //     "Only the first node was added, the others were discarded.\n", _entities[p].dim,
            //     _entities[p].tagGmsh);
            // }
            _entities[p].cncID = "Point0D";
            _entities[p].nNodePerElem = 1;
            _entities[p].nEdgePerElem = 0;
            _entities[p].connecNodes.resize(1);
            _entities[p].connecElem.resize(1);
            break;
          }
          default:
            return feErrorMsg(FE_STATUS_ERROR, "Unsupported Gmsh element type.");
        } // switch(elemType)

        std::map<int, int>::const_iterator it;
        for(int iElm = 0; iElm < numElementsInBlock; ++iElm) {
          // int nElemNodes = nodes_of_gmsh_element[elemType-1];
          int nElemNodes = _entities[p].nNodePerElem;
          std::vector<int> elemNodesGmsh(nElemNodes);
          std::vector<int> elemNodes(nElemNodes);
          // elementTag(size_t) nodeTag(size_t) ...
          input >> serialNumber; // Unused

          // Point entity has vector connecElem of size 1 : do not modify it here
          if(elemType != 15) {
            _entities[p].connecElem[iElm] = countElems++;
          }

          for(int j = 0; j < nElemNodes; ++j) {
            input >> ph1;
            // Verification
            it = _verticesMap.find(ph1);
            if(it == _verticesMap.end()) {
              return feErrorMsg(FE_STATUS_ERROR, "Element node %d does not exist /-:", ph1);
            }
            elemNodesGmsh[j] = ph1; // The Gmsh number stored in ph1 is used to create the edges
            elemNodes[j] =
              it->second; // The node number (0...nNode) is used to create the connectivity
          }

          switch(elemType) {
            case 1:
              [[gnu::fallthrough]]; //  2-node line
            case 8:
              [[gnu::fallthrough]]; // TODO : Decider de la numerotation pour les P2+ : {0,2,1} ou
                                    // {0,1,2} //  3-node line (2nd order)
            case 26:
              [[gnu::fallthrough]]; //  4-node line (3rd order)
            case 27:
              [[gnu::fallthrough]]; //  5-node line (4th order)
            case 28:
              [[gnu::fallthrough]]; //  6-node line (5th order)
            case 62:
              [[gnu::fallthrough]]; //  7-node line (6th order)
            case 63:
              [[gnu::fallthrough]]; //  8-node line (7th order)
            case 64:
              [[gnu::fallthrough]]; //  9-node line (8th order)
            case 65:
              [[gnu::fallthrough]]; // 10-node line (9th order)
            case 66: // 11-node line (10th order)
            {
              // Keep Gmsh numbering for the high order nodes ?
              for(int j = 0; j < nElemNodes; ++j) {
                _entities[p].connecNodes[nElemNodes * iElm + j] = elemNodes[j];
              }
              break;
            }
            case 2:
              [[gnu::fallthrough]]; //  3-node triangle
            case 9:
              [[gnu::fallthrough]]; //  6-node triangle (2nd order)
            case 21:
              [[gnu::fallthrough]]; // 10-node triangle (3rd order)
            case 23:
              [[gnu::fallthrough]]; // 15-node triangle (4th order)
            case 25:
              [[gnu::fallthrough]]; // 21-node triangle (5th order)
            case 42:
              [[gnu::fallthrough]]; // 28-node triangle (6th order)
            case 43:
              [[gnu::fallthrough]]; // 36-node triangle (7th order)
            case 44:
              [[gnu::fallthrough]]; // 45-node triangle (8th order)
            case 45:
              [[gnu::fallthrough]]; // 55-node triangle (9th order)
            case 46: // 66-node triangle (10th order)
            {
              if(curved) {
                _entities[p].connecNodes[nElemNodes * iElm + 0] = elemNodes[0];
                _entities[p].connecNodes[nElemNodes * iElm + 1] = elemNodes[2];
                _entities[p].connecNodes[nElemNodes * iElm + 2] = elemNodes[1];
                _entities[p].connecNodes[nElemNodes * iElm + 3] = elemNodes[5];
                _entities[p].connecNodes[nElemNodes * iElm + 4] = elemNodes[4];
                _entities[p].connecNodes[nElemNodes * iElm + 5] = elemNodes[3];
              } else { // Curved are inverted at the moment
                for(int j = 0; j < nElemNodes; ++j) {
                  _entities[p].connecNodes[nElemNodes * iElm + j] = elemNodes[j];
                  // _entities[p].connecNodes[nElemNodes * iElm + j] = elemNodesGmsh[j];
                }
              }

              // Construct the triangle edges :
              Vertex *v0, *v1;
              for(int k = 0; k < 3; ++k) {
                if(!curved) {
                  if(k == 2) {
                    v0 = &_vertices[_verticesMap[elemNodesGmsh[2]]];
                    v1 = &_vertices[_verticesMap[elemNodesGmsh[0]]];
                  } else {
                    v0 = &_vertices[_verticesMap[elemNodesGmsh[k]]];
                    v1 = &_vertices[_verticesMap[elemNodesGmsh[k + 1]]];
                  }
                } else { // Curved are inverted at the moment
                  if(k == 0) {
                    v0 = &_vertices[_verticesMap[elemNodesGmsh[0]]];
                    v1 = &_vertices[_verticesMap[elemNodesGmsh[2]]];
                  } else if(k == 1) {
                    v0 = &_vertices[_verticesMap[elemNodesGmsh[2]]];
                    v1 = &_vertices[_verticesMap[elemNodesGmsh[1]]];
                  } else {
                    v0 = &_vertices[_verticesMap[elemNodesGmsh[1]]];
                    v1 = &_vertices[_verticesMap[elemNodesGmsh[0]]];
                  }
                }
                std::pair<std::set<Edge, EdgeLessThan>::iterator, bool> ret;
                Edge e(v0, v1, nEdges, _entities[p].physicalTags[0]);
                ret = _edges.insert(e);
                if(ret.second) {
                  // Edge was added to the set : nEdges is added to connecEdges
                  _entities[p].connecEdges[3 * iElm + k] = nEdges++;
                } else {
                  // Edge is already in the set : the negative is added to connecEdges
                  // Assumes an edge is shared by only two triangles in 2D
                  // More tests are required in 3D, where an edge can be shared by N tets
                  _entities[p].connecEdges[3 * iElm + k] = -ret.first->getTag();
                }
              }
              break;
            }
            case 3:
              [[gnu::fallthrough]]; //   4-node quadrangle
            case 10:
              [[gnu::fallthrough]]; //   9-node quadrangle (2nd order)
            case 36:
              [[gnu::fallthrough]]; //  16-node quadrangle (3rd order)
            case 37:
              [[gnu::fallthrough]]; //  25-node quadrangle (4th order)
            case 38:
              [[gnu::fallthrough]]; //  36-node quadrangle (5th order)
            case 47:
              [[gnu::fallthrough]]; //  49-node quadrangle (6th order)
            case 48:
              [[gnu::fallthrough]]; //  64-node quadrangle (7th order)
            case 49:
              [[gnu::fallthrough]]; //  81-node quadrangle (8th order)
            case 50:
              [[gnu::fallthrough]]; // 100-node quadrangle (9th order)
            case 51: // 121-node quadrangle (10th order)
            {
              return feErrorMsg(FE_STATUS_ERROR, "Interpolant pas pris en charge (quad).\n");
            }
            case 4:
              [[gnu::fallthrough]]; //   4-node tetrahedron
            case 11:
              [[gnu::fallthrough]]; //  10-node tetrahedron (2nd order)
            case 29:
              [[gnu::fallthrough]]; //  20-node tetrahedron (3rd order)
            case 30:
              [[gnu::fallthrough]]; //  35-node tetrahedron (4th order)
            case 31:
              [[gnu::fallthrough]]; //  56-node tetrahedron (5th order)
            case 71:
              [[gnu::fallthrough]]; //  84-node tetrahedron (6th order)
            case 72:
              [[gnu::fallthrough]]; // 120-node tetrahedron (7th order)
            case 73:
              [[gnu::fallthrough]]; // 165-node tetrahedron (8th order)
            case 74:
              [[gnu::fallthrough]]; // 220-node tetrahedron (9th order)
            case 75: // 286-node tetrahedron (10th order)
            {
              return feErrorMsg(FE_STATUS_ERROR, "Interpolant pas pris en charge (tet).\n");
            }
            case 5:
              [[gnu::fallthrough]]; //    8-node hexahedron
            case 12:
              [[gnu::fallthrough]]; //   27-node hexahedron (2nd order)
            case 92:
              [[gnu::fallthrough]]; //   64-node hexahedron (3rd order)
            case 93:
              [[gnu::fallthrough]]; //  125-node hexahedron (4th order)
            case 94:
              [[gnu::fallthrough]]; //  216-node hexahedron (5th order)
            case 95:
              [[gnu::fallthrough]]; //  343-node hexahedron (6th order)
            case 96:
              [[gnu::fallthrough]]; //  512-node hexahedron (7th order)
            case 97:
              [[gnu::fallthrough]]; //  729-node hexahedron (8th order)
            case 98: // 1000-node hexahedron (9th order)
            {
              return feErrorMsg(FE_STATUS_ERROR, "Interpolant pas pris en charge (hex).\n");
            }
            case 15: // 1-node point
            {
              // We only allow a geometric Point to have a single node
              if(_entities[p].nElm == 1) {
                _nNodesWithNoPhysical++;
                if(printNodeWarning) {
                  feWarning(
                    "More than one node are present in geometric entity (dim = %d, gmshTag = %d). "
                    "Only the first node was added, the others were discarded.",
                    _entities[p].dim, _entities[p].tagGmsh);
                  printNodeWarning = false;
                }
              } else {
                _entities[p].nElm = 1;
                _entities[p].connecElem[0] = countElems++;
                _entities[p].connecNodes[0] = elemNodes[0];
              }
              break;
            }
            default: // any other element
              return feErrorMsg(FE_STATUS_ERROR, "Unsupported Gmsh element type.");
          } // switch (type_of_element)

          // If not a high-order element, remaining nodes were not parsed
          getline(input, buffer);

        } // for numElemInEntity
      } // for numEntities
    } // if buffer = "Elements"
  } // while input

  if(_physicalEntities.size() == 0) {
    return feErrorMsg(FE_STATUS_ERROR, "No physical entities defined on the mesh.\n");
  }

  return FE_STATUS_OK;
}

feStatus feMesh2DP1::readGmsh(std::string meshName, bool curved,
                              mapType physicalEntitiesDescription)
{
  _gmshVersion = 0;
  std::filebuf fb;

  feInfoCond(FE_VERBOSE > 0, "Reading mesh file : %s", meshName.c_str());

  std::ifstream f(meshName.c_str());
  if(!f.good()) {
    return feErrorMsg(FE_STATUS_READ_ERROR, "Mesh file does not exist.");
  }

  // Mesh ID is unused for now
  _ID = "myBeautifulMesh";

  if(fb.open(meshName, std::ios::in)) {
    std::istream input(&fb);

    std::string buffer;
    int dsize;
    getline(input, buffer);
    input >> _gmshVersion >> _isBinary >> dsize;

    if(_isBinary) {
      return feErrorMsg(FE_STATUS_ERROR, "Only reading ASCII files.");
    }

    if(_gmshVersion == 2.2) {
      feCheck(readMsh2(input, curved, physicalEntitiesDescription));
    } else if(_gmshVersion >= 4) {
      feCheck(readMsh4(input, curved, physicalEntitiesDescription));
    } else {
      return feErrorMsg(FE_STATUS_ERROR, "Mesh version should be 2.2 or 4.1 or higher.");
    }

    fb.close();
  } // if fb.open

  if(FE_VERBOSE > 1) {
    feInfo("List of geometric entities with their attributes :");
    for(auto const &x : _entities) {
      entity e = x.second;
      feInfo("Entity #%2d : dim = %1d - gmshTag = %3d - nElm = %6d - nNodePerElem = %2d - "
             "geometric interpolant : %10s "
             "- part of physical groups : %2d",
             e.tag, e.dim, e.tagGmsh, e.nElm, e.nNodePerElem, e.cncID.c_str(),
             (e.numPhysicalTags > 0) ? e.physicalTags[0] : -1);
      // std::cout<<"elems : "<<std::endl;
      // for(auto val : e.connecElem)
      //   std::cout<<val<<" ";
      // std::cout<<std::endl;
      // std::cout<<"nodes : "<<std::endl;
      // for(auto val : e.connecNodes)
      //   std::cout<<val<<" ";
      // std::cout<<std::endl;
      // std::cout<<"edges : "<<std::endl;
      // for(auto val : e.connecEdges)
      //   std::cout<<val<<" ";
      // std::cout<<std::endl;
    }
  }

  // Geometric connectivities are defined on the physical entities (named domains) :
  // we need to transfer the information to the physical entities.
  // Sequential tag for all entities, different from the Gmsh entity tag (reset for each dimension)
  for(auto const &x : _entities) {
    entity e = x.second;
    // std::cout<<"Entity "<<e.tag<<" (dim = "<<e.dim<<" tag Gmsh = "<<e.tagGmsh<<")"<<std::endl;
    for(auto &p : _physicalEntities) {
      physicalEntity &pE = p.second;
      // printf("Physical \"%s\" has %d geometric entities :\n", pE.name.c_str(),
      // pE.listEntities.size());
      bool error = false;
      for(auto ent : pE.listEntities) {
        if(ent == e.tag && pE.dim == e.dim) {
          // std::cout<<"    Physical "<<pE.name<<" - entity "<<ent<<" (tag Gmsh = "<<e.tagGmsh<<")
          // of dimension
          // "<<e.dim<<" and "<<e.nElm<<" elements : match"<<std::endl;
          pE.nElm += e.nElm;
          // std::cout<<pE.nElm<<std::endl;
          // Assign and check nNodePerElem
          if(pE.nNodePerElem == -1)
            pE.nNodePerElem = e.nNodePerElem;
          else {
            if(pE.nNodePerElem != e.nNodePerElem) error = true;
          }
          // Assign and check nEdgePerElem
          if(pE.nEdgePerElem == -1)
            pE.nEdgePerElem = e.nEdgePerElem;
          else {
            if(pE.nEdgePerElem != e.nEdgePerElem) error = true;
          }
          // Assign and check geoSpace and cncID
          if(pE.cncID == "")
            pE.cncID = e.cncID;
          else {
            if(pE.cncID != e.cncID) error = true;
          }
        }
      }
      if(error) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Multiple geometric connectivities on physical entity "
                          "%s"
                          ".\n",
                          pE.name.c_str());
      }
    }
  }

  if(FE_VERBOSE > 1) {
    feInfo("List of physical groups with their attributes :");
    for(auto &p : _physicalEntities) {
      physicalEntity &pE = p.second;
      feInfo("Physical %15s : dim = %1d - nEntities = %2d - nElm = %9d - nNodePerElem = %2d - "
             "nEdgePerElem = %1d - "
             "geometric interpolant : %10s",
             pE.name.c_str(), pE.dim, pE.listEntities.size(), pE.nElm, pE.nNodePerElem,
             pE.nEdgePerElem, pE.cncID.c_str());
    }
  }

  // Construct the intrinsic (physical) connectivities
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    pE.connecElem.resize(pE.nElm);
    pE.connecNodes.resize(pE.nElm * pE.nNodePerElem);
    pE.connecEdges.resize(pE.nElm * pE.nEdgePerElem);
    int countElm = 0;
    for(auto const &x : _entities) {
      entity e = x.second;
      for(auto ent : pE.listEntities) {
        if(ent == e.tag) {
          // printf("In physical %s - geometric entity %d : \n", pE.name.c_str(), e.tagGmsh);
          // std::vector<int> tmp = e.connecNodes;
          // std::sort(tmp.begin(), tmp.end());
          // tmp.erase( std::unique(tmp.begin(), tmp.end()), tmp.end() );
          // std::cout<<tmp.size()<<std::endl;
          // printf("Copying info from entity (%d,%d) with %d elements and cnc %s to physical
          // (%d,%d) named %s\n",
          //   e.dim, e.tag, e.nElm, e.cncID.c_str(), pE.dim, pE.tag, pE.name.c_str());
          // printf("Physical has %d elements and connecElem.size() = %d\n", pE.nElm,
          // pE.connecElem.size());
          for(int iElm = 0; iElm < e.nElm; ++iElm) {
            // Connec elem
            pE.connecElem[countElm] = e.connecElem[iElm];
            // Connec nodes
            for(int j = 0; j < e.nNodePerElem; ++j) {
              pE.connecNodes[e.nNodePerElem * countElm + j] =
                e.connecNodes[e.nNodePerElem * iElm + j];
            }
            // Connec edges
            for(int j = 0; j < e.nEdgePerElem; ++j) {
              pE.connecEdges[e.nEdgePerElem * countElm + j] =
                e.connecEdges[e.nEdgePerElem * iElm + j];
            }
            ++countElm;
          }
        }
      }
    }
  }

  // Construct a boundary representation (b-rep) for the physical entities
  // Not used at the moment and only works if the boundary curved of all considered surfaces
  // are part of a physical entity, so commented for now.

  // if(_gmshVersion >= 4) {
  //   for(auto const &x : _entities) {
  //     entity e = x.second;
  //     if(e.dim > 1) { // Not considering b-rep for curves at the moment
  //       if(e.numPhysicalTags > 0) { // if e is part of a physical entity
  //         std::cout<<"Entity "<<e.tag<<" is part of physical : "<<e.physicalTags[0] << std::endl;
  //         std::cout<<"and is bounded by geometric entities"<<std::endl;
  //         for(int boundedByE : e.isBoundedBy) { // loop over geometric entities bounding e...
  //           std::pair<int, int> p = {e.dim - 1, fabs(boundedByE)};
  //           int pp = (_entities[p].numPhysicalTags > 0) ? _entities[p].physicalTags[0] : -1;
  //           std::cout<<boundedByE<<" which is part of physical "<< pp << std::endl;
  //           int tag = e.physicalTags[0] -
  //                     1; // Physical entities are numbered starting at 1, sequentially by
  //                     dimension
  //           int boundedByTag = _entities[p].physicalTags[0] - 1;
  //           // Basically unique() :
  //           if(std::find(_physicalEntities[tag].isBoundedBy.begin(),
  //                        _physicalEntities[tag].isBoundedBy.end(),
  //                        boundedByTag) == _physicalEntities[tag].isBoundedBy.end())
  //             _physicalEntities[tag].isBoundedBy.push_back(boundedByTag);
  //           if(std::find(_physicalEntities[boundedByTag].isBoundaryOf.begin(),
  //                        _physicalEntities[boundedByTag].isBoundaryOf.end(),
  //                        tag) == _physicalEntities[boundedByTag].isBoundaryOf.end())
  //             _physicalEntities[boundedByTag].isBoundaryOf.push_back(tag);
  //         }
  //       }
  //     }
  //   }
  // }

  // Transfer edge/face connectivity to boundary elements
  int maxDim = 0;
  for(auto &pE : _physicalEntities) {
    maxDim = fmax(maxDim, pE.second.dim);
  }
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(maxDim == 2) {
      // Transfer edge connectivity to boundary edges
      if(pE.dim == maxDim - 1) {
        pE.nEdgePerElem = 1;
        pE.connecEdges.resize(pE.nElm);
        // Brute force : build and find edge in the edges set
        Vertex *v0, *v1;
        std::set<Edge, EdgeLessThan>::iterator it;
        for(int i = 0; i < pE.nElm; ++i) {
          v0 = &_vertices[pE.connecNodes[pE.nNodePerElem * i]];
          v1 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 1]];
          Edge e(v0, v1);
          it = _edges.find(e);
          if(it != _edges.end()) {
            pE.connecEdges[i] = it->getTag();
            if(pE.connecEdges[i] < 0) {
              // Boundary edges should be positively oriented
              feWarning("Boundary edge (%d,%d) orientation is negative.", v0->getTag(),
                        v1->getTag());
            }
          } else {
            // Edge should be in the set...
            return feErrorMsg(FE_STATUS_ERROR,
                              "Boundary edge (%d,%d) was not found in the set of edges...",
                              v0->getTag(), v1->getTag());
          }
        }
      }
    } else if(maxDim == 3) {
      // Transfer face and edge connectivities to elements
    }
  }

  // Assign pointer to geometric interpolation space
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(pE.cncID == "Point0D") {
      pE.geoSpace = new feSpace1DP0("xyz");
    } else if(pE.cncID == "LineP1") {
      pE.geoSpace = new feSpace1DP1("xyz");
    } else if(pE.cncID == "LineP2") {
      pE.geoSpace = new feSpace1DP2("xyz");
    } else if(pE.cncID == "LineP3") {
      pE.geoSpace = new feSpace1DP3("xyz");
    } else if(pE.cncID == "TriP1") {
      pE.geoSpace = new feSpaceTriP1("xyz");
    } else if(pE.cncID == "TriP2") {
      pE.geoSpace = new feSpaceTriP2("xyz");
    } else {
      return feErrorMsg(FE_STATUS_ERROR, "Unknown geometric connectivity \"%s\" on domain \"%s\".",
                        pE.cncID.c_str(), pE.name.c_str());
    }
  }

  // Finally, create geometric connectivities
  int nCncGeo = 0;
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;

    feCncGeo *cnc =
      new feCncGeo(nCncGeo, pE.dim, pE.nNodePerElem, pE.nElm, pE.nEdgePerElem, pE.name, pE.cncID,
                   pE.geoSpace, pE.connecNodes, pE.connecElem, pE.connecEdges);

    _cncGeo.push_back(cnc);
    _cncGeoMap[pE.name] = nCncGeo;
    cnc->getFeSpace()->setCncGeoTag(nCncGeo++);
  }

  _nEdg = _edges.size();
  _dim = maxDim;

  _nInteriorElm = 0;
  _nBoundaryElm = 0;
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(pE.dim == _dim) _nInteriorElm += pE.nElm;
    if(pE.dim == _dim - 1) _nBoundaryElm += pE.nElm;
  }
  // Add all interior elements (for now only triangles) to the _elements vector
  // Physical entities are non-overlapping so we just loop over them and elements shouldn't be
  // duplicated
  _elements.resize(_nInteriorElm);
  int cnt = 0;
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(_dim == 2 && pE.dim == _dim) {
      Vertex *v0, *v1, *v2;
      for(int i = 0; i < pE.nElm; ++i) {
        v0 = &_vertices[pE.connecNodes[pE.nNodePerElem * i]];
        v1 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 1]];
        v2 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 2]];
        // The triangle tag is i because it is local to the pE's cncgeo, hence the tag is not
        // unique...
        _elements[cnt++] = new Triangle(v0, v1, v2, i, pE.tag);
      }
    }
  }

  // TODO : Also add boundaryElements when necessary

  // Assign a pointer to this mesh to each of its geometric connectivities and their fespace
  for(feCncGeo *cnc : _cncGeo) {
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }
  _nCncGeo = _cncGeo.size();

  // Create an RTree storing the elements of highest dimension
  for(int i = 0; i < _elements.size(); ++i) {
    Triangle *t = _elements[i];
    SBoundingBox3d bbox;
    for(int i = 0; i < 3; ++i) {
      Vertex *v = t->getVertex(i);
      SPoint3 pt(v->x(), v->y(), v->z());
      bbox += pt;
    }
    _rtree.Insert((double *)(bbox.min()), (double *)(bbox.max()), t->getTag());
  }

  return FE_STATUS_OK;
}
