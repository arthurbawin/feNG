#include "feNG.h"
#include "feMesh.h"
#include "feTriangle.h"

#include "SBoundingBox3d.h"

extern int FE_VERBOSE;
static int VERBOSE_THRESHOLD = 0;

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

feStatus feMesh2DP1::readMsh2(std::istream &input, const bool curved, const bool reversed,
                              const mapType physicalEntitiesDescription)
{
  UNUSED(reversed);

  std::string buffer;
  int ph1; // Placeholder
  // std::map<int, int> _verticesMap; // Gmsh tag (which may include gaps) to sequential tag. Just
  // in case : not sure there are gaps in msh2...
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
      pE.nTriPerElem = -1;
      pE.nTetPerElem = -1;
      pE.interp = geometricInterpolant::NONE;
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
        pE.interp = geometricInterpolant::NONE;
        _physicalEntities[{pE.dim, pE.tag}] = pE;
      }

      /* Save the description of the physical entities : this is used when looping and adapting the
      mesh. Since mmg does not keep track of the name of the physical entities, they have to be
      saved and loaded during the following adaptation cycles. */
      for(auto pair : _physicalEntities) {
        _physicalEntitiesDescription[pair.first] = _physicalEntities[pair.first].name;
      }
    }

    // There is no Entities block in msh2 : entities information is in the elements

    else if(buffer == "$Nodes") { // Read nodes

      if(_nPhysicalEntities == 0) {
        return feErrorMsg(FE_STATUS_ERROR, "No physical entities defined on the mesh.");
      }

      int gmshNumber, countVertex = 0;
      double x, y, z;
      input >> _nVertices;
      _vertices.resize(_nVertices);
      getline(input, buffer);
      for(int iNode = 0; iNode < _nVertices; ++iNode) {
        // node-number x-coord y-coord z-coord
        input >> gmshNumber >> x >> y >> z;
        _verticesMap[gmshNumber] = countVertex;
        _vertices[countVertex++] = Vertex(x, y, z, gmshNumber);
      }

      if((int)_verticesMap.size() != _nVertices) {
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
      /* number-of-elements */
      input >> _nElm;
      _nTotalElm = _nElm;

      int serialNumberElm = 0; // To number the elements when skipping the Points
      for(int iElm = 0; iElm < _nElm; ++iElm) {
        getline(input, buffer);
        /* elm-number elm-type number-of-tags < tag > … node-number-list */
        input >> gmshElemTag >> elemType >> numEntities >> physicalTag;
        /* Physical and geometric entities are not separated, and physical seem to be given first.
        Since we only accept non-overlapping physical domains (i.e. each element belongs to a
        single physical entity), let's just assume the first number is the physical domain and
        the following are the geometric entities. */
        if(numEntities != 2) {
          if(numEntities < 2)
            return feErrorMsg(FE_STATUS_ERROR,
                              "Element with Gmsh tag %d is not part of a physical entity.",
                              gmshElemTag);
          else
            return feErrorMsg(FE_STATUS_ERROR,
                              "Element with Gmsh tag %d is part of more than one physical entity.",
                              gmshElemTag);
        }

        int entityDim = dim_of_gmsh_element[elemType - 1];

        /* Physical tag (number of the physical entity) should be strictly positive.
        Meshes produced by mmg2d include points with physicalTag = 0, which we skip. */
        if(physicalTag <= 0) continue;

        /* numEntities should be 2 since an element belongs to max 1 physical entity
        and should belong to a single geometric entity (?) : no need to loop */
        // for(int i = 1; i < numEntities; ++i){
        input >> geometricTag;

        std::pair<int, int> p = {entityDim, geometricTag};

        /* Meshes produced by mmg2d seem to include unwanted 0D physical entities (?).
        If a description of the original physical entities is provided, we make sure
        the read physical entities is indeed part of these entities. */
        if(physicalEntitiesDescription.size() > 0) {
          bool isPartOfTheOriginalEntities = false;
          for(auto pair : physicalEntitiesDescription) {
            /* pair = { {dim,tag}, name } */
            // printf("Comparing entity dim = %d - tag = %d with stored entity %s - dim = %d - tag =
            // %d\n", entityDim, physicalTag, pair.second.data(), pair.first.first,
            // pair.first.second);
            if(entityDim == pair.first.first && physicalTag == pair.first.second)
              isPartOfTheOriginalEntities = true;
          }
          if(!isPartOfTheOriginalEntities) {
            continue;
          }
        }

        /* If geometric entity doesn't exist, create it : */
        if(_entities.find(p) == _entities.end()) {
          entity e;
          e.dim = entityDim;
          e.tag = tagCount++;
          e.tagGmsh = geometricTag;
          e.numPhysicalTags = 1;
          e.physicalTags.push_back(physicalTag);
          e.nElm = 0;
          std::pair<int, int> pPhys = {e.dim, physicalTag};
          _physicalEntities[pPhys].listEntities.push_back(e.tag);
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
              _entities[p].interp = geometricInterpolant::LINEP2;
            }
            [[gnu::fallthrough]]; //  3-node line (2nd order)
          case 26:
            if(curved) {
              _entities[p].interp = geometricInterpolant::LINEP3;
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
                _entities[p].interp = geometricInterpolant::LINEP1;
              }
              _entities[p].nNodePerElem = (curved) ? nNodePerElem : 2;
              _entities[p].nEdgePerElem = 0;
              break;
            }
          case 2:
            [[gnu::fallthrough]]; //  3-node triangle
          case 9:
            if(curved) {
              _entities[p].interp = geometricInterpolant::TRIP2;
            }
            [[gnu::fallthrough]]; //  6-node triangle (2nd order)
          case 21:
            if(curved) {
              _entities[p].interp = geometricInterpolant::TRIP3;
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
              _entities[p].interp = geometricInterpolant::TRIP1;
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
                              "Interpolant pas (encore) pris en charge pour la géométrie "
                              "de l'entité (dim = %d, tag = %d) (quad).\n",
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
                              "Interpolant pas (encore) pris en charge pour la géométrie "
                              "de l'entité (dim = %d, tag = %d) (tet).\n",
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
                              "Interpolant pas (encore) pris en charge pour la géométrie "
                              "de l'entité (dim = %d, tag = %d) (hex).\n",
                              entityDim, geometricTag);
          }
          case 15: // 1-node point
          {
            _entities[p].interp = geometricInterpolant::POINTP0;
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

//
// FIXME: Give the list of modified attributes in feMesh
//
feStatus feMesh2DP1::readMsh4(std::istream &input, const bool curved, const bool reversed,
                              const mapType physicalEntitiesDescription)
{
  std::string buffer;
  // Placeholders
  int ph1, ph2, ph3;
  double ph1D, ph2D, ph3D, ph4D, ph5D, ph6D;
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
      pE.nTriPerElem = -1;
      pE.nTetPerElem = -1;
      pE.interp = geometricInterpolant::NONE;
      _physicalEntities[pair.first] = pE;
      feInfoCond(FE_VERBOSE > VERBOSE_THRESHOLD, "\t\tCreated Physical entity \"%s\"",
                 pE.name.data());
    }
  }

  while(input >> buffer)
  {
    //
    // Read physical entities (named groups of geometric entities)
    //
    if(buffer == "$PhysicalNames")
    {
      input >> _nPhysicalEntities;

      if(_nPhysicalEntities == 0) {
        return feErrorMsg(FE_STATUS_ERROR, "No physical entities defined on the mesh.");
      }

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
        pE.nTriPerElem = -1;
        pE.nTetPerElem = -1;
        pE.interp = geometricInterpolant::NONE;
        _physicalEntities[{pE.dim, pE.tag}] = pE;
        feInfoCond(FE_VERBOSE > VERBOSE_THRESHOLD,
                   "\t\tCreated Physical entity \"%s\" with dimension %d and gmsh tag %d",
                   pE.name.data(), pE.dim, pE.tag);
      }

      // Save the description of the physical entities : this is used when looping and adapting the
      // mesh. Since mmg does not keep track of the name of the physical entities, they have to be
      // saved and loaded during the following adaptation cycles.
      // Update : this should be dealt with during adaptation by using Gmsh's API to add physicals
      // after adaptation with MMG
      for(auto pair : _physicalEntities) {
        _physicalEntitiesDescription[pair.first] = _physicalEntities[pair.first].name;
      }
    }

    //
    // Read geometric entities (building blocks of the geometry)
    //
    else if(buffer == "$Entities")
    {
      // numPoints(size_t) numCurves(size_t) numSurfaces(size_t) numVolumes(size_t)
      input >> _nPoints >> _nCurves >> _nSurfaces >> _nVolumes;

      int tagCount = 1;
      int tagUnnumberedEntities = 0;

      //
      // Points
      //
      for(int i = 0; i < _nPoints; ++i) {
        entity e;
        e.dim = 0;
        e.tag = tagCount++;
        // We only allow a Point geometric to have a single node
        e.nElm = 0;
        //  pointTag(int) X(double) Y(double) Z(double)
        //  numPhysicalTags(size_t) physicalTag(int)
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> e.numPhysicalTags;

        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          if(_physicalEntities.find({e.dim, ph1}) != _physicalEntities.end()) {
            _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
          } else if(_physicalEntities.find({e.dim, -ph1}) != _physicalEntities.end()) {
            // Physical entity numbered in reverse
            _physicalEntities[{e.dim, -ph1}].listEntities.push_back(e.tag);
            feWarning("Physical entity \"%s\" is negative on geometric entity (%d,%d) and is "
                      "thus numbered in reverse order.",
                      _physicalEntities[{e.dim, -ph1}].name.data(), e.dim, e.tagGmsh);
          }
        }

        if(e.numPhysicalTags > 1) {
          feInfo("Geometric entity with (dim,tagGmsh) = (%d,%d) is part of %d physical entities :",
                 e.dim, e.tagGmsh, e.numPhysicalTags);
          for(int ii = 0; ii < e.numPhysicalTags; ++ii) {
            feInfo("    is part of physical %s",
                   _physicalEntities[{e.dim, e.physicalTags[ii]}].name.data());
          }
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
        }

        _entities.insert({{e.dim, e.tagGmsh}, e});
      }

      //
      // Curves
      //
      for(int i = 0; i < _nCurves; ++i) {
        entity e;
        e.dim = 1;
        e.tag = tagCount++;
        int numBoundingPoints;
        // curveTag(int) minX(double) minY(double) minZ(double)
        //  maxX(double) maxY(double) maxZ(double)
        //  numPhysicalTags(size_t) physicalTag(int) ...
        //  numBoundingPoints(size_t) pointTag(int) ... >>>>>>>>>>>>> Ignorés pour le moment
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> ph4D >> ph5D >> ph6D >> e.numPhysicalTags;

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
            } else if(_physicalEntities.find({e.dim, -ph1}) != _physicalEntities.end()) {
              // Physical entity numbered in reversed order
              _physicalEntities[{e.dim, -ph1}].listEntities.push_back(e.tag);
              feWarning("Physical entity \"%s\" is negative on geometric entity (%d,%d) and is "
                        "thus numbered in reverse order.",
                        _physicalEntities[{e.dim, -ph1}].name.data(), e.dim, e.tagGmsh);
            }
          } else {
            _physicalEntities[{e.dim, tagUnnumberedEntities++}].listEntities.push_back(e.tag);
          }
        }

        if(e.numPhysicalTags > 1) {
          feInfo("Geometric entity with (dim,tagGmsh) = (%d,%d) is part of %d physical entities :",
                 e.dim, e.tagGmsh, e.numPhysicalTags);
          for(int ii = 0; ii < e.numPhysicalTags; ++ii) {
            feInfo("    is part of physical %s",
                   _physicalEntities[{e.dim, e.physicalTags[ii]}].name.data());
          }
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
        }

        input >> numBoundingPoints;
        for(int j = 0; j < numBoundingPoints; ++j) {
          input >> ph1;
          e.isBoundedBy.push_back(ph1);
        }
        _entities.insert({{e.dim, e.tagGmsh}, e});
        getline(input, buffer);
      }

      //
      // Surfaces
      //
      for(int i = 0; i < _nSurfaces; ++i) {
        entity e;
        e.dim = 2;
        e.tag = tagCount++;
        int numBoundingCurves;
        // surfaceTag(int) minX(double) minY(double) minZ(double)
        //   maxX(double) maxY(double) maxZ(double)
        //   numPhysicalTags(size_t) physicalTag(int) ...
        //   numBoundingCurves(size_t) curveTag(int) ... >>>>>>>>>>>>> Ignorés pour le moment
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> ph4D >> ph5D >> ph6D >> e.numPhysicalTags;

        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          if(ph1 != 0) {
            if(_physicalEntities.find({e.dim, ph1}) != _physicalEntities.end()) {
              _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
            } else if(_physicalEntities.find({e.dim, -ph1}) !=
                      _physicalEntities.end()) { // Physical entity numbered in reverse
              _physicalEntities[{e.dim, -ph1}].listEntities.push_back(e.tag);
              feWarning("Physical entity \"%s\" is negative on geometric entity (%d,%d) and is "
                        "thus numbered in "
                        "reverse order.\n",
                        _physicalEntities[{e.dim, -ph1}].name.data(), e.dim, e.tagGmsh);
            }
          } else {
            _physicalEntities[{e.dim, tagUnnumberedEntities++}].listEntities.push_back(e.tag);
          }
        }

        if(e.numPhysicalTags > 1) {
          feInfo("Geometric entity with (dim,tagGmsh) = (%d,%d) is part of %d physical entities :",
                 e.dim, e.tagGmsh, e.numPhysicalTags);
          for(int ii = 0; ii < e.numPhysicalTags; ++ii) {
            feInfo("    is part of physical %s",
                   _physicalEntities[{e.dim, e.physicalTags[ii]}].name.data());
          }
          return feErrorMsg(
            FE_STATUS_ERROR,
            "Geometric entity can only be part of a maximum of one "
            "physical entity (named domain). Overlapping physical domains are not allowed.");
        }

        input >> numBoundingCurves;
        for(int j = 0; j < numBoundingCurves; ++j) {
          input >> ph1;
          e.isBoundedBy.push_back(ph1);
        }
        _entities.insert({{e.dim, e.tagGmsh}, e});
        getline(input, buffer);
      }

      //
      // Volumes
      //
      for(int i = 0; i < _nVolumes; ++i) {
        entity e;
        e.dim = 3;
        e.tag = tagCount++;
        int numBoundingSurfaces;
        // volumeTag(int) minX(double) minY(double) minZ(double)
        //   maxX(double) maxY(double) maxZ(double)
        //   numPhysicalTags(size_t) physicalTag(int) ...
        //   numBoundingSurfaces(size_t) surfaceTag(int) ... >>>>>>>>>>>>> Ignorés pour le moment
        input >> e.tagGmsh >> ph1D >> ph2D >> ph3D >> ph4D >> ph5D >> ph6D >> e.numPhysicalTags;

        for(int j = 0; j < e.numPhysicalTags; ++j) {
          input >> ph1;
          e.physicalTags.push_back(ph1);
          if(ph1 != 0) {
            _physicalEntities[{e.dim, ph1}].listEntities.push_back(e.tag);
          } else {
            _physicalEntities[{e.dim, tagUnnumberedEntities++}].listEntities.push_back(e.tag);
          }
        }

        if(e.numPhysicalTags > 1) {
          feInfo("Geometric entity with (dim,tagGmsh) = (%d,%d) is part of %d physical entities :",
                 e.dim, e.tagGmsh, e.numPhysicalTags);
          for(int ii = 0; ii < e.numPhysicalTags; ++ii) {
            feInfo("    is part of physical %s",
                   _physicalEntities[{e.dim, e.physicalTags[ii]}].name.data());
          }
          return feErrorMsg(FE_STATUS_ERROR,
                            "Geometric entity can only be part of a maximum of one "
                            "physical entity (named domain). Overlapping domains are not allowed.");
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

    //
    // Read nodes
    //
    else if(buffer == "$Nodes")
    {
      int numEntityBlocks, countVertex = 0;
      double x, y, z;
      // numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
      input >> numEntityBlocks >> _nVertices;
      getline(input, buffer);
      _vertices.resize(_nVertices);
      _sequentialNodeToGmshNode.resize(_nVertices);
      for(int i = 0; i < numEntityBlocks; ++i) {
        int entityDim, numNodes;
        // entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
        input >> entityDim >> ph2 >> ph3 >> numNodes;
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
      if((int)_verticesMap.size() != _nVertices) {
        return feErrorMsg(FE_STATUS_ERROR, "Vertices indices are not unique.");
      }
    }

    //
    // Read elements
    //
    else if(buffer == "$Elements")
    {
      int numEntityBlocks, serialNumber, maxElementTag;
      // numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
      input >> numEntityBlocks >> _nElm >> ph1 >> maxElementTag;

      // Total number of elements of all dimensions, unused for now
      _nTotalElm = _nElm;

      int numEdges = 0;
      int numTriangles = 0;
      int numTetrahedra = 0;

      // Starting facets numbering at 1 to distinguish + or - in the connectivity
      int tagEdges = 1; 
      int tagEdgesDebug = 1;
      int tagTrianglesDebug = 1;

      // Sets of elements, transformed in vectors after all elements are added
      std::set<Edge, EdgeLessThan> allEdges;
      // std::set<Triangle, TriangleLessThan> allTriangles;
      // std::set<Tetrahedron, TetrahedronLessThan> _allTetrahedra;

      int countElems = 0;
      _nVerticesWithNoPhysical = 0;

      getline(input, buffer);

      for(int i = 0; i < numEntityBlocks; ++i)
      {
        //
        // Step 1 : read first line of entity (elementType) and allocate the connectivity tables
        //
        int numElementsInBlock, elemType, entityTag, entityDim;
        // entityDim(int) entityTag(int) elementType(int) numElementsInBlock(size_t)
        input >> entityDim >> entityTag >> elemType >> numElementsInBlock;

        std::pair<int, int> p = {entityDim, entityTag};
        bool printNodeWarning = false;
        _entities[p].nElm = numElementsInBlock;

        // Element connectivity:
        // Each element in the mesh is assigned a unique tag,
        // regardless of its dimension. There may be gaps in Gmsh's element numbering,
        // so a new sequential tag (countElems) is used instead.
        _entities[p].connecElem.resize(numElementsInBlock);

        // These are overwritten afterwards if necessary
        _entities[p].nNodePerElem = 0;
        _entities[p].nEdgePerElem = 0;
        _entities[p].nTriPerElem = 0;
        _entities[p].nTetPerElem = 0;

        // Determine geometric feSpace and allocate nodes connectivity based on elemType
        switch(elemType) {
          case 1:
            if(curved) {
              // We assume the mesh is made of the same elements,
              // hence there should not be 2-node edges if the
              // mesh is curved.
              return feErrorMsg(
                FE_STATUS_ERROR,
                "Curved mesh is enabled, but there are 2-node edges"
                " in the mesh. For now curved meshes can only have high-order edges.");
            }
            [[gnu::fallthrough]]; //  2-node line
          case 8:
            if(curved && _entities[p].interp == geometricInterpolant::NONE) {
              _entities[p].interp = geometricInterpolant::LINEP2;
            }
            [[gnu::fallthrough]]; //  3-node line (2nd order)
          case 26:
            if(curved && _entities[p].interp == geometricInterpolant::NONE) {
              _entities[p].interp = geometricInterpolant::LINEP3;
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
                _entities[p].interp = geometricInterpolant::LINEP1;
              }
              _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType - 1] : 2;
              _entities[p].nEdgePerElem = 1;
              _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
              _entities[p].connecEdges.resize(numElementsInBlock);
              break;
            }
          case 2: //  3-node triangle

            if(curved) {
              // We assume the mesh is made of the same elements,
              // hence there should not be 3-node triangles if the
              // mesh is curved.
              return feErrorMsg(
                FE_STATUS_ERROR,
                "Curved mesh is enabled, but there are 3-node triangles"
                " in the mesh. For now curved meshes can only have high-order triangles.");
            }

            [[gnu::fallthrough]];
          case 9:
            if(curved && _entities[p].interp == geometricInterpolant::NONE) {
              _entities[p].interp = geometricInterpolant::TRIP2;
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
            // For now we assume the mesh is only P1 or only P2
            if(!curved && elemType != 2) {
              return feErrorMsg(
                FE_STATUS_ERROR,
                "Curved mesh is disabled, but there are 6-node triangles"
                " in the mesh. For now high-order meshes must be read with curve = true.");
            }

            if(!curved) {
              _entities[p].interp = geometricInterpolant::TRIP1;
            }
            _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType - 1] : 3;
            _entities[p].nEdgePerElem = 3;
            _entities[p].nTriPerElem = 1;
            _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
            _entities[p].connecEdges.resize(numElementsInBlock * 3);
            _entities[p].connecEdgesDebug.resize(numElementsInBlock * 3);
            _entities[p].connecTri.resize(numElementsInBlock);
            _entities[p].connecTrianglesDebug.resize(numElementsInBlock);
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
            // _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType - 1] : 4;
            // _entities[p].nEdgePerElem = 4;
            // _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
            // _entities[p].connecEdges.resize(numElementsInBlock * _entities[p].nEdgePerElem);
            return feErrorMsg(FE_STATUS_ERROR,
                              "Geometric entity %d : quadrangles are not supported.\n",
                              entityTag);
            break;
          }
          case 4: //   4-node tetrahedron
            _entities[p].interp = geometricInterpolant::TETP1;
            _entities[p].nNodePerElem = nodes_of_gmsh_element[elemType - 1];
            _entities[p].nEdgePerElem = 6;
            _entities[p].nTriPerElem = 4;
            // Number of *boundary* tet per element is 0, and attribute should be removed in general
            // _entities[p].nTetPerElem = 0;
            _entities[p].connecNodes.resize(numElementsInBlock * _entities[p].nNodePerElem);
            _entities[p].connecEdges.resize(numElementsInBlock * 6);
            _entities[p].connecTri.resize(numElementsInBlock * 4);
            _entities[p].connecTrianglesDebug.resize(numElementsInBlock * _entities[p].nTriPerElem);
            // _entities[p].connecTet.resize(numElementsInBlock);
            break;
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
                              "Geometric entity %d : high-order tetrahedra are not supported.\n",
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
                              "Geometric entity %d : hexahedra are not supported.\n",
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
            _entities[p].interp = geometricInterpolant::POINTP0;
            _entities[p].nNodePerElem = 1;
            _entities[p].nEdgePerElem = 0;
            _entities[p].connecNodes.resize(1);
            _entities[p].connecElem.resize(1);
            break;
          }
          default:
            return feErrorMsg(FE_STATUS_ERROR, "Unsupported Gmsh element type.");
        } // switch(elemType)

        //
        // Step 2 : read remaining lines of entity and create elements
        //
        std::map<int, int>::const_iterator it;
        for(int iElm = 0; iElm < numElementsInBlock; ++iElm) {
          // int nElemNodes = nodes_of_gmsh_element[elemType-1];
          int nElemNodes = _entities[p].nNodePerElem;
          std::vector<int> elemNodesGmsh(nElemNodes);
          std::vector<int> elemNodes(nElemNodes);
          // elementTag(size_t) nodeTag(size_t) ...
          input >> serialNumber; // Unused

          // Element connectivity, i.e., unique global tag assigned to
          // all elements regardless of their dimension.
          // "Point" entities have connecElem of size 1, assigned a bit further
          if(elemType != 15) {
            _entities[p].connecElem[iElm] = countElems++;
          }

          for(int j = 0; j < nElemNodes; ++j) {
            input >> ph1;

            // FIXME: This is slow! Each vertex of each element is searched :/
            // Check that vertex exists in the map
            it = _verticesMap.find(ph1);
            if(it == _verticesMap.end()) {
              return feErrorMsg(FE_STATUS_ERROR, "Element node %d does not exist /-:", ph1);
            }

            // The Gmsh number stored in ph1 is used to create the edges
            elemNodesGmsh[j] = ph1;
            // The node number (0...nNode) is used to create the connectivity
            elemNodes[j] = it->second;
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
                if(reversed || _entities[p].physicalTags[0] < 0) {
                  // Physical entity is reversed
                  _entities[p].connecNodes[nElemNodes * iElm + j] = elemNodes[(nElemNodes - 1) - j];
                } else {
                  // Physical entity is well oriented
                  _entities[p].connecNodes[nElemNodes * iElm + j] = elemNodes[j];
                }
              }

              // FIXME: We should only do this when maxDim = 1 actually
              // This will in general be overwritten when transferring cnc from higher-dimensional entities
              _entities[p].connecEdges[iElm] = numEdges;
              numEdges++;
              break;
            }
            case 2: //  3-node triangle

              // FIXME: We should only do this when maxDim = 2 actually
              // In general, we construct the mesh elements starting from the highest-dimensional elements.
              // Then lower dimensional elements are identified from among the existing set of triangles or edges.
              // Cleaner version:
              {
                Vertex *v0 = &_vertices[_verticesMap[elemNodesGmsh[0]]];
                Vertex *v1 = &_vertices[_verticesMap[elemNodesGmsh[1]]];
                Vertex *v2 = &_vertices[_verticesMap[elemNodesGmsh[2]]];

                Triangle *tri;
                if(reversed || _entities[p].physicalTags[0] < 0) {
                  tri = new Triangle(v0, v2, v1, numTriangles);
                } else {
                  tri = new Triangle(v0, v1, v2, numTriangles);
                }

                // Create the boundary edges of tri
                // They are added to the set if they (or their reverse) do not exist yet
                tri->createBoundary(allEdges, tagEdgesDebug);
                _allTriangles.insert(*tri);

                //
                // Update connectivities on this geometric entity:
                //

                // Vertices
                // TODO, see below

                // Edges
                // Add the edges (with orientation) to the edge connectivity of this entity
                for(int k = 0; k < 3; ++k) {
                  _entities[p].connecEdgesDebug[3 * iElm + k] = tri->getEdgeOrientation(k) * tri->getEdge(k)->getTag();
                }

                // Triangles (trivial) (needed?)
                _entities[p].connecTrianglesDebug[iElm] = numTriangles;
                numTriangles++;

                // Do not break while debugging (-:
                // feWarning("Creating both data structures for 3-node triangles for debug");
                // break;
              }
              [[gnu::fallthrough]];
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
                // Reverse if required or if Physical Entity is negative
                if(reversed || _entities[p].physicalTags[0] < 0) {
                  _entities[p].connecNodes[nElemNodes * iElm + 0] = elemNodes[0];
                  _entities[p].connecNodes[nElemNodes * iElm + 1] = elemNodes[2];
                  _entities[p].connecNodes[nElemNodes * iElm + 2] = elemNodes[1];
                  _entities[p].connecNodes[nElemNodes * iElm + 3] = elemNodes[5];
                  _entities[p].connecNodes[nElemNodes * iElm + 4] = elemNodes[4];
                  _entities[p].connecNodes[nElemNodes * iElm + 5] = elemNodes[3];
                } else {
                  _entities[p].connecNodes[nElemNodes * iElm + 0] = elemNodes[0];
                  _entities[p].connecNodes[nElemNodes * iElm + 1] = elemNodes[1];
                  _entities[p].connecNodes[nElemNodes * iElm + 2] = elemNodes[2];
                  _entities[p].connecNodes[nElemNodes * iElm + 3] = elemNodes[3];
                  _entities[p].connecNodes[nElemNodes * iElm + 4] = elemNodes[4];
                  _entities[p].connecNodes[nElemNodes * iElm + 5] = elemNodes[5];
                }
              } else {
                // Reverse if required or if Physical Entity is negative
                if(reversed || _entities[p].physicalTags[0] < 0) {
                  _entities[p].connecNodes[nElemNodes * iElm + 0] = elemNodes[0];
                  _entities[p].connecNodes[nElemNodes * iElm + 1] = elemNodes[2];
                  _entities[p].connecNodes[nElemNodes * iElm + 2] = elemNodes[1];
                } else {
                  for(int j = 0; j < nElemNodes; ++j) {
                    _entities[p].connecNodes[nElemNodes * iElm + j] = elemNodes[j];
                  }
                }
              }

              // Construct the triangle edges :
              Vertex *v0, *v1, *vMid;
              for(int k = 0; k < 3; ++k) {
                if(reversed || _entities[p].physicalTags[0] < 0) {
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
                } else {
                  v0 = &_vertices[_verticesMap[elemNodesGmsh[k]]];
                  v1 = &_vertices[_verticesMap[elemNodesGmsh[(k + 1) % 3]]];
                }
                std::pair<std::set<Edge, EdgeLessThan>::iterator, bool> ret;
                Edge e(v0, v1, tagEdges, _entities[p].physicalTags[0]);
                ret = _edges.insert(e);
                if(ret.second) {
                  // Edge was added to the set : tagEdges is added to connecEdges
                  _entities[p].connecEdges[3 * iElm + k] = tagEdges++;

                  if(curved) {
                    // Compute displacement of midnode (alpha)
                    vMid = &_vertices[_verticesMap[elemNodesGmsh[k + 3]]];
                    double xMid = (v0->x() + v1->x()) / 2.;
                    double yMid = (v0->y() + v1->y()) / 2.;
                    // double normAlpha = sqrt((vMid->x() - xMid) * (vMid->x() - xMid) +
                    //                         (vMid->y() - yMid) * (vMid->y() - yMid));
                    auto edgeIt = _edges.find(e);
                    _edge2midnode[&(*edgeIt)] = vMid;
                    _edge2alpha[&(*edgeIt)][0] = vMid->x() - xMid;
                    _edge2alpha[&(*edgeIt)][1] = vMid->y() - yMid;
                  }

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
              return feErrorMsg(FE_STATUS_ERROR, "Quadrangles are not supported.\n");
            }
            case 4: //   4-node tetrahedron
            {
              Vertex *v0 = &_vertices[_verticesMap[elemNodesGmsh[0]]];
              Vertex *v1 = &_vertices[_verticesMap[elemNodesGmsh[1]]];
              Vertex *v2 = &_vertices[_verticesMap[elemNodesGmsh[2]]];
              Vertex *v3 = &_vertices[_verticesMap[elemNodesGmsh[3]]];

              Tetrahedron *tet;
              if(reversed || _entities[p].physicalTags[0] < 0) {
                tet = new Tetrahedron(v0, v3, v2, v1, numTetrahedra);
              } else {
                tet = new Tetrahedron(v0, v1, v2, v3, numTetrahedra);
              }

              // Create the boundary facets of tet
              // They are added to the set if they (or their reverse) do not exist yet
              // This also creates the mesh edges, as boundaries of the facets.
              tet->createBoundary(_allTriangles, tagTrianglesDebug, allEdges, tagEdgesDebug);

              // Create the boundary edges of the facets of tet
              _allTetrahedra.insert(*tet);

              //
              // Update connectivities on this geometric entity:
              //

              // Vertices
              if(reversed || _entities[p].physicalTags[0] < 0) {
                _entities[p].connecNodes[4 * iElm + 0] = elemNodes[0];
                _entities[p].connecNodes[4 * iElm + 1] = elemNodes[3];
                _entities[p].connecNodes[4 * iElm + 2] = elemNodes[2];
                _entities[p].connecNodes[4 * iElm + 3] = elemNodes[1];
              } else {
                for(int j = 0; j < 4; ++j) {
                  _entities[p].connecNodes[4 * iElm + j] = elemNodes[j];
                }
              }

              // Edges
              // Add the facets (with orientation) to the triangle connectivity of this entity
              // TODO
              feWarning("Add edge connectivity for 3D entities if needed");

              // Triangles
              // Add the facets (with orientation) to the triangle connectivity of this entity
              for(int k = 0; k < 4; ++k) {
                _entities[p].connecTri           [4 * iElm + k] = tet->getFacetOrientation(k) * tet->getFacet(k)->getTag();
                _entities[p].connecTrianglesDebug[4 * iElm + k] = tet->getFacetOrientation(k) * tet->getFacet(k)->getTag();
              }

              // Tet (trivial) (needed? redundant with connecElem?)
              // _entities[p].connecTet[iElm] = numTetrahedra;
              numTetrahedra++;
              
              feWarning("Creating data structures for 4-node tetrahedron for debug");
              break;
            }
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
              return feErrorMsg(FE_STATUS_ERROR, "High-order tetrahedra are not supported.\n");
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
              return feErrorMsg(FE_STATUS_ERROR, "Hexahedra are not supported.\n");
            }
            case 15: // 1-node point
            {
              // We only allow a geometric Point to have a single node
              if(!printNodeWarning) {
                _entities[p].nElm = 1;
                _entities[p].connecElem[0] = countElems++;
                _entities[p].connecNodes[0] = elemNodes[0];
                printNodeWarning = true;
              } else {
                _nVerticesWithNoPhysical++;
                feWarning(
                  "More than one node are present in geometric entity (dim = %d, gmshTag = %d). "
                  "Only the first node was added, the others were discarded.",
                  _entities[p].dim, _entities[p].tagGmsh);
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

      feInfo("There are %d mesh edges", allEdges.size());
      feInfo("There are %d mesh triangles", _allTriangles.size());
      feInfo("There are %d mesh tetrahedra", _allTetrahedra.size());

      _tetrahedra.assign(_allTetrahedra.begin(), _allTetrahedra.end());

    } // if buffer = "Elements"
  } // while input

  return FE_STATUS_OK;
}

//
// Read a mesh in Gmsh's .msh format.
// The elements and connectivity tables are created here:
// The highest dimensional elements are created and lower dimensional
// elements are created as their boundary.
//
// For instance, in a 2D triangulation, the mesh edges are created as
// boundaries of the triangles. If there are 1D physical entities in the mesh,
// for example to apply boundary conditions, their 1D elements are identified
// to the edges of existing triangles.
//
// This assumes that there are no "stray" 1D physical entities (i.e., the mesh
// is the triangulation of a manifold).
//
// For a 3D mesh, the tetrahedra are created, then triangles on 2D physical entities
// are identified to existing facets of tetrahedra, then 1D physical entities (if any)
// are identified to the edges of the facets.
//
feStatus feMesh2DP1::readGmsh(const std::string meshName, const bool curved, const bool reversed,
                              const mapType physicalEntitiesDescription)
{
  _gmshVersion = 0;
  std::filebuf fb;

  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "MESH:");
  feInfoCond(FE_VERBOSE > -1, "\t\tReading mesh file: %s", meshName.data());

  std::ifstream f(meshName.data());
  if(!f.good()) {
    return feErrorMsg(FE_STATUS_READ_ERROR, "Mesh file does not exist.");
  }

  _ID = meshName;

  if(fb.open(meshName, std::ios::in)) {
    std::istream input(&fb);

    std::string buffer;
    int dsize;
    getline(input, buffer);
    input >> _gmshVersion >> _isBinary >> dsize;

    if(_isBinary) {
      return feErrorMsg(FE_STATUS_ERROR, "Only reading ASCII files.");
    }

    //
    // Read mesh file, create elements and connectivity tables of the geometric entities
    //
    if(_gmshVersion == 2.2) {
      feCheck(readMsh2(input, curved, reversed, physicalEntitiesDescription));
    } else if(_gmshVersion >= 4) {
      feCheck(readMsh4(input, curved, reversed, physicalEntitiesDescription));
    } else {
      return feErrorMsg(FE_STATUS_ERROR, "Mesh version should be 2.2 or 4.1 or higher.");
    }

    fb.close();
  } // if fb.open

  feInfoCond(FE_VERBOSE > VERBOSE_THRESHOLD,
             "\t\tList of geometric entities (belonging to physical) with their attributes:");
  for(auto const &x : _entities) {
    entity e = x.second;
    if(e.numPhysicalTags > 0) {
      feInfoCond(FE_VERBOSE > VERBOSE_THRESHOLD,
                 "\t\t\tEntity #%2d : dim = %1d - gmshTag = %3d - nElm = %6d - nNodePerElem = %2d - "
                 "geometric interpolant : %10s "
                 "- part of physical groups : %2d",
                 e.tag, e.dim, e.tagGmsh, e.nElm, e.nNodePerElem, toString(e.interp).data(),
                 (e.numPhysicalTags > 0) ? e.physicalTags[0] : -1);
    }
  }

  // Check that the geometric entities defined on a physical entity are compatible.
  // Each geometric entity has a sequential tag in [0, nEntitiesOfThisDimension].
  // The tag is reset for each dimension and is thus different from the Gmsh entity tag.
  for(auto const &x : _entities) {
    entity e = x.second;
    for(auto &p : _physicalEntities) {
      physicalEntity &pE = p.second;
      bool error = false;
      for(auto ent : pE.listEntities) {

        // If geometric entity is in tested physical entity, transfer the information.
        if(ent == e.tag && pE.dim == e.dim)
        {
          pE.nElm += e.nElm;

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

          // Assign and check nTriPerElem
          if(pE.nTriPerElem == -1)
            pE.nTriPerElem = e.nTriPerElem;
          else {
            if(pE.nTriPerElem != e.nTriPerElem) error = true;
          }

          // Assign and check nTetPerElem
          if(pE.nTetPerElem == -1)
            pE.nTetPerElem = e.nTetPerElem;
          else {
            if(pE.nTetPerElem != e.nTetPerElem) error = true;
          }

          // Assign and check interpolant for the geometry
          if(pE.interp == geometricInterpolant::NONE) {
            pE.interp = e.interp;
          } else {
            if(pE.interp != e.interp) error = true;
          }
        }
      }
      if(error) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Geometric entities with conflicting properties are defined on physical entity "
                          "%s.\n",
                          pE.name.data());
      }
    }
  }

  // Transfer the connectivity tables from the geometric entities
  // to the physical entities (named domains).
  for(auto &p : _physicalEntities)
  {
    physicalEntity &pE = p.second;
    pE.connecElem.resize(pE.nElm);
    pE.connecNodes.resize(pE.nElm * pE.nNodePerElem);
    pE.connecEdges.resize(pE.nElm * pE.nEdgePerElem);
    pE.connecTri.resize(pE.nElm * pE.nTriPerElem);

    int countElm = 0;
    for(auto const &x : _entities) {
      entity e = x.second;
      for(auto ent : pE.listEntities) {
        if(ent == e.tag) {
          for(int iElm = 0; iElm < e.nElm; ++iElm)
          {
            // Element connectivity
            pE.connecElem[countElm] = e.connecElem[iElm];

            // Connectivity of vertices
            for(int j = 0; j < e.nNodePerElem; ++j) {
              pE.connecNodes[e.nNodePerElem * countElm + j] =
                e.connecNodes[e.nNodePerElem * iElm + j];
            }

            // Connectivity of boundary edges
            for(int j = 0; j < e.nEdgePerElem; ++j) {
              pE.connecEdges[e.nEdgePerElem * countElm + j] =
                e.connecEdges[e.nEdgePerElem * iElm + j];
              // feInfo("V1: %s : connecEdges[%d] = %d", pE.name.data(),e.nEdgePerElem * countElm + j, pE.connecEdges[e.nEdgePerElem * countElm + j]);
            }

            // Connectivity of boundary triangles (facets)
            for(int j = 0; j < e.nTriPerElem; ++j) {
              pE.connecTri[e.nTriPerElem * countElm + j] =
                e.connecTri[e.nTriPerElem * iElm + j];
            }

            ++countElm;
          }
        }
      }
    }
  }

  // Transfer connectivity to boundary elements.
  // This overwrites the existing connectivity tables
  // of n-dimensional physical entities with n < maxDim.
  //
  // Maybe we could *only* create the highest-dimensional elements,
  // then transfer the connectivities to boundary entities?
  //
  // For now, we loop over lower dimensional entities and
  // recreate the Edges and Triangles then find them in the sets,
  // which does not seem optimal...
  int maxDim = 0;
  for(auto &pE : _physicalEntities) {
    maxDim = fmax(maxDim, pE.second.dim);
  }

  for(auto &p : _physicalEntities)
  {
    physicalEntity &pE = p.second;

    //
    // Transfer edge connectivity to boundary edges, i.e.,
    // set the edge connectivity of 1D phys. entities with the
    // edge tags created as boundaries of 2D phys. entities.
    //
    if(pE.dim == 1 && pE.dim < maxDim)
    {
      // Create and find edge in the set of all edges
      Vertex *v0, *v1;
      for(int i = 0; i < pE.nElm; ++i) {
        v0 = &_vertices[pE.connecNodes[pE.nNodePerElem * i]];
        v1 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 1]];
        Edge e(v0, v1);
        auto it = _edges.find(e);
        if(it != _edges.end()) {
          pE.connecEdges[i] = it->getTag();
          if(v0->getTag() != it->getTag(0)) {
            // Boundary edges is reversed w.r.t. the stored edge,
            // which should not happen if we have correctly respected
            // Gmsh's ordering.
            return feErrorMsg(FE_STATUS_ERROR,
              "Error when transferring connectivity tables to lower-dimensional physical entity:\n"
              "Boundary edge (%d,%d) on entity %s with dimension %d is negatively oriented.\n"
              "Check that the element ordering on the boundary is consistent with the higher-dimensional entity.",
              pE.name.data(), pE.dim, v0->getTag(), v1->getTag());
          }
        } else {
          // Edge is not in the set: this should not happen
          // unless this entity is *not* the boundary of a higher-dimensional entity
          // (i.e., if the mesh is not a manifold).
          return feErrorMsg(FE_STATUS_ERROR,
            "Error when transferring connectivity tables to lower-dimensional physical entity:\n"
            "Boundary edge (%d,%d) on entity %s with dimension %d was not found in the set of mesh edges.\n"
            "Check that this entity is indeed the boundary of higher-dimensional physical entity.\n"
            "(i.e., check that the mesh is a manifold)",
            pE.name.data(), pE.dim, v0->getTag(), v1->getTag());
        }
      }

      // for(auto val : pE.connecEdges) {
      //   feInfo("V2: %s : connecEdges[] = %d", pE.name.data(), val);
      // } 
    }

    //
    // Transfer triangle connectivity to boundary edges, i.e.,
    // set the tri connectivity of 2D phys. entities with the
    // triangles tags created as facets of 3D phys. entities.
    //
    if(pE.dim == 2 && pE.dim < maxDim)
    {
      // Create and find triangle in the set of all edges
      Vertex *v0, *v1, *v2;
      for(int i = 0; i < pE.nElm; ++i) {
        v0 = &_vertices[pE.connecNodes[pE.nNodePerElem * i]];
        v1 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 1]];
        v2 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 2]];
        Triangle t(v0, v1, v2);
        auto it = _allTriangles.find(t);
        if(it != _allTriangles.end()) {
          pE.connecEdges[i] = it->getTag();
          bool sameOrder = (v0->getTag() == it->getTag(0) && v1->getTag() == it->getTag(1) && v2->getTag() == it->getTag(2));
          if(!sameOrder) {
            // Boundary triangle is reversed w.r.t. the stored triangle,
            // which should not happen if we have correctly respected
            // Gmsh's ordering.
            return feErrorMsg(FE_STATUS_ERROR,
              "Error when transferring connectivity tables to lower-dimensional physical entity:\n"
              "Boundary triangle (%d,%d,%d) on entity %s with dimension %d is negatively oriented.\n"
              "Check that the element ordering on the boundary is consistent with the higher-dimensional entity.",
              pE.name.data(), pE.dim, v0->getTag(), v1->getTag(), v2->getTag());
          }
        } else {
          // Triangle is not in the set: this should not happen
          // unless this entity is *not* the boundary of a higher-dimensional entity
          // (i.e., if the mesh is not a manifold).
          return feErrorMsg(FE_STATUS_ERROR,
            "Error when transferring connectivity tables to lower-dimensional physical entity:\n"
            "Boundary edge (%d,%d,%d) on entity %s with dimension %d was not found in the set of mesh edges.\n"
            "Check that this entity is indeed the boundary of higher-dimensional physical entity.\n"
            "(i.e., check that the mesh is a manifold)",
            pE.name.data(), pE.dim, v0->getTag(), v1->getTag(), v2->getTag());
        }
      }
    }
  }

  // Assign pointer to geometric interpolation space
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(pE.interp == geometricInterpolant::POINTP0) {
      pE.geometry = geometryType::POINT;
      pE.geoSpace = new feSpace1DP0("xyz");
    } else if(pE.interp == geometricInterpolant::LINEP1) {
      pE.geometry = geometryType::LINE;
      pE.geoSpace = new feSpace1DP1("xyz");
    } else if(pE.interp == geometricInterpolant::LINEP2) {
      pE.geometry = geometryType::LINE;
      pE.geoSpace = new feSpace1DP2("xyz");
    } else if(pE.interp == geometricInterpolant::LINEP3) {
      pE.geometry = geometryType::LINE;
      pE.geoSpace = new feSpace1DP3("xyz");
    } else if(pE.interp == geometricInterpolant::TRIP1) {
      pE.geometry = geometryType::TRI;
      pE.geoSpace = new feSpaceTriP1("xyz");
    } else if(pE.interp == geometricInterpolant::TRIP2) {
      pE.geometry = geometryType::TRI;
      pE.geoSpace = new feSpaceTriP2("xyz");
    } else if(pE.interp == geometricInterpolant::TRIP3) {
      pE.geometry = geometryType::TRI;
      pE.geoSpace = new feSpaceTriP3("xyz");
    } else if(pE.interp == geometricInterpolant::TETP1) {
      pE.geometry = geometryType::TET;
      // pE.geoSpace = new feSpaceTetP1("xyz");
      pE.geoSpace = new feSpaceTetPn(1, "xyz");
    } else {
      return feErrorMsg(FE_STATUS_ERROR,
                        "Unknown geometric connectivity interpolant \"%s\" on domain \"%s\".",
                        toString(pE.interp).data(), pE.name.data());
    }
  }

  feInfoCond(FE_VERBOSE > VERBOSE_THRESHOLD, "\t\tList of physical groups with their attributes :");
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    feInfoCond(
      FE_VERBOSE > VERBOSE_THRESHOLD,
      "\t\t\tPhysical %15s : dim = %1d - nEntities = %2d - nElm = %9d - nNodePerElem = %2d - "
      "nEdgePerElem = %1d - geometric interpolant : %10s",
      pE.name.data(), pE.dim, pE.listEntities.size(), pE.nElm, pE.nNodePerElem, pE.nEdgePerElem,
      toString(pE.interp).data());
  }

  // CONTINUER ICI
  //
  // Finally, create geometric connectivities
  //
  // FIXME: Name is misleading as these are connectivities on Physical entities,
  // not on geometric entities. Physical entities should be replaced by feCncgeo entirely.
  int nCncGeo = 0;
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;

    feCncGeo *cnc =
      new feCncGeo(nCncGeo, pE.dim, maxDim, pE.nNodePerElem, pE.nElm, pE.nEdgePerElem, pE.name, pE.geometry,
                   pE.interp, pE.geoSpace, pE.connecNodes, pE.connecElem, pE.connecEdges, pE.connecTri);

    _cncGeo.push_back(cnc);
    _cncGeoMap[pE.name] = nCncGeo;
    cnc->getFeSpace()->setCncGeoTag(nCncGeo++);

    // // Create a vector of vertices for each Physical Entity
    // std::vector<Vertex> entityVertices(pE.connecNodes.size());
    // for(size_t i = 0; i < pE.connecNodes.size(); ++i) {
    //   entityVertices[i] = _vertices[pE.connecNodes[i]];
    // }
    // _verticesPerCnc[cnc] = entityVertices;
  }

  _nEdges = _edges.size();
  _dim = maxDim;

  _nInteriorElm = 0;
  _nBoundaryElm = 0;
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(pE.dim == _dim) _nInteriorElm += pE.nElm;
    if(pE.dim == _dim - 1) _nBoundaryElm += pE.nElm;
  }

  // Add all interior elements (for now only P1 or P2 triangles) to the _elements vector
  // Physical entities are non-overlapping so we just loop over them and elements shouldn't be
  // duplicated
  _elements.resize(_nInteriorElm);
  int cnt = 0;
  for(auto &p : _physicalEntities) {
    physicalEntity &pE = p.second;
    if(_dim == 2 && pE.dim == _dim) {
      Vertex *v0, *v1, *v2, *v3, *v4, *v5;
      for(int i = 0; i < pE.nElm; ++i) {
        if(pE.interp == geometricInterpolant::TRIP1) {
          v0 = &_vertices[pE.connecNodes[pE.nNodePerElem * i]];
          v1 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 1]];
          v2 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 2]];
          // The triangle tag is i because it is local to the pE's cncgeo, hence the tag is not
          // unique...
          _elements[cnt] = new Triangle(v0, v1, v2, cnt, i, pE.tag);
        } else if(pE.interp == geometricInterpolant::TRIP2) {
          v0 = &_vertices[pE.connecNodes[pE.nNodePerElem * i]];
          v1 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 1]];
          v2 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 2]];
          v3 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 3]];
          v4 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 4]];
          v5 = &_vertices[pE.connecNodes[pE.nNodePerElem * i + 5]];
          _elements[cnt] = new TriangleP2(v0, v1, v2, v3, v4, v5, cnt, i, pE.tag);
        } else {
          return feErrorMsg(FE_STATUS_ERROR,
                            "Could not create a Triangle for physical entity \"%s\" with geometric "
                            "interpolant \"%s\"",
                            pE.name.data(), toString(pE.interp).data());
        }

        cnt++;
      }
    }
  }

  // TODO : Also add boundaryElements when necessary

  // Create an RTree storing the elements of highest dimension
  // Elements are hardcoded "Triangle" for now
  if(maxDim == 2)
  {
    double min2[2], max2[2];
    for(size_t iElm = 0; iElm < _elements.size(); ++iElm) {
      Triangle *t = _elements[iElm];
      SBoundingBox3d bbox;
      for(int j = 0; j < t->getNumVertices(); ++j) {
        Vertex *v = t->getVertex(j);
        SPoint3 pt(v->x(), v->y(), v->z());
        bbox += pt;
      }
      // _rtree.Insert((double *)(bbox.min()), (double *)(bbox.max()), t->getTag());
      min2[0] = bbox.min()[0];
      min2[1] = bbox.min()[1];
      max2[0] = bbox.max()[0];
      max2[1] = bbox.max()[1];
      _rtree2d.Insert(min2, max2, t->getTag());
    }
  }

  if(maxDim == 3)
  {
    double min3[3], max3[3];
    for(auto &tet : _tetrahedra) {
      SBoundingBox3d bbox;
      for(int j = 0; j < tet.getNumVertices(); ++j) {
        Vertex *v = tet.getVertex(j);
        SPoint3 pt(v->x(), v->y(), v->z());
        bbox += pt;
      }
      // _rtree.Insert((double *)(bbox.min()), (double *)(bbox.max()), t->getTag());
      min3[0] = bbox.min()[0];
      min3[1] = bbox.min()[1];
      min3[2] = bbox.min()[2];
      max3[0] = bbox.max()[0];
      max3[1] = bbox.max()[1];
      max3[2] = bbox.max()[2];
      _rtree.Insert(min3, max3, tet.getTag());
    }
  }

  // Assign a pointer to this mesh to each of its geometric connectivities and their fespace
  for(feCncGeo *cnc : _cncGeo) {
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }

  if(_cncGeo.size() > 0) {
    _cncGeo[1]->printColoring("coloring.pos");
    _cncGeo[1]->printColoring2("coloring2.pos");
  }

  // Assign to each geometric FE space a pointer to its connectivity
  for(feCncGeo *cnc : _cncGeo) {
    cnc->getFeSpace()->setCncPtr(cnc);
  }
  _nCncGeo = _cncGeo.size();

  // Add edges pointer to a vector for faster access (is it true tho?)
  _edgesVec.resize(_edges.size());
  cnt = 0;
  for(std::set<Edge, EdgeLessThan>::iterator it = _edges.begin(); it != _edges.end(); ++it) {
    _edgesVec[cnt++] = &(*it);
  }

  // Ordered map of edges for visualization lookup
  for(auto *edge : _edgesVec) {
    _edgesMap[edge->getTag()] = edge;
  }

  return FE_STATUS_OK;
}