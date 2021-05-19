#include "feMesh.h"
#include "feElement.h"
#include "feLine.h"
#include "feSpaceTriangle.h"

#include <iostream>
#include <fstream>

// number of nodes for each type of Gmsh elements, type is the index of the array + 1
 int nodes_of_gmsh_element[] =
 {
    2,   // 2-node line.
    3,   // 3-node triangle.
    4,   // 4-node quadrangle.
    4,   // 4-node tetrahedron.
    8,   // 8-node hexahedron.
    6,   // 6-node prism.
    5,   // 5-node pyramid.
    3,   /* 3-node second order line (2 nodes associated with the
            vertices and 1 with the edge). */
    6,   /* 6-node second order triangle (3 nodes associated with the
            vertices and 3 with the edges). */
    9,   /* 9-node second order quadrangle (4 nodes associated with the
            vertices, 4 with the edges and 1 with the face). */
    10,  /* 10-node second order tetrahedron (4 nodes associated with
            the vertices and 6 with the edges). */
    27,  /* 27-node second order hexahedron (8 nodes associated with the
            vertices, 12 with the edges, 6 with the faces and 1 with the
            volume). */
    18,  /* 18-node second order prism (6 nodes associated with the
            vertices, 9 with the edges and 3 with the quadrangular
            faces). */
    14,  /* 14-node second order pyramid (5 nodes associated with the
            vertices, 8 with the edges and 1 with the quadrangular
            face). */
    1,   // 1-node point.
    8,   /* 8-node second order quadrangle (4 nodes associated with the
            vertices and 4 with the edges). */
    20,  /* 20-node second order hexahedron (8 nodes associated with the
            vertices and 12 with the edges). */
    15,  /* 15-node second order prism (6 nodes associated with the
            vertices and 9 with the edges). */
    13,  /* 13-node second order pyramid (5 nodes associated with the
            vertices and 8 with the edges). */
    9,   /* 9-node third order incomplete triangle (3 nodes associated
            with the vertices, 6 with the edges) */
    10,  /* 10-node third order triangle (3 nodes associated with the
            vertices, 6 with the edges, 1 with the face) */
    12,  /* 12-node fourth order incomplete triangle (3 nodes associated
            with the vertices, 9 with the edges) */
    15,  /* 15-node fourth order triangle (3 nodes associated with the
            vertices, 9 with the edges, 3 with the face) */
    15,  /* 15-node fifth order incomplete triangle (3 nodes associated
            with the vertices, 12 with the edges) */
    21,  /* 21-node fifth order complete triangle (3 nodes associated
            with the vertices, 12 with the edges, 6 with the face) */
    4,   /* 4-node third order edge (2 nodes associated with the
            vertices, 2 internal to the edge) */
    5,   /* 5-node fourth order edge (2 nodes associated with the
            vertices, 3 internal to the edge) */
    6,   /* 6-node fifth order edge (2 nodes associated with the
            vertices, 4 internal to the edge) */
    20,  /* 20-node third order tetrahedron (4 nodes associated with the
            vertices, 12 with the edges, 4 with the faces) */
    35,  /* 35-node fourth order tetrahedron (4 nodes associated with
            the vertices, 18 with the edges, 12 with the faces, and 1
            with the volume) */
    56,  /* 56-node fifth order tetrahedron (4 nodes associated with the
            vertices, 24 with the edges, 24 with the faces, and 4 with
            the volume) */
    -1,-1, /* unsupported tetrahedral types */
    -1,-1, /* unsupported polygonal and polyhedral types */
    16,  /* 16-node third order quadrilateral (4 nodes associated with
            the vertices, 8 with the edges, 4 wth the face) */
    25,  /* 25-node fourth order quadrilateral (4 nodes associated with
            the vertices, 12 with the edges, 9 wth the face) */
    36,  /* 36-node fifth order quadrilateral (4 nodes associated with
            the vertices, 16 with the edges, 16 wth the face) */
    -1,-1,-1, /* unsupported quadrilateral types */
    28,  /* 28-node sixth order complete triangle (3 nodes associated
            with the vertices, 15 with the edges, 10 with the face) */
    36,  /* 36-node seventh order complete triangle (3 nodes associated
            with the vertices, 18 with the edges, 15 with the face) */
    45,  /* 45-node eighth order complete triangle (3 nodes associated
            with the vertices, 21 with the edges, 21 with the face) */
    55,  /* 55-node ninth order complete triangle (3 nodes associated
            with the vertices, 24 with the edges, 28 with the face) */
    66,  /* 66-node tenth order complete triangle (3 nodes associated
            with the vertices, 27 with the edges, 36 with the face) */
    49,  /* 49-node sixth order quadrilateral (4 nodes associated with
            the vertices, 20 with the edges, 25 wth the face) */
    64,  /* 64-node seventh order quadrilateral (4 nodes associated with
            the vertices, 24 with the edges, 36 wth the face) */
    81,  /* 81-node eighth order quadrilateral (4 nodes associated with
            the vertices, 28 with the edges, 49 wth the face) */
    100, /* 100-node ninth order quadrilateral (4 nodes associated with
            the vertices, 32 with the edges, 64 wth the face) */
    121, /* 121-node tenth order quadrilateral (4 nodes associated with
            the vertices, 36 with the edges, 81 wth the face) */
    -1,-1,-1,-1,-1, /* unsupported triangular types */
    -1,-1,-1,-1,-1, /* unsupported quadrilateral types */
    7,   /* 7-node sixth order edge (2 nodes associated with the
            vertices, 5 internal to the edge) */
    8,   /* 8-node seventh order edge (2 nodes associated with the
            vertices, 6 internal to the edge) */
    9,   /* 9-node eighth order edge (2 nodes associated with the
            vertices, 7 internal to the edge) */
    10,  /* 10-node ninth order edge (2 nodes associated with the
            vertices, 8 internal to the edge) */
    11,  /* 11-node tenth order edge (2 nodes associated with the
            vertices, 9 internal to the edge) */
    -1,  /* unsupported linear types */
    -1,-1,-1, /* unsupported types */
    84,  /* 84-node sixth order tetrahedron (4 nodes associated with the
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
    -1,-1,-1,       /* undefined types */
    -1,-1,-1,-1,-1, /* unsupported tetrahedral types */
    -1,-1,-1,-1,-1,-1, /* unsupported types */
    40,  /* 40-node third order prism (6 nodes associated with the
            vertices, 18 with the edges, 14 with the faces, and 2 with
            the volume) */
    75,  /* 75-node fourth order prism (6 nodes associated with the
            vertices, 27 with the edges, 33 with the faces, and 9 with
            the volume) */
    64,  /* 64-node third order hexahedron (8 nodes associated with the
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
    1000,/* 1000-node ninth order hexahedron (8 nodes associated with
            the vertices, 96 with the edges, 384 with the faces and 512
            with the volume).*/
    -1,-1,-1,-1,-1,-1,-1, /* unsupported hexahedron types */
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
    -1,-1,-1,-1,-1,-1,-1, /* unsupported prism types */
    30,  /* 30-node third order pyramid (5 nodes associated with the
            vertices, 16 with the edges and 8 with the faces, and 1 with
            the volume). */
    55,  /* 55-node fourth order pyramid (5 nodes associated with the
            vertices, 24 with the edges and 21 with the faces, and 5
            with the volume). */
    91,  /* 91-node fifth order pyramid (5 nodes associated with the
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
    385  /* 385-node ninth order pyramid (5 nodes associated with the
            vertices, 64 with the edges and 176 with the faces, and 140
            with the volume). */
 };

int feMesh2DP1::readGmsh(std::string meshName, bool curved){
  std::filebuf fb;
  if(fb.open(meshName,std::ios::in)){
    std::istream input(&fb);

    std::string buffer, foo;
    double version;
    int binary, dsize;
    getline(input, buffer);
    input >> version >> binary >> dsize;
    if(version < 4){
      printf("In readGmsh : Error - Mesh version should be 4.1 or higher.\n");
      return 1;
    }
    if(binary){
      printf("In readGmsh : Error - Only reading ASCII files.\n");
      return 1;
    }
    
    int ph1, ph2, ph3, ph4, ph5, ph6; // Placeholders

    std::map<int, int> verticesMap;
    std::vector<int> numNodesInBlock;

    while(input >> buffer){

      if(buffer == "$PhysicalNames"){ // Read physical entities
        input >> _nPhysicalEntities;
        getline(input, buffer);
        for(int i = 0; i < _nPhysicalEntities; ++i){
          physicalEntity pE;
          // dimension(ASCII int) - physicalTag(ASCII int) - "name"(127 characters max)
          input >> pE.dim >> pE.tag >> pE.name;
          if(pE.name == ""){
            printf("In readGmsh : Error - Empty name for physical entity.\n");
            return 1;
          }
          pE.name.pop_back();
          pE.name.erase(pE.name.begin());
          pE.nElm = 0;
          pE.nNodePerElem = -1;
          pE.nEdgePerElem = -1;
          pE.cncID = "";
          _physicalEntities.push_back(pE);
        }
      }

      else if(buffer == "$Entities"){ // Read geometric entities
        // numPoints(size_t) numCurves(size_t) numSurfaces(size_t) numVolumes(size_t)
        input >> _nPoints >> _nCurves >> _nSurfaces >> _nVolumes;

        // Pour le moment : pas de cncGeo sur les Points
        numNodesInBlock.resize(_nCurves + _nSurfaces + _nVolumes);

        int tagCount = 1;
        // Points
        for(int i = 0; i < _nPoints; ++i){
          entity e;
          e.dim = 0;
          e.tag = tagCount++;
          //  pointTag(int) X(double) Y(double) Z(double) 
          //  numPhysicalTags(size_t) physicalTag(int)
          input >> e.tagGmsh >> ph1 >> ph2 >> ph3 >> e.numPhysicalTags;
          for(int j = 0; j < e.numPhysicalTags; ++j){
            input >> ph1;
            e.physicalTags.push_back(ph1);
            _physicalEntities[ph1-1].listEntities.push_back(e.tag);
          }
          _entities.insert({{e.dim,e.tagGmsh},e});
        }
        // Curves
        for(int i = 0; i < _nCurves; ++i){
          entity e;
          e.dim = 1;
          e.tag = tagCount++;
          // curveTag(int) minX(double) minY(double) minZ(double)
          //  maxX(double) maxY(double) maxZ(double)
          //  numPhysicalTags(size_t) physicalTag(int) ...
          //  numBoundingPoints(size_t) pointTag(int) ... >>>>>>>>>>>>> Ignorés pour le moment
          input >> e.tagGmsh >> ph1 >> ph2 >> ph3 >> ph4 >> ph5 >> ph6 >> e.numPhysicalTags;
          for(int j = 0; j < e.numPhysicalTags; ++j){
            input >> ph1;
            e.physicalTags.push_back(ph1);
            _physicalEntities[ph1-1].listEntities.push_back(e.tag);
          }
          _entities.insert({{e.dim,e.tagGmsh},e});
          getline(input, buffer);
        }
        // Surfaces
        for(int i = 0; i < _nSurfaces; ++i){
          entity e;
          e.dim = 2;
          e.tag = tagCount++;
          // surfaceTag(int) minX(double) minY(double) minZ(double)
          //   maxX(double) maxY(double) maxZ(double)
          //   numPhysicalTags(size_t) physicalTag(int) ...
          //   numBoundingCurves(size_t) curveTag(int) ... >>>>>>>>>>>>> Ignorés pour le moment
          input >> e.tagGmsh >> ph1 >> ph2 >> ph3 >> ph4 >> ph5 >> ph6 >> e.numPhysicalTags;
          for(int j = 0; j < e.numPhysicalTags; ++j){
            input >> ph1;
            e.physicalTags.push_back(ph1);
            _physicalEntities[ph1-1].listEntities.push_back(e.tag);
          }
          _entities.insert({{e.dim,e.tagGmsh},e});
          getline(input, buffer);
        }
        // Volumes
        for(int i = 0; i < _nVolumes; ++i){
          entity e;
          e.dim = 3;
          e.tag = tagCount++;
          // volumeTag(int) minX(double) minY(double) minZ(double)
          //   maxX(double) maxY(double) maxZ(double)
          //   numPhysicalTags(size_t) physicalTag(int) ...
          //   numBoundingSurfaces(size_t) surfaceTag(int) ... >>>>>>>>>>>>> Ignorés pour le moment
          input >> e.tagGmsh >> ph1 >> ph2 >> ph3 >> ph4 >> ph5 >> ph6 >> e.numPhysicalTags;
          for(int j = 0; j < e.numPhysicalTags; ++j){
            input >> ph1;
            e.physicalTags.push_back(ph1);
            _physicalEntities[ph1-1].listEntities.push_back(e.tag);
          }
          _entities.insert({{e.dim,e.tagGmsh},e});
          getline(input, buffer);
        }
      }

      else if(buffer == "$Nodes"){ // Read nodes
        int numEntityBlocks, countVertex = 0;
        double x, y, z;
        // numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
        input >> numEntityBlocks >> _nNod; // != nVertices, mais alloue _vertices : ambigu (:
        getline(input, buffer);
        _vertices.resize(_nNod);
        // A counter for the dim > 0 entities
        int entityCount = 0;
        for(int i = 0; i < numEntityBlocks; ++i){
          int entityDim, numNodes;
          // entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
          input >> entityDim >> ph2 >> ph3 >> numNodes;
          if(entityDim > 0) numNodesInBlock[entityCount++] = numNodes;
          std::vector<int> serialNumber(numNodes);
          for(int iNode = 0; iNode < numNodes; ++iNode){
            input >> ph1;
            serialNumber[iNode] = ph1;
            // std::cout<<serialNumber<<std::endl;
            verticesMap[serialNumber[iNode]] = countVertex++;
          }
          countVertex -= numNodes;
          for(int iNode = 0; iNode < numNodes; ++iNode){
            input >> x >> y >> z;
            _vertices[countVertex++] = Vertex(x,y,z,serialNumber[iNode]);
          }
        }
        if(verticesMap.size() != _nNod){
          printf("In readGmsh : Error - vertices indices are not unique.\n");
          return 1;
        }
      }

      else if(buffer == "$Elements"){ // Read elements
        int numEntityBlocks, serialNumber, maxElementTag;
        // numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
        input >> numEntityBlocks >> _nElm >> ph1 >> maxElementTag;
        _nTotalElm = _nElm;

        // Check if we can use the Gmsh element tag in our element connectivity :
        if(maxElementTag != _nElm){
          printf("In readGmsh : Error - maxElementTag differs from the number of elem.\n");
          printf("Possible gap in element numbering.\n");
          return 1;
        }

        // Inspired by MFEM :
        // std::vector<feElement*> elements0D, elements1D, elements2D, elements3D;
        // elements0D.reserve(_nElm);
        // elements1D.reserve(_nElm);
        // elements2D.reserve(_nElm);
        // elements3D.reserve(_nElm);
        getline(input, buffer);
        int nEdges = 1; // Starting edge numbering at 0 to distinguish + or - in the connectivity

        for(int i = 0; i < numEntityBlocks; ++i){

          int numElementsInBlock, elemType, entityTag, entityDim;
          // entityDim(int) entityTag(int) elementType(int) numElementsInBlock(size_t)
          input >> entityDim >> entityTag >> elemType >> numElementsInBlock;

          std::pair<int,int> p = {entityDim, entityTag};
          _entities[p].nElm = numElementsInBlock;
          // Element connectivity : use Gmsh element tag ? Not sure if sequential, otherwise create a map
          // Also, it starts at 1 and not 0.
          _entities[p].connecElem.resize(numElementsInBlock);

          // Determine geometric feSpace and allocate nodes connectivity based on elemType
          switch(elemType){
            case  1: //  2-node line
            case  8: if(curved){  _entities[p].cncID = "LineP2"; }//  3-node line (2nd order)
            case 26: if(curved){  _entities[p].cncID = "LineP3"; }//  4-node line (3rd order)
            case 27: if(curved){ /* Not supported */                          }//  5-node line (4th order)
            case 28: if(curved){ /* Not supported */                          }//  6-node line (5th order)
            case 62: if(curved){ /* Not supported */                          }//  7-node line (6th order)
            case 63: if(curved){ /* Not supported */                          }//  8-node line (7th order)
            case 64: if(curved){ /* Not supported */                          }//  9-node line (8th order)
            case 65: if(curved){ /* Not supported */                          }// 10-node line (9th order)
            case 66: if(curved){ /* Not supported */                          }// 11-node line (10th order)
            {
              // Default is linear interpolation for the geometry
              if(!curved){
                _entities[p].cncID = "LineP1";
              }
              _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType-1] : 2;
              _entities[p].nEdgePerElem = 0;
              _entities[p].connecNodes.resize(numElementsInBlock*_entities[p].nNodePerElem);
              break;
            }
            case  2: //  3-node triangle
            case  9: // if(curved){ geoSpace = new feSpaceTriP2("TriP2"); connecNodes.resize(numElementsInBlock*6); break; }//  6-node triangle (2nd order)
            case 21: ; // 10-node triangle (3rd order)
            case 23: ; // 15-node triangle (4th order)
            case 25: ; // 21-node triangle (5th order)
            case 42: ; // 28-node triangle (6th order)
            case 43: ; // 36-node triangle (7th order)
            case 44: ; // 45-node triangle (8th order)
            case 45: ; // 55-node triangle (9th order)
            case 46: ; // 66-node triangle (10th order)
            {
              if(!curved){
                _entities[p].cncID = "TriP1";
              }
              _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType-1] : 3;
              _entities[p].nEdgePerElem = 3;
              _entities[p].connecNodes.resize(numElementsInBlock*_entities[p].nNodePerElem);
              _entities[p].connecEdges.resize(numElementsInBlock*_entities[p].nEdgePerElem);
              break;
            }
            case  3: ; //   4-node quadrangle
            case 10: ; //   9-node quadrangle (2nd order)
            case 36: ; //  16-node quadrangle (3rd order)
            case 37: ; //  25-node quadrangle (4th order)
            case 38: ; //  36-node quadrangle (5th order)
            case 47: ; //  49-node quadrangle (6th order)
            case 48: ; //  64-node quadrangle (7th order)
            case 49: ; //  81-node quadrangle (8th order)
            case 50: ; // 100-node quadrangle (9th order)
            case 51: ; // 121-node quadrangle (10th order)
            {
              _entities[p].nNodePerElem = (curved) ? nodes_of_gmsh_element[elemType-1] : 4;
              _entities[p].nEdgePerElem = 4;
              _entities[p].connecNodes.resize(numElementsInBlock*_entities[p].nNodePerElem);
              _entities[p].connecEdges.resize(numElementsInBlock*_entities[p].nEdgePerElem);
              printf("In readMesh : Error - Interpolant pas pris en charge pour la géométrie de l'entité %d (quad).\n", entityTag);
              return 1;
              break;
            }
            case  4: ; //   4-node tetrahedron
            case 11: ; //  10-node tetrahedron (2nd order)
            case 29: ; //  20-node tetrahedron (3rd order)
            case 30: ; //  35-node tetrahedron (4th order)
            case 31: ; //  56-node tetrahedron (5th order)
            case 71: ; //  84-node tetrahedron (6th order)
            case 72: ; // 120-node tetrahedron (7th order)
            case 73: ; // 165-node tetrahedron (8th order)
            case 74: ; // 220-node tetrahedron (9th order)
            case 75: ; // 286-node tetrahedron (10th order)
            {
              printf("In readMesh : Error - Interpolant pas pris en charge pour la géométrie de l'entité %d (tet).\n", entityTag);
              return 1;
              break;
            }
            case  5: ; //    8-node hexahedron
            case 12: ; //   27-node hexahedron (2nd order)
            case 92: ; //   64-node hexahedron (3rd order)
            case 93: ; //  125-node hexahedron (4th order)
            case 94: ; //  216-node hexahedron (5th order)
            case 95: ; //  343-node hexahedron (6th order)
            case 96: ; //  512-node hexahedron (7th order)
            case 97: ; //  729-node hexahedron (8th order)
            case 98: ; // 1000-node hexahedron (9th order)
            {
              printf("In readMesh : Error - Interpolant pas pris en charge pour la géométrie de l'entité %d (hex).\n", entityTag);
              return 1;
              break;
            }
            case 15: // 1-node point
            {
              _entities[p].cncID = "Point0D";
              _entities[p].nNodePerElem = 1;
              _entities[p].nEdgePerElem = 0;
              _entities[p].connecNodes.resize(numElementsInBlock);
              break;
            }
              default: printf("Unsupported Gmsh element type."); 
              return 1;
              break;
          } // switch(elemType)

          // std::cout<<"Interpolant assigné : "<<_entities[p].cncID<<std::endl;

          std::map<int, int>::const_iterator it;

          for(int iElm = 0; iElm < numElementsInBlock; ++iElm){
            int nElemNodes = nodes_of_gmsh_element[elemType-1];
            std::vector<int> elemNodesGmsh(nElemNodes);
            std::vector<int> elemNodes(nElemNodes);
            // elementTag(size_t) nodeTag(size_t) ...
            input >> serialNumber;
            _entities[p].connecElem[iElm] = serialNumber-1; // Gmsh elem tag starts at 1

            for(int j = 0; j < nElemNodes; ++j){
              input >> ph1;
              // Verification
              it = verticesMap.find(ph1);
              if(it == verticesMap.end()){
                printf("In readGmsh : Error - element node does not exist /-:\n");
                return 1;
              }
              elemNodesGmsh[j] = ph1;    // The Gmsh number stored in ph1 is used to create the edges
              elemNodes[j] = it->second; // The node number (0...nNode) is used to create the connectivity
            }

            int el_order = 11;
            switch (elemType){
              case  1: //  2-node line
              case  8: // TODO : Decider de la numerotation pour les P2+ : {0,2,1} ou {0,1,2} //  3-node line (2nd order)
              case 26: //  4-node line (3rd order)
              case 27: //  5-node line (4th order)
              case 28: //  6-node line (5th order)
              case 62: //  7-node line (6th order)
              case 63: //  8-node line (7th order)
              case 64: //  9-node line (8th order)
              case 65: // 10-node line (9th order)
              case 66: // 11-node line (10th order)
              {
                // !!! ConnecNodes valide juste pour des P1 géométriques pour le moment (décider de la numérotation)
                _entities[p].connecNodes[2*iElm  ] = elemNodes[0];
                _entities[p].connecNodes[2*iElm+1] = elemNodes[1];
                 // elements_1D.push_back(
                 //    new Segment(&vert_indices[0], phys_domain));
                 // if (type_of_element != 1)
                 // {
                 //    Array<int> * hov = new Array<int>;
                 //    hov->Append(&vert_indices[0], n_elem_nodes);
                 //    ho_verts_1D.push_back(hov);
                 //    el_order = n_elem_nodes - 1;
                 //    ho_el_order_1D.push_back(el_order);
                 // }
                 break;
              }
              case  2: el_order--; //  3-node triangle
              case  9: el_order--; //  6-node triangle (2nd order)
              case 21: el_order--; // 10-node triangle (3rd order)
              case 23: el_order--; // 15-node triangle (4th order)
              case 25: el_order--; // 21-node triangle (5th order)
              case 42: el_order--; // 28-node triangle (6th order)
              case 43: el_order--; // 36-node triangle (7th order)
              case 44: el_order--; // 45-node triangle (8th order)
              case 45: el_order--; // 55-node triangle (9th order)
              case 46: el_order--; // 66-node triangle (10th order)
              {
                _entities[p].connecNodes[3*iElm  ] = elemNodes[0];
                _entities[p].connecNodes[3*iElm+1] = elemNodes[1];
                _entities[p].connecNodes[3*iElm+2] = elemNodes[2];
                // Construct the triangle edges :
                Vertex *v0, *v1;
                for(int k = 0; k < 3; ++k){
                  if(k == 2){
                    v0 = &_vertices[verticesMap[elemNodesGmsh[2]]];
                    v1 = &_vertices[verticesMap[elemNodesGmsh[0]]];
                  } else{
                    v0 = &_vertices[verticesMap[elemNodesGmsh[k]]];
                    v1 = &_vertices[verticesMap[elemNodesGmsh[k+1]]];
                  }
                  std::pair<std::set<Edge, EdgeLessThan>::iterator,bool> ret;
                  Edge e(v0,v1, nEdges);
                  ret = _edges.insert(e);
                  if(ret.second){
                    // Edge was added to the set : nEdges is added to connecEdges
                    _entities[p].connecEdges[3*iElm+k] = nEdges++;
                  } else{
                    // Edge is already in the set : the negative is added to connecEdges
                    // Assumes an edge is shared by only two triangles in 2D
                    // More tests required in 3D, where an edge can be shared by N tets
                    _entities[p].connecEdges[3*iElm+k] = -ret.first->getTag();
                  }
                }
                // elements_2D.push_back(
                //    new Triangle(&vert_indices[0], phys_domain));
                // if (el_order > 1)
                // {
                //    Array<int> * hov = new Array<int>;
                //    hov->Append(&vert_indices[0], n_elem_nodes);
                //    ho_verts_2D.push_back(hov);
                //    ho_el_order_2D.push_back(el_order);
                // }
                break;
              }
              case  3: el_order--; //   4-node quadrangle
              case 10: el_order--; //   9-node quadrangle (2nd order)
              case 36: el_order--; //  16-node quadrangle (3rd order)
              case 37: el_order--; //  25-node quadrangle (4th order)
              case 38: el_order--; //  36-node quadrangle (5th order)
              case 47: el_order--; //  49-node quadrangle (6th order)
              case 48: el_order--; //  64-node quadrangle (7th order)
              case 49: el_order--; //  81-node quadrangle (8th order)
              case 50: el_order--; // 100-node quadrangle (9th order)
              case 51: el_order--; // 121-node quadrangle (10th order)
              {
                printf("In readMesh : Error - Interpolant pas pris en charge (quad).\n");
                // elements_2D.push_back(
                //    new Quadrilateral(&vert_indices[0], phys_domain));
                // if (el_order > 1)
                // {
                //    Array<int> * hov = new Array<int>;
                //    hov->Append(&vert_indices[0], n_elem_nodes);
                //    ho_verts_2D.push_back(hov);
                //    ho_el_order_2D.push_back(el_order);
                // }
                break;
              }
              case  4: el_order--; //   4-node tetrahedron
              case 11: el_order--; //  10-node tetrahedron (2nd order)
              case 29: el_order--; //  20-node tetrahedron (3rd order)
              case 30: el_order--; //  35-node tetrahedron (4th order)
              case 31: el_order--; //  56-node tetrahedron (5th order)
              case 71: el_order--; //  84-node tetrahedron (6th order)
              case 72: el_order--; // 120-node tetrahedron (7th order)
              case 73: el_order--; // 165-node tetrahedron (8th order)
              case 74: el_order--; // 220-node tetrahedron (9th order)
              case 75: el_order--; // 286-node tetrahedron (10th order)
              {
                printf("In readMesh : Error - Interpolant pas pris en charge (tet).\n");
                return 1;
                // elements_3D.push_back(
                //    new Tetrahedron(&vert_indices[0], phys_domain));
                // if (el_order > 1)
                // {
                //    Array<int> * hov = new Array<int>;
                //    hov->Append(&vert_indices[0], n_elem_nodes);
                //    ho_verts_3D.push_back(hov);
                //    ho_el_order_3D.push_back(el_order);
                // }
                break;
              }
              case  5: el_order--; //    8-node hexahedron
              case 12: el_order--; //   27-node hexahedron (2nd order)
              case 92: el_order--; //   64-node hexahedron (3rd order)
              case 93: el_order--; //  125-node hexahedron (4th order)
              case 94: el_order--; //  216-node hexahedron (5th order)
              case 95: el_order--; //  343-node hexahedron (6th order)
              case 96: el_order--; //  512-node hexahedron (7th order)
              case 97: el_order--; //  729-node hexahedron (8th order)
              case 98: el_order--; // 1000-node hexahedron (9th order)
              {
                printf("In readMesh : Error - Interpolant pas pris en charge (hex).\n");
                return 1;
                // el_order--;
                // elements_3D.push_back(
                //    new Hexahedron(&vert_indices[0], phys_domain));
                // if (el_order > 1)
                // {
                //    Array<int> * hov = new Array<int>;
                //    hov->Append(&vert_indices[0], n_elem_nodes);
                //    ho_verts_3D.push_back(hov);
                //    ho_el_order_3D.push_back(el_order);
                // }
                break;
              }
              case 15: // 1-node point
              {
                _entities[p].connecNodes[iElm] = elemNodes[0];
                // elements_0D.push_back(
                //    new Point(&vert_indices[0], phys_domain));
                break;
              }
              default: // any other element
                printf("Unsupported Gmsh element type.");
                return 1;
                break;
            } // switch (type_of_element)
          } //for numElemInEntity

          // std::cout<<"Num nodes in entity "<<i<<" = "<<numNodesInBlock[i]<<std::endl;
        } // for numEntities
      } // if buffer = "Elements"
    } // while input
    fb.close();
  } // if fb.open

  // Geometric connectivities are defined on the physical entities (named domains) :
  // we need to transfer the information to the physical entities.
  // Sequential tag for all entities, different from the Gmsh entity tag (reset for each dimension)
  for(auto const& x : _entities){
    entity e = x.second;
    // std::cout<<"Entity "<<e.tag<<std::endl;
    for(auto& pE : _physicalEntities){
      bool error = false;
      for(auto ent : pE.listEntities){
        if(ent == e.tag){
          // std::cout<<"Physical "<<pE.name<<" - entity "<<ent<<" : match"<<std::endl;
          pE.nElm += e.nElm;
          // Assign and check nNodePerElem
          if(pE.nNodePerElem == -1) pE.nNodePerElem = e.nNodePerElem;
          else{
            if(pE.nNodePerElem != e.nNodePerElem) error = true;
          }
          // Assign and check nEdgePerElem
          if(pE.nEdgePerElem == -1) pE.nEdgePerElem = e.nEdgePerElem;
          else{
            if(pE.nEdgePerElem != e.nEdgePerElem) error = true;
          }
          // Assign and check geoSpace and cncID
          if(pE.cncID == "") pE.cncID = e.cncID;
          else{
            if(pE.cncID != e.cncID) error = true;
          }
        }
      }
      if(error){
        printf("In readMesh : Error - Multiple geometric connectivities on physical entity ""%s"".\n", pE.name.c_str());
        return 1;
      }
    }
  }

  // Assign pointer to geometric space
  for(auto& pE : _physicalEntities){
    if     (pE.cncID == "Point0D"){ pE.geoSpace = new feSpace1DP0("xyz");  }
    else if(pE.cncID == "LineP1" ){ pE.geoSpace = new feSpace1DP1("xyz");  }
    else if(pE.cncID == "LineP2" ){ pE.geoSpace = new feSpace1DP2("xyz");  }
    else if(pE.cncID == "LineP3" ){ pE.geoSpace = new feSpace1DP3("xyz");  }
    else if(pE.cncID == "TriP1"  ){ pE.geoSpace = new feSpaceTriP1("xyz"); }
     else{
        printf("In readMesh : Error - Unknown geometric connectivity on domain \"%s\".\n", pE.name.c_str());
        return 1;
    }
  }

  // Connectivities
  int nCncGeo = 0;
  int maxDim = 0;
  for(auto& pE : _physicalEntities){
    maxDim = fmax(maxDim, pE.dim);
    pE.connecElem.resize(pE.nElm);
    pE.connecNodes.resize(pE.nElm*pE.nNodePerElem);
    pE.connecEdges.resize(pE.nElm*pE.nEdgePerElem);
    int countElm = 0;
    for(auto const& x : _entities){
      entity e = x.second;
      for(auto ent : pE.listEntities){
        if(ent == e.tag){
          for(int iElm = 0; iElm < e.nElm; ++iElm){
            // Connec elem
            pE.connecElem[countElm] = e.connecElem[iElm];
            // Connec nodes
            for(int j = 0; j < e.nNodePerElem; ++j){
              pE.connecNodes[e.nNodePerElem*countElm+j] = e.connecNodes[e.nNodePerElem*iElm+j];
            }
            // Connec edges
            for(int j = 0; j < e.nEdgePerElem; ++j){
              pE.connecEdges[e.nEdgePerElem*countElm+j] = e.connecEdges[e.nEdgePerElem*iElm+j];
            }
            ++countElm;
          }
        }
      }
    }
    // Finally, create geometric connectivity
    feCncGeo *cnc = new feCncGeo(pE.nNodePerElem, pE.nElm, pE.nEdgePerElem, pE.name, pE.cncID, pE.geoSpace, pE.connecNodes, pE.connecElem, pE.connecEdges);
    // printf("Created cnc on domain %s with nNod = %d - nElm = %d - connecNodes.size() = %d - connecElem.size() = %d\n",
    //   pE.name.c_str(),
    //   cnc->getNbNodePerElem(),
    //   cnc->getNbElm(),
    //   pE.connecNodes.size(),
    //   pE.connecElem.size());
    _cncGeo.push_back(cnc);
    // std::cout<<pE.name<<std::endl;
    _cncGeoMap[pE.name] = nCncGeo;
    cnc->getFeSpace()->setCncGeoTag(nCncGeo++);
  }

  _nEdg = _edges.size();
  
  _dim = maxDim;
  _nInteriorElm = 0;
  _nBoundaryElm = 0;
  for(auto& pE : _physicalEntities){
    if(pE.dim == _dim)   _nInteriorElm += pE.nElm;
    if(pE.dim == _dim-1) _nBoundaryElm += pE.nElm;
  }

  // TODO : Check for overlapping physical entities

  return 0;
}

