#include "Mesh.h"

Mesh::Mesh(const InputReader& inputReader)
    : inputReader_(inputReader) {
    // You may want to initialize other member variables here
}

void Mesh::getHexahedralElements() {
    gmsh::initialize();
    gmsh::open(inputReader_.getMeshFileName());

    // get element type for Hexahedron of order 1
    int order=1;
    if (numPointsInHex>8){order=2;};
    int elementType = gmsh::model::mesh::getElementType("Hexahedron", order);

    // get tags for all hexahedrons of order 1
    std::vector<double>  parametricCoord;
    gmsh::model::mesh::getElementsByType(elementType, elementTags, elementNodeTags);  // get all Element of given type 
    gmsh::model::mesh::getNodesByElementType(elementType,nodeTagselem,coord,parametricCoord); // get all nodes from Element type (related to size of system and dofs)
    std::vector<double> nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, coord, nodeParams, -1, -1); // get all nodes in the mesh

    // Check if no tags were found and print a warning
    if (elementTags.empty()) {
        std::cout << "!!!WARNING!!! No hexahedral order 1 element tags found." << std::endl;
    } else {
        // Output the number of nodes and elements
        numelem= elementTags.size();
        numnodes= nodeTagselem.size();
        numallnodes =nodeTags.size();
        std::cout << "Number of Hexahedral Order 1 Elements: " << numelem << std::endl;
        std::cout << "Total Number of Nodes: " << numallnodes << std::endl;
        std::cout << "Number of Nodes: " << numnodes << std::endl;
    }

    gmsh::finalize();
}

