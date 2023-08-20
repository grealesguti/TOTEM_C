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
    std::vector<double> nodeParams,coord1,coord2;
    gmsh::model::mesh::getElementsByType(elementType, elementTags, elementNodeTags);  // get all Element of given type 
    gmsh::model::mesh::getNodesByElementType(elementType,nodeTagselem,coord1,parametricCoord); // get all nodes from Element type (related to size of system and dofs)
    gmsh::model::mesh::getNodes(nodeTags, coord2, nodeParams, -1, -1); // get all nodes in the mesh
    
    // THE nodeTags are not in order!!!
        std::ofstream outputFile("nodeTags.txt"); // Open the file for writing
        if (outputFile.is_open()) { // Check if the file was successfully opened
            for (int tag : nodeTags) {
            outputFile << tag << std::endl; // Write each tag to a new line in the file
            }
            outputFile.close(); // Close the file
            std::cout << "Contents of nodeTags have been written to output.txt" << std::endl;
        } else {
            std::cerr << "Unable to open the file for writing" << std::endl;
        }
        std::ofstream outputFile2("coord.txt"); // Open the file for writing
        if (outputFile2.is_open()) { // Check if the file was successfully opened
            for (int i = 0; i <  nodeTags.size(); ++i) {
            outputFile2 << coord2[i*3]<<" "<<coord2[i*3+1]<<" "<<coord2[i*3+2]<<" " << std::endl; // Write each tag to a new line in the file
            }
            outputFile2.close(); // Close the file
            std::cout << "Contents of nodeTags have been written to output.txt" << std::endl;
        } else {
            std::cerr << "Unable to open the file for writing" << std::endl;
        }
    numallnodes =nodeTags.size();
    int maxNodeTag = *std::max_element(nodeTags.begin(), nodeTags.end());
    coord.assign(maxNodeTag * 3, 0.0);
        std::ofstream outputFile4("debug.txt"); // Open the file for writing
    for (int i = 0; i <  numallnodes; ++i) {
        coord[(nodeTags[i]-1)*3]=coord2[i*3];
        coord[(nodeTags[i]-1)*3+1]=coord2[i*3+1];
        coord[(nodeTags[i]-1)*3+2]=coord2[i*3+2];
        outputFile4 <<i<<" "<<nodeTags[i]<<" "<< (nodeTags[i]-1)*3<<" "<<(nodeTags[i]-1)*3+1<<" "<<(nodeTags[i]-1)*3+2<<" " <<coord2[(nodeTags[i]-1)*3]<<" "<<coord2[(nodeTags[i]-1)*3+1]<<" "<<coord2[(nodeTags[i]-1)*3+2]<< std::endl; // Write each tag to a new line in the file
    }     
            outputFile4.close(); // Close the file

    std::ofstream outputFile3("coord2.txt"); // Open the file for writing
        if (outputFile3.is_open()) { // Check if the file was successfully opened
            for (int i = 0; i <  nodeTags.size(); ++i) {
            outputFile3 << coord[i*3]<<" "<<coord[i*3+1]<<" "<<coord[i*3+2]<<" " << std::endl; // Write each tag to a new line in the file
            }
            outputFile3.close(); // Close the file
            std::cout << "Contents of nodeTags have been written to output.txt" << std::endl;
        } else {
            std::cerr << "Unable to open the file for writing" << std::endl;
        }
    // Check if no tags were found and print a warning
    if (elementTags.empty()) {
        std::cout << "!!!WARNING!!! No hexahedral order 1 element tags found." << std::endl;
    } else {
        // Output the number of nodes and elements
        numelem= elementTags.size();
        numnodes= nodeTagselem.size();

        int numelementNodeTags =elementNodeTags.size();

        std::cout << "Number of Hexahedral Order 1 Elements: " << numelem << std::endl;
        std::cout << "Size of elementNodeTags: " << numelementNodeTags << std::endl;
        std::cout << "for elem "<<elementTags[0] <<", nodes ";
        for (int i = 0; i < 8; ++i) {
            std::cout << elementNodeTags[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "with cords "<< std::endl;
        for (int i = 0; i < 8; ++i) {
           std::cout << coord[(elementNodeTags[i]-1)*3] << " "<<coord[(elementNodeTags[i]-1)*3+1] << " "<<coord[(elementNodeTags[i]-1)*3+2]<< std::endl;
        }
        std::cout << "Total Number of Nodes: " << numallnodes << std::endl;
        std::cout << "Number of Nodes: " << numnodes << std::endl;
    }

    gmsh::finalize();
}

