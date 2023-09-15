#include "Mesh.h"
using namespace arma;

Mesh::Mesh(const InputReader& inputReader)
    : inputReader_(inputReader) {
    // You may want to initialize other member variables here
    gmsh::initialize();
    gmsh::open(inputReader_.getMeshFileName());
    
    // set all nodes coordinates and order according to the element tag
    std::vector<double>  parametricCoord;
    std::vector<double> nodeParams,coord2;
    gmsh::model::mesh::getNodes(nodeTags, coord2, nodeParams, -1, -1); // get all nodes in the mesh
    
    numallnodes =nodeTags.size();
    int maxNodeTag = *std::max_element(nodeTags.begin(), nodeTags.end());
    coord.assign(maxNodeTag * 3, 0.0);
    for (int i = 0; i <  numallnodes; ++i) {
        coord[(nodeTags[i]-1)*3]=coord2[i*3];
        coord[(nodeTags[i]-1)*3+1]=coord2[i*3+1];
        coord[(nodeTags[i]-1)*3+2]=coord2[i*3+2];
    }     

    InitMeshEntityElements();
    applyElementsMaterials();
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
void Mesh::getHexahedralElements() {

    // get element type for Hexahedron of order 1
    int order=1;
    if (numPointsInHex>8){order=2;};
    int elementType = gmsh::model::mesh::getElementType("Hexahedron", order);

    // get tags for all hexahedrons of order 1
    std::vector<double>  parametricCoord;
    std::vector<double> nodeParams,coord1,coord2;
    gmsh::model::mesh::getElementsByType(elementType, elementTags, elementNodeTags);  // get all Element of given type 
    gmsh::model::mesh::getNodesByElementType(elementType,nodeTagselem,coord1,parametricCoord); // get all nodes from Element type (related to size of system and dofs)
    
        // THE nodeTags are not in order!!!
        /* Testing nodetags and initial coords written to file
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
        }*/

    /*  Testing coords written to file after order modification
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
        */
    // Check if no tags were found and print a warning
    if (elementTags.empty()) {
        std::cout << "!!!WARNING!!! No hexahedral order 1 element tags found." << std::endl;
    } else {
        // Output the number of nodes and elements
        numelem= elementTags.size();
        element_materials.assign(numelem, 0);
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
        freedofs_.resize(2 * numallnodes, 0);
        for (int i = 0; i < elementNodeTags.size(); ++i) {
                freedofs_[elementNodeTags[i]*2]=1;
                freedofs_[elementNodeTags[i]*2+1]=1;
        }

    }


}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<long long unsigned int> Mesh::getNodesForPhysicalGroup(const std::string& desiredGroupName) {
    int matchedDimension;

    std::vector<std::pair<int, int>> dimTags;

    // Get all the physical groups
    gmsh::model::getPhysicalGroups(dimTags, -1);

    int matchingTag = -1;
    for (const auto& pair : dimTags) {
        int entityDimension = pair.first;
        int entityTag = pair.second;

        // Get the name of the physical group
        std::string groupName;
        gmsh::model::getPhysicalName(entityDimension, entityTag, groupName);

        // Compare the group name with the desired name
        if (groupName == desiredGroupName) {
            matchingTag = entityTag;
            matchedDimension = entityDimension; // Store the matched dimension
            break; // Exit the loop when a match is found
        }
    }
    std::vector<long long unsigned int> nodeTagsgroup; // Declare a non-reference vector

    if (matchingTag != -1) {
        // Get the node tags and coordinates associated with the physical entity tag
        std::vector<double> coord;
        gmsh::model::mesh::getNodesForPhysicalGroup(matchedDimension, matchingTag, nodeTagsgroup, coord);
        if (nodeTagsgroup.empty()) {
            std::cout << "Warning!!! No nodes found for the physical group with name: " << desiredGroupName << std::endl;
        }else {std::cout << "Found nodes for group: "<< desiredGroupName<< std::endl;}

    } else {
        std::cout << "Warning!!! No matching physical group found with name: " << desiredGroupName << std::endl;
    }
        return nodeTagsgroup; // Return an empty vector to indicate no match
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::size_t> Mesh::getElementsForPhysicalGroup(const std::string& desiredGroupName) {
    // Get the desired entity name from InputReader
    std::cout << "Looking for element types in  "<< desiredGroupName<< std::endl;
    std::vector<std::size_t> elementsingroup;
    int matchedDimension = -1;
    std::vector<std::pair<int, int>> dimTags;

    // Get all the physical groups
    gmsh::model::getPhysicalGroups(dimTags, -1);

    int matchingTag = -1;
    for (const auto& pair : dimTags) {
        int entityDimension = pair.first;
        int entityTag = pair.second;

        // Get the name of the physical group
        std::string groupName;
        gmsh::model::getPhysicalName(entityDimension, entityTag, groupName);

        // Compare the group name with the desired name
        if (groupName == desiredGroupName) {
            matchingTag = entityTag;
            matchedDimension = entityDimension; // Store the matched dimension
            break; // Exit the loop when a match is found
        }
    }

    if (matchingTag == -1) {
        std::cout << "Physical group with name '" << desiredGroupName << "' not found." << std::endl;
    }else{
        std::cout << "Physical group with name '" << desiredGroupName << "' found with tag "<< matchingTag << " and dimension "<< matchedDimension << std::endl;
    }
    std::vector<int> entityTags; // Declare a non-reference vector
    gmsh::model::getEntitiesForPhysicalGroup(matchedDimension, matchingTag, entityTags);

    if (entityTags.empty()) {
        std::cout << "No entities found in physical group with name '" << desiredGroupName << "'" << std::endl;
    }
    for (int entityTag : entityTags) {
        std::vector<std::vector<std::size_t>> entity_NodeTags, elementTagsInEntity;
        std::vector<int> elemTypes;
        // Modify the function call to read the dimension as well
        gmsh::model::mesh::getElements(elemTypes, elementTagsInEntity, entity_NodeTags, matchedDimension, entityTag); 
        std::cout << "Get elements and nodes for entity " << entityTag << "." << std::endl;
        for (const auto& e : elementTagsInEntity) {
            for (const auto& element : e) {
                elementsingroup.push_back(element);
            }
        }
    }
    return elementsingroup;
}

void Mesh::finalizeGmsh() {
    // Finalize Gmsh and clean up
    gmsh::finalize();
}

void Mesh::reinitializeGmsh() {
    // Initialize Gmsh
    gmsh::initialize();
    gmsh::open(inputReader_.getMeshFileName());
}

void Mesh::setFixedof(int dof) {
            freedofs_[dof] = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
Mesh::getElementsAndNodeTagsForPhysicalGroup(const std::string& desiredGroupName) {
    // Get the desired entity name from InputReader
    std::cout << "Looking for element types in  "<< desiredGroupName<< std::endl;
    std::vector<std::size_t> elementsingroup;
    int matchedDimension = -1;
    std::vector<std::pair<int, int>> dimTags;

    // Get all the physical groups
    gmsh::model::getPhysicalGroups(dimTags, -1);

    int matchingTag = -1;
    for (const auto& pair : dimTags) {
        int entityDimension = pair.first;
        int entityTag = pair.second;

        // Get the name of the physical group
        std::string groupName;
        gmsh::model::getPhysicalName(entityDimension, entityTag, groupName);

        // Compare the group name with the desired name
        if (groupName == desiredGroupName) {
            matchingTag = entityTag;
            matchedDimension = entityDimension; // Store the matched dimension
            break; // Exit the loop when a match is found
        }
    }

    if (matchingTag == -1) {
        std::cout << "Physical group with name '" << desiredGroupName << "' not found." << std::endl;
    }else{
        std::cout << "Physical group with name '" << desiredGroupName << "' found with tag "<< matchingTag << " and dimension "<< matchedDimension << std::endl;
    }
    std::vector<int> entityTags; // Declare a non-reference vector
    gmsh::model::getEntitiesForPhysicalGroup(matchedDimension, matchingTag, entityTags);

    if (entityTags.empty()) {
        std::cout << "No entities found in physical group with name '" << desiredGroupName << "'" << std::endl;
    }
    std::vector<std::size_t> nodesingroup; // Vector to store all entity_NodeTags

    for (int entityTag : entityTags) {
        std::vector<std::vector<std::size_t>> entity_NodeTags, elementTagsInEntity;
        std::vector<int> elemTypes;

        // Modify the function call to read the dimension as well
        gmsh::model::mesh::getElements(elemTypes, elementTagsInEntity, entity_NodeTags, matchedDimension, entityTag);

        std::cout << "Get elements and nodes for entity " << entityTag << "." << std::endl;

        for (const auto& e : entity_NodeTags) {
            nodesingroup.insert(nodesingroup.end(), e.begin(), e.end()); // Flatten the node tags
        }

        for (const auto& e : elementTagsInEntity) {
            for (const auto& element : e) {
                elementsingroup.push_back(element);
            }
        }
    }

    // Return a pair of vectors: elementsingroup and nodesingroup
    return std::make_pair(elementsingroup, nodesingroup);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// Print element types in a physical group with a desired name
std::vector<int> Mesh::printElementTypesInPhysicalGroupByName(std::string desiredGroupName) {
    std::cout << "Looking for element types in  "<< desiredGroupName<< std::endl;

    int matchedDimension = -1;
    std::vector<std::pair<int, int>> dimTags;

    // Get all the physical groups
    gmsh::model::getPhysicalGroups(dimTags, -1);

    int matchingTag = -1;
    for (const auto& pair : dimTags) {
        int entityDimension = pair.first;
        int entityTag = pair.second;

        // Get the name of the physical group
        std::string groupName;
        gmsh::model::getPhysicalName(entityDimension, entityTag, groupName);

        // Compare the group name with the desired name
        if (groupName == desiredGroupName) {
            matchingTag = entityTag;
            matchedDimension = entityDimension; // Store the matched dimension
            break; // Exit the loop when a match is found
        }
    }

    if (matchingTag == -1) {
        std::cout << "Physical group with name '" << desiredGroupName << "' not found." << std::endl;
    }else{
        std::cout << "Physical group with name '" << desiredGroupName << "' found with tag "<< matchingTag << " and dimension "<< matchedDimension << std::endl;
    }
    std::vector<int> entityTags; // Declare a non-reference vector
    gmsh::model::getEntitiesForPhysicalGroup(matchedDimension, matchingTag, entityTags);

    if (entityTags.empty()) {
        std::cout << "No entities found in physical group with name '" << desiredGroupName << "'" << std::endl;
    }

    std::cout << "Element types in physical group with name '" << desiredGroupName << "':" << std::endl;

    for (int entityTag : entityTags) {
        std::vector<int> elementTypes;
        std::vector<int> entityTags;
        gmsh::model::mesh::getElementTypes(elementTypes, matchedDimension, entityTag);

        if (elementTypes.empty()) {
            std::cout << "No elements found for entity with tag " << entityTag << std::endl;
        } else {
            std::cout << "Element types for entity with tag " << entityTag << ":" << std::endl;
            for (int elementType : elementTypes) {
                // Print element type number
                std::cout << "Element type " << elementType << ": ";

                // Get element properties
                std::string elementName;
                int dim, order, numNodes, numPrimaryNodes;
                std::vector<double> localNodeCoord;

                gmsh::model::mesh::getElementProperties(elementType, elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);

                // Print element name next to the element type number
                std::cout << elementName << std::endl;

                // Print other properties if needed
                std::cout << "Dimension: " << dim << std::endl;
                std::cout << "Order: " << order << std::endl;
                std::cout << "Number of Nodes: " << numNodes << std::endl;
                std::cout << "Number of Primary Nodes: " << numPrimaryNodes << std::endl;
                
                // Print local node coordinates if needed
                std::cout << "Local Node Coordinates: ";
                for (double coord : localNodeCoord) {
                    std::cout << coord << " ";
                }
                std::cout << std::endl;
            }
        }
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
void Mesh::InitMeshEntityElements() {
    std::cout << "#MESH::InitMeshEntityElements" << std::endl;

    // Resize the element_materials vector to the size of numelem
    element_materials.resize(numelem);

    // Get the desired entity name from InputReader
    std::string desiredEntityName = inputReader_.getMeshEntityName();    
    std::cout << "Looking for element types in  "<< desiredEntityName<< std::endl;
    
    int matchedDimension = -1;
    std::vector<std::pair<int, int>> dimTags;

    // Get all the physical groups
    gmsh::model::getPhysicalGroups(dimTags, -1);

    int matchingTag = -1;
    for (const auto& pair : dimTags) {
        int entityDimension = pair.first;
        int entityTag = pair.second;

        // Get the name of the physical group
        std::string groupName;
        gmsh::model::getPhysicalName(entityDimension, entityTag, groupName);

        // Compare the group name with the desired name
        if (groupName == desiredEntityName) {
            matchingTag = entityTag;
            matchedDimension = entityDimension; // Store the matched dimension
            break; // Exit the loop when a match is found
        }
    }
    // Write to screen if the entity is found and its tags
    if (matchingTag == -1) {
        std::cout << "Physical group with name '" << desiredEntityName << "' not found." << std::endl;
    }else{
        std::cout << "Physical group with name '" << desiredEntityName << "' found with tag "<< matchingTag << " and dimension "<< matchedDimension << std::endl;
    }
    // Get the elements of the found entity
    std::vector<int> entityTags; // Declare a non-reference vector
    gmsh::model::getEntitiesForPhysicalGroup(matchedDimension, matchingTag, entityTags);

    if (entityTags.empty()) {
        std::cout << "No entities found in physical group with name '" << desiredEntityName << "'" << std::endl;
    }
    for (int entityTag : entityTags) {
        std::vector<std::vector<std::size_t>> entity_NodeTags, elementTagsInEntity;
        std::vector<int> elemTypes;
        // TODO: Modify the function call to read the dimension as well
        gmsh::model::mesh::getElements(elemTypes, elementTagsInEntity, entity_NodeTags, matchedDimension, entityTag); 
        std::cout << "Get elements and nodes for entity " << entityTag << "." << std::endl;
        for (const auto& e : elementTagsInEntity) {
            for (const auto& element : e) {
                elementTags.push_back(element);
            }
        }
        for (const auto& e : entity_NodeTags) {
            for (const auto& element : e) {
                elementNodeTags.push_back(element);
            }
        }
    }

    // Check if no tags were found and print a warning
    if (elementTags.empty()) {
        std::cout << "!!!WARNING!!! No elements stored in elementTags." << std::endl;
    } else {
        // Output the number of nodes and elements
        numelem= elementTags.size();

        int numelementNodeTags =elementNodeTags.size();

        std::cout << "Number of Elements: " << numelem << std::endl;
        std::cout << "Size of elementNodeTags: " << numelementNodeTags << std::endl;
        std::cout << "Total Number of Nodes: " << numallnodes << std::endl;
        std::cout << "Number of Nodes: " << numnodes << std::endl;
        freedofs_.resize(2 * numallnodes, 0);
        for (int i = 0; i < elementNodeTags.size(); ++i) {
                freedofs_[elementNodeTags[i]*2]=1;
                freedofs_[elementNodeTags[i]*2+1]=1;
        }

    }

}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
    // Method to get coordinates for a list of nodeTags
    mat Mesh::getCoordinates(const std::vector<int>& nodeTags) {
        // Determine the number of coordinates per node (assuming it's 3)
        const int numCoordinatesPerNode = 3;

        // Determine the number of nodes to retrieve
        const int numNodesToRetrieve = nodeTags.size();

        // Initialize the arma::mat to store the coordinates
        arma::mat coordinates(numCoordinatesPerNode, numNodesToRetrieve);

        // Populate the coordinates matrix with the requested nodes' coordinates
        for (int i = 0; i < numNodesToRetrieve; ++i) {
            int nodeTag = nodeTags[i]-1;
            if (nodeTag >= 0 && nodeTag < coord.size() / numCoordinatesPerNode) {
                // Calculate the starting index of the node's coordinates in the coord vector
                int startIndex = nodeTag * numCoordinatesPerNode;

                // Extract the coordinates for the node and store them in the arma::mat
                for (int j = 0; j < numCoordinatesPerNode; ++j) {
                    coordinates(j, i) = coord[startIndex + j];
                }
            } else {
                // Handle invalid nodeTag, for example, by setting the coordinates to NaN
                coordinates.col(i).fill(arma::datum::nan);
            }
        }

        return coordinates;
    }


    void Mesh::getElementInfo(int elementTag, std::vector<int> & nodeTags_el) {
        int dim, entityTag, elementType;
        std::vector<std::size_t> nodeTags_size_t;
        gmsh::model::mesh::getElement(elementTag, elementType, nodeTags_size_t, dim, entityTag);

        // Convert the std::vector<std::size_t> to std::vector<int>
        nodeTags_el.assign(nodeTags_size_t.begin(), nodeTags_size_t.end());
    }

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
void Mesh::applyElementsMaterials() {
    std::cout << "#MESH::ApplyElementMaterials" << std::endl;

    // Get all material properties from the InputReader
    std::map<std::string, std::map<std::string, double>> materialProperties = inputReader_.getAllMaterialProperties();
    std::size_t currentIndex = 1; // Start indexing from 1

    // Assign indices to unique material names
    for (const auto& entry : materialProperties) {
        const std::string& materialName = entry.first;

        // Check if the material name is not already assigned an index
        if (materialIndices.find(materialName) == materialIndices.end()) {
            materialIndices[materialName] = currentIndex;
            currentIndex++;
        }
    }

    // Initialize the element_materials vector with zeros equal to the largest element index
    auto maxElementIter = std::max_element(elementTags.begin(), elementTags.end());
    element_materials.assign(*maxElementIter+1, 0);
    std::cerr << "Maximum element index:" << *maxElementIter << std::endl;

    // Iterate through the material names and assign material indices to elements
    for (const auto& entry : materialProperties) {
        const std::string& materialName = entry.first;
        const std::vector<std::size_t> elementsForMaterial = getElementsForPhysicalGroup(materialName);
        // Retrieve the material index for this material name
        std::size_t materialIndex = materialIndices[materialName];
        std::cout << "MaterialName: " << materialName << ", MaterialIndex: " << materialIndex << ", Elements: " << elementsForMaterial.size() << std::endl;

        // Assign the material index to the corresponding elements
        for (std::size_t j : elementsForMaterial) {
            if (j <= *maxElementIter+1) {
                element_materials[j] = static_cast<int>(materialIndex);
            } else {
                std::cerr << "Warning: Element index out of range, "<< j << " ." << std::endl;
            }
        }
    }
    // Check that no element in element_materials is equal to zero
    std::vector<std::size_t> elementsWithoutMaterial;

    for (std::size_t i = 0; i < elementTags.size(); i++) {
        std::size_t elementIndex = elementTags[i];
        if (elementIndex < element_materials.size()) {
            int materialIndex = element_materials[elementIndex];
            if (materialIndex == 0) {
                elementsWithoutMaterial.push_back(elementIndex);
            }
        } else {
            std::cerr << "Warning: Element tag " << elementIndex << " exceeds the size of element_materials." << std::endl;
        }
    }

    // Print the list of elements without material properties
    if (!elementsWithoutMaterial.empty()) {
        std::cerr << "Warning: The following elements have no material properties:" << std::endl;
        for (std::size_t elementIndex : elementsWithoutMaterial) {
            std::cerr << elementIndex << " ";
        }
        std::cerr << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
double Mesh::getMaterialPropertyForElement(std::size_t elementIndex, const std::string& propertyName) const {
    double propertyValue = 0.0; // Initialize to a default value

    // Check if the element index is within the bounds of element_materials
    if (elementIndex < numelem) {
        // Retrieve the material index for the given element
        int materialIndex = element_materials[elementIndex];

        // Find the corresponding material name using materialIndices
        std::string materialName;

        for (const auto& entry : materialIndices) {
            if (entry.second == static_cast<std::size_t>(materialIndex)) {
                materialName = entry.first;
                break;
            }
        }
        std::map<std::string, std::map<std::string, double>> materialProperties = inputReader_.getAllMaterialProperties();

        // Find the material properties for the materialName in materialProperties
        auto materialPropsIterator = materialProperties.find(materialName);

        if (materialPropsIterator != materialProperties.end()) {
            // Try to find the specified property in the material properties map
            const std::map<std::string, double>& materialProps = materialPropsIterator->second;
            auto propertyIterator = materialProps.find(propertyName);

            if (propertyIterator != materialProps.end()) {
                // Retrieve the property value
                propertyValue = propertyIterator->second;
            } else {
                std::cerr << "Warning: Property '" << propertyName << "' not found for element index " << elementIndex << std::endl;
            }
        } else {
            std::cerr << "Warning: Material properties not found for element index " << elementIndex << std::endl;
        }
    } else {
        std::cerr << "Error: Element index " << elementIndex << " is out of bounds." << std::endl;
    }

    return propertyValue;
}
