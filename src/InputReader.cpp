#include "InputReader.h"

InputReader::InputReader(const std::string& filename) : filename_(filename) {
    readFile();
    std::cout << "Read Input." << std::endl;
}

bool InputReader::readFile() {
    std::cout << "#INPUTREADER::readFile" << std::endl;
    std::ifstream inputFile(filename_);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return false;
    }

    std::string line;
    setMeshFileName("defaultmesh.gmsh");

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string keyword;
        iss >> keyword;

        if (keyword == "mesh_file") {
            // Read mesh file and associated entity name
            std::string meshFileName;
            std::string entityName;
            
            iss >> meshFileName >> entityName;
            
            if (!iss) {
                std::cerr << "Warning: Invalid mesh file and entity name." << std::endl;
                continue;
            }else{
                std::cout << "New Mesh File and Mesh entity:"<<meshFileName<<" "<< entityName<< std::endl;
            }
            
            setMeshFileName(meshFileName);       // Store the mesh file name in the class
            setMeshEntityName(entityName);           // Set the full mesh entity name
            // MATERIAL PROPERTIES
        }else if (keyword == "material") {
            std::string volumeName;
            iss >> volumeName;

            if (!iss) {
                std::cerr << "Warning: Invalid material ID." << std::endl;
                continue;
            }

            std::map<std::string, double> materialProperties;
            int matIndex = 1;
            while (std::getline(inputFile, line)) {
                if (line.empty() || line.find("##") == 0) {
                    break;
                }

                std::istringstream lineIss(line);
                std::string propertyName;
                double propertyValue;
                lineIss >> propertyName >> propertyValue;

                if (lineIss) {
                    materialProperties[propertyName] = propertyValue;
                } else {
                    std::cerr << "Warning: Invalid property format." << std::endl;
                    std::cerr << "Line Content: " << line << std::endl;
                }
            }

            MaterialProperties_[volumeName] = materialProperties;

            Material_index[matIndex] = volumeName; // Store the index and name in Material_index map
            std::cout<< "Index and material:" << matIndex<<" " <<Material_index[matIndex]<<std::endl;
            ++matIndex;  // Increment the index for the next material
        }else if (keyword == "bc") {
            std::string boundaryName, surfaceName;
            double value;
            iss >> boundaryName >> surfaceName >> value;

            if (!iss) {
                std::cerr << "Warning: Invalid boundary condition format." << std::endl;
                std::cerr << "Line Content: " << line << std::endl;
                continue;
            }

            boundaryConditions_[std::make_pair(boundaryName, surfaceName)] = value;
        }
    }

    inputFile.close();
    return true;
}


// PRINTING OPTIONS
void InputReader::printMaterialProperties() const {
    for (const auto& entry : MaterialProperties_) {
        std::string volumeName = entry.first;
        const std::map<std::string, double>& materialProperties = entry.second;

        std::cout << "Material Volume: " << volumeName << std::endl;
        for (const auto& propertyEntry : materialProperties) {
            std::cout << "Property: " << propertyEntry.first
                      << ", Value: " << propertyEntry.second << std::endl;
        }
    }
}

void InputReader::printBoundaryConditions() const {
    for (const auto& entry : boundaryConditions_) {
        std::pair<std::string, std::string> boundarySurfacePair = entry.first;
        double value = entry.second;

        std::cout << "Boundary Name: " << boundarySurfacePair.first
                  << ", Surface Name: " << boundarySurfacePair.second
                  << ", Value: " << value << std::endl;
    }
}

// MESH READING
void InputReader::setMeshFileName(const std::string& meshFileName) {
    meshFileName_ = meshFileName;
    std::cout << "New mesh file name set: " << meshFileName_ << std::endl;
}

void InputReader::setMeshEntityName(const std::string& MeshEntityName) {
    MeshEntityName_ = MeshEntityName;
    std::cout << "New mesh group name set: " << MeshEntityName_ << std::endl;
}


bool InputReader::readGmshMesh() {
    return true;
}


// Function to get a material property value by index and property namedouble InputReader::getMaterialPropertyValue(int matidx, const std::string& propertyName) const {
double InputReader::getMaterialPropertyValue(int matidx, const std::string& propertyName) const {
    // Find the material name corresponding to matidx
    std::string materialName; // Declare materialName without & to avoid an uninitialized reference
    for (const auto& entry : Material_index) {
        int idx = entry.first;
        if (idx == matidx) {
            materialName = entry.second;
            //std::cout << "Found material name: " << materialName << std::endl;
            break; // No need to continue searching once we've found the materialName
        }
    }

    // Check if materialName was found
    if (materialName.empty()) {
        // Handle the case where matidx doesn't correspond to any material
        // You can throw an exception or return a default value as needed
        // For now, let's return a special value like -1.0 to indicate an error
        std::cout << "Material name not found for matidx: " << matidx << std::endl;
        return -1.0;
    }

    // Search for the property value using materialName
    for (const auto& entry : MaterialProperties_) {
        std::string volumeName = entry.first;
        if (volumeName == materialName) {
            const std::map<std::string, double>& materialProperties = entry.second;
            for (const auto& propertyEntry : materialProperties) {
                if (propertyEntry.first == propertyName) {
                    //std::cout << "Found property value: " << propertyEntry.second << std::endl;
                    return propertyEntry.second;
                }
            }
        }
    }

    // Handle the case where propertyName doesn't correspond to any property
    // You can throw an exception or return a default value as needed
    // For now, let's return a special value like -1.0 to indicate an error
    std::cout << "Property not found for propertyName: " << propertyName << std::endl;
    return -1.0;
}
