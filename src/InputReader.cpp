#include "InputReader.h"

InputReader::InputReader(const std::string& filename) : filename_(filename) {}

bool InputReader::readFile() {
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


