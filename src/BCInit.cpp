#include "BCInit.h"
#include "gmsh.h"

BCInit::BCInit(const InputReader& inputReader)
    : inputReader_(inputReader) {
    // You may want to initialize other member variables here
}

void BCInit::boundaryConditions() const {
    // Get the boundary conditions from the InputReader
    const auto& boundaryConditions = inputReader_.getBoundaryConditions();

    // Print each boundary condition
    for (const auto& bc : boundaryConditions) {
        std::string boundaryName = bc.first.first;
        std::string surfaceName = bc.first.second;
        double value = bc.second;

        std::cout << "Boundary Name: " << boundaryName
                  << ", Surface Name: " << surfaceName
                  << ", Value: " << value << std::endl;
    }
}

void BCInit::createLoadVector() {
    // Read the mesh file using gmsh library
    gmsh::initialize();
    gmsh::open(inputReader_.getMeshFileName());
    // Get the boundary conditions from the InputReader
    const auto& boundaryConditions = inputReader_.getBoundaryConditions();

    // Iterate over the boundary conditions and create the load vector
    for (const auto& bc : boundaryConditions) {
        // Get the boundary name, surface name, and value
        std::string boundaryName = bc.first.first;
        std::string surfaceName = bc.first.second;
        double value = bc.second;
        std::cout << value << std::endl;

        // Use the gmsh library to get the nodes on the surface
        std::vector<std::size_t> nodeTags;
        /*gmsh::model::mesh::getNodesForPhysicalGroup(2, surfaceName, nodeTags);

        // Add the value to the load vector for each node on the surface
        for (std::size_t nodeTag : nodeTags) {
            // Assuming loadVector_ is a std::vector<double> with size equal to the number of nodes in the mesh
            loadVector_[nodeTag] += value;
        }*/
    }
}
