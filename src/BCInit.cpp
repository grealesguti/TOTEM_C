#include "BCInit.h"
using namespace arma;

BCInit::BCInit(const InputReader& inputReader, Mesh& mesh)
    : inputReader_(inputReader), mesh_(mesh), elements_() {
    // Get the number of nodes from the mesh
    int numNodes = mesh_.getNumAllNodes();

    // Initialize the loadVector_ member with a size double the number of nodes
    loadVector_.resize(2 * numNodes);
    initialdofs_.resize(2 * numNodes);

    // You may want to initialize other member variables here
}

void BCInit::boundaryConditions() {
    // Get the boundary conditions from the InputReader
    const auto& boundaryConditions = inputReader_.getBoundaryConditions();

// Iterate through each boundary condition
    for (const auto& bc : boundaryConditions) {
        const std::string& boundaryName = bc.first.first;
        const std::string& surfaceName = bc.first.second;
        double value = bc.second;

        std::cout << "Boundary Name: " << boundaryName
                << ", Surface Name: " << surfaceName
                << ", Value: " << value << std::endl;

        if (boundaryName == "Temperature") {
            // Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
            const auto& nodes = mesh_.getNodesForPhysicalGroup(surfaceName);
            // Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
            for (unsigned long long node : nodes) {
                    // Make sure node is within the bounds of initialdofs_
                    if (node*2 < initialdofs_.size()) {
                        initialdofs_[node*2] = value;
                        mesh_.setFixedof(node*2);
                    } else {
                        // Handle the case where the node index is out of bounds.
                        // This could be an error condition depending on your application.
                        std::cerr << "Node index out of bounds: " << node << std::endl;
                    }
            }
        }else if (boundaryName == "Voltage") {
            // Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
            const auto& nodes = mesh_.getNodesForPhysicalGroup(surfaceName);
            // Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
            for (unsigned long long node : nodes) {
                    // Make sure node is within the bounds of initialdofs_
                    if (node*2+1 < initialdofs_.size()) {
                        initialdofs_[node*2+1] = value;
                        mesh_.setFixedof(node*2+1);
                    } else {
                        // Handle the case where the node index is out of bounds.
                        // This could be an error condition depending on your application.
                        std::cerr << "Node index out of bounds: " << node << std::endl;
                    }
            }
        }else if (boundaryName == "heat_x") {
            // Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
            const auto& nodes = mesh_.getNodesForPhysicalGroup(surfaceName);
            // Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
            for (unsigned long long node : nodes) {
                    // Make sure node is within the bounds of initialdofs_
                    if (node*2+1 < initialdofs_.size()) {
                        initialdofs_[node*2+1] = value;
                        mesh_.setFixedof(node*2+1);
                    } else {
                        // Handle the case where the node index is out of bounds.
                        // This could be an error condition depending on your application.
                        std::cerr << "Node index out of bounds: " << node << std::endl;
                    }
            }
    }

}
}


mat BCInit::CteSurfBC(const mat& natcoords,const mat& coords, double heatvalue) {
    // Define variables
    vec shapeFunctions(4);             // Shape functions
    mat shapeFunctionDerivatives(4, 2); // Shape function derivatives
    mat JM(2, 2);                       // Jacobian matrix
    mat F_q(4, 1, fill::zeros);         // Initialize F_q as a 3x1 zero matrix for heat flow
    double xi = natcoords(0, 0);  // Extracts the first element (a)
    double eta = natcoords(1, 0); // Extracts the second element (b)
    // Calculate shape functions and their derivatives
    elements_.EvaluateLinearQuadrilateralShapeFunctions(xi, eta, shapeFunctions);
    elements_.EvaluateLinearQuadrilateralShapeFunctionDerivatives(xi, eta, shapeFunctionDerivatives);

    // Calculate Jacobian matrix JM
    for (uword i = 0; i < 2; ++i) {
        for (uword j = 0; j < 2; ++j) {
            JM(i, j) = dot(shapeFunctionDerivatives.col(i), coords.row(j));
        }
    }

    // Calculate the determinant of the Jacobian
    double detJ = det(JM);
    // Calculate F_q (heat flow) using your equations
    F_q = detJ * (shapeFunctions * heatvalue);

    // Return the heat flow as a 3x1 Armadillo matrix
    return F_q;
}