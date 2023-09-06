#include "BCInit.h"
using namespace arma;

BCInit::BCInit(const InputReader& inputReader, Mesh& mesh)
    : inputReader_(inputReader), mesh_(mesh), elements_(), utils_() {
    // Get the number of nodes from the mesh
    int numNodes = mesh_.getNumAllNodes();

    // Initialize the loadVector_ member with a size double the number of nodes
    loadVector_.zeros(2 * numNodes, 1);
    initialdofs_.resize(2 * numNodes);

    // Initialize the integrationFunction_ using the assignment operator
    integrationFunction_ = [&](const arma::mat& natcoords, const arma::mat& coords, double value) {
        return CteSurfBC(natcoords, coords, value);
    };

    boundaryConditions();
}


void BCInit::boundaryConditions() {
    std::cout << "#BCInit::boundaryConditions" << std::endl;
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
        }else if (boundaryName == "heat_n") {
            // Assuming 'getNodesForPhysicalGroup' returns a vector of unsigned long long integers
            const auto& elements = mesh_.getElementsForPhysicalGroup(surfaceName);
            // Loop through the nodes and set the corresponding values in initialdofs_ (loadVector_)
            for (unsigned long long element : elements) {

                utils_.gaussIntegrationBC(2, 3, element, mesh_, value, integrationFunction_, loadVector_);

            }
    }

}
}


mat BCInit::CteSurfBC(const mat& natcoords, const mat& coords, double value) {
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
    F_q = detJ * (shapeFunctions * value);

    // Return the heat flow as a 4x1 Armadillo matrix
    return F_q;
}

