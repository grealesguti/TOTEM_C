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
            auto result = mesh_.getElementsAndNodeTagsForPhysicalGroup(surfaceName);
            std::vector<std::size_t> elements = result.first;
            std::vector<std::size_t> element_nodes = result.second;
            std::vector<int> elementindexVector(elements.size());

            // Initialize elementindexVector
            std::iota(elementindexVector.begin(), elementindexVector.end(), 0);

            // Debugging: Print the sizes of elements and element_nodes vectors
            std::cout << "Size of 'elements' vector: " << elements.size() << std::endl;
            std::cout << "Size of 'element_nodes' vector: " << element_nodes.size() << std::endl;


            std::cout << "HEAT INTEGRATION." << std::endl;
            arma::mat element_load_vector(4, 1, arma::fill::zeros);
            // Get and print the number of threads
            int num_threads = omp_get_max_threads();
            std::cout << "Number of threads to be used: " << num_threads << std::endl;

            #pragma omp parallel for shared(loadVector_) private(element_load_vector)
            for (std::size_t elementindex : elementindexVector) {
                // Reset element_load_vector to zeros before each element
                std::size_t element = elements[elementindex];
                element_load_vector.set_size(4, 1); // Resize element_load_vector to 4x1

                element_load_vector.zeros();

                // Calculate element_load_vector for the current element
                utils_.gaussIntegrationBC(2, 3, element, mesh_, value, integrationFunction_, element_load_vector);

                // Debugging: Print the current element being processed
                // std::cout << "Processing element: " << element << std::endl;

                // Loop through element nodes and update loadVector_
                for (std::size_t i = elementindex * 4; i < (elementindex + 1) * 4; ++i) {
                    if (i < element_nodes.size()) {
                        std::size_t node = element_nodes[i];
                        if (node * 2 + 1 < loadVector_.size()) {
                            // Use atomic operation to avoid race condition
                            #pragma omp atomic
                            loadVector_(node * 2 + 1) += element_load_vector(i - elementindex * 4, 0);
                        } else {
                            std::cerr << "Error: Node index out of bounds!" << std::endl;
                        }
                    } else {
                        std::cerr << "Error: Element node index out of bounds!" << std::endl;
                    }
                }
            }


        
        std::cout << "HEAT INTEGRATION FINISHED." << std::endl;

    }

    }
    mesh_.SetFreedofsIdx(); // store the index of the dofs that are not fixed in order
}


mat BCInit::CteSurfBC(const mat& natcoords, const mat& coords, double value) {
    // Define variables
    //std::cout << "Initialize shape functions and derivatives. " << std::endl;
    arma::vec shapeFunctions(4,1);          // Shape functions as a 4x1 vector
    mat shapeFunctionDerivatives(4, 2); // Shape function derivatives
    mat JM(2, 2);                       // Jacobian matrix
    mat F_q(4, 1, fill::zeros);         // Initialize F_q as a 3x1 zero matrix for heat flow
    //std::cout << "Extract natural coordinates. " << std::endl;
    double xi=0, eta=0;
    if (natcoords.n_rows >= 2 && natcoords.n_cols >= 1) {
        xi = natcoords(0, 0);  // Extracts the first element (a)
        eta = natcoords(1, 0); // Extracts the second element (b)
        // Rest of your code...
    } else {
        // Handle the case when natcoords doesn't have the expected dimensions.
        std::cout << "Wrong natural coordinates dimension. " << std::endl;
        // You may want to print an error message or take appropriate action.
    }
    // Calculate shape functions and their derivatives
    //std::cout << "Get shape functions and derivatives. " << std::endl;
    shapeFunctions = elements_.EvaluateLinearQuadrilateralShapeFunctions(xi, eta);
    elements_.EvaluateLinearQuadrilateralShapeFunctionDerivatives(xi, eta, shapeFunctionDerivatives);
    //std::cout << "Calculate jacobian. " << std::endl;
    // Calculate Jacobian matrix JM
    for (uword i = 0; i < 2; ++i) {
        for (uword j = 0; j < 2; ++j) {
            JM(i, j) = dot(shapeFunctionDerivatives.col(i), coords.row(j));
        }
    }
    //std::cout << "Calculate jacobian determinant. " << std::endl;
    // Calculate the determinant of the Jacobian
    double detJ = arma::det(JM);
    //std::cout << "Calculate integrand. " << std::endl;
    // Calculate F_q (heat flow) using your equations
    F_q = detJ * (shapeFunctions * value);
    //std::cout << "integrand calculated " << std::endl;

    // Return the heat flow as a 4x1 Armadillo matrix
    return F_q;
}

