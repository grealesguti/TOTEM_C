#include "BCInit.h"
using namespace arma;

BCInit::BCInit(const InputReader& inputReader, Mesh& mesh)
    : inputReader_(inputReader), mesh_(mesh), elements_(), utils_(), loadVector_(2 * mesh.getNumAllNodes(), 1) {
    // Initialize Utils based on inputReader
    utils_ = inputReader.getDesiredOutput() == "all" ? Utils(true) : Utils(false);

    // Initialize initialdofs_ based on the number of nodes
    initialdofs_.resize(2 * mesh.getNumAllNodes());

    // Initialize integrationFunction_ using a lambda function
    integrationFunction_ = [&](const arma::mat& natcoords, const arma::mat& coords, double value, int element) {
        return CteSurfBC(natcoords, coords, value, element);
    };

    // Call boundaryConditions to perform any necessary setup
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
            std::size_t elementTag = elements[elementindexVector[0]];    std::vector<int> nodeTags_el;
            int etype = mesh_.getElementInfo(elementTag, nodeTags_el);
            // https://docs.juliahub.com/GmshTools/9rYp5/0.4.2/element_types/
            int nodes_per_element=8;// DEFAULT
            if(etype==3){// 4-node quadrangle
                nodes_per_element = 4; // Total number of nodes per element
            }else if(etype==16){// 8-node second order quadrangle
                nodes_per_element = 8; // Total number of nodes per element
            }
            std::cout << "Element of type: " << etype<< " nodes_per_element " <<nodes_per_element<< std::endl;

            arma::mat element_load_vector(nodes_per_element, 1, arma::fill::zeros);
            // Get and print the number of threads
            int num_threads = omp_get_max_threads();
            std::cout << "Number of threads to be used: " << num_threads << std::endl;

            //#pragma omp parallel for shared(loadVector_) private(element_load_vector)
            for (std::size_t elementindex : elementindexVector) {
                // Reset element_load_vector to zeros before each element
                std::size_t element = elements[elementindex];
                element_load_vector.set_size(nodes_per_element, 1); // Resize element_load_vector to 4x1

                element_load_vector.zeros();

                // Calculate element_load_vector for the current element
                utils_.gaussIntegrationBC(2, 3, element, mesh_, value, integrationFunction_, element_load_vector);

                // Debugging: Print the current element being processed
                 std::cout << "Processing element: " << element << std::endl;

                // Loop through element nodes and update loadVector_
                for (std::size_t i = elementindex * nodes_per_element; i < (elementindex + 1) * nodes_per_element; ++i) {
                    if (i < element_nodes.size()) {
                        std::size_t node = element_nodes[i];
                        if (node * 2 + 1 < loadVector_.size()) {
                            // Use atomic operation to avoid race condition
                            //#pragma omp atomic
                            loadVector_(node * 2) += element_load_vector(i - elementindex * nodes_per_element, 0);
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


mat BCInit::CteSurfBC(const mat& natcoords, const mat& coords, double value, int element) {
    // Define variables
    std::size_t elementTag = element;    std::vector<int> nodeTags_el;
    int etype = mesh_.getElementInfo(elementTag, nodeTags_el);
    // https://docs.juliahub.com/GmshTools/9rYp5/0.4.2/element_types/
    int nodes_per_element=8;// DEFAULT
    if(etype==3){// 4-node quadrangle
        nodes_per_element = 4; // Total number of nodes per element
    }else if(etype==16){// 8-node second order quadrangle
        nodes_per_element = 8; // Total number of nodes per element
    }

    arma::vec shapeFunctions(nodes_per_element,1);          // Shape functions as a 4x1 vector
    mat shapeFunctionDerivatives(nodes_per_element, 2); // Shape function derivatives
    //std::cout << "Initialize shape functions and derivatives. " << std::endl;
    mat F_q(nodes_per_element, 1, fill::zeros);         // Initialize F_q as a 3x1 zero matrix for heat flow
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
    if(etype==3){// 4-node quadrangle
        shapeFunctions = elements_.EvaluateLinearQuadrilateralShapeFunctions(xi, eta);
        elements_.EvaluateLinearQuadrilateralShapeFunctionDerivatives(xi, eta, shapeFunctionDerivatives);
    }else if(etype==16){// 8-node second order quadrangle
        elements_.EvaluateQuadraticQuadrilateralShapeFunctions(xi, eta,shapeFunctions);
        elements_.CalculateQuadraticQuadrilateralShapeFunctionDerivatives(xi, eta, shapeFunctionDerivatives);
    }

    //std::cout << "Calculate jacobian. " << std::endl;
    // Calculate Jacobian matrix JM
    mat JM = shapeFunctionDerivatives.t() * coords.t();
    //"Calculate jacobian determinant. " << std::endl;
    // Calculate the determinant of the Jacobian
    double detJ = arma::det(JM);
    //std::cout << "Calculate integrand. " << std::endl;
    // Calculate F_q (heat flow) using your equations
    F_q = detJ * (shapeFunctions * value);
    //std::cout << "integrand calculated " << std::endl;
    if (inputReader_.getDesiredOutput()=="all"){
        utils_.writeDataToFile(shapeFunctions,"Outputs/HeatIntegrationShapeFunctions_"+std::to_string(element)+".txt",true);
        utils_.writeDataToFile(shapeFunctionDerivatives,"Outputs/HeatIntegrationShapeFunctionsDerivatives_"+std::to_string(element)+".txt",true);
        utils_.writeDataToFile(JM,"Outputs/HeatIntegrationJM_"+std::to_string(element)+".txt",true);
        utils_.writeDataToFile(coords,"Outputs/HeatIntegrationElcoords_"+std::to_string(element)+".txt",true);
        utils_.writeDataToFile(natcoords,"Outputs/HeatIntegrationNatcoords_"+std::to_string(element)+".txt",true);
        utils_.writeDataToFile(F_q,"Outputs/HeatIntegrationFq_"+std::to_string(element)+".txt",true);
    }
    // Return the heat flow as a 4x1 Armadillo matrix
    return F_q;
}

/////////////////////////////////////////////////////////////////////////
