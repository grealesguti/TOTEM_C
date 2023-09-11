#include "Solver.h"
using namespace arma;

Solver::Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit)
    : inputReader_(inputReader), mesh_(mesh), elements_(), bcinit_(bcinit), utils_() {
    // Get the number of nodes from the mesh
    int numNodes = mesh_.getNumAllNodes();

    // Initialize the loadVector_ member with a size double the number of nodes
    loadVector_= bcinit_.getloadVector();
    soldofs_.resize(2 * numNodes);
    freedofidxs_ = mesh_.GetFreedofsIdx();

    // Initialize the thermoelectricityintegrationFunction_ using the assignment operator
    thermoelectricityintegrationFunction_ = [&](const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs) -> Utils::IntegrationResult {
        return thermoelectricityintegration(natcoords, coords, dofs);
    };

}

std::pair<arma::mat, arma::mat> Solver::Assembly() {
    int nodes_per_element = 8; // Total number of nodes per element
    int dof_per_node = 2; // Assuming 2 degrees of freedom per node
    const std::vector<std::size_t>& elementTags = mesh_.getElementTags();

    // Create matrices to store global assembly
    arma::mat Ra = arma::zeros<arma::mat>(dof_per_node * nodes_per_element, 1);
    arma::mat KJ = arma::zeros<arma::mat>(dof_per_node * nodes_per_element, dof_per_node * nodes_per_element);

    // Temporary matrices for the current element
    arma::mat element_Ra = arma::zeros<arma::mat>(16, 1);
    arma::mat element_KJ = arma::zeros<arma::mat>(16, 16);
    arma::mat element_dofs = arma::zeros<arma::mat>(16, 1);

    // Iterate over each element tag
    for (std::size_t elementTag : elementTags) {
        // Call the getElementInfo function to retrieve information about the element
        std::vector<int> nodeTags_el;
        mesh_.getElementInfo(elementTag, nodeTags_el);
        
        int cc = 0;
        for (int nodeTag : nodeTags_el) {
            element_dofs[cc] = nodeTag * dof_per_node;
            element_dofs[cc + nodes_per_element] = nodeTag * dof_per_node + 1;
            cc += 1;
        }

        // Resize temporary matrices for the current element
        element_KJ.set_size(16, 16);
        element_KJ.zeros();
        element_Ra.set_size(16, 1);
        element_Ra.zeros();

        // Calculate element_load_vector for the current element
        /*
        utils_.gaussIntegrationAssembly(2, 3, elementTag, mesh_, nodeTags_el, thermoelectricityintegrationFunction_, element_KJ, element_Ra);

        // Assembly in global residual and jacobian matrix
        Ra.submat(element_dofs, 0) += element_Ra;
        KJ.submat(element_dofs, element_dofs) += element_KJ;*/
    }

    /*
    // Reduced System
    arma::mat R_b = Ra.rows(freedofidxs_) - loadVector_(freedofidxs_);
    arma::mat KJ_b = KJ.submat(freedofidxs_, freedofidxs_);*/
    arma::mat R_b = {1};
    arma::mat KJ_b = {1};
    // Return the results
    return std::make_pair(R_b, KJ_b);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


Utils::IntegrationResult Solver::thermoelectricityintegration(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs){
        Utils::IntegrationResult result; // Create a struct to hold KV and R

    // Define variables
    //std::cout << "Initialize shape functions and derivatives. " << std::endl;
    arma::vec shapeFunctions(8,1);          // Shape functions as a 4x1 vector
    mat shapeFunctionDerivatives(8, 3); // Shape function derivatives
    mat JM(3, 3);                       // Jacobian matrix
    // Define the integration result matrices
    arma::mat RT(8, 8, arma::fill::zeros);
    arma::mat RV(8, 1, arma::fill::zeros);
    arma::mat K11(8, 8, arma::fill::zeros);
    arma::mat K12(8, 1, arma::fill::zeros);
    arma::mat K21(1, 8, arma::fill::zeros);
    arma::mat K22(1, 1, arma::fill::zeros);

    //std::cout << "Extract natural coordinates. " << std::endl;
    double xi=0, eta=0, zeta=0;
    if (natcoords.n_rows >= 3 && natcoords.n_cols >= 1) {
        xi = natcoords(0, 0);  // Extracts the first element (a)
        eta = natcoords(1, 0); // Extracts the second element (b)
        zeta = natcoords(2, 0); // Extracts the second element (b)
    } else {
        // Handle the case when natcoords doesn't have the expected dimensions.
        std::cout << "Wrong natural coordinates dimension. " << std::endl;
        // You may want to print an error message or take appropriate action.
    }
    // Calculate shape functions and their derivatives
    //std::cout << "Get shape functions and derivatives. " << std::endl;
    elements_.EvaluateHexahedralLinearShapeFunctions(xi, eta, zeta, shapeFunctions);
    elements_.CalculateHexahedralLinearShapeFunctionDerivatives(xi, eta, zeta, shapeFunctionDerivatives);
    //std::cout << "Calculate jacobian. " << std::endl;
    // Calculate Jacobian matrix JM
    for (uword i = 0; i < 3; ++i) {
        for (uword j = 0; j < 3; ++j) {
            JM(i, j) = dot(shapeFunctionDerivatives.col(i), coords.col(j)); // Fixed the loop indexing
        }
    }

    //std::cout << "Calculate jacobian determinant. " << std::endl;
    // Calculate the determinant of the Jacobian
    double detJ = arma::det(JM);

    // Extract material properties
    double De=1, Da=1, Dk=1, Dde=0, Dda=0, Ddk=0;

    // Assuming 'dofs' is an Armadillo vector
    arma::mat Tee(8,1); // Vector for odd-indexed elements
    arma::mat Vee(8,1); // Vector for even-indexed elements
    // Extract odd-indexed elements into Tee
    for (uword i = 1; i < dofs.n_elem; i += 2) {
        Tee((i-1) / 2,1) = dofs(i);
    }

    // Extract even-indexed elements into Vee
    for (uword i = 0; i < 16; i += 2) {
        Vee(i / 2,1) = dofs(i);
    }
    double Th = arma::dot(shapeFunctions, dofs);

    // Calculate je and qe
    arma::mat transposed = shapeFunctionDerivatives.t();
    arma::mat negTransposed = -transposed;

    arma::mat test = negTransposed * Vee;

    arma::mat je = -De * shapeFunctionDerivatives.t() * Vee - Da * De * shapeFunctionDerivatives.t() * Tee;
    arma::mat qe = Da * (shapeFunctions.t() * Tee) * je - Dk * shapeFunctionDerivatives.t() * Tee;
    arma::mat test = shapeFunctions * Tee;

    // Calculate djdt, djdv, dqdt, and dqdv
    arma::mat djdt, djdv, dqdt, dqdv;
    // Perform the necessary calculations here to compute djdt, djdv, dqdt, dqdv
    arma::mat djdt = -Da * De * shapeFunctionDerivatives.t()
                    - Dda * De * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t()
                    - Dde * (shapeFunctionDerivatives.t() * Vee.t() + Da * shapeFunctionDerivatives.t() * Tee.t()) * shapeFunctions.t();

    arma::mat djdv = -De * shapeFunctionDerivatives;

    arma::mat dqdt = Da * Th * djdt + Da * je * shapeFunctions.t()
                    - Dk * shapeFunctionDerivatives.t()
                    + Dda * Th * je * shapeFunctions.t()
                    - Ddk * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t();

    arma::mat dqdv = -Da * De * Th * shapeFunctionDerivatives;
    // Perform the integration for RT, RV, K11, K12, K21, K22
    RT += (-shapeFunctionDerivatives * qe + shapeFunctions * je.t() * shapeFunctionDerivatives.t() * Vee);
    RV += (-shapeFunctionDerivatives * je);
    K11 +=(shapeFunctionDerivatives * dqdt - shapeFunctions * (djdt.t() * shapeFunctionDerivatives.t() * Vee).t());
    K12 +=(shapeFunctionDerivatives * dqdv - shapeFunctions * (djdv.t() * shapeFunctionDerivatives.t() * Vee).t() - N * (je.t() * shapeFunctionDerivatives));
    K21 +=(shapeFunctionDerivatives * djdt);
    K22 +=(shapeFunctionDerivatives * djdv);


    //std::cout << "Calculate integrand. " << std::endl;
    // Define the size of KV and RT matrices based on the sizes of K11, K12, K21, K22, RT, and RV
    int numRowsKJ = K11.n_rows + K21.n_rows;
    int numColsKJ = K11.n_cols + K12.n_cols;
    int numRowsR = RT.n_rows + RV.n_rows;

    // Initialize KV and RT matrices with zeros
    arma::mat KV(numRowsKJ, numColsKJ, arma::fill::zeros);
    arma::mat R(numRowsR, 1, arma::fill::zeros);

    // Fill in the KV matrix with K11, K12, K21, and K22
    KV.submat(0, 0, K11.n_rows - 1, K11.n_cols - 1) = K11;
    KV.submat(0, K11.n_cols, K11.n_rows - 1, numColsKJ - 1) = K12;
    KV.submat(K11.n_rows, 0, numRowsKJ - 1, K11.n_cols - 1) = K21;
    KV.submat(K11.n_rows, K11.n_cols, numRowsKJ - 1, numColsKJ - 1) = K22;

    // Fill in the R matrix with RT and RV
    R.submat(0, 0, RT.n_rows - 1, 0) = RT;
    R.submat(RT.n_rows, 0, numRowsR - 1, 0) = RV;

    //std::cout << "integrand calculated " << std::endl;
    result.R=R;
    result.KT=KV;
    // Return the heat flow as a 4x1 Armadillo matrix
    return result;
}


/*
// make another function only as solver that takes the assembly and can choose btw modified NR, or normal and returns successive errors etc...
// This is the function to call from main after initialization of the class. This function also identifies the physics to use


// Te, Ve defined in Solver, need function to extract them from node tags
// matvs defined in Mesh
// Introduce natural and overall coords
// xx stored in Mesh
// penalty values stored in Input
// numnodes initially 8 to linear, possible to read from gmsh, type of element
std::tuple<mat, mat> CreateElementIntegrationFunction(const mat& cooro, const mat& sysv)  {
    return [](const mat& naturalCoords, const mat& elementCoords) -> mat {
        // Extract natural coordinates
        double xi1 = naturalCoords(0);
        double xi2 = naturalCoords(1);
        double xi3 = naturalCoords(2);

        // Extract element coordinates
        mat cooro = elementCoords;

        // Initialize matrices based on numNodesPerElement
        mat N = zeros<mat>(numNodesPerElement, 1);
        mat dShape = zeros<mat>(numNodesPerElement, 3);
        mat matp = zeros<mat>(numNodesPerElement, numNodesPerElement); // Adjust dimensions based on your problem
        // Initialize other material-related matrices based on your problem

        // Initialize other matrices (K11, K12, K21, K22, RT, RV) based on your problem's requirements
        mat K11 = zeros<mat>(numNodesPerElement, numNodesPerElement);
        mat K12 = zeros<mat>(numNodesPerElement, numNodesPerElement);
        mat K21 = zeros<mat>(numNodesPerElement, numNodesPerElement);
        mat K22 = zeros<mat>(numNodesPerElement, numNodesPerElement);
        mat RT = zeros<mat>(numNodesPerElement, 1);
        mat RV = zeros<mat>(numNodesPerElement, 1);

        // Get shape functions N and their derivatives dShape
        mat N; // Calculate N based on xi1, xi2, and xi3
        mat dShape; // Calculate dShape based on xi1, xi2, and xi3

        // Calculate Jacobian matrix JM and its inverse Jacinv
        mat JM = dShape * cooro;
        mat Jacinv = inv(JM);

        // Calculate DN and DN0
        mat DN = Jacinv * dShape;
        mat DN0 = DN.t();

        // Extract values from elementCoords
        mat Tee = cooro; // Assuming Tee represents some temperature values
        mat Vee = cooro; // Assuming Vee represents some other values

        // Calculate Th
        double Th = dot(N, Tee);

        // Define material properties
        mat matp, matv, localsys, sysv, xx, p, seebp, rhop; // Define these matrices
        mat De, Da, Dk, Dde, Dda, Ddk; // Calculate De, Da, Dk, Dde, Dda, and Ddk based on material properties

        // Initialize Wp and detJ
        double Wp = 0.0;
        double detJ = det(JM);

        // Calculate je
        mat je = -De * DN0.t() * Vee - Da * De * DN0.t() * Tee;

        // Calculate qe
        mat qe = Da * (N.t() * Tee) * je - Dk * DN0.t() * Tee;

        // Calculate djdt, djdv, dqdt, dqdv
        mat djdt = -Da * De * DN0.t() - Dda * De * DN0.t() * Tee * N.t() - Dde * (DN0.t() * Vee + Da * DN0.t() * Tee) * N.t();
        mat djdv = -De * DN0.t();
        mat dqdt = Da * Th * djdt + Da * je * N.t() - Dk * DN0.t() + Dda * Th * je * N.t() - Ddk * DN0.t() * Tee * N.t();
        mat dqdv = -Da * De * Th * DN0.t();

        // Calculate RT and RV
        RT = detJ * Wp * (-(DN0 * qe) + (N * je) * (DN0.t() * Vee));
        RV = detJ * Wp * (-DN0 * je);

        // Calculate K11, K12, K21, K22
        K11 = detJ * Wp * (DN0 * dqdt - N * (djdt.t() * DN0.t() * Vee));
        K12 = detJ * Wp * (DN0 * dqdv - N * (djdv.t() * DN0.t() * Vee) - N * (je.t() * DN0));
        K21 = detJ * Wp * (DN0 * djdt);
        K22 = detJ * Wp * (DN0 * djdv);

        // Construct the elemental Jacobian KJ and Residues R
        mat KJ = join_horiz(join_vert(K11, K21), join_vert(K12, K22));
        mat R = join_vert(RT.col(0), RV.col(0));

        return std::make_tuple(KJ, R);
    };
}
*/