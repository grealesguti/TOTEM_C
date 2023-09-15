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
    thermoelectricityintegrationFunction_ = [&](const arma::mat& natcoords, const arma::mat& coords, const arma::uvec& dofs, const int elementTag) -> Utils::IntegrationResult {
        return thermoelectricityintegration(natcoords, coords, dofs,elementTag);
    };
    std::cout << "### SOLVER Initialized." << std::endl;

}

Solver::SparseSystem Solver::Assembly() {
    std::cout << "### START ASSEMBLY." << std::endl;
    int nodes_per_element = 8; // Total number of nodes per element
    int dof_per_node = 2; // Assuming 2 degrees of freedom per node
    const std::vector<std::size_t>& elementTags = mesh_.getElementTags();
    std::size_t num_of_elements = elementTags.size();

    // Initialize larger matrices for KJ and Ra
    arma::mat larger_KJ(16, num_of_elements, arma::fill::zeros);
    arma::mat larger_Ra(16, num_of_elements, arma::fill::zeros);
    arma::umat larger_element_dofs(8, num_of_elements, arma::fill::zeros);

    // Full Residual vector
    int totalmeshnodes = mesh_.getNumAllNodes();
    arma::uvec Full_Ra = arma::uvec(totalmeshnodes*dof_per_node, arma::fill::zeros);

    // Temporary matrices for the current element
    arma::mat element_Ra = arma::zeros<arma::mat>(16, 1);
    arma::mat element_KJ = arma::zeros<arma::mat>(16, 16);
    arma::uvec element_dofs = arma::uvec(16, arma::fill::zeros);
    arma::uvec element_dof_values = arma::uvec(16, arma::fill::zeros);
    arma::uvec element_zeros = arma::uvec(16, arma::fill::zeros);

        //#pragma omp parallel num_threads(num_threads)
        //{
            int elementIndex; // Declare this variable within the parallel section

            // Each thread will execute a subset of the loop iterations
        //    #pragma omp for
            for (elementIndex = 0; elementIndex < elementTags.size(); elementIndex++) {
                std::size_t elementTag = elementTags[elementIndex];

                // Rest of your loop code remains the same
                std::vector<int> nodeTags_el;
                mesh_.getElementInfo(elementTag, nodeTags_el);

                int cc = 0;
                for (int nodeTag : nodeTags_el) {
                    element_dofs[cc] = nodeTag * dof_per_node;
                    element_dofs[cc + nodes_per_element] = nodeTag * dof_per_node + 1;
                    element_dof_values[cc] = nodeTag * dof_per_node;
                    element_dof_values[cc + nodes_per_element]= nodeTag * dof_per_node + 1;
                    cc += 1;
                }
                std::cout<<"Desired Output:"<<std::endl;
                std::cout<<inputReader_.getDesiredOutput()<<std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(element_dofs,"Outputs/KTelement_dofs.txt",true);
                    utils_.writeDataToFile(element_dof_values,"Outputs/KTelement_dof_values.txt",true);
                }
                Utils::IntegrationResult elementresult;
                // if physics == 
                std::cout << "element integration "<< elementTag << std::endl;
                elementresult = utils_.gaussIntegrationK(3, 5, elementTag, mesh_, element_dof_values, thermoelectricityintegrationFunction_);
                std::cout << "element "<< elementTag<< " integrated." << std::endl;

                arma::vec vector_KJ = arma::vectorise(elementresult.KT);
                arma::vec vector_Ra = elementresult.R;
                std::cout << "vector format." << std::endl;

                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(elementresult.KT,"Outputs/KTintegration_elKT.txt",true);
                    utils_.writeDataToFile(elementresult.R,"Outputs/KTintegration_elR.txt",true);
                }
                larger_KJ.row(elementIndex) = vector_KJ.t();
                std::cout << "filled KJlarger." << std::endl;

                // Use atomic addition for updating Full_Ra
                //#pragma omp atomic
                for (int i = 0; i < element_dofs.n_elem; i++) {
                    Full_Ra[element_dofs[i]] += vector_Ra[i];
                }
                std::cout << "filled Ra." << std::endl;

                arma::uvec vector_element_dofs = element_dofs.t();
                larger_element_dofs.col(elementIndex) = vector_element_dofs;
                std::cout << "filled dofs." << std::endl;

            }
        //}

    // Create CSR sparse matrix from larger_KJ and larger_element_dofs
    arma::mat sparseMatrixData(larger_KJ.n_elem, 1);
    arma::uvec sparseMatrixRowIndices(larger_KJ.n_elem);
    arma::uvec sparseMatrixColPtrs(num_of_elements + 1);

    // Initialize variables for CSR construction
    int nnz = 0;
    for (int i = 0; i < larger_KJ.n_elem; i++) {
        for (int j = 0; j < larger_element_dofs.n_rows; j++) {
            int colIndex = larger_element_dofs(j, i);
            if (colIndex != 0) {
                sparseMatrixData(nnz) = larger_KJ(i);
                sparseMatrixRowIndices(nnz) = j;
                nnz++;
            }
        }
        sparseMatrixColPtrs(i + 1) = nnz;
    }
    if (inputReader_.getDesiredOutput()=="all"){
        std::cout<<"print all"<<std::endl;
        utils_.writeDataToFile(sparseMatrixRowIndices,"Outputs/KTsparseMatrixRowIndices.txt",true);
        utils_.writeDataToFile(sparseMatrixRowIndices,"Outputs/KTsparseMatrixRowIndices.txt",true);
        utils_.writeDataToFile(sparseMatrixColPtrs,"Outputs/KTsparseMatrixColPtrs.txt",true);
        utils_.writeDataToFile(sparseMatrixData,"Outputs/KTsparseMatrixData.txt",true);
    }
    // Create CSR sparse matrix
    arma::sp_mat KJ_sparse(sparseMatrixRowIndices, sparseMatrixColPtrs, sparseMatrixData, larger_KJ.n_rows, larger_KJ.n_cols);
    
    // reduce the system and store in the return structure
    SparseSystem result;
    // Convert std::vector<int> to arma::uvec
    arma::uvec uVector = arma::conv_to<arma::uvec>::from(freedofidxs_);
    result.KT_sparse_reduced=utils_.spmat_submat(KJ_sparse,freedofidxs_,freedofidxs_);
    for (int i = 0; i < freedofidxs_.size(); i++) {
                    result.R_reduced[i] += Full_Ra[freedofidxs_[i]];
                }
    return result;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


Utils::IntegrationResult Solver::thermoelectricityintegration(const arma::mat& natcoords, const arma::mat& coords, const arma::uvec& dofs, const int elementTag){
        Utils::IntegrationResult result; // Create a struct to hold KV and R

    // Define variables
    //std::cout << "Initialize shape functions and derivatives. " << std::endl;
    arma::vec shapeFunctions(8,1);          // Shape functions as a 4x1 vector
    mat shapeFunctionDerivatives(8, 3); // Shape function derivatives
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

    mat JM = shapeFunctionDerivatives.t() * coords.t(); // Fixed the loop indexing
                std::cout << "JM." << std::endl;
    if (inputReader_.getDesiredOutput()=="all"){
        utils_.writeDataToFile(JM,"Outputs/KTJM.txt",true);
        utils_.writeDataToFile(shapeFunctions,"Outputs/KTshapeFunctions.txt",true);
        utils_.writeDataToFile(shapeFunctionDerivatives,"Outputs/KTshapeFunctionDerivatives.txt",true);
        utils_.writeDataToFile(natcoords,"Outputs/KTnatcoords.txt",true);
        utils_.writeDataToFile(coords,"Outputs/KTcoords.txt",true);
        utils_.writeDataToFile(dofs,"Outputs/KTdofs.txt",true);

    }
    //std::cout << "Calculate jacobian determinant. " << std::endl;
    // Calculate the determinant of the Jacobian
    double detJ = arma::det(JM);

    // Extract material properties
    int materialindex = mesh_.getElementMaterial(elementTag);
    bool flag;
    std::cout << "material index "<<materialindex << std::endl;

    double De=inputReader_.getMaterialPropertyValue(materialindex,"ElectricalConductivity");
    double Da=inputReader_.getMaterialPropertyValue(materialindex,"Seebeck");
    double Dk=inputReader_.getMaterialPropertyValue(materialindex,"ThermalConductivity");

    double Dde=0, Dda=0, Ddk=0;
    std::cout << "materials " << De << " "<< Da<<" "<< Dk << std::endl;

    // Assuming 'dofs' is an Armadillo vector
    arma::mat Tee(8,1); // Vector for odd-indexed elements
    arma::mat Vee(8,1); // Vector for even-indexed elements
    // Extract odd-indexed elements into Tee
    for (uword i = 1; i < dofs.n_elem; i += 2) {
        Tee((i-1) / 2,1) = dofs(i);
    }
    std::cout << "Tee" << std::endl;

    // Extract even-indexed elements into Vee
    for (uword i = 0; i < dofs.n_elem; i += 2) {
        Vee(i / 2,1) = dofs(i);
    }
    std::cout << "Vee" << std::endl;

    double Th = arma::dot(shapeFunctions, dofs);
    std::cout << "Th" << std::endl;

    // Calculate je and qe
    arma::mat je = -De * shapeFunctionDerivatives.t() * Vee - Da * De * shapeFunctionDerivatives.t() * Tee;
    arma::mat qe = Da * (shapeFunctions.t() * Tee) * je - Dk * shapeFunctionDerivatives.t() * Tee;
    std::cout << "qe je." << std::endl;

    // Perform the necessary calculations here to compute djdt, djdv, dqdt, dqdv
    arma::mat djdt = -Da * De * shapeFunctionDerivatives.t()
                    - Dda * De * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t()
                    - Dde * (shapeFunctionDerivatives.t() * Vee.t() + Da * shapeFunctionDerivatives.t() * Tee.t()) * shapeFunctions.t();
    std::cout << "djdt" << std::endl;

    arma::mat djdv = -De * shapeFunctionDerivatives;
    std::cout << "djdv" << std::endl;

    arma::mat dqdt = Da * Th * djdt + Da * je * shapeFunctions.t()
                    - Dk * shapeFunctionDerivatives.t()
                    + Dda * Th * je * shapeFunctions.t()
                    - Ddk * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t();
    std::cout << "dqdt" << std::endl;

    arma::mat dqdv = -Da * De * Th * shapeFunctionDerivatives;
    std::cout << "dqdv" << std::endl;

    // Perform the integration for RT, RV, K11, K12, K21, K22
    RT += (-shapeFunctionDerivatives * qe + shapeFunctions * je.t() * shapeFunctionDerivatives.t() * Vee);
    RV += (-shapeFunctionDerivatives * je);
    K11 +=(shapeFunctionDerivatives * dqdt - shapeFunctions * (djdt.t() * shapeFunctionDerivatives.t() * Vee).t());
    K12 +=(shapeFunctionDerivatives * dqdv - shapeFunctions * (djdv.t() * shapeFunctionDerivatives.t() * Vee).t() - shapeFunctions * (je.t() * shapeFunctionDerivatives));
    K21 +=(shapeFunctionDerivatives * djdt);
    K22 +=(shapeFunctionDerivatives * djdv);
    std::cout << "matrixes" << std::endl;


    //std::cout << "Calculate integrand. " << std::endl;
    // Define the size of KV and RT matrices based on the sizes of K11, K12, K21, K22, RT, and RV
    int numRowsKJ = K11.n_rows + K21.n_rows;
    int numColsKJ = K11.n_cols + K12.n_cols;
    int numRowsR = RT.n_rows + RV.n_rows;

    // Initialize KV and RT matrices with zeros
    arma::mat KV(numRowsKJ, numColsKJ, arma::fill::zeros);
    arma::mat R(numRowsR, 1, arma::fill::zeros);
    std::cout << "matrixes larger" << std::endl;

    // Fill in the KV matrix with K11, K12, K21, and K22
    KV.submat(0, 0, K11.n_rows - 1, K11.n_cols - 1) = K11;
    KV.submat(0, K11.n_cols, K11.n_rows - 1, numColsKJ - 1) = K12;
    KV.submat(K11.n_rows, 0, numRowsKJ - 1, K11.n_cols - 1) = K21;
    KV.submat(K11.n_rows, K11.n_cols, numRowsKJ - 1, numColsKJ - 1) = K22;

    // Fill in the R matrix with RT and RV
    R.submat(0, 0, RT.n_rows - 1, 0) = RT;
    R.submat(RT.n_rows, 0, numRowsR - 1, 0) = RV;
    std::cout << "matrixes larger filled" << std::endl;

    //std::cout << "integrand calculated " << std::endl;
    result.R=R*detJ;
    result.KT=KV*detJ;
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