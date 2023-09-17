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
    thermoelectricityintegrationFunction_ = [&](const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag) -> Utils::IntegrationResult {
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
    arma::mat larger_KJ(16*16, num_of_elements, arma::fill::zeros);
    arma::mat larger_Ra(16, num_of_elements, arma::fill::zeros);
    arma::umat larger_element_dofs(16, num_of_elements, arma::fill::zeros);

    // Full Residual vector
    int totalmeshnodes = mesh_.getNumAllNodes();
    arma::vec Full_Ra = arma::vec(totalmeshnodes*dof_per_node, arma::fill::zeros);

    // Temporary matrices for the current element
    arma::mat element_Ra = arma::zeros<arma::mat>(16, 1);
    arma::mat element_KJ = arma::zeros<arma::mat>(16, 16);
    arma::uvec element_dofs = arma::uvec(16, arma::fill::zeros);
    arma::vec element_dof_values = arma::vec(16, arma::fill::zeros);
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
                    element_dof_values[cc] = bcinit_.getInitialDof( nodeTag * dof_per_node);
                    element_dof_values[cc + nodes_per_element]= bcinit_.getInitialDof( nodeTag * dof_per_node + 1);
                    cc += 1;
                }
                std::cout<<"Desired Output:"<<std::endl;
                std::cout<<inputReader_.getDesiredOutput()<<std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(nodeTags_el,"Outputs/KTnodeTags_el.txt",true);
                    utils_.writeDataToFile(element_dofs,"Outputs/KTelement_dofs.txt",true);
                    utils_.writeDataToFile(element_dof_values,"Outputs/KTelement_dof_values.txt",true);
                }
                // if physics == 
                std::cout << "element integration "<< elementTag << std::endl;
                Utils::IntegrationResult elementresult = utils_.gaussIntegrationK(3, 3, elementTag, mesh_, element_dof_values, thermoelectricityintegrationFunction_);
                std::cout << "element "<< elementTag<< " integrated." << std::endl;

                arma::vec vector_KJ = arma::vectorise(elementresult.KT);
                arma::vec vector_Ra = elementresult.R;
                std::cout << "vector format." << std::endl;

                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(vector_KJ,"Outputs/KTintegration_elKT.txt",true);
                    utils_.writeDataToFile(elementresult.R,"Outputs/KTintegration_elR.txt",true);
                }
                larger_KJ.col(elementIndex) = vector_KJ;
                std::cout << "filled KJlarger." << std::endl;

                // Use atomic addition for updating Full_Ra
                //#pragma omp atomic
                for (int i = 0; i < element_dofs.n_elem; i++) {
                    Full_Ra[element_dofs[i]] += vector_Ra[i];
                }
                std::cout << "filled Ra." << std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(larger_KJ,"Outputs/KTintegration_elKTlarger.txt",true);
                    utils_.writeDataToFile(Full_Ra,"Outputs/KTintegration_elRa.txt",true);
                }
                larger_element_dofs.col(elementIndex) = element_dofs;
                std::cout << "filled dofs." << std::endl;

            }
              std::cout << "finished loop." << std::endl;

        //}
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(larger_element_dofs,"Outputs/KTintegration_larger_element_dofs.txt",true);
                    utils_.writeDataToFile(larger_KJ,"Outputs/KTintegration_elKTlargerFinal.txt",true);
                    utils_.writeDataToFile(Full_Ra,"Outputs/KTintegration_elRaFinal.txt",true);
                }
    // Create a vector to store Eigen triplets
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(16 * 16 * num_of_elements); // Reserve space to avoid reallocation

    // Iterate over the 2x2 matrices and add their entries to the Triplets vector
        for (int i = 0; i < num_of_elements; i++) {
            for (int row = 0; row < 16; row++) {
                for (int col = 0; col < 16; col++) {
                        int dof_row = larger_element_dofs(row, i);
                        int dof_col = larger_element_dofs(col, i);
                    try {
                        // Check bounds before accessing elements
                        triplets.emplace_back(dof_row-2, dof_col-2, larger_KJ( row * 16 + col,i)); // -1 due to start index of 0

                    } catch (const std::out_of_range& e) {
                        // Handle the out-of-range exception here (e.g., print an error message).
                        std::cerr << "Error: " << e.what() << " at i=" << i << ", row=" << row << ", col=" << col << std::endl;
                        // You can also choose to do something else, like logging or terminating the program.
                    }
                }
            }
        }

                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(triplets,"Outputs/KTintegration_triplets.txt",true);
                }
  std::cout << "triplets made." << std::endl;

  // Create the sparse matrix and set its values from the Triplets vector
  Eigen::SparseMatrix<double> KsparseMatrix(mesh_.getNumAllNodes()*2, mesh_.getNumAllNodes()*2);
  std::cout << "init sparse. " <<mesh_.getNumAllNodes()*2 << std::endl;
  KsparseMatrix.setFromTriplets(triplets.begin(), triplets.end());
  std::cout << "made Sparse from eigen." << std::endl;
    // Create CSR sparse matrix

    
    // reduce the system and store in the return structure
    SparseSystem result;
    // Convert std::vector  <int> to arma::uvec
    result.KsubMatrix =  Solver::reduceSystem(KsparseMatrix);
    std::cout << " K Submatrix retrieved." << std::endl;
    // Resize the vector to the desired size and initialize with zeros
    result.R_reduced.resize(freedofidxs_.size(), 0.0);
    for (int i = 0; i < freedofidxs_.size(); i++) {
                    result.R_reduced[i] += Full_Ra[freedofidxs_[i]];
    }

    std::cout << " R Submatrix retrieved." << std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(result.KsubMatrix,"Outputs/KTintegration_reducedK.txt",true);
                    utils_.writeDataToFile(result.R_reduced,"Outputs/KTintegration_reducedRa.txt",true);
                }
    return result;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


Utils::IntegrationResult Solver::thermoelectricityintegration(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag){
        Utils::IntegrationResult result; // Create a struct to hold KV and R

    // Define variables
    //std::cout << "Initialize shape functions and derivatives. " << std::endl;
    arma::vec shapeFunctions(8,1);          // Shape functions as a 4x1 vector
    arma::mat shapeFunctionDerivatives(8, 3); // Shape function derivatives
    // Define the integration result matrices
    arma::mat RT(8, 1, arma::fill::zeros);
    arma::mat RV(8, 1, arma::fill::zeros);
    arma::mat K11(8, 8, arma::fill::zeros);
    arma::mat K12(8, 8, arma::fill::zeros);
    arma::mat K21(8, 8, arma::fill::zeros);
    arma::mat K22(8, 8, arma::fill::zeros);

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

    arma::mat JM = shapeFunctionDerivatives.t() * coords.t(); // Fixed the loop indexing
                //std::cout << "JM." << std::endl;
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
    //std::cout <<"detJ"<< detJ << std::endl;

    // Extract material properties
    int materialindex = mesh_.getElementMaterial(elementTag);
    bool flag;
    //std::cout << "material index "<<materialindex << std::endl;

    double De=inputReader_.getMaterialPropertyValue(materialindex,"ElectricalConductivity");
    double Da=inputReader_.getMaterialPropertyValue(materialindex,"Seebeck");
    double Dk=inputReader_.getMaterialPropertyValue(materialindex,"ThermalConductivity");

    double Dde=0, Dda=0, Ddk=0;
    //std::cout << "materials " << De << " "<< Da<<" "<< Dk << std::endl;

    // Assuming 'dofs' is an Armadillo vector
    arma::mat Tee(8,1); // Vector for odd-indexed elements
    arma::mat Vee(8,1); // Vector for even-indexed elements
    // Extract odd-indexed elements into Tee
    for (uword i = 0; i < 8; i ++) {
        Tee(i) = dofs(i);
        Vee(i) = dofs(i+8);
    }
    //std::cout << "Tee,Vee" << std::endl;
    //std::cout << "Vee" << std::endl;
    // Check dimensions before dot product

    arma::mat Thmat = shapeFunctions.t()*Tee;
    double Th=Thmat(0,0);
    //std::cout << "Th: " << Th << std::endl;

    // Calculate je and qe
    arma::mat je = -De * shapeFunctionDerivatives.t() * Vee - Da * De * shapeFunctionDerivatives.t() * Tee;
    //std::cout << " je." << std::endl;
    arma::mat qe = Da * Th * je - Dk * shapeFunctionDerivatives.t() * Tee;
    //std::cout << "qe ." << std::endl;

    // Perform the necessary calculations here to compute djdt, djdv, dqdt, dqdv
    arma::mat djdt = -Da * De * shapeFunctionDerivatives.t()
                    - Dda * De * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t()
                    - Dde * (shapeFunctionDerivatives.t() * Vee + Da * shapeFunctionDerivatives.t() * Tee) * shapeFunctions.t();
    //std::cout << "djdt" << std::endl;

    arma::mat djdv = -De * shapeFunctionDerivatives.t();
    //std::cout << "djdv" << std::endl;

    arma::mat dqdt = Da * Th * djdt + Da * je * shapeFunctions.t()
                    - Dk * shapeFunctionDerivatives.t()
                    + Dda * Th * je * shapeFunctions.t()
                    - Ddk * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t();
    //std::cout << "dqdt" << std::endl;

    arma::mat dqdv = -Da * De * Th * shapeFunctionDerivatives.t();
    //std::cout << "dqdv" << std::endl;

    // Perform the integration for RT, RV, K11, K12, K21, K22
    RT += (-shapeFunctionDerivatives * qe + shapeFunctions * je.t() * shapeFunctionDerivatives.t() * Vee);
    //std::cout << "RT" << std::endl;
    RV += (-shapeFunctionDerivatives * je);
    //std::cout << "RV" << std::endl;
    K11 +=(shapeFunctionDerivatives * dqdt - shapeFunctions * (djdt.t() * shapeFunctionDerivatives.t() * Vee).t());
    //std::cout << "K11" << std::endl;
    K12 +=(shapeFunctionDerivatives * dqdv - shapeFunctions * (djdv.t() * shapeFunctionDerivatives.t() * Vee).t() - shapeFunctions * (je.t() * shapeFunctionDerivatives.t()));
    //std::cout << "K12" << std::endl;
    K21 +=(shapeFunctionDerivatives * djdt);
    //std::cout << "K21" << std::endl;
    K22 +=(shapeFunctionDerivatives * djdv);
    //std::cout << "matrixes" << std::endl;

    if (inputReader_.getDesiredOutput()=="all"){
        utils_.writeDataToFile(je,"Outputs/KTqe.txt",true);
        utils_.writeDataToFile(qe,"Outputs/KTje.txt",true);
        utils_.writeDataToFile(djdt,"Outputs/KTdjdt.txt",true);
        utils_.writeDataToFile(djdv,"Outputs/KTdjdv.txt",true);
        utils_.writeDataToFile(dqdt,"Outputs/KTdqdt.txt",true);
        utils_.writeDataToFile(dqdv,"Outputs/KTdqdv.txt",true);
        utils_.writeDataToFile(RT,"Outputs/KTRT.txt",true);
        utils_.writeDataToFile(RV,"Outputs/KTRV.txt",true);

        utils_.writeDataToFile(K11,"Outputs/KTK11.txt",true);
        utils_.writeDataToFile(K12,"Outputs/KTK12.txt",true);
        utils_.writeDataToFile(K21,"Outputs/KTK21.txt",true);
        utils_.writeDataToFile(K22,"Outputs/KTK22.txt",true);

    }
    //std::cout << "Calculate integrand. " << std::endl;
    // Define the size of KV and RT matrices based on the sizes of K11, K12, K21, K22, RT, and RV
    int numRowsKJ = K11.n_rows + K21.n_rows;
    int numColsKJ = K11.n_cols + K12.n_cols;
    int numRowsR = RT.n_rows + RV.n_rows;

    // Initialize KV and RT matrices with zeros
    arma::mat KV(numRowsKJ, numColsKJ, arma::fill::zeros);
    arma::mat R(numRowsR, 1, arma::fill::zeros);
    //std::cout << "matrixes larger" << std::endl;

    // Fill in the KV matrix with K11, K12, K21, and K22
    KV.submat(0, 0, K11.n_rows - 1, K11.n_cols - 1) = K11;
    KV.submat(0, K11.n_cols, K11.n_rows - 1, numColsKJ - 1) = K12;
    KV.submat(K11.n_rows, 0, numRowsKJ - 1, K11.n_cols - 1) = K21;
    KV.submat(K11.n_rows, K11.n_cols, numRowsKJ - 1, numColsKJ - 1) = K22;

    // Fill in the R matrix with RT and RV
    R.submat(0, 0, RT.n_rows - 1, 0) = RT;
    R.submat(RT.n_rows, 0, numRowsR - 1, 0) = RV;
    //std::cout << "matrixes larger filled" << std::endl;

    //std::cout << "integrand calculated " << std::endl;
    result.R=R*detJ;
    result.KT=KV*detJ;
    if (inputReader_.getDesiredOutput()=="all"){
        utils_.writeDataToFile(result.KT,"Outputs/KTel_rstKT.txt",true);
        utils_.writeDataToFile(result.R,"Outputs/KTel_rstR.txt",true);
        utils_.writeDataToFile(KV,"Outputs/KTel_KT.txt",true);
        utils_.writeDataToFile(R,"Outputs/KTel_R.txt",true);
    }
    // Return the heat flow as a 4x1 Armadillo matrix
    return result;
}
///////////////////////////////////////////////////////////////////

Eigen::SparseMatrix<double> Solver::reduceSystem(const Eigen::SparseMatrix<double>& K) {
            std::cout << "Reducing system " << std::endl;

    int numDofs = K.rows();
    int numFixedDofs = freedofidxs_.size(); // Assuming freedofidx_ is a private member variable
            std::cout << "freedofs size:  " << numFixedDofs << std::endl;
            std::cout << "Krows  " << numDofs << std::endl;

    // Create the mapping from original degrees of freedom to reduced degrees of freedom
    std::vector<int> dofMap(numDofs, -1);
    for (int i = 0; i < numFixedDofs; i++) {
            std::cout << freedofidxs_[i] << " "<<i<< std::endl;
        dofMap[freedofidxs_[i]] = i;
    }
            std::cout << "dofmap " << std::endl;

    // Create the reduced sparse matrix
    Eigen::SparseMatrix<double> reducedK(numFixedDofs, numFixedDofs);
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            if (dofMap[row] != -1 && dofMap[col] != -1) {
                reducedK.insert(dofMap[row], dofMap[col]) = it.value();
            }
        }
    }
            std::cout << "loop " << std::endl;

    reducedK.makeCompressed();

    return reducedK;
}
/////////////////////////////////////////////////////////////
Eigen::VectorXd Solver::solveSparseSystem(const SparseSystem& system) {
    // Extract the submatrix and reduced right-hand side
    Eigen::SparseMatrix<double> Ksub = system.KsubMatrix;
    Eigen::VectorXd Rreduced(system.R_reduced.size());
    for (size_t i = 0; i < system.R_reduced.size(); ++i) {
        Rreduced(i) = system.R_reduced[i];
    }

    // Create a sparse LU solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(Ksub);

    if (solver.info() != Eigen::Success) {
        // Decomposition failed
        std::cout << "Warning: Sparse LU decomposition failed." << std::endl;
        // Handle the error here, possibly by returning an error code or throwing an exception
        // Example: throw std::runtime_error("Sparse LU decomposition failed.");
    }

    // Solve the system
    Eigen::VectorXd solution = solver.solve(Rreduced);

    if (solver.info() != Eigen::Success) {
        // Solve failed
        std::cout << "Warning: Solving the sparse system failed." << std::endl;
        // Handle the error here, possibly by returning an error code or throwing an exception
        // Example: throw std::runtime_error("Solving the sparse system failed.");
    }

    return solution;
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