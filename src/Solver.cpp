#include "Solver.hpp"
using namespace arma;

Solver::Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit)
    : inputReader_(inputReader), mesh_(mesh), elements_(), bcinit_(bcinit), utils_() {
    // Get the number of nodes from the mesh
    int numNodes = mesh_.getNumAllNodes();

    // Initialize the loadVector_ member with a size double the number of nodes
    loadVector_= bcinit_.getloadVector();
    soldofs_.resize(2 * numNodes);
    soldofs_= bcinit_.getAllInitialDof();
    freedofidxs_ = mesh_.GetFreedofsIdx();

    // Initialize the thermoelectricityintegrationFunction_ using the assignment operator
    thermoelectricityintegrationFunction_ = [&](const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag) -> Utils::IntegrationResult {
        return thermoelectricityintegration(natcoords, coords, dofs,elementTag);
    };
    std::cout << "### SOLVER Initialized." << std::endl;

}

    void Solver::Assembly(Eigen::SparseMatrix<double> KsubMatrix, std::vector<double> R_reduced) {
    std::cout << "### START ASSEMBLY." << std::endl;
    int dof_per_node;
    if(inputReader_.getPhysics()=="thermoelectricity"){
        dof_per_node = 2;
    }else{
        dof_per_node = 2;// by default
    } 
    const std::vector<std::size_t> elementTags = mesh_.getElementTags();
    std::size_t num_of_elements = elementTags.size();

    std::pair<int, int> nodesperelement_etype = mesh_.getNumNodesForElement(elementTags[0]);
    int nodes_per_element = nodesperelement_etype.first; // Number of nodes
    int etype = nodesperelement_etype.second; // Element type    
    int num_dofs_per_element = nodes_per_element*dof_per_node;
    //std::cout<<"element type: "<< etype<<"nodes per element "<<nodes_per_element<<" num_dofs_per_element "<<num_dofs_per_element<<std::endl;

    // Initialize larger matrices for KJ and Ra
    arma::mat larger_KJ(num_dofs_per_element*num_dofs_per_element, num_of_elements, arma::fill::zeros);
    arma::mat larger_Ra(num_dofs_per_element, num_of_elements, arma::fill::zeros);
    arma::umat larger_element_dofs(num_dofs_per_element, num_of_elements, arma::fill::zeros);

    // Full Residual vector
    int totalmeshnodes = mesh_.getNumAllNodes();
    arma::vec Full_Ra = arma::vec(totalmeshnodes*dof_per_node, arma::fill::zeros);

    // Temporary matrices for the current element
    arma::mat element_Ra = arma::zeros<arma::mat>(num_dofs_per_element, 1);
    arma::mat element_KJ = arma::zeros<arma::mat>(num_dofs_per_element, num_dofs_per_element);
    arma::uvec element_dofs = arma::uvec(num_dofs_per_element, arma::fill::zeros);
    arma::vec element_dof_values = arma::vec(num_dofs_per_element, arma::fill::zeros);
    arma::uvec element_zeros = arma::uvec(num_dofs_per_element, arma::fill::zeros);

        //#pragma omp parallel num_threads(num_threads)
        //{
            int elementIndex; // Declare this variable within the parallel section

            // Each thread will execute a subset of the loop iterations
        //    #pragma omp for
            for (elementIndex = 0; elementIndex < elementTags.size(); elementIndex++) {
                std::size_t elementTag = elementTags[elementIndex];

                // Rest of your loop code remains the same
                std::vector<int> nodeTags_el;
                int etype = mesh_.getElementInfo(elementTag, nodeTags_el);
              //  std::cout<<"element type: "<< etype<<"nodes per element "<<nodes_per_element<<std::endl;
                std::pair<int, int> nodesperelement_etype = mesh_.getNumNodesForElement(elementTags[elementIndex]);
                int nodes_per_element = nodesperelement_etype.first; // Number of nodes
                int etype1 = nodesperelement_etype.second; // Element type    
                int num_dofs_per_element = nodes_per_element*dof_per_node;
               // std::cout<<"element type: "<< etype1<<"nodes per element "<<nodes_per_element<<" num_dofs_per_element "<<num_dofs_per_element<<std::endl;

                int cc = 0;
                for (int nodeTag : nodeTags_el) {
                    element_dofs[cc] = nodeTag * dof_per_node;
                    element_dofs[cc + nodes_per_element] = nodeTag * dof_per_node + 1;
                    element_dof_values[cc] = soldofs_[ nodeTag * dof_per_node];
                    element_dof_values[cc + nodes_per_element]= soldofs_[ nodeTag * dof_per_node + 1];
                    cc += 1;
                }
               // std::cout<<"-Desired Output:"<<std::endl;
                //std::cout<<inputReader_.getDesiredOutput()<<std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(nodeTags_el,"Outputs/KTnodeTags_el.txt",true);
                    utils_.writeDataToFile(element_dofs,"Outputs/KTelement_dofs.txt",true);
                    utils_.writeDataToFile(element_dof_values,"Outputs/KTelement_dof_values.txt",true);
                }
                // if physics == 
                std::cout << "-element integration "<< elementTag << std::endl;
                Utils::IntegrationResult elementresult = utils_.gaussIntegrationK(3, 3, elementTag, mesh_, element_dof_values, thermoelectricityintegrationFunction_);
                //std::cout << "-element "<< elementTag<< " integrated." << std::endl;

                arma::vec vector_KJ = arma::vectorise(elementresult.KT);
                arma::vec vector_Ra = elementresult.R;
               // std::cout << "-vector format." << std::endl;

                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(elementresult.KT,"Outputs/KTintegration_elKTnovector.txt",true);
                    utils_.writeDataToFile(vector_KJ,"Outputs/KTintegration_elKT.txt",true);
                    utils_.writeDataToFile(elementresult.R,"Outputs/KTintegration_elR.txt",true);
                }
                larger_KJ.col(elementIndex) = vector_KJ;
                //std::cout << "-filled KJlarger." << std::endl;

                // Use atomic addition for updating Full_Ra
                //#pragma omp atomic
                for (int i = 0; i < element_dofs.n_elem; i++) {
                    // FIXME (Parche temporal) Find out why some values are out of range and use methods like insert, at, an so on
                    if (element_dofs[i]<Full_Ra.size())
                        Full_Ra[element_dofs[i]] += vector_Ra[i];
                }
                //std::cout << "-filled Ra." << std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(larger_KJ,"Outputs/KTintegration_elKTlarger.txt",true);
                    utils_.writeDataToFile(Full_Ra,"Outputs/KTintegration_elRa.txt",true);
                }
                larger_element_dofs.col(elementIndex) = element_dofs;
               // std::cout << "-filled dofs." << std::endl;

            }
              //std::cout << "-finished loop." << std::endl;

        //}
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(larger_element_dofs,"Outputs/KTintegration_larger_element_dofs.txt",true);
                    utils_.writeDataToFile(larger_KJ,"Outputs/KTintegration_elKTlargerFinal.txt",true);
                    utils_.writeDataToFile(Full_Ra,"Outputs/KTintegration_elRaFinal.txt",true);
                }
    // Create a vector to store Eigen triplets
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(num_dofs_per_element * num_dofs_per_element * num_of_elements); // Reserve space to avoid reallocation

    // Iterate over the 2x2 matrices and add their entries to the Triplets vector
        for (int i = 0; i < num_of_elements; i++) {
            for (int row = 0; row < num_dofs_per_element; row++) {
                for (int col = 0; col < num_dofs_per_element; col++) {
                        int dof_row = larger_element_dofs(row, i);
                        int dof_col = larger_element_dofs(col, i);
                    try {
                        // Check bounds before accessing elements
                        triplets.emplace_back(dof_row-2, dof_col-2, larger_KJ( row * num_dofs_per_element + col,i)); // -1 due to start index of 0

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
 // std::cout << "-triplets made." << std::endl;

  // Create the sparse matrix and set its values from the Triplets vector
  Eigen::SparseMatrix<double> KsparseMatrix(mesh_.getNumAllNodes()*2, mesh_.getNumAllNodes()*2);
  //std::cout << "-init sparse. " <<mesh_.getNumAllNodes()*2 << std::endl;
  KsparseMatrix.setFromTriplets(triplets.begin(), triplets.end());
  //std::cout << "-made Sparse from eigen." << std::endl;
    
    // reduce the system and store in the return structure
    // Convert std::vector  <int> to arma::uvec
    Solver::reduceSystem(KsparseMatrix,KsubMatrix);
    //std::cout << " -K Submatrix retrieved." << std::endl;
    // Resize the vector to the desired size and initialize with zeros
    R_reduced.resize(freedofidxs_.size(), 0.0);
    for (int i = 0; i < freedofidxs_.size(); i++) {
                if (abs(Full_Ra[freedofidxs_[i]])>1e-20){
                    R_reduced[i] = Full_Ra[freedofidxs_[i]];
                }
    }           

    //std::cout << " -R Submatrix retrieved." << std::endl;
                if (inputReader_.getDesiredOutput()=="all"){
                    utils_.writeDataToFile(KsubMatrix,"Outputs/KTintegration_reducedK.txt",true);
                    utils_.writeDataToFile(R_reduced,"Outputs/KTintegration_reducedRa.txt",true);
                }
    std::cout << " assembly_ return." << std::endl;
    return;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


Utils::IntegrationResult Solver::thermoelectricityintegration(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag){
        Utils::IntegrationResult result; // Create a struct to hold KV and R

    std::pair<int, int> nodesperelement_etype = mesh_.getNumNodesForElement(elementTag);
    int nodes_per_element = nodesperelement_etype.first; // Number of nodes
    int etype = nodesperelement_etype.second; // Element type    
    //std::cout<< "nodes per element "<<nodes_per_element<<std::endl;
    // Define variables
    //std::cout << "Initialize shape functions and derivatives. " << std::endl;
    arma::vec shapeFunctions(nodes_per_element,1);          // Shape functions as a 4x1 vector
    arma::mat shapeFunctionDerivatives(nodes_per_element, 3); // Shape function derivatives
    // Define the integration result matrices
    arma::mat RT(nodes_per_element, 1, arma::fill::zeros);
    arma::mat RV(nodes_per_element, 1, arma::fill::zeros);
    arma::mat K11(nodes_per_element, nodes_per_element, arma::fill::zeros);
    arma::mat K12(nodes_per_element, nodes_per_element, arma::fill::zeros);
    arma::mat K21(nodes_per_element, nodes_per_element, arma::fill::zeros);
    arma::mat K22(nodes_per_element, nodes_per_element, arma::fill::zeros);

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
    //std::cout << "Get shape functions and derivatives. " << etype << std::endl;
    mesh_.selectShapeFunctionsAndDerivatives(etype,xi,eta,zeta,shapeFunctions,shapeFunctionDerivatives);
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
    //std::cout << "material index "<<materialindex << std::endl;

    double De=inputReader_.getMaterialPropertyValue(materialindex,"ElectricalConductivity");
    double Da=inputReader_.getMaterialPropertyValue(materialindex,"Seebeck");
    double Dk=inputReader_.getMaterialPropertyValue(materialindex,"ThermalConductivity");

    double Dde=0, Dda=0, Ddk=0;
    //std::cout << "materials " << De << " "<< Da<<" "<< Dk << std::endl;

    // Assuming 'dofs' is an Armadillo vector
    arma::mat Tee(nodes_per_element,1); // Vector for odd-indexed elements
    arma::mat Vee(nodes_per_element,1); // Vector for even-indexed elements
    // Extract odd-indexed elements into Tee
    for (uword i = 0; i < nodes_per_element; i ++) {
        Tee(i) = dofs(i);
        Vee(i) = dofs(i+nodes_per_element);
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
   // std::cout << "djdv" << std::endl;

    arma::mat dqdt = Da * Th * djdt + Da * je * shapeFunctions.t()
                    - Dk * shapeFunctionDerivatives.t()
                    + Dda * Th * je * shapeFunctions.t()
                    - Ddk * shapeFunctionDerivatives.t() * Tee * shapeFunctions.t();
    //std::cout << "dqdt" << std::endl;

    arma::mat dqdv = -Da * De * Th * shapeFunctionDerivatives.t();
    //std::cout << "dqdv" << std::endl;
    std::cout << "thermoel_matrix addition" << std::endl;

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
    std::cout << "thermoel_matrixes larger" << std::endl;

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
    //std::cout << "matrixes larger filled" << std::endl;

    result.R=R*detJ;
    result.KT=KV*detJ;
    //std::cout << "integrand calculated " << std::endl;

    if (inputReader_.getDesiredOutput()=="all"){
        utils_.writeDataToFile(result.KT,"Outputs/KTel_rstKT.txt",true);
        utils_.writeDataToFile(result.R,"Outputs/KTel_rstR.txt",true);
        utils_.writeDataToFile(KV,"Outputs/KTel_KT.txt",true);
        utils_.writeDataToFile(R,"Outputs/KTel_R.txt",true);
    }
    // Return the heat flow as a 4x1 Armadillo matrix
    std::cout << "thermoel_return" << std::endl;
    return result;
}
///////////////////////////////////////////////////////////////////

void Solver::reduceSystem(const Eigen::SparseMatrix<double>& K, Eigen::SparseMatrix<double>& reducedK) {
     //       std::cout << "-Reducing system " << std::endl;

    int numDofs = K.rows();
    int numFixedDofs = freedofidxs_.size(); // Assuming freedofidx_ is a private member variable
    //       std::cout << "-freedofs size:  " << numFixedDofs << std::endl;
     //       std::cout << "-Krows  " << numDofs << std::endl;

    // Create the mapping from original degrees of freedom to reduced degrees of freedom
    std::vector<int> dofMap(numDofs, -1);
    for (int i = 0; i < numFixedDofs; i++) {
        dofMap[freedofidxs_[i]] = i;
    }
    //        std::cout << "-dofmap " << std::endl;

    // Create the reduced sparse matrix
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            if (dofMap[row] != -1 && dofMap[col] != -1) {
                if(abs(it.value())>1e-20){
                    reducedK.insert(dofMap[row], dofMap[col]) = it.value();
                }
            }
        }
    }
    //        std::cout << "-loop " << std::endl;

    reducedK.makeCompressed();
}
/////////////////////////////////////////////////////////////
void Solver::solveSparseSystem(
    const Eigen::SparseMatrix<double, 0, int>& KsubMatrix,
    const std::vector<double>& R_reduced,
    Eigen::VectorXd& solution){    // Extract the submatrix and reduced right-hand side
        std::cout << "Solving system" << std::endl;

    Eigen::SparseMatrix<double> Ksub = KsubMatrix;
    Eigen::VectorXd Rreduced(R_reduced.size());
    for (size_t i = 0; i < R_reduced.size(); ++i) {
        Rreduced(i) = R_reduced[i];
    }

    // Set the desired number of threads for Eigen
    Eigen::setNbThreads(4);

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
    solution = solver.solve(Rreduced);

    if (solver.info() != Eigen::Success) {
        // Solve failed
        std::cout << "Warning: Solving the sparse system failed." << std::endl;
        // Handle the error here, possibly by returning an error code or throwing an exception
        // Example: throw std::runtime_error("Solving the sparse system failed.");
    }

    // Reset the number of threads to its default value (optional)
    Eigen::setNbThreads(0);

}
////////////////////////////////////////////////////////////////////////////////

/*
// make another function only as solver that takes the assembly and can choose btw modified NR, or normal and returns successive errors etc...
// This is the function to call from main after initialization of the class. This function also identifies the physics to use
*/
/*
double Solver::runNewtonRaphson() {
    // Initialize some parameters and initial guess
    double tolerance = 1e-6;
    int maxIterations = 5;
    if (inputReader_.getDesiredOutput()=="all"){
            utils_.writeDataToFile(soldofs_,"Outputs/NR_soldofs_"+std::to_string(0)+".txt",true);
    }
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Assembly the system matrix
        SparseSystem system = Assembly();

        // Solve the system using the tangential matrix and residual: KT*dU=R->dU
        Eigen::VectorXd delta_degreesoffreedom = solveSparseSystem(system);

        // Update solution
        for(int i=0; i<freedofidxs_.size();i++){
            soldofs_[freedofidxs_[i]]+=delta_degreesoffreedom[i];
        }

        // Calculate the residual (error);
        Eigen::VectorXd eigenR_reduced = Eigen::Map<Eigen::VectorXd>(system.R_reduced.data(), system.R_reduced.size());        // Convert std::vector to Eigen::VectorXd
        double residual = eigenR_reduced.norm();// Calculate the norm
        std::cout<< "### NR. Iteration "<< iter <<" residual "<< residual<< std::endl;
        // Check for convergence
        if (inputReader_.getDesiredOutput()=="all"){
            utils_.writeDataToFile(soldofs_,"Outputs/NR_soldofs_"+std::to_string(iter+1)+".txt",true);
            utils_.writeDataToFile(delta_degreesoffreedom,"Outputs/NR_delta_degreesoffreedom_"+std::to_string(iter)+".txt",true);
        }
        if (residual < tolerance) {
            // Converged, return both the solution and the final residual
            std::cout<< "### NR. CONVERGED "<< std::endl;
            return residual;
        }

    }

    // If we reach here, the Newton-Raphson method did not converge
    throw std::runtime_error("Newton-Raphson did not converge.");
}*/
////////////////////////////////////////////////////////////////////////////////
     /*
double Solver::runArcLengthSolver() {
    // Initialize some parameters
    double tolerance = 1e-6;
    int maxIterations = 5;
    double lambda = 0.0; // Initial guess for arc-length parameter
    double arcLengthTolerance = 1e-6; // Tolerance for arc-length convergence
    double maxArcLengthIncrement = 0.1; // Maximum arc-length increment

    if (inputReader_.getDesiredOutput() == "all") {
        utils_.writeDataToFile(soldofs_, "Outputs/AL_soldofs_" + std::to_string(0) + ".txt", true);
    }

    for (int iter = 0; iter < maxIterations; ++iter) {
        // Assembly the system matrix
        SparseSystem system = Assembly();
       
        // Modify the right-hand side (load vector) to include the arc-length constraint
        Eigen::VectorXd R_total = system.R_reduced;
        R_total -= lambda * delta_degreesoffreedom;

        // Solve the system using the modified load vector: KT*dU = R_total
        Eigen::VectorXd delta_degreesoffreedom = solveSparseSystem(system, R_total);

        // Update solution
        for (int i = 0; i < freedofidxs_.size(); i++) {
            soldofs_[freedofidxs_[i]] += delta_degreesoffreedom[i];
        }

        // Calculate the normal residual
        Eigen::VectorXd eigenR_reduced = Eigen::Map<Eigen::VectorXd>(system.R_reduced.data(), system.R_reduced.size());
        double residual = eigenR_reduced.norm();

        // Calculate the arc-length residual
        double arcLengthResidual = delta_degreesoffreedom.norm() - maxArcLengthIncrement;

        std::cout << "### AL. Iteration " << iter << " normal residual " << residual << " arc-length residual " << arcLengthResidual << std::endl;

        // Check for convergence
        if (inputReader_.getDesiredOutput() == "all") {
            utils_.writeDataToFile(soldofs_, "Outputs/AL_soldofs_" + std::to_string(iter + 1) + ".txt", true);
            utils_.writeDataToFile(delta_degreesoffreedom, "Outputs/AL_delta_degreesoffreedom_" + std::to_string(iter) + ".txt", true);
        }

        if (residual < tolerance && std::abs(arcLengthResidual) < arcLengthTolerance) {
            // Converged, return both the solution and the final residuals
            std::cout << "### AL. CONVERGED " << std::endl;
            return residual;
        }

        // Update arc-length parameter using a predictor-corrector scheme
        lambda += arcLengthResidual / delta_degreesoffreedom.norm();
    }

    // If we reach here, the arc-length solver did not converge
    throw std::runtime_error("Arc-Length solver did not converge.");

    return 0;
}*/
////////////////////////////////////////////////////////////////////////////////
/*
double Solver::runModifiedNewtonRaphsonSolver(bool applyLoadIncrements) {
    // Initialize some parameters
    double tolerance = 1e-6;
    int maxIterations = 5;

    if (inputReader_.getDesiredOutput() == "all") {
        utils_.writeDataToFile(soldofs_, "Outputs/MNR_soldofs_" + std::to_string(0) + ".txt", true);
    }

    for (int iter = 0; iter < maxIterations; ++iter) {
        // Assembly the system matrix
        SparseSystem system = Assembly();

        // Solve the system using the tangential matrix and residual: KT*dU = R->dU
        Eigen::VectorXd delta_degreesoffreedom = solveSparseSystem(system);

        if (applyLoadIncrements) {
            // Apply load increments to the solution
            for (int i = 0; i < freedofidxs_.size(); i++) {
                soldofs_[freedofidxs_[i]] += delta_degreesoffreedom[i];
            }
        }

        // Calculate the residual (error);
        Eigen::VectorXd eigenR_reduced = Eigen::Map<Eigen::VectorXd>(system.R_reduced.data(), system.R_reduced.size());
        double residual = eigenR_reduced.norm(); // Calculate the norm

        std::cout << "### MNR. Iteration " << iter << " residual " << residual << std::endl;

        // Check for convergence
        if (inputReader_.getDesiredOutput() == "all") {
            utils_.writeDataToFile(soldofs_, "Outputs/MNR_soldofs_" + std::to_string(iter + 1) + ".txt", true);
            utils_.writeDataToFile(delta_degreesoffreedom, "Outputs/MNR_delta_degreesoffreedom_" + std::to_string(iter) + ".txt", true);
        }

        if (residual < tolerance) {
            // Converged, return the final residual
            std::cout << "### MNR. CONVERGED " << std::endl;
            return residual;
        }
    }

    // If we reach here, the Modified Newton-Raphson method did not converge
    throw std::runtime_error("Modified Newton-Raphson did not converge.");
}*/
////////////////////////////////////////////////////////////////////////////////
/*
double Solver::runNewtonRaphsonWithUniformIncrements(const Eigen::VectorXd& totalLoadVector, int numUniformIncrements) {
    // Calculate uniform load increments
    Eigen::VectorXd uniformLoadIncrement = totalLoadVector / numUniformIncrements;
    Eigen::VectorXd currentLoad = Eigen::VectorXd::Zero(totalLoadVector.size());

    double totalResidual = 0.0; // Track the total residual over all increments
    double tolerance=1e-6;
    // Loop over the specified number of increments
    for (int increment = 0; increment < numUniformIncrements; ++increment) {
        // Run the standard Newton-Raphson for each increment
        double residual = runNewtonRaphson();

        // Accumulate the total residual
        totalResidual += residual;

        // Update the current load with the uniform load increment
        currentLoad += uniformLoadIncrement;

        // Print information about the current increment
        std::cout << "### Increment " << increment + 1 << "/" << numUniformIncrements
                  << " completed with residual " << residual << std::endl;

        if (residual < tolerance) {
            // Converged for this increment, continue to the next one
            continue;
        } else {
            // If not converged, you can decide how to handle it, e.g., break or throw an error.
            std::cerr << "### Increment " << increment + 1 << "/" << numUniformIncrements
                      << " did not converge. Terminating." << std::endl;
            break;
        }
    }

    // Return the total residual over all increments
    return totalResidual;
}*/

/*
// Function to decide and run the solver
double Solver::decideAndRunSolver() {
    // Check the settings from inputReader_ to determine which solver to use
    std::string desiredOutput = inputReader_.getDesiredOutput();
    bool useArcLengthSolver = inputReader_.useArcLengthSolver(); // Assuming you have a getter for this setting

    if (useArcLengthSolver) {
        // Use the Arc-Length solver
        ArcLengthSolver arcLengthSolver(inputReader_, utils_);
        return arcLengthSolver.runArcLengthSolver();
    } else if (desiredOutput == "uniformIncrements") {
        // Use the Newton-Raphson solver with uniform increments
        Eigen::VectorXd totalLoadVector = inputReader_.getTotalLoadVector(); // Assuming you have a getter for this
        int numUniformIncrements = inputReader_.getNumUniformIncrements(); // Assuming you have a getter for this
        return runNewtonRaphsonWithUniformIncrements(totalLoadVector, numUniformIncrements);
    } else {
        // Use the standard Newton-Raphson solver
        return runNewtonRaphson();
    }
}
*/