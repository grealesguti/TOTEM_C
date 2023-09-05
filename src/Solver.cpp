#include "Solver.h"
using namespace arma;
/*
Solver::Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit)
    : inputReader_(inputReader), mesh_(mesh), elements_(), bcinit_(bcinit) {
    // Get the number of nodes from the mesh
    int numNodes = mesh_.getNumAllNodes();

    // Initialize the loadVector_ member with a size double the number of nodes
    loadVector_.resize(2 * numNodes);
    initialdofs_.resize(2 * numNodes);

    // You may want to initialize other member variables here
}

arma::mat Solver::Assembly(const arma::mat& U0, const arma::mat& coords, double heatvalue) {
    int ITMAX = 100; // Maximum number of iterations (adjust as needed)
    double tol = 1e-6; // Tolerance for convergence (adjust as needed)
    
    int ntot = coords.n_rows;
    int dof = 2; // Assuming 2 degrees of freedom per node
    
    int nele = order.n_rows;
    
    arma::mat U = U0;
    arma::mat dU = arma::zeros<arma::mat>(dof * ntot, 1);
    
    arma::vec sigma = matp.row(0);
    arma::vec sigmainf = matp.row(1);
    arma::vec a = matp.row(2);
    arma::vec ainf = matp.row(3);
    arma::vec kx = matp.row(4);
    
    arma::mat Ra = arma::zeros<arma::mat>(dof * ntot, 1);
    arma::mat KJs = arma::zeros<arma::mat>(dof * ntot, dof * ntot);
    
    int nnv = 20 * 2 * 20 * 2;
    
    int it = 0;
    double err = 1.0;
    
    while (it < ITMAX && err > tol) {
        it++;
        
        // Assembly
        Ra.zeros();
        KJs.zeros();
        
        for (int ii = 0; ii < nele; ii++) {
            arma::mat Rs = arma::zeros<arma::mat>(dof * ntot, 1);
            
            // Material properties of element ii (you'll need to adapt this)
            // arma::mat De, Da, Dk;
            // GetDeDaDk(sigma(ii), sigmainf(ii), a(ii), ainf(ii), kx(ii), De, Da, Dk);
            
            // Degrees of freedom of element ii (you'll need to adapt this)
            arma::uvec doforderT = (order.row(ii) - 1) * dof;
            arma::uvec doforderV = (order.row(ii) - 1) * dof + 1;
            
            arma::uvec doforder = arma::join_cols(doforderT, doforderV);
            
            arma::mat Te = U(doforderT);
            arma::mat Ve = U(doforderV);
            
            arma::mat KJ, R;
            
            // Call your GaussKAS14A_TO function to get KJ and R (you'll need to adapt this)
            // GaussKAS14A_TO(coords, order.row(ii), Te, Ve, matp, matv, ii, sysv, localval, xx(ii), p, seebp, rhop, KJ, R);
            
            // Assembly in global residual and jacobian matrix
            Rs(doforder, 0) = R;
            Ra += Rs;
            KJs.submat(doforder, doforder) += KJ;
        }
        
        // Solve system
        arma::mat R_b = Ra.rows(freedofs) - F(freedofs) + K_conv.submat(freedofs, freedofs) * U(freedofs);
        arma::mat KJ_b = KJs.submat(freedofs, freedofs) - K_conv.submat(freedofs, freedofs);
        
        arma::mat dUf = arma::solve(KJ_b, R_b); // Calculation of step
        
        U(freedofs) += dUf; // Update current dof values
        
        err = arma::norm(R_b);
        
        ev(it) = err;
    }
    
    return U;
} 
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