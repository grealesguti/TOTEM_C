#ifndef SOLVER_H
#define SOLVER_H

#include "InputReader.h"
#include "Elements.h"
#include "Utils.h"
#include <armadillo>
#include "BCInit.h"
#include "Mesh.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file
#include "C:\msys64\mingw64\include\eigen3\Eigen\Dense"
#include "C:\msys64\mingw64\include\eigen3\Eigen\Sparse"

class Solver {
public:

    struct SparseSystem {
        Eigen::SparseMatrix<double> KsubMatrix;
        std::vector<double> R_reduced;
    };

    Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit);
    
    // SOLVER FUNCTIONS
    // Function to perform assembly and return R_b and KJ_b
    SparseSystem Assembly();
    Eigen::SparseMatrix<double> reduceSystem(const Eigen::SparseMatrix<double>& K);
    Eigen::VectorXd solveSparseSystem(const SparseSystem& system);

    //arma::mat DirectSolver(const arma::mat&,const arma::mat& coords, double heatvalue);
    //arma::mat NewtonRaphson(const arma::mat&,const arma::mat& coords, double heatvalue);

    // PHYSICS
    //arma::mat CoupledThermoelectricity(const arma::mat&,const arma::mat& coords, double heatvalue);
    //arma::mat CoupledThermoelectromechanics(const arma::mat&,const arma::mat& coords, double heatvalue);
    //arma::mat Mechanics(const arma::mat&,const arma::mat& coords, double heatvalue);
    //arma::mat Decoupled_Thermoelectricity_Mechanics(const arma::mat&,const arma::mat& coords, double heatvalue);

private:
    const InputReader& inputReader_;
    Mesh& mesh_;
    Elements elements_; // Add a member variable of type Elements
    BCInit bcinit_; // Add a member variable of type Elements
    Utils utils_;

    std::string meshFileName_;
    std::vector<double> soldofs_;
    arma::mat loadVector_;
    std::vector<int> freedofidxs_;
    Utils::IntegrationResult thermoelectricityintegration(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag);
    std::function<Utils::IntegrationResult(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag)> thermoelectricityintegrationFunction_;

};

#endif // SOLVER_H
