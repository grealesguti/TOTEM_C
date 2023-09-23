#ifndef SOLVER_H
#define SOLVER_H

// Local
#include "InputReader.hpp"
#include "Utils.hpp"
#include "BCInit.hpp"
#include "Mesh.hpp"
#include "utils/data.hpp"

// Armadillo
#include <armadillo>

// Gmsh
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

// Eigen
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

class Solver {
public:

    struct SparseSystem {
        Eigen::SparseMatrix<double> KsubMatrix;
        std::vector<double> R_reduced;
    };

    Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit);
    
    // SOLVER FUNCTIONS
    // Function to perform assembly and return R_b and KJ_b
    void Assembly(EigenSparseMatrix KsubMatrix, Vector<double> R_reduced);
    void reduceSystem(const Eigen::SparseMatrix<double>& K, Eigen::SparseMatrix<double>& reducedK);
    void solveSparseSystem(
        const Eigen::SparseMatrix<double>& KsubMatrix,
        const std::vector<double>& R_reduced,
        Eigen::VectorXd& solution);
    double runNewtonRaphson(); // Return both the solution and final residual
    inline const Vector<double> getAllSolDofs() const {return soldofs_;};
    double getAllSolDof(int i) const {return soldofs_[i];};
    double runArcLengthSolver();
    double runModifiedNewtonRaphsonSolver(bool applyLoadIncrements);
    double runNewtonRaphsonWithUniformIncrements(const Eigen::VectorXd& totalLoadVector, int numUniformIncrements);

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
    BCInit bcinit_; // Add a member variable of type Elements
    Utils utils_;

    std::string meshFileName_;
    std::vector<double> soldofs_;
    arma::mat loadVector_;
    std::vector<int> freedofidxs_;
    Utils::IntegrationResult thermoelectricityintegration(const Armadillo<arma::mat>& natcoords, const Armadillo<arma::mat>& coords, const Armadillo<arma::vec>& dofs, const int elementTag);
    std::function<Utils::IntegrationResult(const arma::mat& natcoords, const Armadillo<arma::mat>& coords, const Armadillo<arma::vec>& dofs, const int elementTag)> thermoelectricityintegrationFunction_;

};

#endif // SOLVER_H
