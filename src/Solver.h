#ifndef SOLVER_H
#define SOLVER_H

#include "InputReader.h"
#include "Elements.h"
#include "Utils.h"
#include <armadillo>
#include "BCInit.h"
#include "Mesh.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class Solver {
public:
    Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit);
    arma::mat Assembly(const arma::mat&,const arma::mat& coords, double heatvalue);
    arma::mat Thermoelectricity(const arma::mat&,const arma::mat& coords, double heatvalue);
    arma::mat NewtonRaphson(const arma::mat&,const arma::mat& coords, double heatvalue);

private:
    const InputReader& inputReader_;
    Mesh& mesh_;
    Elements elements_; // Add a member variable of type Elements
    BCInit bcinit_; // Add a member variable of type Elements

    std::string meshFileName_;
    std::vector<double> loadVector_,initialdofs_;
};

#endif // SOLVER_H
