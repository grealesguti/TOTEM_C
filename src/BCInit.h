#ifndef BCINIT_H
#define BCINIT_H

#include "InputReader.h"
#include "Elements.h"
#include "Utils.h"
#include <armadillo>

#include "Mesh.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class BCInit {
public:
    BCInit(const InputReader& inputReader, Mesh& mesh);
    void boundaryConditions();
    void createLoadVector();

    // Integrators
    arma::mat CteSurfBC(const arma::mat&,const arma::mat& coords, double heatvalue);

private:
    const InputReader& inputReader_;
    Mesh& mesh_;
    Elements elements_; // Add a member variable of type Elements

    std::string meshFileName_;
    std::vector<double> loadVector_,initialdofs_;
};

#endif // BCINIT_H
