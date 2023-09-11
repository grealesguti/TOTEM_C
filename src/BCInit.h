#ifndef BCINIT_H
#define BCINIT_H

#include "InputReader.h"
#include "Elements.h"
#include "Utils.h"
#include <armadillo>
#include <omp.h>

#include "Mesh.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class BCInit {
public:
    BCInit(const InputReader& inputReader, Mesh& mesh);
    void boundaryConditions();
    void createLoadVector();
    const double getInitialDof(int i) const {return initialdofs_[i];}
    arma::mat getloadVector() const {return loadVector_;}

    // Integrators
    arma::mat CteSurfBC(const arma::mat&,const arma::mat& coords, double heatvalue);

private:
    const InputReader& inputReader_;
    Mesh& mesh_;
    Utils utils_;
    Elements elements_; // Add a member variable of type Elements

    std::string meshFileName_;
    std::vector<double> initialdofs_;
    arma::mat loadVector_;
    std::function<arma::mat(const arma::mat&, const arma::mat&, double)> integrationFunction_;

};

#endif // BCINIT_H
