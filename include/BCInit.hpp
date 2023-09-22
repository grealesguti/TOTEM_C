#ifndef BCINIT_H
#define BCINIT_H

#include "InputReader.hpp"
#include "Elements.hpp"
#include "Utils.hpp"
#include <armadillo>
#include <omp.h>

#include "Mesh.hpp"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class BCInit {
public:
    BCInit(const InputReader& inputReader, Mesh& mesh);
    void boundaryConditions();
    void createLoadVector();
    const double getInitialDof(int i) const {return initialdofs_[i];}
    std::vector<double> getAllInitialDof() const {return initialdofs_;}

    arma::mat getloadVector() const {return loadVector_;}

    // Integrators
    arma::mat CteSurfBC(const arma::mat&,const arma::mat& coords, double heatvalue, int element);


private:
    const InputReader& inputReader_;
    Mesh& mesh_;
    Utils utils_;
    Elements elements_; // Add a member variable of type Elements

    std::string meshFileName_;
    std::vector<double> initialdofs_;
    arma::mat loadVector_;
    std::function<arma::mat(const arma::mat&, const arma::mat&, double, int)> integrationFunction_;

};

#endif // BCINIT_H
