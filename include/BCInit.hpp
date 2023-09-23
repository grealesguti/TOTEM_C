#ifndef BCINIT_H
#define BCINIT_H

// Local
#include "InputReader.hpp"
#include "Mesh.hpp"
#include "Utils.hpp"
#include "utils/data.hpp"

// Armadillo
#include <armadillo>

// OpenMP
#include <omp.h>

// Gmsh
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class BCInit {
public:
    BCInit(const InputReader& inputReader, Mesh& mesh);
    void boundaryConditions();
    void createLoadVector();
    const double getInitialDof(int i) const {return initialdofs_[i];}
    inline const Vector<double> getAllInitialDof() const {return initialdofs_;}

    inline const Armadillo<arma::mat> getloadVector() const {return loadVector_;}

    // Integrators
    arma::mat CteSurfBC(const arma::mat&, const Armadillo<arma::mat>& coords, double heatvalue, int element);


private:
    const InputReader& inputReader_;
    Mesh mesh_;
    Utils utils_;

    std::string meshFileName_;
    std::vector<double> initialdofs_;
    arma::mat loadVector_;
    std::function<arma::mat(const arma::mat&, const Armadillo<arma::mat>&, const double, const int)> integrationFunction_;

};

#endif // BCINIT_H
