#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <armadillo>
#include <functional>
#include "Mesh.h"

class Utils {
public:
    Utils(); // Constructor declaration
    struct IntegrationResult {
        arma::mat KT;
        arma::mat R;
    };
    void getGaussWeightsAndPoints(int order, arma::mat& weights, arma::mat& gaussPoints);
    void gaussIntegrationBC(int dimension, int order, int elementTag, Mesh mesh, double bcvalue, std::function<arma::mat(const arma::mat&,const arma::mat&, double)> func, arma::mat& result);
    // Function to perform Gaussian integration
    IntegrationResult gaussIntegrationK(
        int dimension,
        int order,
        int elementTag,
        Mesh mesh,
        std::function<IntegrationResult(const arma::mat& natcoords, const arma::mat& coords, const arma::mat& dofs)> func
    );
    arma::mat TransformCoordinates(const arma::mat& cooro);

};

#endif // UTILS_H
