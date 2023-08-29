#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <armadillo>
#include <functional>
#include "Mesh.h"

class Utils {
public:
    Utils(); // Constructor declaration

    void getGaussWeightsAndPoints(int order, arma::mat& weights, arma::mat& gaussPoints);
    void gaussIntegrationBC(int dimension, int order, int elementTag, Mesh mesh, double bcvalue, std::function<arma::mat(const arma::mat&,const arma::mat&, double)> func, arma::mat& result);
    arma::mat TransformCoordinates(const arma::mat& cooro);

};

#endif // UTILS_H
