#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <armadillo>
#include <functional>

class Utils {
public:
    Utils(); // Constructor declaration

    void getGaussWeightsAndPoints(int order, arma::mat& weights, arma::mat& gaussPoints);
    void gaussIntegration(int dimension, int order, std::function<arma::mat(const arma::mat&)> func, arma::mat& result);
};

#endif // UTILS_H
