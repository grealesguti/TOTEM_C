#include "Utils.h"

using namespace arma;

Utils::Utils() {
    // Initialize a member variable (if needed)
}



void getGaussWeightsAndPoints(int order, mat& weights, mat& gaussPoints) {
    // Check if the order is valid
    if (order < 1 || order > 14) {
        std::cerr << "Invalid order for Gauss integration." << std::endl;
        return;
    }

    // Initialize the matrices with zeros
    weights = mat(order, 1, fill::zeros);
    gaussPoints = mat(order, 1, fill::zeros);

    if (order == 1) {
        weights(0, 0) = 2.0;
        gaussPoints(0, 0) = 0.0;
    } else if (order == 2) {
        weights(0, 0) = 1.0;
        weights(1, 0) = 1.0;
        gaussPoints(0, 0) = -0.577350269189626;
        gaussPoints(1, 0) = 0.577350269189626;
    } else if (order == 3) {
        weights(0, 0) = 0.555555555555556;
        weights(1, 0) = 0.888888888888889;
        weights(2, 0) = 0.555555555555556;
        gaussPoints(0, 0) = -0.774596669241483;
        gaussPoints(1, 0) = 0.0;
        gaussPoints(2, 0) = 0.774596669241483;
    } else if (order == 14) {
        // 14-point Gauss integration rule.
        double cor = 0.335180055401662;
        double cen = 0.886426592797784;
        
        mat G14 = {
            {-cor, -cor, -cor, -cor, +cor, +cor, +cor, +cor, -cen, cen, 0.0, 0.0, 0.0, 0.0},
            {-cor, -cor, +cor, +cor, -cor, -cor, +cor, +cor, 0.0, 0.0, -cen, cen, 0.0, 0.0},
            {-cor, +cor, -cor, +cor, -cor, +cor, -cor, +cor, 0.0, 0.0, 0.0, 0.0, -cen, cen}
        };
        
        vec W14(order);
        W14.subvec(0, 7).fill(cor);
        W14.subvec(8, 13).fill(cen);

        weights = W14;
        gaussPoints = G14;
    }
}


void gaussIntegration(int dimension, int order, std::function<mat(const mat&)> func, mat& result) {
    if (dimension < 1 || order < 1) {
        std::cerr << "Invalid dimension or order for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    mat weights;
    mat gaussPoints;

    // Get the Gauss weights and points for the specified order.
    getGaussWeightsAndPoints(order, weights, gaussPoints);

    if (weights.is_empty() || gaussPoints.is_empty()) {
        std::cerr << "Invalid order for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    if (weights.n_rows != gaussPoints.n_rows) {
        std::cerr << "Weights and Gauss points have mismatched dimensions." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    // Initialize result to zero matrix of appropriate dimensions.
    result = zeros<mat>(dimension, dimension);

    if (dimension == 1) {
        // 1D integration using a single loop.
        for (uword i = 0; i < weights.n_rows; ++i) {
            result += weights(i, 0) * func({gaussPoints(i, 0)});
        }
    } else if (dimension == 2) {
        // 2D integration using a double loop.
        for (uword i = 0; i < weights.n_rows; ++i) {
            for (uword j = 0; j < weights.n_rows; ++j) {
                result += weights(i, 0) * weights(j, 0) * func({gaussPoints(i, 0), gaussPoints(j, 0)});
            }
        }
    } else if (dimension == 3) {
        if (order == 14) {
            // Special case for 3D integration with order 14.
            // Directly use the given points and weights without looping.
            result = weights * weights.t() % weights * weights.t() % weights * weights.t() % func(gaussPoints);
        } else {
            // Generic 3D integration using a triple loop.
            for (uword i = 0; i < weights.n_rows; ++i) {
                for (uword j = 0; j < weights.n_rows; ++j) {
                    for (uword k = 0; k < weights.n_rows; ++k) {
                        result += weights(i, 0) * weights(j, 0) * weights(k, 0) *
                                  func({gaussPoints(i, 0), gaussPoints(j, 0), gaussPoints(k, 0)});
                    }
                }
            }
        }
    } else {
        std::cerr << "Invalid dimension for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
    }
}
