#include "Utils.h"

using namespace arma;

Utils::Utils() {

}



void Utils::getGaussWeightsAndPoints(int order, mat& weights, mat& gaussPoints) {
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

mat Utils::TransformCoordinates(const mat& cooro) {
    // Extract v1 and v2 vectors
    arma::rowvec v1 = cooro.row(1) - cooro.row(0);
    arma::rowvec v2 = cooro.row(2) - cooro.row(0);

    // Calculate thetay1 in radians and degrees
    double thetay1 = -std::acos(arma::dot(v1, v2) / (arma::norm(v1) * arma::norm(v2)));
    double thetay1deg = thetay1 * 180 / arma::datum::pi;

    // Calculate nxy
    arma::rowvec nxy = arma::cross(v1, v2) / arma::norm(arma::cross(v1, v2));

    // Calculate thetaz0 in radians and degrees
    double thetaz0 = std::acos(nxy(0));
    double thetaz0deg = thetaz0 * 180 / arma::datum::pi;

    // Check if thetaz0 is greater than pi/2 and perform necessary adjustments
    if (thetaz0 > arma::datum::pi / 2) {
        double thetaz0old = thetaz0;
        double thetaz0degold = thetaz0deg;
        thetaz0 = arma::datum::pi - thetaz0;
        thetaz0deg = thetaz0 * 180 / arma::datum::pi;
        thetay1 = -thetay1;
    }

    // Create transformation matrices Gz0 and Gy1
    arma::mat Gz0 = {{std::cos(thetaz0), -std::sin(thetaz0), 0},
                     {std::sin(thetaz0), std::cos(thetaz0), 0},
                     {0, 0, 1}};

    arma::mat Gy1 = {{std::cos(thetay1), 0, -std::sin(thetay1)},
                     {0, 1, 0},
                     {std::sin(thetay1), 0, std::cos(thetay1)}};

    // Calculate the transformation matrix G
    arma::mat G = Gz0 * Gy1;

    // Set elements of G close to zero to zero
    G(arma::abs(G) < 1e-12).zeros();

    // Perform the final coordinate transformation
    arma::mat cooro1 = G.t() * cooro.t();

    return cooro1.t();
}

void Utils::gaussIntegrationBC(int dimension, int order, int elementTag, Mesh mesh, double bcvalue, std::function<mat(const mat&,const mat&, double)> func, mat& result) {
    if (dimension < 1 || order < 1) {
        std::cerr << "Invalid dimension or order for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    mat weights;
    mat gaussPoints;
    std::vector<int> nodeTags_el;
    mesh.getElementInfo(elementTag, nodeTags_el);
    arma::mat coordinates(3, size(nodeTags_el));
    arma::mat coordinates_tr(3, size(nodeTags_el));

    coordinates=mesh.getCoordinates(nodeTags_el);
    coordinates_tr=TransformCoordinates(coordinates);

    // Get the Gauss weights and points for the specified order.
    Utils::getGaussWeightsAndPoints(order, weights, gaussPoints);

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
            result += weights(i, 0) * func({gaussPoints(i, 0)},coordinates_tr,bcvalue);
        }
    } else if (dimension == 2) {
        // 2D integration using a double loop.
        for (uword i = 0; i < weights.n_rows; ++i) {
            for (uword j = 0; j < weights.n_rows; ++j) {
                result += weights(i, 0) * weights(j, 0) * func({gaussPoints(i, 0), gaussPoints(j, 0)},coordinates_tr,bcvalue);
            }
        }
    } else if (dimension == 3) {
        if (order == 14) {
            // Special case for 3D integration with order 14.
            // Directly use the given points and weights without looping.
            result = weights * weights.t() % weights * weights.t() % weights * weights.t() % func(gaussPoints,coordinates_tr,bcvalue);
        } else {
            // Generic 3D integration using a triple loop.
            for (uword i = 0; i < weights.n_rows; ++i) {
                for (uword j = 0; j < weights.n_rows; ++j) {
                    for (uword k = 0; k < weights.n_rows; ++k) {
                        result += weights(i, 0) * weights(j, 0) * weights(k, 0) *
                                  func({gaussPoints(i, 0), gaussPoints(j, 0), gaussPoints(k, 0)},coordinates_tr,bcvalue);
                    }
                }
            }
        }
    } else {
        std::cerr << "Invalid dimension for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
    }
}
