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
#include <iostream>
#include <armadillo>

arma::mat Utils::TransformCoordinates(const arma::mat& cooro) {
    // Extract columns v1 and v2 vectors
    arma::colvec v1 = cooro.col(1) - cooro.col(0);
    arma::colvec v2 = cooro.col(2) - cooro.col(0);

    // Calculate thetay1 in radians and degrees
    double thetay1 = -std::acos(arma::dot(v1, v2) / (arma::norm(v1) * arma::norm(v2)));
    double thetay1deg = thetay1 * 180 / arma::datum::pi;

    arma::colvec nxy;

    if (v1.n_elem == 3 && v2.n_elem == 3) {
        // Calculate nxy
        nxy = arma::cross(v1, v2) / arma::norm(arma::cross(v1, v2));
    } else {
        // Handle the case where the input vectors are not valid
        std::cerr << "Utils::TransformCoordinates: Invalid input vectors for cross product." << std::endl;
        // Print the input vectors and their sizes
        std::cout << "Vector v1: " << v1 << std::endl;
        std::cout << "Vector v2: " << v2 << std::endl;
        std::cout << "v1 size: " << v1.n_elem << std::endl;
        std::cout << "v2 size: " << v2.n_elem << std::endl;
        return arma::mat(); // Return an empty matrix to indicate an error
    }
      //  std::cout << "Cross product done" << std::endl;

    // Calculate thetaz0 in radians and degrees
    double thetaz0 = std::acos(nxy(0));
    double thetaz0deg = thetaz0 * 180 / arma::datum::pi;
    //std::cout << "angle calculated" << std::endl;

    // Check if thetaz0 is greater than pi/2 and perform necessary adjustments
    if (thetaz0 > arma::datum::pi / 2) {
      //  std::cout << "angle adjustment" << std::endl;
        double thetaz0old = thetaz0;
        double thetaz0degold = thetaz0deg;
        thetaz0 = arma::datum::pi - thetaz0;
        thetaz0deg = thetaz0 * 180 / arma::datum::pi;
        thetay1 = -thetay1;
    }
    
    //std::cout << "create transformation matrixes" << std::endl;
    // Create transformation matrices Gz0 and Gy1
    arma::mat Gz0 = {{std::cos(thetaz0), -std::sin(thetaz0), 0},
                     {std::sin(thetaz0), std::cos(thetaz0), 0},
                     {0, 0, 1}};

    arma::mat Gy1 = {{std::cos(thetay1), 0, -std::sin(thetay1)},
                     {0, 1, 0},
                     {std::sin(thetay1), 0, std::cos(thetay1)}};

    // Calculate the transformation matrix G
    //std::cout << "create full transformation" << std::endl;
    arma::mat G = Gz0 * Gy1;
    //std::cout << "set to zeros" << std::endl;
    // Set elements of G close to zero to zer    std::cout << "create full transformation" << std::endl;o
    for (arma::uword i = 0; i < G.n_rows; ++i) {
        for (arma::uword j = 0; j < G.n_cols; ++j) {
            if (std::abs(G(i, j)) < 1e-12) {
                G(i, j) = 0.0;
            }
        }
    }
    //std::cout << "do the transformation" << std::endl;
    // Perform the final coordinate transformation
    arma::mat cooro1 = G * cooro;
    //std::cout << "Coordinate transformation done" << std::endl;

    return cooro1;
}


void Utils::gaussIntegrationBC(int dimension, int order, int elementTag, Mesh mesh, double bcvalue, std::function<mat(const mat& natcoords,const mat& coords, double value)> func, mat& result) {
    
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
    //std::cout << "Retrieve gauss integration points" << std::endl;

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
    // result = zeros<mat>(dimension, dimension); Should already be of the right dimension??, in this case 8,1?? make a check that it is the right dimension
    //std::cout << "Start the integration " ;

    if (dimension == 1) {
        // 1D integration using a single loop.
        //std::cout << "of dimension 1" << std::endl;
            arma::mat natcoords(1, 1);
            arma::mat f(2, 1);
        for (uword i = 0; i < weights.n_rows; ++i) {
                // Explicitly use the arma::operator* function for multiplication
                f = arma::operator*(func(natcoords, coordinates_tr, bcvalue), weights(i, 0));
                result += f;                }
    } else if (dimension == 2) {
        //std::cout << "of dimension 2" << std::endl;
            arma::mat natcoords(2, 1);
            arma::mat f(4, 1);
        // 2D integration using a double loop.
        for (uword i = 0; i < weights.n_rows; ++i) {
            for (uword j = 0; j < weights.n_rows; ++j) {
                        natcoords(0, 0) = gaussPoints(i, 0);
                        natcoords(1, 0) = gaussPoints(j, 0);
                        // Create a 3x1 matrix with gaussPoints(i, 0), gaussPoints(j, 0), and gaussPoints(k, 0)
                        // Explicitly use the arma::operator* function for multiplication
                        f = arma::operator*(func(natcoords, coordinates_tr, bcvalue), weights(i, 0) * weights(j, 0));
                        result += f;            
             }
        }
    } else if (dimension == 3) {
        if (order == 14) {
            //std::cout << "of dimension 3 and reduced order" << std::endl;
            // Special case for 3D integration with order 14.
            // Directly use the given points and weights without looping.
            result = weights * weights.t() % weights * weights.t() % weights * weights.t() % func(gaussPoints,coordinates_tr,bcvalue);
        } else {
            arma::mat natcoords(3, 1);
            arma::mat f(8, 1);
            //std::cout << "of dimension 3" << std::endl;
            // Generic 3D integration using a triple loop.
            for (uword i = 0; i < weights.n_rows; ++i) {
                for (uword j = 0; j < weights.n_rows; ++j) {
                    for (uword k = 0; k < weights.n_rows; ++k) {
                        natcoords(0, 0) = gaussPoints(i, 0);
                        natcoords(1, 0) = gaussPoints(j, 0);
                        natcoords(2, 0) = gaussPoints(k, 0);
                        // Create a 3x1 matrix with gaussPoints(i, 0), gaussPoints(j, 0), and gaussPoints(k, 0)
                        // Explicitly use the arma::operator* function for multiplication
                        f = arma::operator*(func(natcoords, coordinates_tr, bcvalue), weights(i, 0) * weights(j, 0) * weights(k, 0));
                        result += f;
                    }
                }
            }
        }
    } else {
        std::cerr << "Invalid dimension for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
    }
}


//////////////////////////////////////////////////////////////////////////////////
Utils::IntegrationResult Utils::gaussIntegrationK(
    int dimension,
    int order,
    int elementTag,
    Mesh mesh,
    std::function<Utils::IntegrationResult(const arma::mat& natcoords, const arma::mat& coords, const arma::mat& dofs)> func
) {    
    Utils::IntegrationResult result; // Create a struct to hold KT and R

    if (dimension < 1 || order < 1) {
        std::cerr << "Invalid dimension or order for Gauss integration." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
        return result;
    }
    int nodes_per_element = 8; // Total number of nodes per element
    int dof_per_node = 2; // Assuming 2 degrees of freedom per node
    arma::mat weights;
    arma::mat gaussPoints;
    std::vector<int> nodeTags_el;
    mesh.getElementInfo(elementTag, nodeTags_el);
    arma::mat coordinates(3, nodeTags_el.size());
    arma::mat coordinates_tr(3, nodeTags_el.size());
    arma::mat element_dofs = arma::zeros<arma::mat>(16, 1);

        int cc = 0;
        for (int nodeTag : nodeTags_el) {
            element_dofs[cc] = nodeTag * dof_per_node;
            element_dofs[cc + nodes_per_element] = nodeTag * dof_per_node + 1;
            cc += 1;
        }
    // get current dofs from elementdofs index, this is the input!!!

    this->getGaussWeightsAndPoints(order, weights, gaussPoints);

    if (weights.is_empty() || gaussPoints.is_empty()) {
        std::cerr << "Invalid order for Gauss integration." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
        return result;
    }

    if (weights.n_rows != gaussPoints.n_rows) {
        std::cerr << "Weights and Gauss points have mismatched dimensions." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
        return result;
    }

    // Initialize result.KT and result.R to zero matrices of appropriate dimensions.
    if (dimension == 1) {
        // 1D integration using a single loop.
        arma::mat natcoords(1, 1);
        arma::mat Re(4, 1, arma::fill::zeros);
        arma::mat KTe(4, 4, arma::fill::zeros);

        for (uword i = 0; i < weights.n_rows; ++i) {
            natcoords(0, 0) = gaussPoints(i, 0);
            Utils::IntegrationResult localResult = func(natcoords, coordinates_tr,element_dofs );
            Re += localResult.R * weights(i, 0);
            KTe += localResult.KT * weights(i, 0);
        }

            result.KT = KTe;
            result.R = Re;

    } else if (dimension == 2) {
        // 2D integration using a double loop.
        arma::mat natcoords(2, 1);
        arma::mat R(8, 1, arma::fill::zeros);
        arma::mat KT(8, 8, arma::fill::zeros);

        for (uword i = 0; i < weights.n_rows; ++i) {
            for (uword j = 0; j < weights.n_rows; ++j) {
                natcoords(0, 0) = gaussPoints(i, 0);
                natcoords(1, 0) = gaussPoints(j, 0);
                Utils::IntegrationResult localResult = func(natcoords, coordinates_tr,element_dofs);
                R += localResult.R * weights(i, 0) * weights(j, 0);
                KT += localResult.KT * weights(i, 0) * weights(j, 0);
            }
        }

            result.KT = KT;
            result.R = R;

    } else if (dimension == 3) {
        if (order == 14) {
            // Special case for 3D integration with order 14.
            // Directly use the given points and weights without looping.
            //result.KT = weights * weights.t() % weights * weights.t() % weights * weights.t() * func(gaussPoints, coordinates_tr,element_dofs);
            //result.R = arma::zeros<arma::mat>(16, 1); // Initialize R to a zero matrix.
        } else {
            arma::mat natcoords(3, 1);
            arma::mat R(16, 1, arma::fill::zeros);
            arma::mat KT(16, 16, arma::fill::zeros);

            for (uword i = 0; i < weights.n_rows; ++i) {
                for (uword j = 0; j < weights.n_rows; ++j) {
                    for (uword k = 0; k < weights.n_rows; ++k) {
                        natcoords(0, 0) = gaussPoints(i, 0);
                        natcoords(1, 0) = gaussPoints(j, 0);
                        natcoords(2, 0) = gaussPoints(k, 0);
                        Utils::IntegrationResult localResult = func(natcoords, coordinates_tr,element_dofs);
                        R += localResult.R * weights(i, 0) * weights(j, 0) * weights(k, 0);
                        KT += localResult.KT * weights(i, 0) * weights(j, 0) * weights(k, 0);
                    }
                }
            }

            result.KT = KT;
            result.R = R;
        }
    } else {
        std::cerr << "Invalid dimension for Gauss integration." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
    }

    return result; // Return the struct containing KT and R.
}
