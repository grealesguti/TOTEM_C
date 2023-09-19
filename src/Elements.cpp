// elements.cpp

#include "Elements.h"

// Constructor
Elements::Elements() {}

// Destructor
Elements::~Elements() {}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TRI 3 nodes /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to evaluate shape functions for a linear triangular element with 3 nodes
void Elements::EvaluateLinearTriangularShapeFunctions(double xi, double eta, arma::vec& shapeFunctions) {
    // Ensure the shapeFunctions vector is of the correct size (3 nodes)
    shapeFunctions.set_size(3);

    // Calculate the shape functions for a linear triangular element
    double xi1 = xi;
    double xi2 = eta;
    double xi3 = 1.0 - xi1 - xi2;

    // Shape functions for a linear triangular element
    shapeFunctions(0) = xi1;
    shapeFunctions(1) = xi2;
    shapeFunctions(2) = xi3;
}

// Function to calculate shape function derivatives for a linear triangular element with 3 nodes
void Elements::CalculateLinearTriangularShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives) {
    // Ensure the shapeFunctionDerivatives matrix is of the correct size (3 nodes x 2 derivatives)
    shapeFunctionDerivatives.set_size(3, 2);

    // Calculate the shape function derivatives for a linear triangular element
    double xi1 = xi;
    double xi2 = eta;
    double xi3 = 1.0 - xi1 - xi2;

    // Derivatives with respect to xi and eta
    shapeFunctionDerivatives(0, 0) = 1.0;  // dN1/dxi
    shapeFunctionDerivatives(0, 1) = 0.0;  // dN1/deta

    shapeFunctionDerivatives(1, 0) = 0.0;  // dN2/dxi
    shapeFunctionDerivatives(1, 1) = 1.0;  // dN2/deta

    shapeFunctionDerivatives(2, 0) = -1.0; // dN3/dxi
    shapeFunctionDerivatives(2, 1) = -1.0; // dN3/deta
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// QUAD 4 nodes /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate shape functions for a linear quadrilateral element with 4 nodes
arma::mat Elements::EvaluateLinearQuadrilateralShapeFunctions(double xi, double eta) {
    // Define the shape functions as a 4x1 matrix
    arma::mat shapeFunctions(4, 1, arma::fill::zeros);

    // Compute the shape functions for a linear quadrilateral element
    shapeFunctions(0, 0) = 0.25 * (1 - xi) * (1 - eta); // N1
    shapeFunctions(1, 0) = 0.25 * (1 + xi) * (1 - eta); // N2
    shapeFunctions(2, 0) = 0.25 * (1 + xi) * (1 + eta); // N3
    shapeFunctions(3, 0) = 0.25 * (1 - xi) * (1 + eta); // N4

    return shapeFunctions;
}

// Evaluate derivatives of shape functions for a linear quadrilateral element with 4 nodes
void Elements::EvaluateLinearQuadrilateralShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives) {
    // Ensure the shapeFunctionDerivatives matrix is of the correct size (4 nodes x 2 derivatives)
    shapeFunctionDerivatives.set_size(4, 2);

    // Define the derivatives of shape functions for a linear quadrilateral element
    shapeFunctionDerivatives(0, 0) = -0.25 * (1 - eta);  // dN1/dxi
    shapeFunctionDerivatives(0, 1) = -0.25 * (1 - xi);   // dN1/deta

    shapeFunctionDerivatives(1, 0) = 0.25 * (1 - eta);   // dN2/dxi
    shapeFunctionDerivatives(1, 1) = -0.25 * (1 + xi);  // dN2/deta

    shapeFunctionDerivatives(2, 0) = 0.25 * (1 + eta);   // dN3/dxi
    shapeFunctionDerivatives(2, 1) = 0.25 * (1 + xi);   // dN3/deta

    shapeFunctionDerivatives(3, 0) = -0.25 * (1 + eta);  // dN4/dxi
    shapeFunctionDerivatives(3, 1) = 0.25 * (1 - xi);   // dN4/deta
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// QUAD 8 nodes /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate shape functions for a quadratic quadrilateral element with 8 nodes
void Elements::EvaluateQuadraticQuadrilateralShapeFunctions(double xi, double eta, arma::vec& shapeFunctions) {
    // Ensure the shapeFunctions vector is of the correct size (8 nodes)
    shapeFunctions.set_size(8);

    // Define the shape functions for a quadratic quadrilateral element
    double xi1 = -1.0;
    double xi2 = 1.0;
    double eta1 = -1.0;
    double eta2 = 1.0;

    shapeFunctions(0) = 0.25 * (xi - xi1) * (eta - eta1);
    shapeFunctions(1) = 0.25 * (xi - xi2) * (eta - eta1);
    shapeFunctions(2) = 0.25 * (xi - xi2) * (eta - eta2);
    shapeFunctions(3) = 0.25 * (xi - xi1) * (eta - eta2);
    shapeFunctions(4) = 0.5 * (1 - xi * xi) * 0.25 * (eta - eta1);
    shapeFunctions(5) = 0.5 * (1 - xi * xi) * 0.25 * (eta - eta2);
    shapeFunctions(6) = 0.5 * (1 - eta * eta) * 0.25 * (xi - xi2);
    shapeFunctions(7) = 0.5 * (1 - eta * eta) * 0.25 * (xi - xi1);
}

// Function to calculate derivatives of shape functions for a quadratic quadrilateral element with 8 nodes
void Elements::CalculateQuadraticQuadrilateralShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives) {
    // Ensure the shapeFunctionDerivatives matrix is of the correct size (8 nodes x 2 derivatives)
    shapeFunctionDerivatives.set_size(8, 2);

    // Define the derivatives of shape functions for a quadratic quadrilateral element
    double xi1 = -1.0;
    double xi2 = 1.0;
    double eta1 = -1.0;
    double eta2 = 1.0;

    shapeFunctionDerivatives(0, 0) = -0.25 * (eta - eta1);  // dN1/dxi
    shapeFunctionDerivatives(0, 1) = -0.25 * (xi - xi1);   // dN1/deta
    shapeFunctionDerivatives(1, 0) = 0.25 * (eta - eta1);   // dN2/dxi
    shapeFunctionDerivatives(1, 1) = -0.25 * (xi - xi2);   // dN2/deta
    shapeFunctionDerivatives(2, 0) = 0.25 * (eta - eta2);   // dN3/dxi
    shapeFunctionDerivatives(2, 1) = 0.25 * (xi - xi2);    // dN3/deta
    shapeFunctionDerivatives(3, 0) = -0.25 * (eta - eta2);  // dN4/dxi
    shapeFunctionDerivatives(3, 1) = 0.25 * (xi - xi1);    // dN4/deta
    shapeFunctionDerivatives(4, 0) = -0.5 * xi * (eta - eta1);  // dN5/dxi
    shapeFunctionDerivatives(4, 1) = -0.5 * (1 - xi * xi1);    // dN5/deta
    shapeFunctionDerivatives(5, 0) = 0.5 * xi * (eta - eta2);   // dN6/dxi
    shapeFunctionDerivatives(5, 1) = -0.5 * (1 - xi * xi2);    // dN6/deta
    shapeFunctionDerivatives(6, 0) = 0.5 * (1 - eta * eta1);    // dN7/dxi
    shapeFunctionDerivatives(6, 1) = -0.5 * eta * (xi - xi2);   // dN7/deta
    shapeFunctionDerivatives(7, 0) = -0.5 * (1 - eta * eta2);   // dN8/dxi
    shapeFunctionDerivatives(7, 1) = -0.5 * eta * (xi - xi1);   // dN8/deta
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TETRA 4 nodes /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to evaluate shape functions for a linear tetrahedral element with 4 nodes
void Elements::EvaluateLinearTetrahedraShapeFunctions(double xi, double eta, double zeta, arma::vec& shapeFunctions) {
    // Ensure the shapeFunctions vector is of the correct size (4 nodes)
    shapeFunctions.set_size(4);

    // Calculate the shape functions for a linear tetrahedral element
    double xi1 = xi;
    double xi2 = eta;
    double xi3 = zeta;
    double xi4 = 1.0 - xi1 - xi2 - xi3;

    // Shape functions for a linear tetrahedral element
    shapeFunctions(0) = xi1;
    shapeFunctions(1) = xi2;
    shapeFunctions(2) = xi3;
    shapeFunctions(3) = xi4;
}
// Function to calculate shape function derivatives for a linear tetrahedral element with 4 nodes
void Elements::CalculateLinearTetrahedraShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives) {
    // Ensure the shapeFunctionDerivatives matrix is of the correct size (4 nodes x 3 derivatives)
    shapeFunctionDerivatives.set_size(4, 3);

    // Calculate the shape function derivatives for a linear tetrahedral element
    double xi1 = xi;
    double xi2 = eta;
    double xi3 = zeta;
    double xi4 = 1.0 - xi1 - xi2 - xi3;

    // Derivatives with respect to xi, eta, and zeta
    shapeFunctionDerivatives(0, 0) = 1.0;  // dN1/dxi
    shapeFunctionDerivatives(0, 1) = 0.0;  // dN1/deta
    shapeFunctionDerivatives(0, 2) = 0.0;  // dN1/dzeta

    shapeFunctionDerivatives(1, 0) = 0.0;  // dN2/dxi
    shapeFunctionDerivatives(1, 1) = 1.0;  // dN2/deta
    shapeFunctionDerivatives(1, 2) = 0.0;  // dN2/dzeta

    shapeFunctionDerivatives(2, 0) = 0.0;  // dN3/dxi
    shapeFunctionDerivatives(2, 1) = 0.0;  // dN3/deta
    shapeFunctionDerivatives(2, 2) = 1.0;  // dN3/dzeta

    shapeFunctionDerivatives(3, 0) = -1.0; // dN4/dxi
    shapeFunctionDerivatives(3, 1) = -1.0; // dN4/deta
    shapeFunctionDerivatives(3, 2) = -1.0; // dN4/dzeta
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HEXA 8 nodes /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate shape functions for a linear hexahedral element with 8 nodes
void Elements::EvaluateHexahedralLinearShapeFunctions(double xi, double eta, double zeta, arma::vec& shapeFunctions) {
    // Ensure the shapeFunctions vector is of the correct size (8 nodes)
    shapeFunctions.set_size(8);

    // Define the shape functions for a linear hexahedral element
    shapeFunctions(0) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
    shapeFunctions(1) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
    shapeFunctions(2) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
    shapeFunctions(3) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
    shapeFunctions(4) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
    shapeFunctions(5) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
    shapeFunctions(6) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
    shapeFunctions(7) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);
}

// Function to calculate derivatives of shape functions for a linear hexahedral element with 8 nodes
void Elements::CalculateHexahedralLinearShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives) {
    // Ensure the shapeFunctionDerivatives matrix is of the correct size (8 nodes x 3 derivatives)
    shapeFunctionDerivatives.set_size(8, 3);

    // Define the derivatives of shape functions for a linear hexahedral element
    shapeFunctionDerivatives(0, 0) = -0.125 * (1 - eta) * (1 - zeta);  // dN1/dxi
    shapeFunctionDerivatives(0, 1) = -0.125 * (1 - xi) * (1 - zeta);   // dN1/deta
    shapeFunctionDerivatives(0, 2) = -0.125 * (1 - xi) * (1 - eta);    // dN1/dzeta

    shapeFunctionDerivatives(1, 0) = 0.125 * (1 - eta) * (1 - zeta);   // dN2/dxi
    shapeFunctionDerivatives(1, 1) = -0.125 * (1 + xi) * (1 - zeta);  // dN2/deta
    shapeFunctionDerivatives(1, 2) = -0.125 * (1 + xi) * (1 - eta);   // dN2/dzeta

    shapeFunctionDerivatives(2, 0) = 0.125 * (1 + eta) * (1 - zeta);   // dN3/dxi
    shapeFunctionDerivatives(2, 1) = 0.125 * (1 + xi) * (1 - zeta);   // dN3/deta
    shapeFunctionDerivatives(2, 2) = -0.125 * (1 + xi) * (1 + eta);  // dN3/dzeta

    shapeFunctionDerivatives(3, 0) = -0.125 * (1 + eta) * (1 - zeta);  // dN4/dxi
    shapeFunctionDerivatives(3, 1) = 0.125 * (1 - xi) * (1 - zeta);   // dN4/deta
    shapeFunctionDerivatives(3, 2) = -0.125 * (1 - xi) * (1 + eta);  // dN4/dzeta

    shapeFunctionDerivatives(4, 0) = -0.125 * (1 - eta) * (1 + zeta);  // dN5/dxi
    shapeFunctionDerivatives(4, 1) = -0.125 * (1 - xi) * (1 + zeta);   // dN5/deta
    shapeFunctionDerivatives(4, 2) = 0.125 * (1 - xi) * (1 - eta);    // dN5/dzeta

    shapeFunctionDerivatives(5, 0) = 0.125 * (1 - eta) * (1 + zeta);   // dN6/dxi
    shapeFunctionDerivatives(5, 1) = -0.125 * (1 + xi) * (1 + zeta);  // dN6/deta
    shapeFunctionDerivatives(5, 2) = 0.125 * (1 + xi) * (1 - eta);   // dN6/dzeta

    shapeFunctionDerivatives(6, 0) = 0.125 * (1 + eta) * (1 + zeta);   // dN7/dxi
    shapeFunctionDerivatives(6, 1) = 0.125 * (1 + xi) * (1 + zeta);   // dN7/deta
    shapeFunctionDerivatives(6, 2) = 0.125 * (1 + xi) * (1 + eta);   // dN7/dzeta

    shapeFunctionDerivatives(7, 0) = -0.125 * (1 + eta) * (1 + zeta);  // dN8/dxi
    shapeFunctionDerivatives(7, 1) = 0.125 * (1 - xi) * (1 + zeta);   // dN8/deta
    shapeFunctionDerivatives(7, 2) = 0.125 * (1 - xi) * (1 + eta);   // dN8/dzeta
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HEXA 20 nodes /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate shape functions for a hexahedral serendipity element with 20 nodes
void Elements::CalculateHexahedralSerendipityShapeFunctions(double xi, double eta, double zeta, arma::vec& shapeFunctions) {
   // Ensure the shapeFunctions vector is of the correct size (20 nodes)
    shapeFunctions.set_size(20);

    // Define the shape functions for a hexahedral serendipity element
    double xi1 = xi;
    double xi2 = eta;
    double xi3 = zeta;

    shapeFunctions(0) = (1. - xi2) * (1. - xi3) * (1. - xi1) * (-xi2 - xi3 - xi1 - 2) / 8;
    shapeFunctions(1) = (1. + xi2) * (1. - xi3) * (1. - xi1) * (xi2 - xi3 - xi1 - 2) / 8;
    shapeFunctions(2) = (1. + xi2) * (1. + xi3) * (1. - xi1) * (xi2 + xi3 - xi1 - 2) / 8;
    shapeFunctions(3) = (1. - xi2) * (1. + xi3) * (1. - xi1) * (-xi2 + xi3 - xi1 - 2) / 8;

    shapeFunctions(4) = (1. - xi2) * (1. - xi3) * (1. + xi1) * (-xi2 - xi3 + xi1 - 2) / 8;
    shapeFunctions(5) = (1. + xi2) * (1. - xi3) * (1. + xi1) * (xi2 - xi3 + xi1 - 2) / 8;
    shapeFunctions(6) = (1. + xi2) * (1. + xi3) * (1. + xi1) * (xi2 + xi3 + xi1 - 2) / 8;
    shapeFunctions(7) = (1. - xi2) * (1. + xi3) * (1. + xi1) * (-xi2 + xi3 + xi1 - 2) / 8;

    shapeFunctions(8) = (1 - xi2 * xi2) * (1 - xi3) * (1 - xi1) / 4;
    shapeFunctions(9) = (1 + xi2) * (1 - xi3 * xi3) * (1 - xi1) / 4;
    shapeFunctions(10) = (1 - xi2 * xi2) * (1 + xi3) * (1 - xi1) / 4;
    shapeFunctions(11) = (1 - xi2) * (1 - xi3 * xi3) * (1 - xi1) / 4;

    shapeFunctions(12) = (1 - xi2 * xi2) * (1 - xi3) * (1 + xi1) / 4;
    shapeFunctions(13) = (1 + xi2) * (1 - xi3 * xi3) * (1 + xi1) / 4;
    shapeFunctions(14) = (1 - xi2 * xi2) * (1 + xi3) * (1 + xi1) / 4;
    shapeFunctions(15) = (1 - xi2) * (1 - xi3 * xi3) * (1 + xi1) / 4;

    shapeFunctions(16) = (1 - xi2) * (1 - xi3) * (1 - xi1 * xi1) / 4;
    shapeFunctions(17) = (1 + xi2) * (1 - xi3) * (1 - xi1 * xi1) / 4;
    shapeFunctions(18) = (1 + xi2) * (1 + xi3) * (1 - xi1 * xi1) / 4;
    shapeFunctions(19) = (1 - xi2) * (1 + xi3) * (1 - xi1 * xi1) / 4;
}

// Function to calculate derivatives of shape functions for a hexahedral serendipity element with 20 nodes
void Elements::CalculateHexahedralSerendipityShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives0) {
    // Ensure the shapeFunctionDerivatives matrix is of the correct size (20 nodes x 3 derivatives)
    shapeFunctionDerivatives0.set_size(20, 3);
    arma::mat shapeFunctionDerivatives;
    shapeFunctionDerivatives.set_size(3, 20);
    // Define the shape function derivatives for a hexahedral serendipity element
    double xi1 = xi;
    double xi2 = eta;
    double xi3 = zeta;

// Define the provided shape function derivatives
    shapeFunctionDerivatives(0, 0) = ((xi2 - 1) * (xi3 - 1) * (xi1 + xi2 + xi3 + 2)) / 8 + ((xi1 - 1) * (xi2 - 1) * (xi3 - 1)) / 8;
    shapeFunctionDerivatives(0, 1) = -((xi2 + 1) * (xi3 - 1) * (xi1 - xi2 + xi3 + 2)) / 8 - ((xi1 - 1) * (xi2 + 1) * (xi3 - 1)) / 8;
    shapeFunctionDerivatives(0, 2) = ((xi2 + 1) * (xi3 + 1) * (xi1 - xi2 - xi3 + 2)) / 8 + ((xi1 - 1) * (xi2 + 1) * (xi3 + 1)) / 8;
    shapeFunctionDerivatives(0, 3) = -((xi2 - 1) * (xi3 + 1) * (xi1 + xi2 - xi3 + 2)) / 8 - ((xi1 - 1) * (xi2 - 1) * (xi3 + 1)) / 8;
    shapeFunctionDerivatives(0, 4) = ((xi1 + 1) * (xi2 - 1) * (xi3 - 1)) / 8 - ((xi2 - 1) * (xi3 - 1) * (xi2 - xi1 + xi3 + 2)) / 8;
    shapeFunctionDerivatives(0, 5) = -((xi2 + 1) * (xi3 - 1) * (xi1 + xi2 - xi3 - 2)) / 8 - ((xi1 + 1) * (xi2 + 1) * (xi3 - 1)) / 8;
    shapeFunctionDerivatives(0, 6) = ((xi2 + 1) * (xi3 + 1) * (xi1 + xi2 + xi3 - 2)) / 8 + ((xi1 + 1) * (xi2 + 1) * (xi3 + 1)) / 8;
    shapeFunctionDerivatives(0, 7) = -((xi2 - 1) * (xi3 + 1) * (xi1 - xi2 + xi3 - 2)) / 8 - ((xi1 + 1) * (xi2 - 1) * (xi3 + 1)) / 8;
    shapeFunctionDerivatives(0, 8) = -((xi2 * xi2 - 1) * (xi3 - 1)) / 4;
    shapeFunctionDerivatives(0, 9) = ((xi3 * xi3 - 1) * (xi2 + 1)) / 4;
    shapeFunctionDerivatives(0, 10) = ((xi2 * xi2 - 1) * (xi3 + 1)) / 4;
    shapeFunctionDerivatives(0, 11) = -((xi3 * xi3 - 1) * (xi2 - 1)) / 4;
    shapeFunctionDerivatives(0, 12) = ((xi2 * xi2 - 1) * (xi3 - 1)) / 4;
    shapeFunctionDerivatives(0, 13) = -((xi3 * xi3 - 1) * (xi2 + 1)) / 4;
    shapeFunctionDerivatives(0, 14) = -((xi2 * xi2 - 1) * (xi3 + 1)) / 4;
    shapeFunctionDerivatives(0, 15) = ((xi3 * xi3 - 1) * (xi2 - 1)) / 4;
    shapeFunctionDerivatives(0, 16) = -(xi1 * (xi2 - 1) * (xi3 - 1)) / 2;
    shapeFunctionDerivatives(0, 17) = (xi1 * (xi2 + 1) * (xi3 - 1)) / 2;
    shapeFunctionDerivatives(0, 18) = -(xi1*(xi2 + 1)*(xi3 + 1))/2;
    shapeFunctionDerivatives(0, 19) = (xi1*(xi2 - 1)*(xi3 + 1))/2;

    shapeFunctionDerivatives(1, 0) =((xi1 - 1)*(xi3 - 1)*(xi1 + xi2 + xi3 + 2))/8 + ((xi1 - 1)*(xi2 - 1)*(xi3 - 1))/8;
    shapeFunctionDerivatives(1, 1) =((xi1 - 1)*(xi2 + 1)*(xi3 - 1))/8 - ((xi1 - 1)*(xi3 - 1)*(xi1 - xi2 + xi3 + 2))/8;
    shapeFunctionDerivatives(1, 2) =((xi1 - 1)*(xi3 + 1)*(xi1 - xi2 - xi3 + 2))/8 - ((xi1 - 1)*(xi2 + 1)*(xi3 + 1))/8;
    shapeFunctionDerivatives(1, 3) =- ((xi1 - 1)*(xi3 + 1)*(xi1 + xi2 - xi3 + 2))/8 - ((xi1 - 1)*(xi2 - 1)*(xi3 + 1))/8;
    shapeFunctionDerivatives(1, 4) =- ((xi1 + 1)*(xi3 - 1)*(xi2 - xi1 + xi3 + 2))/8 - ((xi1 + 1)*(xi2 - 1)*(xi3 - 1))/8;
    shapeFunctionDerivatives(1, 5) =- ((xi1 + 1)*(xi3 - 1)*(xi1 + xi2 - xi3 - 2))/8 - ((xi1 + 1)*(xi2 + 1)*(xi3 - 1))/8;
    shapeFunctionDerivatives(1, 6) = ((xi1 + 1)*(xi3 + 1)*(xi1 + xi2 + xi3 - 2))/8 + ((xi1 + 1)*(xi2 + 1)*(xi3 + 1))/8;
    shapeFunctionDerivatives(1, 7) = ((xi1 + 1)*(xi2 - 1)*(xi3 + 1))/8 - ((xi1 + 1)*(xi3 + 1)*(xi1 - xi2 + xi3 - 2))/8;
    shapeFunctionDerivatives(1, 8) = -(xi2*(xi1 - 1)*(xi3 - 1))/2;
    shapeFunctionDerivatives(1, 9) = ((xi3*xi3 - 1)*(xi1 - 1))/4;
    shapeFunctionDerivatives(1, 10) = (xi2*(xi1 - 1)*(xi3 + 1))/2;
    shapeFunctionDerivatives(1, 11) = -((xi3*xi3 - 1)*(xi1 - 1))/4;
    shapeFunctionDerivatives(1, 12) = (xi2*(xi1 + 1)*(xi3 - 1))/2;
    shapeFunctionDerivatives(1, 13) = -((xi3*2 - 1)*(xi1 + 1))/4;
    shapeFunctionDerivatives(1, 14) = -(xi2*(xi1 + 1)*(xi3 + 1))/2;
    shapeFunctionDerivatives(1, 15) = ((xi3*xi3 - 1)*(xi1 + 1))/4;
    shapeFunctionDerivatives(1, 16) = -((xi1*xi1 - 1)*(xi3 - 1))/4;
    shapeFunctionDerivatives(1, 17) = ((xi1*xi1 - 1)*(xi3 - 1))/4;
    shapeFunctionDerivatives(1, 18) = -((xi1*xi1 - 1)*(xi3 + 1))/4;
    shapeFunctionDerivatives(1, 19) = ((xi1*xi1 - 1)*(xi3 + 1))/4;

    shapeFunctionDerivatives(2, 0) =((xi1 - 1)*(xi2 - 1)*(xi1 + xi2 + xi3 + 2))/8 + ((xi1 - 1)*(xi2 - 1)*(xi3 - 1))/8;
    shapeFunctionDerivatives(2, 1) =- ((xi1 - 1)*(xi2 + 1)*(xi1 - xi2 + xi3 + 2))/8 - ((xi1 - 1)*(xi2 + 1)*(xi3 - 1))/8;
    shapeFunctionDerivatives(2, 2) =((xi1 - 1)*(xi2 + 1)*(xi1 - xi2 - xi3 + 2))/8 - ((xi1 - 1)*(xi2 + 1)*(xi3 + 1))/8;
    shapeFunctionDerivatives(2, 3) =((xi1 - 1)*(xi2 - 1)*(xi3 + 1))/8 - ((xi1 - 1)*(xi2 - 1)*(xi1 + xi2 - xi3 + 2))/8;
    shapeFunctionDerivatives(2, 4) =- ((xi1 + 1)*(xi2 - 1)*(xi2 - xi1 + xi3 + 2))/8 - ((xi1 + 1)*(xi2 - 1)*(xi3 - 1))/8;
    shapeFunctionDerivatives(2, 5) =((xi1 + 1)*(xi2 + 1)*(xi3 - 1))/8 - ((xi1 + 1)*(xi2 + 1)*(xi1 + xi2 - xi3 - 2))/8;
    shapeFunctionDerivatives(2, 6) =((xi1 + 1)*(xi2 + 1)*(xi1 + xi2 + xi3 - 2))/8 + ((xi1 + 1)*(xi2 + 1)*(xi3 + 1))/8;
    shapeFunctionDerivatives(2, 7) =- ((xi1 + 1)*(xi2 - 1)*(xi1 - xi2 + xi3 - 2))/8 - ((xi1 + 1)*(xi2 - 1)*(xi3 + 1))/8;
    shapeFunctionDerivatives(2, 8) =-((xi2*xi2 - 1)*(xi1 - 1))/4;
    shapeFunctionDerivatives(2, 9) =(xi3*(xi1 - 1)*(xi2 + 1))/2;
    shapeFunctionDerivatives(2, 10) =((xi2*xi2 - 1)*(xi1 - 1))/4;
    shapeFunctionDerivatives(2, 11) =-(xi3*(xi1 - 1)*(xi2 - 1))/2;
    shapeFunctionDerivatives(2, 12) =((xi2*xi2 - 1)*(xi1 + 1))/4; 
    shapeFunctionDerivatives(2, 13) =-(xi3*(xi1 + 1)*(xi2 + 1))/2;
    shapeFunctionDerivatives(2, 14) =-((xi2*xi2 - 1)*(xi1 + 1))/4;
    shapeFunctionDerivatives(2, 15) =(xi3*(xi1 + 1)*(xi2 - 1))/2;
    shapeFunctionDerivatives(2, 16) =-((xi1*xi1 - 1)*(xi2 - 1))/4;
    shapeFunctionDerivatives(2, 17) =((xi1*xi1 - 1)*(xi2 + 1))/4;
    shapeFunctionDerivatives(2, 18) =-((xi1*xi1 - 1)*(xi2 + 1))/4;
    shapeFunctionDerivatives(2, 19) =((xi1*xi1 - 1)*(xi2 - 1))/4;

    shapeFunctionDerivatives0=shapeFunctionDerivatives.t();
}
