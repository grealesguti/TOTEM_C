#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <vector>
#include <armadillo>


// Evaluate shape functions for a linear triangular element with 4 nodes
arma::vec EvaluateLinearTriangularShapeFunctions(const double xi,const  double eta);
arma::mat CalculateLinearTriangularShapeFunctionDerivatives();

// Evaluate shape functions for a linear quadrilateral element with 4 nodes
arma::mat EvaluateLinearQuadrilateralShapeFunctions(double xi, double eta);
arma::mat EvaluateLinearQuadrilateralShapeFunctionDerivatives(const double xi, const double eta);

// Evaluate shape functions for a quadratic quadrilateral element with 8 nodes
arma::vec EvaluateQuadraticQuadrilateralShapeFunctions(const double xi, const double eta);
arma::mat CalculateQuadraticQuadrilateralShapeFunctionDerivatives(const double xi, const double eta);

// Evaluate shape functions for a linear triangular element with 4 nodes
arma::vec EvaluateLinearTetrahedraShapeFunctions(const double xi, const double eta, const double zeta);
arma::mat CalculateLinearTetrahedraShapeFunctionDerivatives();

// Evaluate shape functions for a linear hexahedral element with 8 nodes
arma::vec EvaluateHexahedralLinearShapeFunctions(const double xi, const double eta, const double zeta);
arma::mat CalculateHexahedralLinearShapeFunctionDerivatives(const double xi, const double eta, const double zeta);

// Evaluate shape functions for a hexahedral serendipity element with 20 nodes
arma::vec CalculateHexahedralSerendipityShapeFunctions(const double xi, const double eta, const double zetas);
arma::mat CalculateHexahedralSerendipityShapeFunctionDerivatives(const double xi, const double eta, const double zeta);


#endif // ELEMENTS_H
