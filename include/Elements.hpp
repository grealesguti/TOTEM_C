#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <vector>
#include <armadillo>

class Elements {
public:
    // Constructor
    Elements();

    // Destructor
    ~Elements();

    // Evaluate shape functions for a linear triangular element with 4 nodes
    arma::vec EvaluateLinearTriangularShapeFunctions(const double xi,const  double eta);
    void CalculateLinearTriangularShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a linear quadrilateral element with 4 nodes
    arma::mat EvaluateLinearQuadrilateralShapeFunctions(double xi, double eta);
    void EvaluateLinearQuadrilateralShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a quadratic quadrilateral element with 8 nodes
    arma::vec EvaluateQuadraticQuadrilateralShapeFunctions(const double xi, const double eta);
    void CalculateQuadraticQuadrilateralShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a linear triangular element with 4 nodes
    arma::vec EvaluateLinearTetrahedraShapeFunctions(const double xi, const double eta, const double zeta);
    void CalculateLinearTetrahedraShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a linear hexahedral element with 8 nodes
    arma::vec EvaluateHexahedralLinearShapeFunctions(const double xi, const double eta, const double zeta);
    void CalculateHexahedralLinearShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a hexahedral serendipity element with 20 nodes
    arma::vec CalculateHexahedralSerendipityShapeFunctions(const double xi, const double eta, const double zetas);
    void CalculateHexahedralSerendipityShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives);

private:
    // Add any private methods or data members needed for shape function calculations
};

#endif // ELEMENTS_H
