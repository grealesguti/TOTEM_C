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
    void EvaluateLinearTriangularShapeFunctions(double xi, double eta, arma::vec& shapeFunctions);
    void CalculateLinearTriangularShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a linear quadrilateral element with 4 nodes
    void EvaluateLinearQuadrilateralShapeFunctions(double xi, double eta, arma::vec& shapeFunctions);
    void EvaluateLinearQuadrilateralShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a quadratic quadrilateral element with 8 nodes
    void EvaluateQuadraticQuadrilateralShapeFunctions(double xi, double eta, arma::vec& shapeFunctions);
    void CalculateQuadraticQuadrilateralShapeFunctionDerivatives(double xi, double eta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a linear triangular element with 4 nodes
    void EvaluateLinearTetrahedraShapeFunctions(double xi, double eta, double zeta, arma::vec& shapeFunctions);
    void CalculateLinearTetrahedraShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a linear hexahedral element with 8 nodes
    void EvaluateHexahedralLinearShapeFunctions(double xi, double eta, double zeta, arma::vec& shapeFunctions);
    void CalculateHexahedralLinearShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives);

    // Evaluate shape functions for a hexahedral serendipity element with 20 nodes
    void CalculateHexahedralSerendipityShapeFunctions(double xi, double eta, double zeta, arma::vec& shapeFunctions);
    void CalculateHexahedralSerendipityShapeFunctionDerivatives(double xi, double eta, double zeta, arma::mat& shapeFunctionDerivatives);

private:
    // Add any private methods or data members needed for shape function calculations
};

#endif // ELEMENTS_H
