#ifndef UTILS_H
#define UTILS_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <functional>
#include "Mesh.hpp"
#include <filesystem> // Include the filesystem header
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Utils {
public:
    Utils(bool someFlag = false); // Added 'bool someFlag' parameter with a default value of 'false'
    struct IntegrationResult {
        arma::mat KT;
        arma::mat R;
    };
    void getGaussWeightsAndPoints(int order, arma::mat& weights, arma::mat& gaussPoints);
    void gaussIntegrationBC(int dimension, int order, int elementTag, Mesh mesh, double bcvalue, std::function<arma::mat(const arma::mat&,const arma::mat&, double, int)> func, arma::mat& result);
    // Function to perform Gaussian integration
    IntegrationResult gaussIntegrationK(
        int dimension,
        int order,
        int elementTag,
        Mesh mesh,
        arma::vec element_dof_values,
        std::function<IntegrationResult(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag)> func
    );
    arma::mat TransformCoordinates(const arma::mat& cooro);
        // Function to write an Armadillo matrix or vector to a file
    template<typename T>
    static bool writeDataToFile(const T& data, const std::string& filename, bool append = false); // Added 'bool append' parameter
    // Function to delete all files in a given folder path
    void deleteFilesInFolder(const std::string& folderPath);
    static arma::mat calculate_T3(const arma::mat& nodes);
    arma::mat calculate_inverse_T3(const arma::mat& nodes);
    static arma::sp_mat spmat_submat(const arma::sp_mat& spmatrix, const std::vector<int>& row_indices, const std::vector<int>& col_indices);

private:
bool writeflag_;

};

#endif // UTILS_H
