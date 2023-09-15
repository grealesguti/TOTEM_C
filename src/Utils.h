#ifndef UTILS_H
#define UTILS_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <functional>
#include "Mesh.h"
#include <filesystem> // Include the filesystem header

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
        std::function<IntegrationResult(const arma::mat& natcoords, const arma::mat& coords, const arma::mat& dofs)> func
    );
    arma::mat TransformCoordinates(const arma::mat& cooro);
        // Function to write an Armadillo matrix or vector to a file
    template<typename T>
    static bool writeDataToFile(const T& data, const std::string& filename, bool append = false); // Added 'bool append' parameter
    // Function to delete all files in a given folder path
    static bool deleteFilesInFolder(const std::string& folderPath) {
        try {
            std::filesystem::directory_iterator it(folderPath);
            for (const auto& entry : it) {
                if (std::filesystem::is_regular_file(entry)) {
                    std::filesystem::remove(entry.path());
                }
            }
            return true; // Deletion successful
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return false; // Deletion failed
        }
    }
    static arma::mat calculate_T3(const arma::mat& nodes);
    arma::mat calculate_inverse_T3(const arma::mat& nodes);

private:
bool writeflag_;

};

#endif // UTILS_H
