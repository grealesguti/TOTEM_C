#ifndef DATA_HPP
#define DATA_HPP

// ------------------------------------------------------------
// Includes

// STL
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>

// Armadillo
#include <armadillo>

// Eigen
#include <eigen3/Eigen/Sparse>

// ------------------------------------------------------------
// Functions

/**
 * @brief Delete all files in a given folder path
 * 
 * @param[in] folderPath Path to the folder to be deleted.
 */
void deleteFilesInFolder(const std::string& folderPath);

// ------------------------------------------------------------
// Classes

template <class T>
class Data
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        /**
         * @brief Constructor.
         * 
         * @param[in] data
         */
        Data(const T& data);

        // ------------------------------------------------------------
        // Methods

        const T getData() const;
        void setData(const T& data);

        /**
         * @brief Function to write an Armadillo matrix or vector to a file
         * 
         * @param[in] filename File where the data will be stored.
         * @param[in] append If true, the file will be opened in append mode, otherwise in truncation mode.
         * 
         * @return True if the data was written to the file, false if not.
         */
        bool writeToFile(const std::string& filename, const bool append = false) const;

    protected:
        // ------------------------------------------------------------
        // Protected methods

        virtual bool save(std::ofstream& file) const = 0;

        // ------------------------------------------------------------
        // Protected attributes

        T data_; //!< Data
};

// ============================================================

class ArmadilloSpMat : public Data<arma::sp_mat>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        ArmadilloSpMat(const arma::sp_mat& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

template <class T>
class Armadillo : public Data<T>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        Armadillo(const T& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

template <class T>
class Vector : public Data<std::vector<T>>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        Vector(const std::vector<T>& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

class EigenDoubleTripletVector : public Data<std::vector<Eigen::Triplet<double>>>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        EigenDoubleTripletVector(const std::vector<Eigen::Triplet<double>>& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};


class EigenSparseMatrix : public Data<Eigen::SparseMatrix<double>>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        EigenSparseMatrix(const Eigen::SparseMatrix<double>& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

class EigenVectorXd : public Data<Eigen::VectorXd>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        EigenVectorXd(const Eigen::VectorXd& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

#endif // DATA_HPP
