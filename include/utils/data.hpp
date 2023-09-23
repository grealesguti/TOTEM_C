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

        Data() = default;
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
        bool writeDataToFile(const std::string& filename, const bool append = false) const;

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

        ArmadilloSpMat() = default;
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

        Armadillo() = default;
        Armadillo(const T& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

template <class T>
class ArmadilloVector : public Armadillo<arma::Col<T>>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        ArmadilloVector() = default;
        ArmadilloVector(const arma::Col<T>& data);

        // ------------------------------------------------------------
        // Operators

        /**
         * @brief Access the element stored at index.
         *
         * @param[in] index
         *
         * @throw exception if the requested element is out of bounds.
        */
        const T operator()(const std::size_t index) const;
        T & operator()(const std::size_t index);

        // ------------------------------------------------------------
        // Methods

        const std::size_t size() const;
};

template <class T>
class ArmadilloMatrix : public Armadillo<arma::Mat<T>>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        ArmadilloMatrix() = default;
        ArmadilloMatrix(const arma::Mat<T>& data);

        // ------------------------------------------------------------
        // Operators

        /**
         * @brief Access the element/object stored at index i under the assumption of a flat layout, with column-major ordering of data (ie. column by column).
         *
         * @param[in] index
         *
         * @throw exception if the requested element is out of bounds.
        */
        const T operator()(const std::size_t index) const;
        T & operator()(const std::size_t index);

        /**
         * @brief Access the element/object stored at row @p r and column @p c.
         *
         * @param[in] r Row.
         * @param[in] c Column.
         *
         * @throw exception if the requested element is out of bounds.
        */
        const T operator()(const std::size_t r, const std::size_t c) const;
        T & operator()(const std::size_t r, const std::size_t c);


        /**
         * @brief Creation of subview (column vector).
         *
         * @param[in] c Column.
         *
         * @throw exception if the requested element is out of bounds.
        */
        const arma::subview_col<T> col(const std::size_t c) const;
        arma::subview_col<T> col(const std::size_t c);
};

template <class T>
class Vector : public Data<std::vector<T>>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        Vector() = default;
        Vector(const std::vector<T>& data);

        // ------------------------------------------------------------
        // Methods

        const std::size_t size() const;

        void assign(std::vector<std::size_t>::iterator first, std::vector<std::size_t>::iterator last);

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

        EigenDoubleTripletVector() = default;
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

        EigenSparseMatrix() = default;
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

        EigenVectorXd() = default;
        EigenVectorXd(const Eigen::VectorXd& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

#endif // DATA_HPP
