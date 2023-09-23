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

        // ------------------------------------------------------------
        // Methods

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
};

// ============================================================

class ArmadilloSpMat : public Data<arma::sp_mat>, public arma::sp_mat
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        ArmadilloSpMat() = default;

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

template <class T>
class Armadillo : public Data<T>, public T
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        Armadillo() = default;
        Armadillo(const T& d);

        // ------------------------------------------------------------
        // Operators

        template <class D>
        inline T operator*(const D& rhs) const
        {
            return static_cast<const T>(*this) * rhs;
            // return *this * rhs;
        }

        inline Armadillo<T> operator-() const {
            -(static_cast<const T>(*this));

            return *this;
        }

        Armadillo<T>& operator=(const Armadillo<T>& other) = default;
        inline Armadillo<T>& operator=(const T& arma)
        {
            T::operator=(arma);
            return *this;
        }

        // ------------------------------------------------------------
        // Methods

        inline const T get() const {
            return *this;
        }

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

template <class T>
inline T operator*(const double rhs, const Armadillo<T>& lhs)
{
    return rhs * static_cast<const T>(lhs);
}

template <class T>
inline T operator*(const T& rhs, const Armadillo<T>& lhs)
{
    return rhs * static_cast<const T>(lhs);
}

template <class T>
class Vector : public Data<std::vector<T>>, public std::vector<T>
{
    public:
        // ------------------------------------------------------------
        // Constructors and destructor

        Vector() = default;
        Vector(const std::vector<T>& data);

    private:
        // ------------------------------------------------------------
        // Private methods

        bool save(std::ofstream& file) const;
};

class EigenDoubleTripletVector : public Data<std::vector<Eigen::Triplet<double>>>, public std::vector<Eigen::Triplet<double>>
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


class EigenSparseMatrix : public Data<Eigen::SparseMatrix<double>>, public Eigen::SparseMatrix<double>
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

class EigenVectorXd : public Data<Eigen::VectorXd>, public Eigen::VectorXd
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
