// Class header
#include "utils/data.hpp"

// STL
#include <iostream>


// ------------------------------------------------------------

template <class T>
bool Data<T>::writeDataToFile(const std::string& filename, const bool append) const {
    // Open the file for writing
    std::ofstream file;

    if (append) {
        // Open the file in append mode
        file.open(filename, std::ios::app);
    } else {
        // Open the file in truncation mode (default)
        file.open(filename);
    }

    // Check if the file opened successfully
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return false;
    }

    // Print "New Data" to the file
    file << "New Data:" << std::endl;

    // Write the data to the file
    const bool saved = save(file);

    // Close the file
    file.close();

    // Check if the write was successful
    if (!file.good()) {
        std::cerr << "Error: Failed to write data to the file." << std::endl;
        return false;
    }

    return saved;
}

template class Data<arma::mat>;
template class Data<arma::vec>;
template class Data<arma::uvec>;
template class Data<arma::umat>;

template class Data<std::vector<double>>;
template class Data<std::vector<int>>;

template class Data<std::vector<Eigen::Triplet<double>>>;
template class Data<Eigen::SparseMatrix<double>>;
template class Data<Eigen::VectorXd>;

//////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------
// ============================================================

bool ArmadilloSpMat::save(std::ofstream& file) const {
    std::cerr << "Error: Writing sparse matrices in ASCII format is not supported." << std::endl;

    return false;
}

// ------------------------------------------------------------
// ============================================================

template <class T>
Armadillo<T>::Armadillo(const T& arma)
    : T(arma)
{};

template <class T>
bool Armadillo<T>::save(std::ofstream& file) const {
    // For Armadillo dense data (vector, matrix, or uvec)
    T::save(file, arma::raw_ascii);

    return true;
}

template class Armadillo<arma::vec>;
template class Armadillo<arma::uvec>;
template class Armadillo<arma::mat>;
template class Armadillo<arma::umat>;

// ------------------------------------------------------------
// ============================================================

template <class T>
Vector<T>::Vector(const std::vector<T>& data)
    : std::vector<T>(data)
{}

template <class T>
bool Vector<T>::save(std::ofstream& file) const {
    for (const auto& item : *this) {
        file << item << " ";
    }

    return true;
}

template class Vector<double>;
template class Vector<int>;

// ------------------------------------------------------------
// ============================================================

EigenDoubleTripletVector::EigenDoubleTripletVector(const std::vector<Eigen::Triplet<double>>& data)
    : std::vector<Eigen::Triplet<double>>(data)
{}

bool EigenDoubleTripletVector::save(std::ofstream& file) const {
    // For std::vector<Eigen::Triplet<double>> (sparse matrix in triplet format)
    for (const auto& triplet : *this) {
        file << triplet.row() << " " << triplet.col() << " " << triplet.value() << std::endl;
    }

    return true;
}

// ------------------------------------------------------------
// ============================================================

EigenSparseMatrix::EigenSparseMatrix(const Eigen::SparseMatrix<double>& data)
    : Eigen::SparseMatrix<double>(data)
{}

bool EigenSparseMatrix::save(std::ofstream& file) const {
    // For Eigen sparse matrix (you may need to adjust this based on the actual Eigen sparse matrix type you are using)
    // Write Eigen sparse matrix to file in a custom format
    for (int i = 0; i < this->rows(); ++i) {
        for (int j = 0; j < this->cols(); ++j) {
            file << this->coeff(i, j) << " ";
        }
        file << std::endl;
    }

    return true;
}

// ------------------------------------------------------------
// ============================================================

EigenVectorXd::EigenVectorXd(const Eigen::VectorXd& data)
    : Eigen::VectorXd(data)
{}
bool EigenVectorXd::save(std::ofstream& file) const {
    // For Eigen::VectorXd
    for (int i = 0; i < this->size(); ++i) {
        file << (*this)(i) << " ";
    }

    return true;
}

// ------------------------------------------------------------
// ============================================================
///////////////////////////////////////////////////////////////

void deleteFilesInFolder(const std::string& folderPath) {
    try {
        std::filesystem::directory_iterator it(folderPath);

        for (const auto& entry : it) {
            if (std::filesystem::is_regular_file(entry)) {
                std::filesystem::remove(entry.path());
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}