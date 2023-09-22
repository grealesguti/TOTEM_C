// Class header
#include "utils/data.hpp"

// STL
#include <iostream>


// ------------------------------------------------------------

template <class T>
Data<T>::Data(const T& data)
    : data_(data)
{}

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

template <class T>
const T Data<T>::getData() const {
    return data_;
}

template <class T>
void Data<T>::setData(const T& data) {
    data_ = data;
}

template class Data<arma::mat>;
template class Data<arma::vec>;
template class Data<arma::uvec>;
template class Data<arma::umat>;

template class Data<std::vector<double>>;
template class Data<std::vector<int>>;

//////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------
// ============================================================

ArmadilloSpMat::ArmadilloSpMat(const arma::sp_mat& data)
    : Data<arma::sp_mat>(data)
{}

bool ArmadilloSpMat::save(std::ofstream& file) const {
    std::cerr << "Error: Writing sparse matrices in ASCII format is not supported." << std::endl;

    return false;
}

// ------------------------------------------------------------
// ============================================================

template <class T>
Armadillo<T>::Armadillo(const T& data)
    : Data<T>(data)
{}

template <class T>
bool Armadillo<T>::save(std::ofstream& file) const {
    // For Armadillo dense data (vector, matrix, or uvec)
    this->data_.save(file, arma::raw_ascii);

    return true;
}

template class Armadillo<arma::mat>;
template class Armadillo<arma::vec>;
template class Armadillo<arma::uvec>;
template class Armadillo<arma::umat>;

// ------------------------------------------------------------
// ============================================================

template <class T>
Vector<T>::Vector(const std::vector<T>& data)
    : Data<std::vector<T>>(data)
{}

template <class T>
bool Vector<T>::save(std::ofstream& file) const {
    for (const auto& item : this->data_) {
        file << item << " ";
    }

    return true;
}

template <class T>
const std::size_t Vector<T>::size() const {
    return this->data_.size();
}

template <class T>
void Vector<T>::assign(std::vector<std::size_t>::iterator first, std::vector<std::size_t>::iterator last){
    this->data_.assign(first, last);
}

template class Vector<double>;
template class Vector<int>;

// ------------------------------------------------------------
// ============================================================

EigenDoubleTripletVector::EigenDoubleTripletVector(const std::vector<Eigen::Triplet<double>>& data)
    : Data<std::vector<Eigen::Triplet<double>>>(data)
{}

bool EigenDoubleTripletVector::save(std::ofstream& file) const {
    // For std::vector<Eigen::Triplet<double>> (sparse matrix in triplet format)
    for (const auto& triplet : data_) {
        file << triplet.row() << " " << triplet.col() << " " << triplet.value() << std::endl;
    }

    return true;
}

// ------------------------------------------------------------
// ============================================================

EigenSparseMatrix::EigenSparseMatrix(const Eigen::SparseMatrix<double>& data)
    : Data<Eigen::SparseMatrix<double>>(data)
{}

bool EigenSparseMatrix::save(std::ofstream& file) const {
    // For Eigen sparse matrix (you may need to adjust this based on the actual Eigen sparse matrix type you are using)
    // Write Eigen sparse matrix to file in a custom format
    for (int i = 0; i < data_.rows(); ++i) {
        for (int j = 0; j < data_.cols(); ++j) {
            file << data_.coeff(i, j) << " ";
        }
        file << std::endl;
    }

    return true;
}

// ------------------------------------------------------------
// ============================================================

EigenVectorXd::EigenVectorXd(const Eigen::VectorXd& data)
    : Data<Eigen::VectorXd>(data)
{}

bool EigenVectorXd::save(std::ofstream& file) const {
    // For Eigen::VectorXd
    for (int i = 0; i < data_.size(); ++i) {
        file << data_(i) << " ";
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