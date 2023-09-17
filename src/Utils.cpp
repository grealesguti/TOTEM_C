#include "Utils.h"

using namespace arma;

Utils::Utils(bool someFlag): writeflag_(someFlag) {

}

void Utils::getGaussWeightsAndPoints(int order, mat& weights, mat& gaussPoints) {
    // Check if the order is valid
    if (order < 1 || order > 14) {
        std::cerr << "Invalid order for Gauss integration." << std::endl;
        return;
    }

    // Initialize the matrices with zeros
    weights = mat(order, 1, fill::zeros);
    gaussPoints = mat(order, 1, fill::zeros);

    if (order == 1) {
        weights(0, 0) = 2.0;
        gaussPoints(0, 0) = 0.0;
    } else if (order == 2) {
        weights(0, 0) = 1.0;
        weights(1, 0) = 1.0;
        gaussPoints(0, 0) = -0.577350269189626;
        gaussPoints(1, 0) = 0.577350269189626;
    } else if (order == 3) {
        weights(0, 0) = 0.555555555555556;
        weights(1, 0) = 0.888888888888889;
        weights(2, 0) = 0.555555555555556;
        gaussPoints(0, 0) = -0.774596669241483;
        gaussPoints(1, 0) = 0.0;
        gaussPoints(2, 0) = 0.774596669241483;
    } else if (order == 14) {
        // 14-point Gauss integration rule.
        double cor = 0.335180055401662;
        double cen = 0.886426592797784;
        
        mat G14 = {
            {-cor, -cor, -cor, -cor, +cor, +cor, +cor, +cor, -cen, cen, 0.0, 0.0, 0.0, 0.0},
            {-cor, -cor, +cor, +cor, -cor, -cor, +cor, +cor, 0.0, 0.0, -cen, cen, 0.0, 0.0},
            {-cor, +cor, -cor, +cor, -cor, +cor, -cor, +cor, 0.0, 0.0, 0.0, 0.0, -cen, cen}
        };
        
        vec W14(order);
        W14.subvec(0, 7).fill(cor);
        W14.subvec(8, 13).fill(cen);

        weights = W14;
        gaussPoints = G14;
    }
}
#include <iostream>
#include <armadillo>

arma::mat Utils::TransformCoordinates(const arma::mat& cooro) {
    // Extract columns v1 and v2 vectors
    arma::colvec v1 = cooro.col(1) - cooro.col(0);
    arma::colvec v2 = cooro.col(2) - cooro.col(0);

    // Calculate thetay1 in radians and degrees
    double thetay1 = -std::acos(arma::dot(v1, v2) / (arma::norm(v1) * arma::norm(v2)));
    double thetay1deg = thetay1 * 180 / arma::datum::pi;

    arma::colvec nxy;

    if (v1.n_elem == 3 && v2.n_elem == 3) {
        // Calculate nxy
        nxy = arma::cross(v1, v2) / arma::norm(arma::cross(v1, v2));
    } else {
        // Handle the case where the input vectors are not valid
        std::cerr << "Utils::TransformCoordinates: Invalid input vectors for cross product." << std::endl;
        // Print the input vectors and their sizes
        std::cout << "Vector v1: " << v1 << std::endl;
        std::cout << "Vector v2: " << v2 << std::endl;
        std::cout << "v1 size: " << v1.n_elem << std::endl;
        std::cout << "v2 size: " << v2.n_elem << std::endl;
        return arma::mat(); // Return an empty matrix to indicate an error
    }
      //  std::cout << "Cross product done" << std::endl;

    // Calculate thetaz0 in radians and degrees
    double thetaz0 = std::acos(nxy(0));
    double thetaz0deg = thetaz0 * 180 / arma::datum::pi;
    //std::cout << "angle calculated" << std::endl;

    // Check if thetaz0 is greater than pi/2 and perform necessary adjustments
    if (thetaz0 > arma::datum::pi / 2) {
      //  std::cout << "angle adjustment" << std::endl;
        double thetaz0old = thetaz0;
        double thetaz0degold = thetaz0deg;
        thetaz0 = arma::datum::pi - thetaz0;
        thetaz0deg = thetaz0 * 180 / arma::datum::pi;
        thetay1 = -thetay1;
    }
    
    //std::cout << "create transformation matrixes" << std::endl;
    // Create transformation matrices Gz0 and Gy1
    arma::mat Gz0 = {{std::cos(thetaz0), -std::sin(thetaz0), 0},
                     {std::sin(thetaz0), std::cos(thetaz0), 0},
                     {0, 0, 1}};

    arma::mat Gy1 = {{std::cos(thetay1), 0, -std::sin(thetay1)},
                     {0, 1, 0},
                     {std::sin(thetay1), 0, std::cos(thetay1)}};

    // Calculate the transformation matrix G
    //std::cout << "create full transformation" << std::endl;
    arma::mat G = Gz0 * Gy1;
    //std::cout << "set to zeros" << std::endl;
    // Set elements of G close to zero to zer    std::cout << "create full transformation" << std::endl;o
    for (arma::uword i = 0; i < G.n_rows; ++i) {
        for (arma::uword j = 0; j < G.n_cols; ++j) {
            if (std::abs(G(i, j)) < 1e-12) {
                G(i, j) = 0.0;
            }
        }
    }
    //std::cout << "do the transformation" << std::endl;
    // Perform the final coordinate transformation
    arma::mat cooro1 = G * cooro;
    //std::cout << "Coordinate transformation done" << std::endl;

    return cooro1;
}


void Utils::gaussIntegrationBC(int dimension, int order, int elementTag, Mesh mesh, double bcvalue, std::function<mat(const mat& natcoords,const mat& coords, double value, int element)> func, mat& result) {
    
    if (dimension < 1 || order < 1) {
        std::cerr << "Invalid dimension or order for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    mat weights;
    mat gaussPoints;
    std::vector<int> nodeTags_el;
    mesh.getElementInfo(elementTag, nodeTags_el);
    arma::mat coordinates(3, size(nodeTags_el));
    arma::mat coordinates_tr(3, size(nodeTags_el));


    coordinates=mesh.getCoordinates(nodeTags_el);
    arma::mat iT3= Utils::calculate_inverse_T3(coordinates);
    coordinates_tr=iT3*coordinates;
    // Create a new matrix to store the result without the last row
    arma::mat coordinates_tr_XY;
    double tolerance = 1e-6; // Define your tolerance here


    // Iterate through the rows of coordinates_tr and copy rows that do not meet the tolerance condition
    for (arma::uword i = 0; i < coordinates_tr.n_rows; ++i) {
        bool removeRow = true; // Assume we want to remove the row by default
        double firstValue = coordinates_tr(i, 0); // Store the first value in the row

        // Check if all values in the row are within the tolerance
        for (arma::uword j = 1; j < coordinates_tr.n_cols; ++j) {
            if (std::abs(coordinates_tr(i, j) - firstValue) > tolerance) {
                removeRow = false; // Row contains a different value outside the tolerance, so we keep it
                break;
            }
        }

        if (!removeRow) {
            // Copy the row to coordinates_tr_XY
            if (coordinates_tr_XY.is_empty()) {
                coordinates_tr_XY = coordinates_tr.row(i);
            } else {
                coordinates_tr_XY = arma::join_vert(coordinates_tr_XY, coordinates_tr.row(i));
            }
        }
    }

    if (coordinates_tr_XY.is_empty()) {
        // Handle the case when all rows were removed, resulting in an empty matrix.
        std::cout << "Warning: All rows removed. coordinates_tr_XY is empty." << std::endl;
        // You may want to print a warning message or take appropriate action.
    }

    // Get the Gauss weights and points for the specified order.
    Utils::getGaussWeightsAndPoints(order, weights, gaussPoints);
    //std::cout << "Retrieve gauss integration points" << std::endl;
    if(writeflag_==true){
        Utils::writeDataToFile(nodeTags_el,"Outputs/BCNodeTags_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(coordinates,"Outputs/BCGaussCoords_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(coordinates_tr,"Outputs/BCGaussCoordsTr_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(coordinates_tr_XY,"Outputs/BCGaussCoordsTrXY_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(gaussPoints,"Outputs/BCGaussPoints_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(weights,"Outputs/BCGaussWeights_"+std::to_string(elementTag)+".txt",true);
    }

    if (weights.is_empty() || gaussPoints.is_empty()) {
        std::cerr << "Invalid order for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    if (weights.n_rows != gaussPoints.n_rows) {
        std::cerr << "Weights and Gauss points have mismatched dimensions." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
        return;
    }

    // Initialize result to zero matrix of appropriate dimensions.
    // result = zeros<mat>(dimension, dimension); Should already be of the right dimension??, in this case 8,1?? make a check that it is the right dimension
    //std::cout << "Start the integration " ;

    if (dimension == 1) {
        // 1D integration using a single loop.
        //std::cout << "of dimension 1" << std::endl;
            arma::mat natcoords(1, 1);
            arma::mat f(2, 1);
        for (uword i = 0; i < weights.n_rows; ++i) {
                // Explicitly use the arma::operator* function for multiplication
                f = arma::operator*(func(natcoords, coordinates_tr_XY, bcvalue,elementTag), weights(i, 0));
                result += f;                }
    } else if (dimension == 2) {
        //std::cout << "of dimension 2" << std::endl;
            arma::mat natcoords(2, 1);
            arma::mat f(4, 1);
        // 2D integration using a double loop.
        for (uword i = 0; i < weights.n_rows; ++i) {
            for (uword j = 0; j < weights.n_rows; ++j) {
                        natcoords(0, 0) = gaussPoints(i, 0);
                        natcoords(1, 0) = gaussPoints(j, 0);
                        // Create a 3x1 matrix with gaussPoints(i, 0), gaussPoints(j, 0), and gaussPoints(k, 0)
                        // Explicitly use the arma::operator* function for multiplication
                        f = arma::operator*(func(natcoords, coordinates_tr_XY, bcvalue,elementTag), weights(i, 0) * weights(j, 0));
                        result += f;            
             }
        }
    } else if (dimension == 3) {
        if (order == 14) {
            //std::cout << "of dimension 3 and reduced order" << std::endl;
            // Special case for 3D integration with order 14.
            // Directly use the given points and weights without looping.
            result = weights * weights.t() % weights * weights.t() % weights * weights.t() % func(gaussPoints,coordinates_tr_XY,bcvalue,elementTag);
        } else {
            arma::mat natcoords(3, 1);
            arma::mat f(8, 1);
            //std::cout << "of dimension 3" << std::endl;
            // Generic 3D integration using a triple loop.
            for (uword i = 0; i < weights.n_rows; ++i) {
                for (uword j = 0; j < weights.n_rows; ++j) {
                    for (uword k = 0; k < weights.n_rows; ++k) {
                        natcoords(0, 0) = gaussPoints(i, 0);
                        natcoords(1, 0) = gaussPoints(j, 0);
                        natcoords(2, 0) = gaussPoints(k, 0);
                        // Create a 3x1 matrix with gaussPoints(i, 0), gaussPoints(j, 0), and gaussPoints(k, 0)
                        // Explicitly use the arma::operator* function for multiplication
                        f = arma::operator*(func(natcoords, coordinates_tr_XY, bcvalue,elementTag), weights(i, 0) * weights(j, 0) * weights(k, 0));
                        result += f;
                    }
                }
            }
        }
    } else {
        std::cerr << "Invalid dimension for Gauss integration." << std::endl;
        result = zeros<mat>(1, 1); // Initialize result to a 1x1 matrix with zero value.
    }
}


//////////////////////////////////////////////////////////////////////////////////
Utils::IntegrationResult Utils::gaussIntegrationK(
    int dimension,
    int order,
    int elementTag,
    Mesh mesh,
    arma::vec element_dof_values,
    std::function<Utils::IntegrationResult(const arma::mat& natcoords, const arma::mat& coords, const arma::vec& dofs, const int elementTag)> func
) {    
    Utils::IntegrationResult result; // Create a struct to hold KT and R

    if (dimension < 1 || order < 1) {
        std::cerr << "Invalid dimension or order for Gauss integration." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
        return result;
    }
    int nodes_per_element = 8; // Total number of nodes per element
    int dof_per_node = 2; // Assuming 2 degrees of freedom per node
    arma::mat weights;
    arma::mat gaussPoints;
    std::vector<int> nodeTags_el;
    mesh.getElementInfo(elementTag, nodeTags_el);
    arma::mat coordinates(3, nodeTags_el.size());
    coordinates=mesh.getCoordinates(nodeTags_el);

    // get current dofs from elementdofs index, this is the input!!!
    this->getGaussWeightsAndPoints(order, weights, gaussPoints);

    if (weights.is_empty() || gaussPoints.is_empty()) {
        std::cerr << "Invalid order for Gauss integration." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
        return result;
    }

    if (weights.n_rows != gaussPoints.n_rows) {
        std::cerr << "Weights and Gauss points have mismatched dimensions." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
        return result;
    }
    if(writeflag_==true){
        Utils::writeDataToFile(nodeTags_el,"Outputs/KTNodeTags_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(coordinates,"Outputs/KTGaussCoords_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(gaussPoints,"Outputs/KTGaussPoints_"+std::to_string(elementTag)+".txt",true);
        Utils::writeDataToFile(weights,"Outputs/KTGaussWeights_"+std::to_string(elementTag)+".txt",true);
    }
    // Initialize result.KT and result.R to zero matrices of appropriate dimensions.
    if (dimension == 1) {
        // 1D integration using a single loop.
        arma::mat natcoords(1, 1);
        arma::mat Re(4, 1, arma::fill::zeros);
        arma::mat KTe(4, 4, arma::fill::zeros);

        for (uword i = 0; i < weights.n_rows; ++i) {
            natcoords(0, 0) = gaussPoints(i, 0);
            Utils::IntegrationResult localResult = func(natcoords, coordinates,element_dof_values,elementTag );
            Re += localResult.R * weights(i, 0);
            KTe += localResult.KT * weights(i, 0);
        }

            result.KT = KTe;
            result.R = Re;

    } else if (dimension == 2) {
        // 2D integration using a double loop.
        arma::mat natcoords(2, 1);
        arma::mat R(8, 1, arma::fill::zeros);
        arma::mat KT(8, 8, arma::fill::zeros);

        for (uword i = 0; i < weights.n_rows; ++i) {
            for (uword j = 0; j < weights.n_rows; ++j) {
                natcoords(0, 0) = gaussPoints(i, 0);
                natcoords(1, 0) = gaussPoints(j, 0);
                Utils::IntegrationResult localResult = func(natcoords, coordinates,element_dof_values,elementTag);
                R += localResult.R * weights(i, 0) * weights(j, 0);
                KT += localResult.KT * weights(i, 0) * weights(j, 0);
            }
        }

            result.KT = KT;
            result.R = R;

    } else if (dimension == 3) {
        if (order == 14) {
            // Special case for 3D integration with order 14.
            // Directly use the given points and weights without looping.
            //result.KT = weights * weights.t() % weights * weights.t() % weights * weights.t() * func(gaussPoints, coordinates_tr,element_dofs);
            //result.R = arma::zeros<arma::mat>(16, 1); // Initialize R to a zero matrix.
        } else {
            arma::mat natcoords(3, 1);
            arma::mat R(16, 1, arma::fill::zeros);
            arma::mat KT(16, 16, arma::fill::zeros);

            for (uword i = 0; i < weights.n_rows; ++i) {
                for (uword j = 0; j < weights.n_rows; ++j) {
                    for (uword k = 0; k < weights.n_rows; ++k) {
                        natcoords(0, 0) = gaussPoints(i, 0);
                        natcoords(1, 0) = gaussPoints(j, 0);
                        natcoords(2, 0) = gaussPoints(k, 0);
                        Utils::IntegrationResult localResult = func(natcoords, coordinates,element_dof_values,elementTag);
                        R += localResult.R * weights(i, 0) * weights(j, 0) * weights(k, 0);
                        KT += localResult.KT * weights(i, 0) * weights(j, 0) * weights(k, 0);
                    }
                }
            }

            result.KT = KT;
            result.R = R;
        }
    } else {
        std::cerr << "Invalid dimension for Gauss integration." << std::endl;
        result.KT = arma::zeros<arma::mat>(1, 1); // Initialize KT to a 1x1 matrix with zero value.
        result.R = arma::zeros<arma::mat>(1, 1);  // Initialize R to a 1x1 matrix with zero value.
    }

    return result; // Return the struct containing KT and R.
}

//////////////////////////////////////////////////////////////////////////////////
template <typename T>
bool Utils::writeDataToFile(const T& data, const std::string& filename, bool append) {
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
    // Print "new line" to the file
    file << "New Data:" << std::endl;
    // Write the data to the file
    if constexpr (std::is_same_v<T, arma::mat> || std::is_same_v<T, arma::vec> || std::is_same_v<T, arma::uvec>|| std::is_same_v<T, arma::umat>) {
        // For Armadillo dense data (vector, matrix, or uvec)
        data.save(file, arma::raw_ascii);
    } else if constexpr (std::is_same_v<T, arma::sp_mat>) {
        // For Armadillo sparse matrix (CSR format)
        std::cerr << "Error: Writing sparse matrices in ASCII format is not supported." << std::endl;
        return false;
    } else if constexpr (std::is_same_v<T, std::vector<double>> || std::is_same_v<T, std::vector<int>>) {
        // For double or int vectors
        for (const auto& item : data) {
            file << item << " ";
        }
    } else if constexpr (std::is_same_v<T, std::vector<Eigen::Triplet<double>>>) {
            // For std::vector<Eigen::Triplet<double>> (sparse matrix in triplet format)
            for (const auto& triplet : data) {
                file << triplet.row() << " " << triplet.col() << " " << triplet.value() << std::endl;
            }
    } else if constexpr (std::is_same_v<T, Eigen::SparseMatrix<double>>) {
        // For Eigen sparse matrix (you may need to adjust this based on the actual Eigen sparse matrix type you are using)
        // Write Eigen sparse matrix to file in a custom format
        for (int i = 0; i < data.rows(); ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                file << data.coeff(i, j) << " ";
            }
            file << std::endl;
        }
    }    else {
        // Handle unsupported types here
        std::cerr << "Error: Unsupported data type." << std::endl;
        return false;
    }

    // Close the file
    file.close();

    // Check if the write was successful
    if (!file.good()) {
        std::cerr << "Error: Failed to write data to the file." << std::endl;
        return false;
    }

    return true;
}

// Explicit template specialization for arma::mat
template bool Utils::writeDataToFile(const arma::mat& data, const std::string& filename, bool append);

// Explicit template specialization for arma::vec
template bool Utils::writeDataToFile(const arma::vec& data, const std::string& filename, bool append);

// Explicit template specialization for arma::sp_mat (sparse matrix)
template bool Utils::writeDataToFile(const arma::sp_mat& data, const std::string& filename, bool append);

// Explicit template specialization for std::vector<double>
template bool Utils::writeDataToFile(const std::vector<double>& data, const std::string& filename, bool append);

// Explicit template specialization for std::vector<int>
template bool Utils::writeDataToFile(const std::vector<int>& data, const std::string& filename, bool append);

// Explicit template specialization for arma::uvec
template bool Utils::writeDataToFile(const arma::uvec& data, const std::string& filename, bool append);

// Explicit template specialization for arma::umat
template bool Utils::writeDataToFile(const arma::umat& data, const std::string& filename, bool append);

// Explicit template specialization for arma::umat
template bool Utils::writeDataToFile<Eigen::SparseMatrix<double>>(const Eigen::SparseMatrix<double>& data, const std::string& filename, bool append);

// Explicit template specialization for std::vector<Eigen::Triplet<double>>
template bool Utils::writeDataToFile(const std::vector<Eigen::Triplet<double>>& data, const std::string& filename, bool append);

/////////////////////////////////////////////////////////////////////////////////////
arma::mat Utils::calculate_T3(const arma::mat& nodes) { // the input is a 3x4 matrix containing the coords of each node in the columns.
// LINK: http://what-when-how.com/the-finite-element-method/fem-for-frames-finite-element-method-part-2/
    if (nodes.n_rows != 3 || nodes.n_cols != 4) {
        std::cerr << "Input matrix must be 3x4." << std::endl;
        return arma::mat();  // Return an empty matrix to indicate an error
    }

    // Extract node coordinates
    arma::vec V1 = nodes.col(0);
    arma::vec V2 = nodes.col(1);
    arma::vec V3 = nodes.col(2);

    // Calculate the unit normal vector XN
    arma::vec V12 = V2 - V1;
    arma::vec V31 = V3 - V1;
    arma::vec XN = V12 / arma::norm(V12);
    double lx = XN(0);
    double mx = XN(1);
    double nx = XN(2);

    // Calculate the unit normal vector ZN
    arma::vec ZN = arma::cross(V12, V31) / arma::norm(arma::cross(V12, V31));
    double lz = ZN(0);
    double mz = ZN(1);
    double nz = ZN(2);

    // Calculate ly, my, and ny
    double ly = mz * nx - nz * mx;
    double my = nz * lx - lz * nx;
    double ny = lz * mx - mz * lx;

    // Create and return the transformation matrix T3
    arma::mat T3(3, 3);
    T3 << lx << mx << nx << arma::endr
       << ly << my << ny << arma::endr
       << lz << mz << nz << arma::endr;

    return T3;
}

arma::mat Utils::calculate_inverse_T3(const arma::mat& nodes) {
    // Calculate T3
    arma::mat T3 = calculate_T3(nodes);

    if (T3.is_empty()) {
        std::cerr << "Error in calculating T3. Cannot compute the inverse." << std::endl;
        return arma::mat();  // Return an empty matrix to indicate an error
    }

    // Calculate the inverse of T3
    arma::mat inverse_T3 = arma::inv(T3);

    return inverse_T3;
}
///////////////////////////////////////////////////////////////////////////////////////
arma::sp_mat Utils::spmat_submat(const arma::sp_mat& spmatrix, const std::vector<int>& row_indices, const std::vector<int>& col_indices) {
    // Validate input dimensions
    if (row_indices.size() != col_indices.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    // Get the number of rows and columns
    const int num_rows = spmatrix.n_rows;
    const int num_cols = spmatrix.n_cols;

    // Create vectors to store CSR format data
    std::vector<double> csr_values;
    std::vector<int> csr_column_indices;
    std::vector<int> csr_row_pointers;

    // Initialize the first row pointer
    csr_row_pointers.push_back(0);

    // Loop through the rows
    for (int i = 0; i < row_indices.size(); ++i) {
        int row = row_indices[i];
        int col = col_indices[i];

        // Check if the row and column indices are within bounds
        if (row < 0 || row >= num_rows || col < 0 || col >= num_cols) {
            throw std::out_of_range("Row or column index out of bounds.");
        }

        // Find the corresponding element in the sparse matrix
        double value = spmatrix(row, col);

        // Skip zero values
        if (value != 0.0) {
            // Store the value and column index
            csr_values.push_back(value);
            csr_column_indices.push_back(col);
        }

        // Update the row pointer if we are moving to a new row
        if (i + 1 < row_indices.size() && row_indices[i + 1] != row) {
            csr_row_pointers.push_back(csr_values.size());
        }
    }

    // Create the CSR format sparse matrix
    arma::sp_mat csr_matrix(
        arma::conv_to<arma::uvec>::from(csr_row_pointers),
        arma::conv_to<arma::uvec>::from(csr_column_indices),
        arma::vec(csr_values),
        num_rows,
        num_cols
    );

    return csr_matrix;
}