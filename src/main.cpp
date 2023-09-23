#include "InputReader.hpp"
#include "Mesh.hpp"
#include "PostProcessing.hpp"
#include "BCInit.hpp"
#include "Utils.hpp"
#include "Solver.hpp"
#include "utils/data.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <gmsh.h>
#include <armadillo>



int main(int argc, char* argv[]) {
    //#pragma omp parallel
    //{
    //    printf("Hello World... from thread = %d\n", omp_get_thread_num());
    //}

      // Set the number of threads programmatically
   // omp_set_num_threads(4); // Set to use 4 threads

    Utils utils;
    utils.deleteFilesInFolder("Outputs");

    std::string inputFileName = "input.txt"; // Default input file name

    // Parse command-line arguments
    if (argc >= 3 && std::string(argv[1]) == "-i") {  // -i <Input file name.txt>
        inputFileName = argv[2];
    } else { // we must at least provide an input file name for the code to run
        std::cout << "Usage: " << argv[0] << " -i <input_file.txt>" << std::endl;
        return 1;
    }

    // Read the input file 
    InputReader reader(inputFileName); 
    std::string filename = reader.getMeshFileName();
    // Remove the last 4 characters using string slicing
    filename = filename.substr(0, filename.length() - 4);    

    Mesh mesh(reader);
    mesh.getAllElementMaterials().writeDataToFile("Outputs/mesh_ElementMaterialIndex.txt",false);
    mesh.GetFreedofsIdx().writeDataToFile("Outputs/mesh_GetFreedofsIdx.txt",false);

    BCInit BC(reader,mesh);
    //std::cout<<"Desired Output:"<<std::endl;
    //std::cout<<reader.getDesiredOutput()<<std::endl;
    BC.getAllInitialDof().writeDataToFile("Outputs/Out_DegreesOfFreedom.txt",false);
    std::vector<int> freedofidxs_ = mesh.GetFreedofsIdx();

    // Initialization of Solver variablesç
    int numFreeDofs = freedofidxs_.size(); // Assuming freedofidx_ is a private member variable
    Eigen::SparseMatrix<double> reducedK(numFreeDofs, numFreeDofs);
    Vector<double> reducedR({numFreeDofs, 0.0}); // Initialize with size and value 0.0
    EigenVectorXd solution;    solution.resize(numFreeDofs);    solution.setZero();    // Assembly
    //std::shared_ptr<Eigen::SparseMatrix<double>> reducedK = std::make_shared<Eigen::SparseMatrix<double>>();
    //std::shared_ptr<std::vector<double>> reducedR = std::make_shared<std::vector<double>>();

    // Initialize solver and assembly
    Solver solver(reader, mesh, BC);
    solver.Assembly(reducedK,reducedR);
    std::cout<<"Assembly finished:"<<std::endl;

    //solver.solveSparseSystem(system);
    //solution.writeDataToFile("Outputs/main_solution.txt",false);
    //std::cout<<"System solved"<<std::endl;
    //Solver::SparseSystem system2 = solver.Assembly();
    //std::cout<<"Assembly finished:"<<std::endl;
    //EigenVectorXd solution2 = solver.solveSparseSystem(system);
    //solution2.writeDataToFile("Outputs/main_solution.txt",false);
    //std::cout<<"System solved"<<std::endl;
    //solver.runNewtonRaphson();
    //solver.getAllSolDofs().writeDataToFile("Outputs/main_NRsolution.txt",false);

    if (reader.getDesiredOutput()=="all"){
        std::cout << "OUTPUT:: ALL" << std::endl;
        mesh.getAllCoordinates().writeDataToFile("Outputs/Out_CoordVector.txt",false);
        //BC.getloadVector().writeDataToFile("Outputs/Out_LoadVector.txt",false);
        PostProcessing psp(reader,mesh);
        psp.WriteUnstructuredMeshToVTK();
    }else if (reader.getDesiredOutput()=="mesh"){
        PostProcessing psp(reader,mesh);
        psp.WriteUnstructuredMeshToVTK();
    }

    mesh.finalizeGmsh();
    return 0;
}
