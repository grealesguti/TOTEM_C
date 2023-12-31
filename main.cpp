#include "src/InputReader.h"
#include "src/Mesh.h"
#include "src/PostProcessing.h"
#include "src/BCInit.h"
#include "src/Utils.h"
#include "src/Solver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <gmsh.h>
#include <armadillo>

#include <cstdlib>
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

int main(int argc, char* argv[]) {
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | -_CRTDBG_LEAK_CHECK_DF);
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
    utils.writeDataToFile(mesh.getAllElementMaterials(),"Outputs/mesh_ElementMaterialIndex.txt",false);
    utils.writeDataToFile(mesh.GetFreedofsIdx(),"Outputs/mesh_GetFreedofsIdx.txt",false);

    BCInit BC(reader,mesh);
    //std::cout<<"Desired Output:"<<std::endl;
    //std::cout<<reader.getDesiredOutput()<<std::endl;
    utils.writeDataToFile(BC.getAllInitialDof(),"Outputs/Out_DegreesOfFreedom.txt",false);
    std::vector<int> freedofidxs_ = mesh.GetFreedofsIdx();
    
    // Initialization of Solver variablesç
    int numFreeDofs = freedofidxs_.size(); // Assuming freedofidx_ is a private member variable
    Eigen::SparseMatrix<double> reducedK(numFreeDofs, numFreeDofs);
    std::vector<double> reducedR(numFreeDofs, 0.0); // Initialize with size and value 0.0
    Eigen::VectorXd solution;    solution.resize(numFreeDofs);    solution.setZero();    // Assembly

    // Initialize solver and assembly
    Solver solver(reader, mesh, BC);
    solver.Assembly(reducedK,reducedR);
    std::cout<<"Assembly finished:"<<std::endl;

    //solver.solveSparseSystem(system);
    //utils.writeDataToFile(solution,"Outputs/main_solution.txt",false);
    //std::cout<<"System solved"<<std::endl;
    //Solver::SparseSystem system2 = solver.Assembly();
    //std::cout<<"Assembly finished:"<<std::endl;
    //Eigen::VectorXd solution2 = solver.solveSparseSystem(system);
    //utils.writeDataToFile(solution,"Outputs/main_solution.txt",false);
    //std::cout<<"System solved"<<std::endl;
    //solver.runNewtonRaphson();
    //utils.writeDataToFile(solver.getAllSolDofs(),"Outputs/main_NRsolution.txt",false);

    if (reader.getDesiredOutput()=="all"){
        std::cout << "OUTPUT:: ALL" << std::endl;
        utils.writeDataToFile(mesh.getAllCoordinates(),"Outputs/Out_CoordVector.txt",false);
        //utils.writeDataToFile(BC.getloadVector(),"Outputs/Out_LoadVector.txt",false);
        PostProcessing psp(reader,mesh);
        psp.WriteUnstructuredMeshToVTK();
    }else if (reader.getDesiredOutput()=="mesh"){
        PostProcessing psp(reader,mesh);
        psp.WriteUnstructuredMeshToVTK();
    }

    _CrtDumpMemoryLeaks();

    return 0;
}
