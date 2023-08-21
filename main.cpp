#include "src/InputReader.h"
#include "src/Mesh.h"
#include "src/PostProcessing.h"
#include "src/BCInit.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <gmsh.h>
#include <armadillo>

int testArmadillo() {
  arma::mat A = {{1, 2}, {3, 4}};
  arma::mat B = {{5, 6}, {7, 8}};
  arma::mat result = A + B;

  // Assert the expected result
  arma::mat expected = {{6, 8}, {10, 12}};
  if (arma::approx_equal(result, expected, "absdiff", 1e-6)) {
    std::cout << "Armadillo integration test passed!" << std::endl;
    return 0;
  } else {
    std::cout << "Armadillo integration test failed!" << std::endl;
  }
}

int main(int argc, char* argv[]) {
  
    testArmadillo();

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

    Mesh mesh(reader);
    mesh.getHexahedralElements();

    BCInit BC(reader,mesh);
    BC.boundaryConditions();

    PostProcessing psp(reader,mesh);
    psp.WriteUnstructuredMeshToVTK();

    return 0;
}
