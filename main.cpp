#include "src/InputReader.h"
#include "src/Mesh.h"
#include "src/PostProcessing.h"
#include "src/BCInit.h"
#include "src/Testing.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <gmsh.h>
#include <armadillo>


int main(int argc, char* argv[]) {
  
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

    BCInit BC(reader,mesh);

    PostProcessing psp(reader,mesh);
    psp.WriteUnstructuredMeshToVTK();

    return 0;
}
