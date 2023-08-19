#include "src/InputReader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <gmsh.h>

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
    if (reader.readFile()) {
        reader.printMaterialProperties();
        return 0;
    } else {
        return 1;
    }
}
