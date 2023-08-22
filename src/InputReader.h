#ifndef INPUTREADER_H
#define INPUTREADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include "gmsh.h"

class InputReader {
public:
    InputReader(const std::string& filename);
    bool readFile();

    // Printing to screen values found for Material and Boundary conditions
    void printMaterialProperties() const;
    void printBoundaryConditions() const;

    // Setter methods //
    void setMeshFileName(const std::string& meshFileName);
    void setMeshEntityName(const std::string& MeshEntityName);

    // Getter methods //
    const std::string& getMeshFileName() const {return meshFileName_;}
    const std::string& getMeshEntityName() const {return MeshEntityName_;}

    const std::map<int, std::map<std::string, double>>& getAllMaterialProperties() const;
    const std::map<std::pair<std::string, std::string>, double>& getBoundaryConditions() const {
        return boundaryConditions_;
    }

    // OTHER METHODS //
    // Method to read Gmsh mesh file
    bool readGmshMesh();
private:
    std::string filename_;
    std::string meshFileName_,MeshEntityName_; // Variable to store the mesh file name
    std::map<std::string, std::map<std::string, double>> MaterialProperties_;
    std::map<std::pair<std::string, std::string>, double> boundaryConditions_;

};

#endif // INPUTREADER_H
