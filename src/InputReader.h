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

    const std::map<std::string, std::map<std::string, double>> getAllMaterialProperties() const {
        return MaterialProperties_;
    }
    const std::map<std::pair<std::string, std::string>, double>& getBoundaryConditions() const {
        return boundaryConditions_;
    }
    // Setter method for DesiredOutput_
    void setDesiredOutput(const std::string& desiredOutput) {
        DesiredOutput_ = desiredOutput;
    }

    // Getter method for DesiredOutput_
    const std::string& getDesiredOutput() const {
        return DesiredOutput_;
    }
    // Getter method for DesiredOutput_
    const std::string& getPhysics() const {
        return physics_;
    }
    double getPropertyValue(int index, const std::string& key, bool& success) const;
    std::string getMaterial_index(int index) const;
    double getMaterialPropertyValue(int materialIndex, const std::string& propertyName) const;

    // OTHER METHODS //
    // Method to read Gmsh mesh file
    bool readGmshMesh();
private:
    std::string filename_;
    std::string meshFileName_,MeshEntityName_,DesiredOutput_,physics_; // Variable to store the mesh file name
    std::map<std::string, std::map<std::string, double>> MaterialProperties_;
    std::map<std::pair<std::string, std::string>, double> boundaryConditions_;
    // Initialize a vector to store material indices for elements
    std::vector<std::size_t> elementMaterialIndices;
    std::map<int, std::string> Material_index;

};

#endif // INPUTREADER_H
