#ifndef BCINIT_H
#define BCINIT_H

#include "InputReader.h"
#include "Mesh.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class BCInit {
public:
    BCInit(const InputReader& inputReader, Mesh& mesh);
    void boundaryConditions();
    void createLoadVector();

private:
    const InputReader& inputReader_;
    Mesh& mesh_;

    std::string meshFileName_;
    std::vector<double> loadVector_,initialdofs_;
};

#endif // BCINIT_H
