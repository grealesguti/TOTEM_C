#ifndef BCINIT_H
#define BCINIT_H

#include "InputReader.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file

class BCInit {
public:
    BCInit(const InputReader& inputReader);
    void boundaryConditions() const;
    void createLoadVector();

private:
    const InputReader& inputReader_;
    std::string meshFileName_;
    std::vector<double> loadVector_;
};

#endif // BCINIT_H
