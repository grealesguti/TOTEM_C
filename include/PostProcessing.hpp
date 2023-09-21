#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include "InputReader.hpp"
#include "Mesh.hpp"

#include <iostream>
#include <vector>

class PostProcessing {
public:
    PostProcessing(const InputReader& inputReader, const Mesh& meshinput);
    void getHexahedralElements();
    //void getHexahedral20Elements();
    void convertHexMshToVtk();
    void WriteUnstructuredMeshToVTK();

private:
    const InputReader& inputReader_;
    const Mesh& meshinput_;

    int numPointsInHex = 8; // Get the number of points in the hexahedron i from the MSH file
};

#endif // POSTPROCESSING_H
