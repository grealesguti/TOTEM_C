#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include "InputReader.h"
#include "Mesh.h"

#include <iostream>
#include <vector>
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkHexahedron.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLUnstructuredGridWriter.h"



class PostProcessing {
public:
    PostProcessing(const InputReader& inputReader, const Mesh& meshinput);
    void getHexahedralElements();
    //void getHexahedral20Elements();
    void convertHexMshToVtk();

private:
    const InputReader& inputReader_;
    const Mesh& meshinput_;

    int numPointsInHex = 8; // Get the number of points in the hexahedron i from the MSH file
};

#endif // POSTPROCESSING_H
