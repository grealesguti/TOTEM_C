#include "PostProcessing.h"
/*
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
*/

PostProcessing::PostProcessing(const InputReader& inputReader, const Mesh& meshinput)
    : inputReader_(inputReader), meshinput_(meshinput) {
    // You may want to initialize other member variables here
}

void PostProcessing::convertHexMshToVtk() {
    // Initialize Gmsh
    gmsh::initialize();
    gmsh::open(inputReader_.getMeshFileName());

    // Create VTK objects
         /*
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Set up point coordinates
    for (int i = 0; i < meshinput_.getNumAllNodes()*3; i++) {
        std::cout << "Setting vtk nodes... " << std::endl;
        double px = meshinput_.getCoordi(i*3+0); // Get x-coordinate of point i from MSH file
        double py = meshinput_.getCoordi(i*3+2);// Get y-coordinate of point i from MSH file
        double pz = meshinput_.getCoordi(i*3+3);// Get z-coordinate of point i from MSH file
        points->InsertNextPoint(px, py, pz);
    }

    // Set up hexahedral elements
    for (int i = 0; i < meshinput_.getNumElements(); i++) {
        std::cout << "Setting vtk elements... " << std::endl;
        hexahedron->GetPointIds()->SetNumberOfIds(numPointsInHex);
        std::vector<size_t> nodeTags_i;
        for (int j = 0; j < 8; j++) {
            int pointIndex = meshinput_.getelementNodeTagi(j*numPointsInHex+j);// Get the point index for the j-th point of hexahedron i from the MSH file
            hexahedron->GetPointIds()->SetId(j, pointIndex);
        }

        cells->InsertNextCell(hexahedron);
    }

    // Set up unstructured grid
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_HEXAHEDRON, cells);

    // Write VTK file
    std::cout << "Writing VTK... " << std::endl;
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("hexmesh.vtk");
    writer->SetInputData(unstructuredGrid);
    writer->Write();*/

    gmsh::finalize();
}

void PostProcessing::WriteUnstructuredMeshToVTK() {
    std::ofstream vtkFile("test.vtk");

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Unstructured Mesh" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write node coordinates
    vtkFile << "POINTS " <<  meshinput_.getNumAllNodes() << " double" << std::endl;
    for (int i = 0; i <  meshinput_.getNumAllNodes(); ++i) {
        vtkFile <<  meshinput_.getCoordi(i * 3) << " " <<  meshinput_.getCoordi(i * 3 + 1) << " " <<  meshinput_.getCoordi(i * 3 + 2) << std::endl;
    }


    vtkFile << "CELLS " << meshinput_.getNumElements() << " " << meshinput_.getNumElements()*(8+1) << std::endl;
    for (int i = 0; i < meshinput_.getNumElements(); ++i) {
        vtkFile << 8;
        for (int j = 0; j < 8; ++j) {
            vtkFile << " " << (meshinput_.getelementNodeTagi(i*8+j)-1); //Starts from index 0 in VTK and 1 in MSH
        }
        vtkFile << std::endl;
    }

    // Write cell types (assuming all elements are hexahedra)
    vtkFile << "CELL_TYPES " << meshinput_.getNumElements() << std::endl;
    for (int i = 0; i < meshinput_.getNumElements(); ++i) {
        vtkFile << "12 "; // VTK_HEXAHEDRON code
    }
    vtkFile << std::endl;

    vtkFile.close();
}