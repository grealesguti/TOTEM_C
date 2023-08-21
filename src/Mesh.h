#ifndef MESH_H
#define MESH_H


#include "InputReader.h"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>

class Mesh {
public:
    Mesh(const InputReader& inputReader);
    void getHexahedralElements();

    // Get Physical entities
    std::vector<long long unsigned int> getNodesForPhysicalGroup(const std::string& desiredGroupName);
    std::pair<size_t, std::vector<int>> getElementsForPhysicalGroup(const std::string& desiredGroupName);

    // Setters
    void setFixedof(int dof); 

    // Getter methods for private variables
    int getNumElements() const {        return numelem;    }
    int getNumNodes() const {        return numnodes;    }
    int getNumAllNodes() const {        return numallnodes;    }
    const std::vector<std::size_t>& getElementTags() const {        return elementTags;    }
    const std::vector<std::size_t>& getNodeTags() const {        return nodeTags;    }
    const std::vector<double>& getCoordinates() const {        return coord;    }
    const double getCoordi(int i) const {        return coord[i];    }
    const int getelementNodeTagi(int i) const {        return elementNodeTags[i];    }
    const int getNodeTagi(int i) const {        return nodeTags[i];    }

private:
    const InputReader& inputReader_;
    std::vector<std::size_t> elementTags,nodeTags,elementNodeTags,nodeTagselem;
    int numPointsInHex = 8; // Get the number of points in the hexahedron i from the MSH file
    int numelem, numnodes, numallnodes;
    std::vector<double> coord;
    std::vector<int> freedofs_;

};

#endif // MESH_H
