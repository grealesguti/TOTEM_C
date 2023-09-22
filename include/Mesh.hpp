#ifndef MESH_H
#define MESH_H


#include "InputReader.hpp"
#include "gmsh.h" // Assuming you have a gmsh library for reading the mesh file
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include "Elements.hpp"

#include "utils/data.hpp"


class Mesh {
public:
    Mesh(const InputReader& inputReader);
    void getHexahedralElements();

    // Finalize Gmsh method
    void finalizeGmsh();
    void reinitializeGmsh();

    // Get Physical entities
    std::vector< long unsigned int> getNodesForPhysicalGroup(const std::string& desiredGroupName);
    std::vector<std::size_t> getElementsForPhysicalGroup(const std::string& desiredGroupName);

    // Setters
    void setFixedof(int dof); 
    // Function to populate freedofsidx based on freedofs_
    void SetFreedofsIdx() {
        freedofsidx_.clear();  // Clear the vector if it contains any previous data

        for (int i = 0; i < freedofs_.size(); ++i) {
            if (freedofs_[i] == 1) {
                freedofsidx_.push_back(i);
            }
        }
    }
    // Function to retrieve the freedofsidx_ vector
    const std::vector<int>& GetFreedofsIdx() const {
        return freedofsidx_;
    }
    std::vector<int> printElementTypesInPhysicalGroupByName(std::string desiredGroupName);
    void InitMeshEntityElements();


    // Getter methods for private variables
    int getNumElements() const {        return numelem;    }
    int getNumNodes() const {        return numnodes;    }
    int getNumAllNodes() const {        return numallnodes;    }
    const std::vector<std::size_t> getElementTags() const {        return elementTags;    }
    const std::vector<std::size_t> getNodeTags() const {        return nodeTags;    }
    const std::vector<double>& getCoordinates() const {        return coord;    }
    const double getCoordi(int i) const {        return coord[i];    }
    const int getelementNodeTagi(int i) const {        return elementNodeTags[i];    }
    const int getNodeTagi(int i) const {        return nodeTags[i];    }
    const int getElementMaterial(int i) const {        return element_materials[i];    }
    std::vector<int> getAllElementMaterials() const {        return element_materials;    }

    arma::mat getCoordinates(const std::vector<int>& nodeTags);
    int getElementInfo(int elementTag, std::vector<int> & nodeTags_el);
    void applyElementsMaterials();
    // Declaration for getMaterialPropertyForElement function
    double getMaterialPropertyForElement(std::size_t elementIndex, const std::string& propertyName) const;
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    getElementsAndNodeTagsForPhysicalGroup(const std::string& desiredGroupName);
    // Function to return coordinates
    const std::vector<double>& getAllCoordinates() const {
        return coord;
    }
    std::pair<int, int> getNumNodesForElement(int elementTag);

    /**
     * @brief Select shape functions and derivatives based on element type.
     *
     * @param[in] etype  Element type.
     * @param[in] xi     .
     * @param[in] eta     .
     * @param[in] zeta     .
     * @param[out] shapeFunctions     .
     * @param[out] shapeFunctionDerivatives     .
     */
    void selectShapeFunctionsAndDerivatives(const int etype, const double xi, const double eta, const double zeta,
                                            Armadillo<arma::vec>& shapeFunctions, Armadillo<arma::mat>& shapeFunctionDerivatives);

private:
    Elements elements_;
    const InputReader& inputReader_;
    std::vector<std::size_t> elementTags,nodeTags,elementNodeTags,nodeTagselem;
    int numPointsInHex = 8; // Get the number of points in the hexahedron i from the MSH file
    int numelem, numnodes, numallnodes;
    std::vector<double> coord;
    std::vector<int> freedofs_,freedofsidx_, element_materials;
    std::map<std::string, std::size_t> materialIndices;
};

#endif // MESH_H
