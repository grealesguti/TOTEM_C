#include "Solver.h"
using namespace arma;

Solver::Solver(const InputReader& inputReader, Mesh& mesh, BCInit& bcinit)
    : inputReader_(inputReader), mesh_(mesh), elements_(), bcinit_(bcinit) {
    // Get the number of nodes from the mesh
    int numNodes = mesh_.getNumAllNodes();

    // Initialize the loadVector_ member with a size double the number of nodes
    loadVector_.resize(2 * numNodes);
    initialdofs_.resize(2 * numNodes);

    // You may want to initialize other member variables here
}