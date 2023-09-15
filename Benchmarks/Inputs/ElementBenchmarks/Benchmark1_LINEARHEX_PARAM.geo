// Define parameters for element size along each dimension
lc_x = 1;
lc_y = 2;
lc_z = 1;

// Define points with parameterized coordinates
Point(1) = {0, 0, 0, lc_x};         // Point 1
Point(2) = {lc_x, 0, 0, lc_x};      // Point 2
Point(3) = {lc_x, lc_y, 0, lc_x};   // Point 3
Point(4) = {0, lc_y, 0, lc_x};      // Point 4
Point(5) = {0, 0, lc_z, lc_z};      // Point 5
Point(6) = {lc_x, 0, lc_z, lc_z};   // Point 6
Point(7) = {lc_x, lc_y, lc_z, lc_z};// Point 7
Point(8) = {0, lc_y, lc_z, lc_z};   // Point 8

// Define lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Define curve loops and plane surfaces
Curve Loop(1) = {5, -10, -1, 9};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Curve Loop(3) = {8, -9, -4, 12};
Plane Surface(3) = {3};

Curve Loop(4) = {1, 2, 3, 4};
Plane Surface(4) = {4};

Curve Loop(5) = {12, -7, -11, 3};
Plane Surface(5) = {5};

Curve Loop(6) = {6, -11, -2, 10};
Plane Surface(6) = {6};

// Define surface loop and volume
Surface Loop(1) = {2, 1, 6, 5, 3, 4};
Volume(1) = {1};

// Define transfinite meshing
Transfinite Volume {1};
Transfinite Curve {8, 6, 2, 4} = 3 Using Progression 1;
Transfinite Surface {2};
Transfinite Surface {6};
Transfinite Surface {4};
Transfinite Surface {3};
Transfinite Surface {5};
Transfinite Surface {1};

// Define physical entities
Physical Volume("Volume", 13) = {1};
Physical Surface("minY", 14) = {1};
Physical Surface("maxY", 15) = {5};
Physical Surface("maxZ", 16) = {2};
Physical Surface("minZ", 17) = {4};
Physical Surface("minX", 18) = {3};
Physical Surface("maxX", 19) = {6};

// Set:
// 2D algorithm-> MeshAdapt
// 3D algorithm-> Delaunay
// 2D recombination algorithm-> Simple
// Tick Recobine all triangular meshes
// Subdivision Algorithm-> None (gives issues with hexs!)
// Smoothing Steps-> 1
// Element Size Factor -> 1
// Min/Max -> 2/1e+22
// Export mesh to msh


Meshing parameters
Mesh.Algorithm = 3;         // 3D algorithm -> Delaunay
Mesh.RecombineAll = 1;      // Recombine all triangular meshes
Mesh.SmoothingSteps = 1;    // Smoothing Steps -> 1
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.ElementSizeFactor = 1; // Element Size Factor -> 1
Mesh.MinSize = 2;
Mesh.MaxSize = 1e22;

Transfinite Curve {9, 5, 10, 1, 12, 7, 11, 3} = 2 Using Progression 1;

// Export mesh to msh
Mesh 3;

// Run using: gmsh script.geo //+
Physical Volume("testmaterial", 20) = {1};
//+
Physical Volume("testmaterial", 20) += {1};
