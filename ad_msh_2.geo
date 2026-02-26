// --- Geometry Parameters ---
L = 1.0; 
H = 1.0; 

// --- Points ---
Point(1) = {-L, 0, 0};
Point(2) = {0,  0, 0};   // Origin
Point(3) = {L,  0, 0};
Point(4) = {L,  H, 0};
Point(5) = {0,  H, 0};   // Top Partition
Point(6) = {-L, H, 0};
Point(7) = {-0.5, 0, 0}; // Refinement Point

// --- Lines ---
Line(1) = {1, 7};   // Part of Inflow
Line(2) = {7, 2};   // Part of Inflow
Line(3) = {2, 3};   // Outflow
Line(4) = {3, 4};   // Right Wall
Line(5) = {4, 5};   // Top Wall Right
Line(6) = {5, 6};   // Top Wall Left
Line(7) = {6, 1};   // Left Wall
Line(8) = {2, 5};   // INTERNAL PARTITION (No Physical Tag)

// --- Surfaces ---
Curve Loop(1) = {1, 2, 8, 6, 7}; 
Plane Surface(1) = {1};
Curve Loop(2) = {3, 4, 5, -8}; 
Plane Surface(2) = {2};

Point{7} In Surface{1};

// --- Physical Groups (Keys) ---
// We group Lines 1 & 2 together as they both sit on the Inflow boundary
Physical Curve("Inflow") = {1, 2};
Physical Curve("Outflow") = {3};
Physical Curve("Wall") = {4, 5, 6, 7};

// You must also define the Surface as a Physical Group, 
// otherwise some solvers will ignore the elements inside.
Physical Surface("Domain") = {1, 2};

// --- Mesh Fields (The "COMSOL" Bleed) ---
Field[1] = Distance;
Field[1].PointsList = {2, 7};

Field[2] = MathEval;
Field[2].F = "0.008 + 0.12 * F1"; 

Background Field = 2;

// --- Mesh Aesthetics ---
Mesh.Algorithm = 6; 
Mesh.Smoothing = 10;
Mesh.MeshSizeFactor = 0.7;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
