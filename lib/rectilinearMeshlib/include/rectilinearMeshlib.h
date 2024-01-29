/* --------------------------------------------------------------------------- *\

					Create a rectilinear mesh for the FGM data

\* --------------------------------------------------------------------------- */

#ifndef RECTILINEARMESHLIB_H
#define RECTILINEARMESHLIB_H

# include "FGMlib.h"
# include "KNNlookuplib.h"
# include <math.h>
# include <omp.h>

// Structure representing the rectinlinear mesh
typedef struct {
    float *reducedMesh;  // Pointer to the reduced mesh
    int *mapMesh;         // Pointer to the mapped mesh
    int validPointCount;  // Number of valid points in the reduced mesh
	float *meshMin;
	float *meshMax;
} Mesh;

// Function to create a rectilinear mesh.
// @param fgm: Pointer to the FGM struct containing point cloud data.
// @param N: Number of grid points in each dimension.
// @return: Pointer to the rectilinear mesh.
Mesh createRectilinearMesh
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K
);

// Create a standard rectilinear mesh from point cloud data.
// @param fgm: Pointer to the FGM struct containing point cloud data.
// @param KDTree: Pointer to the root of the KD Tree constructed from FGM data.
// @param Ngrid: Array specifying the number of grid points in each dimension.
// @param K: The number of nearest neighbors to consider in K-D Tree.
// @return: A Mesh struct representing the standard rectilinear mesh.
Mesh createRectilinearMesh
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K
);

// Create a parallelized standard rectilinear mesh from point cloud data.
// @param fgm: Pointer to the FGM struct containing point cloud data.
// @param KDTree: Pointer to the root of the KD Tree constructed from FGM data.
// @param Ngrid: Array specifying the number of grid points in each dimension.
// @param K: The number of nearest neighbors to consider in K-D Tree.
// @return: A Mesh struct representing the parallelized standard rectilinear mesh.
Mesh createRectilinearMesh_parallel
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K
);

// Create a reduced rectilinear mesh from point cloud data, excluding points outside the cloud.
// @param fgm: Pointer to the FGM struct containing point cloud data.
// @param KDTree: Pointer to the root of the KD Tree constructed from FGM data.
// @param Ngrid: Array specifying the number of grid points in each dimension.
// @param K: The number of nearest neighbors to consider in K-D Tree.
// @param threshold: The maximum distance threshold to include a point in the mesh.
// @return: A Mesh struct representing the reduced rectilinear mesh.
Mesh createRectilinearMesh_reduced
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K, 
	float threshold
);

// Create a parallelized reduced rectilinear mesh from point cloud data, excluding points outside the cloud.
// @param fgm: Pointer to the FGM struct containing point cloud data.
// @param KDTree: Pointer to the root of the KD Tree constructed from FGM data.
// @param Ngrid: Array specifying the number of grid points in each dimension.
// @param K: The number of nearest neighbors to consider in K-D Tree.
// @param threshold: The maximum distance threshold to include a point in the mesh.
// @return: A Mesh struct representing the parallelized reduced rectilinear mesh.
Mesh createRectilinearMesh_reduced_parallel
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K, 
	float threshold
);

#endif // RECTILINEARMESHLIB_H
