/* --------------------------------------------------------------------------- *\

            Perform N nearest neighbor lookup including interpolation

\* --------------------------------------------------------------------------- */

#ifndef KNNLOOKUPLIB_H
#define KNNLOOKUPLIB_H

# include "FGMlib.h"
# include "KdTreelib.h"

// Structure representing a Neighbor
typedef struct {
    Node* node;     // Pointer to a node in the k-d tree
    double distance; // Distance from the query point to this node
} Neighbor;

// Structure representing a NeighborHeap
typedef struct {
    Neighbor* elements; // Array of Neighbor elements
    int size;           // Current number of elements in the heap
    int capacity;       // Maximum capacity of the heap
} NeighborHeap;

// Normalizes a query point using the same min-max values used for the original FGM data.
// @param queryPoint: The control variables of the query point to be normalized.
// @param mins: Array containing the minimum value of each control variable.
// @param maxs: Array containing the maximum value of each control variable.
// @param Ncv: The number of control variables (dimensions).
void normalizeQueryPoint
(
    double *queryPoint, 
	double *mins, 
	double *maxs, 
	int Ncv
);

// Denormalizes a query point back to its original scale using the min-max values.
// @param queryPoint: The normalized control variables of the query point.
// @param mins: Array containing the minimum value of each control variable.
// @param maxs: Array containing the maximum value of each control variable.
// @param Ncv: The number of control variables (dimensions).
void denormalizeQueryPoint
(
    double *queryPoint, 
	double *mins, 
	double *maxs, 
	int Ncv
);

// Performs a nearest neighbor search in the k-d tree.
// @param root: The current node in the k-d tree.
// @param queryPoint: The control variables for which the nearest neighbor is searched.
// @param depth: The current depth in the k-d tree.
// @param Ncv: The number of control variables (dimensions).
// @return: The nearest node to the query point in the k-d tree.
void nearestNeighborSearch
(
    Node* root, 
	Neighbor *nearest, 
	double queryPoint[], 
	int depth, 
	int Ncv
);

// Lookup the nearest neighbor in a KD Tree and return the result in 'f'.
// @param fgm: Pointer to the FGM struct containing data and normalization parameters.
// @param KDTree: Pointer to the root node of the KD Tree.
// @param x: The query point to find the nearest neighbor.
// @param f: Output array to store the nearest neighbor values.
// @return: EXIT_SUCCESS if the lookup is successful, EXIT_FAILURE otherwise.
int NNlookupFGM
(
    FGM *fgm, 
	Node *KDTree, 
	double *x, 
	double *f
);

// Lookup the nearest neighbor in a KD Tree with interpolation and return the result in 'f'.
// @param fgm: Pointer to the FGM struct containing data and normalization parameters.
// @param KDTree: Pointer to the root node of the KD Tree.
// @param x: The query point to find the nearest neighbor.
// @param f: Output array to store the interpolated values.
// @param K: Number of nearest neighbors used for interpolation
// @return: EXIT_SUCCESS if the lookup is successful, EXIT_FAILURE otherwise.
int KNNlookupFGM_Interp
(
    FGM *fgm, 
	Node *KDTree, 
	double *x, 
	double *f, 
	int K
);

#endif // KNNLOOKUPLIB_H

// *************************************************************************** //