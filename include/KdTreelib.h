/* --------------------------------------------------------------------------- *\

                      Create the K-D Tree from the FGM data

\* --------------------------------------------------------------------------- */

#ifndef KDTREELIB_H
#define KDTREELIB_H

# include "FGMlib.h"

// Structure representing a node in the K-D Tree
typedef struct Node {
    float *point; // Pointer to the point in FGM data array
    struct Node *left, *right;
} Node;

// Builds the k-d tree from the FGM struct's data with normalization.
// @param fgm: Pointer to the FGM struct containing the data.
// @return: The root of the built k-d tree.
Node* buildKDTree
(
    FGM *fgm
);

// Prints the k-d tree in a pre-order traversal manner.
// @param root: The root of the tree (or subtree).
// @param depth: The current depth in the tree.
// @param Ncv: The number of control variables (dimensions).
void printKDTree
(
    Node* root, 
	int depth, 
	int Ncv
);

#endif // KDTREELIB_H

// *************************************************************************** //