/* --------------------------------------------------------------------------- *\
							Private functions
\* --------------------------------------------------------------------------- */

# include "KdTreelib.h"

// Normalizes the control variables in the FGM struct to a range of [0, 1].
// @param fgm: Pointer to the FGM struct containing the data.
// @param mins: Array containing the minimum value of each control variable.
// @param maxs: Array containing the maximum value of each control variable.
void normalizeFGM
(
    FGM *fgm, 
	double *mins, 
	double *maxs
)
{
    // Calculate the total number of points for dynamic dimensions
    int totalPoints = 1;
    for (int dim = 0; dim < fgm->Ncv; dim++) {
        totalPoints *= fgm->Ngrid[dim];
    }

    // Normalize each control variable individually
    for (int i = 0; i < totalPoints; i++) {
        for (int d = 0; d < fgm->Ncv; d++) {
            double val = fgm->data[i * fgm->Nvar + d];
            // Apply min-max normalization to scale the variable to [0, 1]
            fgm->data[i * fgm->Nvar + d] = (val - mins[d]) / (maxs[d] - mins[d]);
        }
    }
}

// Creates a new Node for the k-d tree.
// @param point: A pointer to the point in the FGM data array.
// @return: A pointer to the newly created Node.
Node* newNode
(
    double *point
)
{
    Node* temp = (Node *)malloc(sizeof(Node));
    temp->point = point;
    temp->left = temp->right = NULL;
    return temp;
}

// Recursively inserts a point into the k-d tree.
// @param root: The root of the k-d tree (or subtree).
// @param point: The point to be inserted.
// @param depth: The current depth in the tree.
// @param Ncv: The number of control variables (dimensions).
// @return: The new root of the k-d tree (or subtree).
Node* insertRec
(
    Node* root, 
	double *point, 
	int depth, 
	int Ncv
)
{
    if (root == NULL) return newNode(point);

    // Calculate current dimension based on depth
    int cd = depth % Ncv;

    // Recursive case: compare the point with the current node and decide to go left or right
    if (point[cd] < root->point[cd])
        root->left = insertRec(root->left, point, depth + 1, Ncv);
    else
        root->right = insertRec(root->right, point, depth + 1, Ncv);

    return root;
}

// Helper function to print a point
void printPoint
(
    double *point, 
	int Ncv
)
{
    printf("(");
    for (int i = 0; i < Ncv; i++) {
        printf("%lf", point[i]);
        if (i < Ncv - 1) printf(", ");
    }
    printf(")");
}

/* --------------------------------------------------------------------------- *\
							Public functions
\* --------------------------------------------------------------------------- */

Node* buildKDTree
(
    FGM *fgm
)
{
	printf("\n Building K-d Tree...\n\n");
	
    // Normalize the FGM data
    normalizeFGM(fgm, fgm->mins, fgm->maxs);

    // Build the KD Tree
    Node* root = NULL;
    // Calculate the total number of points for dynamic dimensions
    int totalPoints = 1;
    for (int dim = 0; dim < fgm->Ncv; dim++) {
        totalPoints *= fgm->Ngrid[dim];
    }

    // Loop through each point in the FGM data and insert it into the tree
    for (int i = 0; i < totalPoints; ++i) {
        root = insertRec(root, &fgm->data[i * fgm->Nvar], 0, fgm->Ncv);
    }

	printf("\nK-d Tree constructed from FGM data\n\n");

    return root;
}

void printKDTree
(
    Node* root, 
	int depth, 
	int Ncv
)
{
    if (root == NULL) return;

    // Print the depth and the point at the current node
    for (int i = 0; i < depth; i++) printf("   "); // Indentation for better visualization
    printf("Depth %d: ", depth);
    printPoint(root->point, Ncv);
    printf("\n");

    // Recursively print left and right subtrees
    printKDTree(root->left, depth + 1, Ncv);
    printKDTree(root->right, depth + 1, Ncv);
}

// *************************************************************************** //