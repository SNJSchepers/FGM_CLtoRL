/* --------------------------------------------------------------------------- *\
							Private functions
\* --------------------------------------------------------------------------- */

# include "KNNlookuplib.h"

// Calculates the squared Euclidean distance between two points.
// @param point1: The first point (array of floats).
// @param point2: The second point (array of floats).
// @param Ncv: The number of control variables (dimensions of the points).
// @return: The squared Euclidean distance between point1 and point2.
float distanceSquared
(
    float point1[], 
	float point2[], 
	int Ncv
)
{
    float dist = 0.0;
    // Loop over each dimension and add the squared difference
    // between the coordinates to the total distance.
    for (int i = 0; i < Ncv; i++) {
        float diff = point1[i] - point2[i];
        dist += diff * diff;
    }
    return dist;  // Return the total squared distance
}

// Create a new NeighborHeap with the specified capacity.
// @param capacity: The maximum capacity of the heap.
// @return: A pointer to the newly created NeighborHeap, or NULL on failure.
NeighborHeap* createNeighborHeap
(
    int capacity
)
{
    NeighborHeap* heap = (NeighborHeap*)malloc(sizeof(NeighborHeap));
    if (!heap) return NULL;

    heap->elements = (Neighbor*)malloc(capacity * sizeof(Neighbor));
    if (!heap->elements) {
        free(heap);
        return NULL;
    }

    heap->size = 0;
    heap->capacity = capacity;
    return heap;
}

// Heapify upwards to maintain the max-heap property after insertion.
// @param heap: Pointer to the NeighborHeap being modified.
// @param index: The index at which to start the heapify operation.
void heapifyUp
(
    NeighborHeap* heap, 
	int index
)
{
    while (index != 0 && heap->elements[index].distance > heap->elements[(index - 1) / 2].distance) {
        // Swap with parent
        Neighbor temp = heap->elements[index];
        heap->elements[index] = heap->elements[(index - 1) / 2];
        heap->elements[(index - 1) / 2] = temp;

        index = (index - 1) / 2;
    }
}

// Heapify downwards to maintain the max-heap property after insertion or replacement.
// @param heap: Pointer to the NeighborHeap being modified.
// @param index: The index at which to start the heapify operation.
void heapifyDown
(
    NeighborHeap* heap, 
	int index
)
{
    int largest = index;
    int left = 2 * index + 1; 
    int right = 2 * index + 2;

    if (left < heap->size && heap->elements[left].distance > heap->elements[largest].distance) {
        largest = left;
    }

    if (right < heap->size && heap->elements[right].distance > heap->elements[largest].distance) {
        largest = right;
    }

    if (largest != index) {
        // Swap with the larger child
        Neighbor temp = heap->elements[index];
        heap->elements[index] = heap->elements[largest];
        heap->elements[largest] = temp;

        // Recursively heapify the affected subtree
        heapifyDown(heap, largest);
    }
}

// Insert a Neighbor into the NeighborHeap.
// @param heap: Pointer to the NeighborHeap where the Neighbor is inserted.
// @param neighbor: The Neighbor to be inserted.
// @return: true if the insertion was successful, false otherwise.
bool insertNeighborHeap
(
    NeighborHeap* heap, 
	Neighbor neighbor
)
{
    if (heap->size == heap->capacity) {
        // Heap is full, replace the top element if the new neighbor is closer (smaller distance)
        if (neighbor.distance < heap->elements[0].distance) {
            heap->elements[0] = neighbor;
            heapifyDown(heap, 0); // Restore heap property
        }
    } else {
        // If the heap is not full, insert the new neighbor and heapify up
        heap->elements[heap->size] = neighbor;
        heapifyUp(heap, heap->size);
        heap->size++;
    }
    return true;
}

// Perform a N nearest neighbor search in the KD Tree.
// @param root: Pointer to the root node of the KD Tree.
// @param queryPoint: The query point for which neighbors are searched.
// @param depth: Current depth in the KD Tree (usually start with 0).
// @param Ncv: The number of control variables (dimensions).
// @param heap: Pointer to the NeighborHeap where results are stored.
void KNNSearch
(
    Node* root, 
	float queryPoint[], 
	int depth, 
	int Ncv, 
	NeighborHeap* heap
)
{
    if (root == NULL) return;

    int cd = depth % Ncv;  // Calculate current dimension based on depth

    // Determine which branch (left or right) is closer to the query point
    bool goesLeft = queryPoint[cd] < root->point[cd];
    Node* closerBranch = goesLeft ? root->left : root->right;
    Node* fartherBranch = goesLeft ? root->right : root->left;

    // First, visit the closer branch
    KNNSearch(closerBranch, queryPoint, depth + 1, Ncv, heap);
 
    // Calculate distance from the query point to the current node
    float dist = distanceSquared(queryPoint, root->point, Ncv);
    
    // Create a Neighbor struct for the current node
    Neighbor neighbor = { root, dist };

    // Insert or replace in the heap
    insertNeighborHeap(heap, neighbor);

    // Check if we need to visit the farther branch
    float radiusSquared = (queryPoint[cd] - root->point[cd]) * (queryPoint[cd] - root->point[cd]);
    if (heap->size < heap->capacity || radiusSquared < heap->elements[0].distance) {
        // Visit the farther branch if there's potential for closer neighbors
        KNNSearch(fartherBranch, queryPoint, depth + 1, Ncv, heap);
    }
	
	//printf("%f\t%f\n",root->point[0],root->point[1]);
}

/* --------------------------------------------------------------------------- *\
							Public functions
\* --------------------------------------------------------------------------- */

void normalizeQueryPoint
(
    float *queryPoint, 
	FGM *fgm
)
{
    for (int d = 0; d < fgm->Ncv; d++) {
        // Apply min-max normalization to each dimension of the query point
        queryPoint[d] = (queryPoint[d] - fgm->mins[d]) / (fgm->maxs[d] - fgm->mins[d]);
	}
}

void denormalizeQueryPoint
(
    float *queryPoint,
	FGM *fgm
)
{
    for (int d = 0; d < fgm->Ncv; d++) {
        // Apply inverse min-max normalization to each dimension of the query point
        queryPoint[d] = queryPoint[d] * (fgm->maxs[d] - fgm->mins[d]) + fgm->mins[d];
    }
}

void nearestNeighborSearch
(
    Node* root, 
	Neighbor *nearest, 
	float queryPoint[], 
	int depth, 
	int Ncv
)
{
	if (root == NULL) {
		nearest->node = NULL;
		nearest->distance = INFINITY;
		return;
	}

    int cd = depth % Ncv;  // Calculate current dimension based on depth

    Node* nextBranch = (queryPoint[cd] < root->point[cd]) ? root->left : root->right;
    Node* otherBranch = (queryPoint[cd] < root->point[cd]) ? root->right : root->left;

	// Recursively search the next branch
	Neighbor tempNearest;
	nearestNeighborSearch(nextBranch, &tempNearest, queryPoint, depth + 1, Ncv);

	// Calculate the squared distance from the query point to the current node
	float dRoot = distanceSquared(queryPoint, root->point, Ncv);
	float dBest = tempNearest.node ? tempNearest.distance : INFINITY;
	
	// Check if the current node is closer than what we have found so far
	if (dRoot < dBest) {
		nearest->node = root;
		nearest->distance = dRoot;
	} else {
		nearest->node = tempNearest.node;
		nearest->distance = dBest;
	}
		
	// Check if we need to search the other branch
	float radiusSquared = queryPoint[cd] - root->point[cd];
	radiusSquared *= radiusSquared;

	if (radiusSquared < nearest->distance) {
		Neighbor otherNearest;
		nearestNeighborSearch(otherBranch, &otherNearest, queryPoint, depth + 1, Ncv);

		// Update best node if a closer node is found in the other branch
		if (otherNearest.distance < nearest->distance) {
			nearest->node = otherNearest.node;
			nearest->distance = otherNearest.distance;
		}
	}
}

int NNlookupFGM
(
    FGM *fgm, 
	Node *KDTree, 
	float *x, 
	float *f
)
{
    // Normalize the query point using the same parameters used for the FGM data
    normalizeQueryPoint(x, fgm);

	// Create instance of neighbor
    Neighbor* nearest;
	
    // Perform the nearest neighbor search in the k-d tree
    nearestNeighborSearch(KDTree, nearest, x, 0, fgm->Ncv);
    
    if (nearest == NULL) {
        // If no nearest neighbor is found, return failure
        return EXIT_FAILURE;
    }

    // Copy the data from the nearest point to the output array 'f'
    for (int i = 0; i < fgm->Nvar; i++) {
        f[i] = nearest->node->point[i];
    }

    return EXIT_SUCCESS;  // Return success after filling the array 'f'
}

int KNNlookupFGM_Interp
(
    FGM *fgm, 
	Node *KDTree, 
	float *x, 
	float *f, 
	int K
)
{
    // Normalize the query point using the same parameters used for the FGM data
    normalizeQueryPoint(x, fgm);
    
    // Create the neighbor heap with capacity K
    NeighborHeap* neighborHeap = createNeighborHeap(K);
    if (neighborHeap == NULL) {
        // Handle the error in case the heap creation failed
        return EXIT_FAILURE;
    } 
	
	// Perfom N nearest neighbour search and fill the heap
    KNNSearch(KDTree, x, 0, fgm->Ncv, neighborHeap);

    // Variables to store the weighted sum and total weight
    float weightedSum[fgm->Nvar];
    float totalWeight = 0.0;

    // Initialize weightedSum to 0
    for (int i = 0; i < fgm->Nvar; i++) {
        weightedSum[i] = 0.0;
    }

    // Loop through the neighbor heap
    for (int i = 0; i < neighborHeap->size; i++) {
		//printf("%f\t%f\n", neighborHeap->elements[i].node->point[0],neighborHeap->elements[i].node->point[1]);
		//printf("%f\t%f\n", neighborHeap->elements[i].node->point[8],neighborHeap->elements[i].node->point[9]);
        float dist = neighborHeap->elements[i].distance;
		//printf("%f\n",dist);
        float weight = 1.0 / (dist + 1e-10); // Adding a small constant to avoid division by zero

        for (int j = 0; j < fgm->Nvar; j++) {
            weightedSum[j] += neighborHeap->elements[i].node->point[j] * weight;
        }
        totalWeight += weight;
    }

    // Compute the final interpolated values
    for (int i = 0; i < fgm->Nvar; i++) {
        f[i] = weightedSum[i] / totalWeight;
        //printf("Value = %f\n", f[i]);
    }

    // Free the neighbor heap
    free(neighborHeap->elements);
    free(neighborHeap);

    return EXIT_SUCCESS;
}

// *************************************************************************** //