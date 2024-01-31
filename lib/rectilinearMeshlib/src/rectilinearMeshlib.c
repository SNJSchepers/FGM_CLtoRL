#include "rectilinearMeshlib.h"

/* --------------------------------------------------------------------------- *\
							Public functions
\* --------------------------------------------------------------------------- */

Mesh createRectilinearMesh
(
    FGM *fgm, 
	Node* KDTree, 
	int* Ngrid, 
	int K
)
{
    printf("\n Building regular rectilinear mesh...\n\n");

    int dimensions = fgm->Ncv;               // Dimensions of the mesh
    int totalVars = fgm->Nvar;               // Total number of variables
    int totalPoints = 1;         
    for (int j = 0; j < dimensions; j++) {
        totalPoints *= Ngrid[j];             // Total number of points
    }
    int meshSize = totalPoints * totalVars;  // meshSize
    
    float *mesh = (float *)malloc(meshSize * sizeof(float));
	float *meshMin = (float *)malloc(dimensions * sizeof(float));
	float *meshMax = (float *)malloc(dimensions * sizeof(float));
    float variables[totalVars];

    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
		
        // Print progress
        printf("Point %d of %d\n", i + 1, totalPoints);
        
        int idx = i * totalVars; // Index in the mesh array for the ith point
        
		float queryPoint[dimensions]; // Fill this with the control variables for the current mesh point
		
        // Calculate the control variable values for this point
        for (int j = 0; j < dimensions; j++) {
            float range = fgm->maxs[j] - fgm->mins[j];
            meshMin[j] = fgm->mins[j] - 0.05 * range;
            meshMax[j] = fgm->maxs[j] + 0.05 * range;
            float extendedRange = meshMax[j] - meshMin[j];
            float delta = extendedRange / (Ngrid[j] - 1);

            int gridIdx = 0;
            int divider = 1;
            for (int k = 0; k < j; k++) {
                divider *= Ngrid[k];
            }
            gridIdx = (i / divider) % Ngrid[j];

            queryPoint[j] = meshMin[j] + gridIdx * delta;
        }
        
        // Perform the K-Nearest Neighbors lookup with inverse distance weighting interpolation
        KNNlookupFGM_Interp(fgm, KDTree, queryPoint, variables, K);
        
        // Denormalize the query point back to original
        denormalizeQueryPoint(queryPoint, fgm);
        
        // Add this point to the mesh
        for (int k = 0; k < totalVars; k++) {
            if (k < dimensions) {
                mesh[i * totalVars + k] = queryPoint[k];
            } else {
                mesh[i * totalVars + k] = variables[k];
            }
        }
    }

	// Allocate and populate the Mesh struct
	Mesh meshStruct;
    meshStruct.reducedMesh = (float*)malloc(meshSize * sizeof(float));
	for (int i = 0; i < meshSize; i++) {
		meshStruct.reducedMesh[i] = mesh[i];
	}
	meshStruct.meshMin = (float*)malloc(dimensions * sizeof(float));
	meshStruct.meshMax = (float*)malloc(dimensions * sizeof(float));
	for (int i = 0; i < dimensions; i++) {
		meshStruct.meshMin[i] = meshMin[i];
		meshStruct.meshMax[i] = meshMax[i];
	}
	free(mesh);    // Free the mesh
	free(meshMin); // Free the mesh min values
	free(meshMax); // Free the mesh max values
	
    printf("\n Rectilinear mesh constructed\n\n");
    
    return meshStruct;
}

Mesh createRectilinearMesh_parallel
(
    FGM *fgm, 
	Node* KDTree, 
	int* Ngrid, 
	int K
)
{
    printf("\n Building regular rectilinear mesh in parallel...\n\n");
    
    int dimensions = fgm->Ncv;                // Dimensions of the mesh
    int totalVars = fgm->Nvar;                // Total number of variables
    int totalPoints = 1;
    for (int j = 0; j < dimensions; j++) {
        totalPoints *= Ngrid[j];              // Total number of points
    }
    int meshSize = totalPoints * totalVars;   // Mesh size              
    
    float *mesh = (float *)malloc(meshSize * sizeof(float));
	float *meshMin = (float *)malloc(dimensions * sizeof(float));
	float *meshMax = (float *)malloc(dimensions * sizeof(float));
    float variables[totalVars];
	
	float queryPoint[dimensions]; // Fill this with the control variables for the current mesh point
	
	int progress = 0; // Shared variable for progress tracking
	
	// Parallel processing
	#pragma omp parallel for private(variables, queryPoint) shared(mesh, progress)
	
    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
        
        int idx = i * totalVars; // Index in the mesh array for the ith point
        
        // Calculate the control variable values for this point
        for (int j = 0; j < dimensions; j++) {
            float range = fgm->maxs[j] - fgm->mins[j];
            meshMin[j] = fgm->mins[j] - 0.05 * range;
            meshMax[j] = fgm->maxs[j] + 0.05 * range;
            float extendedRange = meshMax[j] - meshMin[j];
            float delta = extendedRange / (Ngrid[j] - 1);

            int gridIdx = 0;
            int divider = 1;
            for (int k = 0; k < j; k++) {
                divider *= Ngrid[k];
            }
            gridIdx = (i / divider) % Ngrid[j];

            queryPoint[j] = meshMin[j] + gridIdx * delta;
        }
        
        // Perform the K-Nearest Neighbors lookup with inverse distance weighting interpolation
        KNNlookupFGM_Interp(fgm, KDTree, queryPoint, variables, K);
        
        // Denormalize the query point back to original
        denormalizeQueryPoint(queryPoint, fgm);
        
        // Add this point to the mesh
        for (int k = 0; k < totalVars; k++) {
            if (k < dimensions) {
                mesh[i * totalVars + k] = queryPoint[k];
            } else {
                mesh[i * totalVars + k] = variables[k];
            }
        }
		
		// Critical section to update and print progress
		#pragma omp critical
		{
			progress++;
			if (progress % 100 == 0) {  // Print progress every 100 points, for example
				printf("Progress: %d of %d points processed\n", progress, totalPoints);
			}
		}
	
    }

	// Allocate and populate the Mesh struct
	Mesh meshStruct;
    meshStruct.reducedMesh = (float*)malloc(meshSize * sizeof(float));
	for (int i = 0; i < meshSize; i++) {
		meshStruct.reducedMesh[i] = mesh[i];
	}
	meshStruct.meshMin = (float*)malloc(dimensions * sizeof(float));
	meshStruct.meshMax = (float*)malloc(dimensions * sizeof(float));
	for (int i = 0; i < dimensions; i++) {
		meshStruct.meshMin[i] = meshMin[i];
		meshStruct.meshMax[i] = meshMax[i];
	}
	free(mesh); // Free the temporary mesh
	
    printf("\n Rectilinear mesh constructed\n\n");
    
    return meshStruct;
}

Mesh createRectilinearMesh_reduced
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K, 
	float threshold
)
{	
	printf("\n Building reduced rectilinear mesh...\n\n");
	
    int dimensions = fgm->Ncv;                         // Dimensions of the mesh
    int totalVars = fgm->Nvar;                         // Total number of variables
    int totalPoints = 1;
	for (int j = 0; j < dimensions; j++) {
		totalPoints *= Ngrid[j];                       // Total number of points
	}
    int tempMeshSize = totalPoints * totalVars;        // Temporary mesh size
	int mapMeshSize = totalPoints * (dimensions + 1);  // Original mesh size
	int validPointCount = 0;                           // Initialise the number of valid points
	
	int *mapMesh     = (int *)malloc(mapMeshSize * sizeof(int));
    float *tempMesh = (float *)malloc(tempMeshSize * sizeof(float));
	
	float *meshMin = (float *)malloc(dimensions * sizeof(float));
	float *meshMax = (float *)malloc(dimensions * sizeof(float));
	
	float variables[totalVars];

    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
		
		// Print progress
		printf("Point %d of %d\n", i+1, totalPoints);
		
        int idx = i * totalVars; // Index in the mesh array for the ith point
		
		float queryPoint[dimensions]; // Fill this with the control variables for the current mesh point
		
		// Calculate the control variable values for this point
		for (int j = 0; j < dimensions; j++) {
			int N = Ngrid[j];
			meshMin[j] = fgm->mins[j] - 0.05 * (fgm->maxs[j] - fgm->mins[j]);
			meshMax[j] = fgm->maxs[j] + 0.05 * (fgm->maxs[j] - fgm->mins[j]);
			float extendedRange = meshMax[j] - meshMin[j];
			float delta = extendedRange / (N - 1);

			int gridIdx = 0;
			int divider = 1;
			for (int k = 0; k < j; k++) {
				divider *= Ngrid[k];
			}
			gridIdx = (i / divider) % N;

			queryPoint[j] = meshMin[j] + gridIdx * delta;
			mapMesh[i * (dimensions + 1) + j] = gridIdx;
		}
		
		if (i == 882){
			printf("%f\t%f\n",queryPoint[0],queryPoint[1]);
			//exit(EXIT_FAILURE);
		}
		
		// Normalize the query point using the same parameters used for the FGM data
		normalizeQueryPoint(queryPoint, fgm);
		
		if (i == 882){
			printf("%f\t%f\n",queryPoint[0],queryPoint[1]);
			//exit(EXIT_FAILURE);
		}
		//printf("%f\t%f\n",queryPoint[0],queryPoint[1]);
		
		// Perform nearest neighbor search
		Neighbor nearest;
		nearestNeighborSearch(KDTree, &nearest, queryPoint, 0, dimensions);
		
		//printf("%f\t%f\n",nearest.node->point[0],nearest.node->point[1]);
		
		//printf("%f\n",nearest.distance);
		
		//exit(EXIT_FAILURE);
		// Denormalize the query point back to original
		denormalizeQueryPoint(queryPoint, fgm);
	
		if (i == 882){
			printf("%f\t%f\n",nearest.distance,threshold);
			//exit(EXIT_FAILURE);
		}
		
		// Check if distance exceeds the threshold
		if (nearest.distance < threshold) {

			// Perform the K-Nearest Neighbors lookup with inverse distance weighting interpolation
			KNNlookupFGM_Interp(fgm, KDTree, queryPoint, variables, K);
			
			if (i == 882){
				printf("%f\t%f\n",variables[0],variables[1]);
				printf("%f\t%f\n",variables[8],variables[9]);
				//exit(EXIT_FAILURE);
			}
		
			
			// Denormalize the query point back to original
			denormalizeQueryPoint(queryPoint, fgm);
			
			// Add this point to tempMesh and increment validPointCount
			for (int k = 0; k < totalVars; k++) {
				if (k < dimensions){
					tempMesh[validPointCount * totalVars + k] = queryPoint[k];
				}
				else {
					tempMesh[validPointCount * totalVars + k] = variables[k];
				}
			}
			validPointCount++;
			
			mapMesh[i * (dimensions + 1) + dimensions] = validPointCount - 1; // -1 to adjust for 0-based indexing
		} else {
			// Store -1 in mapMesh to indicate the point is outside the threshold
			mapMesh[i * (dimensions + 1) + dimensions] = -1;
		}
		
		//printf("mapMesh (%d,%d,%d)\n", mapMesh[i * (dimensions + 1) + 0], mapMesh[i * (dimensions + 1) + 1], mapMesh[i * (dimensions + 1) + dimensions]);
    }
	
	int meshSize = validPointCount * totalVars;

	// Allocate and populate the Mesh struct
    Mesh meshStruct;
    meshStruct.reducedMesh = (float*)malloc(validPointCount * totalVars * sizeof(float));
	for (int i = 0; i < validPointCount * totalVars; i++) {
		meshStruct.reducedMesh[i] = tempMesh[i];
	}
	meshStruct.mapMesh     = (int*)malloc(mapMeshSize * sizeof(int));
	for (int i = 0; i < mapMeshSize; i++) {
		meshStruct.mapMesh[i] = mapMesh[i];
	}
    //meshStruct.mapMesh = mapMesh; // Assuming mapMesh is already allocated
    meshStruct.validPointCount = validPointCount;
	
	meshStruct.meshMin = (float*)malloc(dimensions * sizeof(float));
	meshStruct.meshMax = (float*)malloc(dimensions * sizeof(float));
	for (int i = 0; i < dimensions; i++) {
		meshStruct.meshMin[i] = meshMin[i];
		meshStruct.meshMax[i] = meshMax[i];
	}
	
	free(tempMesh); // Free the temporary mesh
	free(mapMesh);  // Free the map mesh
	free(meshMin);  // Free the mesh minimum values
	free(meshMax);  // Free the mesh maximum values
	
	printf("\n Rectilinear mesh constructed\n\n");
	
    return meshStruct;
}

Mesh createRectilinearMesh_reduced_parallel
(
    FGM *fgm, 
	Node* KDTree, 
	int *Ngrid, 
	int K, 
	float threshold
)
{	
	printf("\n Building reduced rectilinear mesh in parallel...\n\n");
	
    int dimensions = fgm->Ncv;                         // Dimensions of the mesh
    int totalVars = fgm->Nvar;                         // Total number of variables
    int totalPoints = 1;
	for (int j = 0; j < dimensions; j++) {
		totalPoints *= Ngrid[j];                       // Total number of points
	}
    int tempMeshSize = totalPoints * totalVars;        // Temporary mesh size
	int mapMeshSize = totalPoints * (dimensions + 1);  // Original mesh size
	int validPointCount = 0;                           // Initialise the number of valid points
	
	int *mapMesh     = (int *)malloc(mapMeshSize * sizeof(int));
    float *tempMesh = (float *)malloc(tempMeshSize * sizeof(float));
	
	float *meshMin = (float *)malloc(dimensions * sizeof(float));
	float *meshMax = (float *)malloc(dimensions * sizeof(float));
	
	float variables[totalVars];
	float queryPoint[dimensions];
	
	Neighbor nearest; // Create instance of Neighbor struct
	
	int progress = 0; // Shared variable for progress tracking

    // Enable parallel processing
	#pragma omp parallel for private(variables, queryPoint, nearest) shared(tempMesh, mapMesh, progress,  validPointCount)
	
    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
		
        int idx = i * totalVars; // Index in the mesh array for the ith point
		
		// Calculate the control variable values for this point
		for (int j = 0; j < dimensions; j++) {
			int N = Ngrid[j];
			meshMin[j] = fgm->mins[j] - 0.05 * (fgm->maxs[j] - fgm->mins[j]);
			meshMax[j] = fgm->maxs[j] + 0.05 * (fgm->maxs[j] - fgm->mins[j]);
			float extendedRange = meshMax[j] - meshMin[j];
			float delta = extendedRange / (N - 1);

			int gridIdx = 0;
			int divider = 1;
			for (int k = 0; k < j; k++) {
				divider *= Ngrid[k];
			}
			gridIdx = (i / divider) % N;

			queryPoint[j] = meshMin[j] + gridIdx * delta;
			mapMesh[i * (dimensions + 1) + j] = gridIdx;
		}
		
		// Normalize the query point using the same parameters used for the FGM data
		normalizeQueryPoint(queryPoint, fgm);
		
		// Perform nearest neighbor search
		nearestNeighborSearch(KDTree, &nearest, queryPoint, 0, dimensions);
		
		// Denormalize the query point back to original
		denormalizeQueryPoint(queryPoint, fgm);

		// Check if distance exceeds the threshold
		if (nearest.distance < threshold) {

			// Perform the K-Nearest Neighbors lookup with inverse distance weighting interpolation
			KNNlookupFGM_Interp(fgm, KDTree, queryPoint, variables, K);
			
			// Denormalize the query point back to original
			denormalizeQueryPoint(queryPoint, fgm);
			
			// Add this point to tempMesh
			int localValidPointCount; // Local variable to store the current valid point count
			#pragma omp critical
			{
				localValidPointCount = validPointCount;
				validPointCount++;
			}

			for (int k = 0; k < totalVars; k++) {
				if (k < dimensions) {
					tempMesh[localValidPointCount * totalVars + k] = queryPoint[k];
				} else {
					tempMesh[localValidPointCount * totalVars + k] = variables[k];
				}
			}

			mapMesh[i * (dimensions + 1) + dimensions] = localValidPointCount; // Adjusted for 0-based indexing
			
		} else {
			// Store -1 in mapMesh to indicate the point is outside the threshold
			mapMesh[i * (dimensions + 1) + dimensions] = -1;
		}
		
		// Critical section to update and print progress
		#pragma omp critical
		{
			progress++;
			if (progress % 100 == 0) {  // Print progress every 100 points, for example
				printf("Progress: %d of %d points processed\n", progress, totalPoints);
			}
		}
    }
	
	int meshSize = validPointCount * totalVars;

	// Allocate and populate the Mesh struct
    Mesh meshStruct;
    meshStruct.reducedMesh = (float*)malloc(validPointCount * totalVars * sizeof(float));
	for (int i = 0; i < validPointCount * totalVars; i++) {
		meshStruct.reducedMesh[i] = tempMesh[i];
	}
	meshStruct.mapMesh     = (int*)malloc(mapMeshSize * sizeof(int));
	for (int i = 0; i < mapMeshSize; i++) {
		meshStruct.mapMesh[i] = mapMesh[i];
	}
    //meshStruct.mapMesh = mapMesh; // Assuming mapMesh is already allocated
    meshStruct.validPointCount = validPointCount;
	
	meshStruct.meshMin = (float*)malloc(dimensions * sizeof(float));
	meshStruct.meshMax = (float*)malloc(dimensions * sizeof(float));
	for (int i = 0; i < dimensions; i++) {
		meshStruct.meshMin[i] = meshMin[i];
		meshStruct.meshMax[i] = meshMax[i];
	}

	free(tempMesh); // Free the temporary mesh
	free(mapMesh);  // Free the map mesh
	free(meshMin);  // Free the mesh minimum values
	free(meshMax);  // Free the mesh maximum values
	
	printf("\n Rectilinear mesh constructed\n\n");
	
    return meshStruct;
}