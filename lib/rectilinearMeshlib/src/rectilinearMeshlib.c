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
    
    double *mesh = (double *)malloc(meshSize * sizeof(double));
	double *meshMin = (double *)malloc(dimensions * sizeof(double));
	double *meshMax = (double *)malloc(dimensions * sizeof(double));
    double variables[totalVars];

    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
		
        // Print progress
        printf("Point %d of %d\n", i + 1, totalPoints);
        
        int idx = i * totalVars; // Index in the mesh array for the ith point
        
		double queryPoint[dimensions]; // Fill this with the control variables for the current mesh point
		
        // Calculate the control variable values for this point
        for (int j = 0; j < dimensions; j++) {
            double range = fgm->maxs[j] - fgm->mins[j];
            meshMin[j] = fgm->mins[j] - 0.1 * range;
            meshMax[j] = fgm->maxs[j] + 0.1 * range;
            double extendedRange = meshMax[j] - meshMin[j];
            double delta = extendedRange / (Ngrid[j] - 1);

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
        denormalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);
        
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
    meshStruct.reducedMesh = (double*)malloc(meshSize * sizeof(double));
	for (int i = 0; i < meshSize; i++) {
		meshStruct.reducedMesh[i] = mesh[i];
	}
	meshStruct.meshMin = (double*)malloc(dimensions * sizeof(double));
	meshStruct.meshMax = (double*)malloc(dimensions * sizeof(double));
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
    
    double *mesh = (double *)malloc(meshSize * sizeof(double));
	double *meshMin = (double *)malloc(dimensions * sizeof(double));
	double *meshMax = (double *)malloc(dimensions * sizeof(double));
    double variables[totalVars];
	
	double queryPoint[dimensions]; // Fill this with the control variables for the current mesh point
	
	int progress = 0; // Shared variable for progress tracking
	
	// Parallel processing
	#pragma omp parallel for private(variables, queryPoint) shared(mesh, progress)
	
    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
        
        int idx = i * totalVars; // Index in the mesh array for the ith point
        
        // Calculate the control variable values for this point
        for (int j = 0; j < dimensions; j++) {
            double range = fgm->maxs[j] - fgm->mins[j];
            meshMin[j] = fgm->mins[j] - 0.1 * range;
            meshMax[j] = fgm->maxs[j] + 0.1 * range;
            double extendedRange = meshMax[j] - meshMin[j];
            double delta = extendedRange / (Ngrid[j] - 1);

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
        denormalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);
        
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
    meshStruct.reducedMesh = (double*)malloc(meshSize * sizeof(double));
	for (int i = 0; i < meshSize; i++) {
		meshStruct.reducedMesh[i] = mesh[i];
	}
	meshStruct.meshMin = (double*)malloc(dimensions * sizeof(double));
	meshStruct.meshMax = (double*)malloc(dimensions * sizeof(double));
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
	double threshold
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
    double *tempMesh = (double *)malloc(tempMeshSize * sizeof(double));
	
	double *meshMin = (double *)malloc(dimensions * sizeof(double));
	double *meshMax = (double *)malloc(dimensions * sizeof(double));
	
	double variables[totalVars];

    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
		
		// Print progress
		printf("Point %d of %d\n", i+1, totalPoints);
		
        int idx = i * totalVars; // Index in the mesh array for the ith point
		
		double queryPoint[dimensions]; // Fill this with the control variables for the current mesh point
		
		// Calculate the control variable values for this point
		for (int j = 0; j < dimensions; j++) {
			int N = Ngrid[j];
			meshMin[j] = fgm->mins[j] - 0.1 * (fgm->maxs[j] - fgm->mins[j]);
			meshMax[j] = fgm->maxs[j] + 0.1 * (fgm->maxs[j] - fgm->mins[j]);
			double extendedRange = meshMax[j] - meshMin[j];
			double delta = extendedRange / (N - 1);

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
		normalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);
		
		// Perform nearest neighbor search
		Neighbor nearest;
		nearestNeighborSearch(KDTree, &nearest, queryPoint, 0, dimensions);
		
		// Denormalize the query point back to original
		denormalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);

		// Check if distance exceeds the threshold
		if (nearest.distance < threshold) {

			// Perform the K-Nearest Neighbors lookup with inverse distance weighting interpolation
			KNNlookupFGM_Interp(fgm, KDTree, queryPoint, variables, K);
			
			// Denormalize the query point back to original
			denormalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);
			
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
    meshStruct.reducedMesh = (double*)malloc(validPointCount * totalVars * sizeof(double));
	for (int i = 0; i < validPointCount * totalVars; i++) {
		meshStruct.reducedMesh[i] = tempMesh[i];
	}
	meshStruct.mapMesh     = (int*)malloc(mapMeshSize * sizeof(int));
	for (int i = 0; i < mapMeshSize; i++) {
		meshStruct.mapMesh[i] = mapMesh[i];
	}
    //meshStruct.mapMesh = mapMesh; // Assuming mapMesh is already allocated
    meshStruct.validPointCount = validPointCount;
	
	meshStruct.meshMin = (double*)malloc(dimensions * sizeof(double));
	meshStruct.meshMax = (double*)malloc(dimensions * sizeof(double));
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
	double threshold
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
    double *tempMesh = (double *)malloc(tempMeshSize * sizeof(double));
	
	double *meshMin = (double *)malloc(dimensions * sizeof(double));
	double *meshMax = (double *)malloc(dimensions * sizeof(double));
	
	double variables[totalVars];
	double queryPoint[dimensions];
	
	Neighbor nearest; // Create instance of Neighbor struct
	
	int progress = 0; // Shared variable for progress tracking

    // Enable parallel processing
	#pragma omp parallel for private(variables, queryPoint, nearest) shared(tempMesh, mapMesh, progress)
	
    // Populate control variable values for each mesh point
    for (int i = 0; i < totalPoints; i++) {
		
        int idx = i * totalVars; // Index in the mesh array for the ith point
		
		// Calculate the control variable values for this point
		for (int j = 0; j < dimensions; j++) {
			int N = Ngrid[j];
			meshMin[j] = fgm->mins[j] - 0.1 * (fgm->maxs[j] - fgm->mins[j]);
			meshMax[j] = fgm->maxs[j] + 0.1 * (fgm->maxs[j] - fgm->mins[j]);
			double extendedRange = meshMax[j] - meshMin[j];
			double delta = extendedRange / (N - 1);

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
		normalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);
		
		// Perform nearest neighbor search
		
		nearestNeighborSearch(KDTree, &nearest, queryPoint, 0, dimensions);
		
		// Denormalize the query point back to original
		denormalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);

		// Check if distance exceeds the threshold
		if (nearest.distance < threshold) {

			// Perform the K-Nearest Neighbors lookup with inverse distance weighting interpolation
			KNNlookupFGM_Interp(fgm, KDTree, queryPoint, variables, K);
			
			// Denormalize the query point back to original
			denormalizeQueryPoint(queryPoint, fgm->mins, fgm->maxs, fgm->Ncv);
			
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
    meshStruct.reducedMesh = (double*)malloc(validPointCount * totalVars * sizeof(double));
	for (int i = 0; i < validPointCount * totalVars; i++) {
		meshStruct.reducedMesh[i] = tempMesh[i];
	}
	meshStruct.mapMesh     = (int*)malloc(mapMeshSize * sizeof(int));
	for (int i = 0; i < mapMeshSize; i++) {
		meshStruct.mapMesh[i] = mapMesh[i];
	}
    //meshStruct.mapMesh = mapMesh; // Assuming mapMesh is already allocated
    meshStruct.validPointCount = validPointCount;
	
	meshStruct.meshMin = (double*)malloc(dimensions * sizeof(double));
	meshStruct.meshMax = (double*)malloc(dimensions * sizeof(double));
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