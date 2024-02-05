#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "FGMlib.h"
#include "KdTreelib.h"
#include "KNNlookuplib.h"
#include "rectilinearMeshlib.h"

int main(void)
{
	// Settings
	
		// Number of nearest neighbors used for interpolation
		int K = 6;
		
		// Grid points of the rectilinear mesh (must be equal to the number of control variables
		int Ngrid[] = {300,300};
		
		// Boolean to exclude points that lay outside of the curvilinear mesh to reduce memory
		int memory_reduction = 1;
	
		// Distance threshold for which points are included in the rectlininear mesh
		float threshold = 0.001;

		// Boolean for parallel processing
		int parallel = 1;
	
	// Read FGM and create the KD Tree for nearest neighbor interpolation
		
		// FGM object
		FGM *fgm;
		
		// KD Tree object
		Node* KDTree;
	
		// Read the FGM database
		fgm = readFGM("data/input/database.fgm.2D.PV_h");
		
		// Build the KD Tree
		KDTree = buildKDTree(fgm);
		
	// Build the rectilinear mesh
		
		// Check if parellel processing is selected
		if (parallel) {
			
			// Check for memory reduction
			if (memory_reduction) {
				
				// Start the mesh generation timer
				float start = omp_get_wtime();
				
				// Create the rectilinear mesh with memory reduction
				Mesh mesh = createRectilinearMesh_reduced_parallel(fgm, KDTree, Ngrid, K, threshold);
				
				// End the mesh generation timer
				float end = omp_get_wtime();
			
				// Print elapsed time on screen
				printf("Time to create mesh: %.2f seconds\n", end - start);
			
				// Create new FGM object for the reduced rectilinear mesh
				FGM *fgm_RL_recuded;
				
				// Create new FGM instance with the mapping of the original mesh on the reduced mesh
				FGM_MAP *fgm_mapping;
				
				// The number of grid points of the reduced rectilinear mesh
				int Ngrid_reduced[1] = {mesh.validPointCount};
				
				// Fill the FGM structure with the reduced rectlininear mesh 
				fgm_RL_recuded = createFGM(1, Ngrid_reduced, fgm->gridpower, fgm->Nvar, fgm->varname, mesh.reducedMesh);
				
				// Fill the FGM structure with the original mesh mapping
				fgm_mapping = createFGM_MAP(fgm->Ncv, Ngrid, fgm->gridpower, mesh.meshMin, mesh.meshMax, fgm->varname, mesh.mapMesh);
			
				// free the original fgm memory
				freeFGM(fgm);
				
				// Write the reduced FGM structure to a text file
				writeFGM(fgm_RL_recuded, "data/output/database.fgm.2D_RL_reduced");
				
				// Write the original mapping FGM structure to a text file
				writeFGM_MAP(fgm_mapping, "data/output/database.fgm.2D_RL_mapping");
				
				// free the reduced FGM memory
				freeFGM(fgm_RL_recuded);
				
				// free the original mapping FGM memory
				freeFGM_MAP(fgm_mapping);
				
			// If no memory reduction is selected
			} else {
				
				// Start the mesh generation timer
				float start = omp_get_wtime();
				
				// Create the rectilinear mesh without memory reduction
				Mesh mesh = createRectilinearMesh_parallel(fgm, KDTree, Ngrid, K);

				// End the mesh generation timer
				float end = omp_get_wtime();
			
				// Print elapsed time on screen
				printf("Time to create mesh: %.2f seconds\n", end - start);
				
				// Create new FGM instance for the rectilinear mesh
				FGM *fgm_RL;

				// Fill the FGM structure with the rectilinear mesh
				fgm_RL = createFGM(fgm->Ncv, Ngrid, fgm->gridpower, fgm->Nvar, fgm->varname, mesh.reducedMesh);
				
				// free the original fgm memory
				freeFGM(fgm);
				
				// Write the new FGM to a text filelength
				writeFGM(fgm_RL, "data/output/database.fgm.RL_2");
				
				// free the rectilinear fgm memory
				freeFGM(fgm_RL);
			}
			
		// If no parallel processing is selected
		} else {
			
			// Check for memory reduction
			if (memory_reduction) {
				
				// Start the mesh generation timer
				clock_t start = clock();
				
				// Create the rectilinear mesh with memory reduction
				Mesh mesh = createRectilinearMesh_reduced(fgm, KDTree, Ngrid, K, threshold);
				//Mesh mesh = createRectilinearMesh_reduced_parallel(fgm, KDTree, Ngrid, K, threshold);
				
				// End the mesh generation timer
				clock_t end = clock();
			
				// Print elapsed time on screen
				printf("Time to create mesh: %.2f seconds\n", ((float)(end - start)) / CLOCKS_PER_SEC);
			
				// Create new FGM object for the reduced rectilinear mesh
				FGM *fgm_RL_recuded;
				
				// Create new FGM instance with the mapping of the original mesh on the reduced mesh
				FGM_MAP *fgm_mapping;
				
				// The number of grid points of the reduced rectilinear mesh
				int Ngrid_reduced[1] = {mesh.validPointCount};
				
				// Fill the FGM structure with the reduced rectlininear mesh 
				fgm_RL_recuded = createFGM(1, Ngrid_reduced, fgm->gridpower, fgm->Nvar, fgm->varname, mesh.reducedMesh);
				
				// Fill the FGM structure with the original mesh mapping
				fgm_mapping = createFGM_MAP(fgm->Ncv, Ngrid, fgm->gridpower, mesh.meshMin, mesh.meshMax, fgm->varname, mesh.mapMesh);
			
				// free the original fgm memory
				freeFGM(fgm);
				
				// Write the reduced FGM structure to a text file
				writeFGM(fgm_RL_recuded, "data/output/database.fgm.2D_RL_reduced");
				
				// Write the original mapping FGM structure to a text file
				writeFGM_MAP(fgm_mapping, "data/output/database.fgm.2D_RL_mapping");
				
				// free the reduced FGM memory
				freeFGM(fgm_RL_recuded);
				
				// free the original mapping FGM memory
				freeFGM_MAP(fgm_mapping);
				
			// If no memory reduction is selected
			} else {
				
				// Start the mesh generation timer
				clock_t start = clock();
				
				// Create the rectilinear mesh without memory reduction
				Mesh mesh = createRectilinearMesh(fgm, KDTree, Ngrid, K);

				// End the mesh generation timer
				clock_t end = clock();
				
				// Print elapsed time on screen
				printf("Time to create mesh: %.2f seconds\n", ((float)(end - start)) / CLOCKS_PER_SEC);
				
				// Create new FGM instance for the rectilinear mesh
				FGM *fgm_RL;

				// Fill the FGM structure with the rectilinear mesh
				fgm_RL = createFGM(fgm->Ncv, Ngrid, fgm->gridpower, fgm->Nvar, fgm->varname, mesh.reducedMesh);
				
				// free the original fgm memory
				freeFGM(fgm);
				
				// Write the new FGM to a text filelength
				writeFGM(fgm_RL, "data/output/database.fgm.RL");
				
				// free the rectilinear fgm memory
				freeFGM(fgm_RL);
			}
		}
	
	return EXIT_SUCCESS;
}
