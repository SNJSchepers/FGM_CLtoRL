
# include "FGMlib.h"

/* --------------------------------------------------------------------------- *\
							Public functions
\* --------------------------------------------------------------------------- */

FGM * readFGM
(
    const char filename[]
)
{
    char line[MAX_LINE_LENGTH];
    FILE *fid;
  
    int Ncv = 0;
    int *Ngrid = NULL;
    float *gridpower = NULL;
    int Nvar = 0;
    char (*varname)[VAR_NAME_LENGTH];
    float *data = NULL;

    FGM *fgm = NULL;
  
    int i,j;
    float f;

    // Open the file
    fid = fopen(filename, "r");
    if (fid == NULL) {
        fprintf(stderr, "Error: Failed to open FGM database file\n");
        exit(EXIT_FAILURE);
    }
  
    // Check the identifier [FGM]
    fgets(line,MAX_LINE_LENGTH,fid);
    if (strncmp(line,"[FGM]",5) != 0) {
        fprintf(stderr, "Error: Incorrect FGM database format\n");
        exit(EXIT_FAILURE);
    }
  
    // Read the number of cv's
    if (locateKeyWord(fid,"[DIMENSION]")) {
        fprintf(stderr, "Error: Couldn't find keyword DIMENSION\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fid,"%i",&Ncv);
    printf("Ncv = %i\n",Ncv);
  
    // Read space for the data size
    Ngrid = (int *)malloc(sizeof(int) * Ncv);
    if (Ngrid == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for data size\n");
        exit(EXIT_FAILURE);
    }
  
    // Read the data size */
    if (locateKeyWord(fid,"[DATASIZE]")) {
        fprintf(stderr, "Error: Couldn't find keyword DATASIZE\n");
        exit(EXIT_FAILURE);
    }
    int Ntotal = 1;
    for (i = 0; i < Ncv; i++) {
        fscanf(fid,"%i",&j);
        Ngrid[i] = j;
        printf("Ngrid[%i] = %i\n",i,Ngrid[i]);
        Ntotal = Ntotal * Ngrid[i];
    }
 
    // Read the number of variables
    if (locateKeyWord(fid,"[VARIABLES]")) {
        fprintf(stderr, "Error: Couldn't find keyword VARIABLES\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fid,"%i",&Nvar);
    printf("Nvar = %i\n",Nvar);

    // Allocate space for the variable names
    varname = malloc(sizeof *varname * Nvar);
    if (varname == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for variable names\n");
        exit(EXIT_FAILURE);
    }
  
    // Read the variable names 
    for (j = 0; j < Nvar; j++) {
        fscanf(fid,"%s",varname[j]);
        printf("%s\n",varname[j]);
    }
  
    // Allocate space for the data
    data = malloc(sizeof(float) * Ntotal * Nvar);
    if (data == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for data\n");
        exit(EXIT_FAILURE);
    }
  
    // Read the data
    if (locateKeyWord(fid,"[DATA]")) {
        fprintf(stderr, "Error: Couldn't find keyword DATA\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < Ntotal; i++) {
        for (j = 0; j < Nvar; j++) {
            fscanf(fid,"%f",&f);
            data[i*Nvar + j] = f;
        }
    }
  
    // Allocate space for grid powers
    gridpower = (float *)malloc(sizeof(float) * Ncv);
    if (gridpower == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for gridpower\n");
        exit(EXIT_FAILURE);
    }
  
    // Read grid powers for non-linear distribution of mesh points
    if (locateKeyWord(fid,"[GRIDPOWER]")) {
        fprintf(stderr, "Warning: Couldn't find keyword GRIDPOWER\n");
        for (i = 0; i < Ncv; i++) {
            gridpower[i] = -1.0; // Set to negative value
        }
    }
    else {
        for (i = 0; i < Ncv; i++) {
            fscanf(fid,"%f",&f);
            // If the power is equal to 1 then set it to a negative value
            if ( fabs(f-1.0) < 1.e-6 ) {
                gridpower[i] = -1.0;
            }
            else {
                gridpower[i] = 1.0 / f; // Store the reciprocal values
            }
            printf("Gridpower[%i] = %f, %f\n", i, gridpower[i], f);
        }
    }
   
    // Close file
    fclose(fid);
     
    // Allocate memory for FGM
    fgm = (FGM *)malloc(sizeof(FGM));
    if (fgm == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for fgm\n");
        exit(EXIT_FAILURE);
    }


    // Assign values
    fgm->Ncv = Ncv;
    fgm->Ngrid = Ngrid;
    fgm->gridpower = gridpower;
    fgm->Nvar = Nvar;
    fgm->varname = varname;
    fgm->data = data;

    // Allocate memory for mins and maxs arrays
    fgm->mins = malloc(fgm->Ncv * sizeof(float));
    fgm->maxs = malloc(fgm->Ncv * sizeof(float));

    // Check if memory allocation was successful
    if (!fgm->mins || !fgm->maxs) {
        fprintf(stderr, "Error: Unable to allocate space for mins or maxs\n");
        exit(EXIT_FAILURE);
    }

    // Initialize mins and maxs
    for (int i = 0; i < fgm->Ncv; i++) {
        fgm->mins[i] = INFINITY;
        fgm->maxs[i] = -INFINITY;
    }

    // Calculate the total number of points
    int totalPoints = 1;
    for (int i = 0; i < fgm->Ncv; i++) {
        totalPoints *= fgm->Ngrid[i];
    }

    // Iterate over the data to find min and max values for each control variable
    for (int i = 0; i < totalPoints; i++) {
        for (int j = 0; j < fgm->Ncv; j++) {
            float val = fgm->data[i * fgm->Nvar + j];
            if (val < fgm->mins[j]) fgm->mins[j] = val;
            if (val > fgm->maxs[j]) fgm->maxs[j] = val;
        }
    }

    return fgm;
};

int freeFGM
(
    FGM *fgm
)
{
    if (fgm != NULL) {
        if (fgm->Ngrid != NULL) { free(fgm->Ngrid); }
        if (fgm->gridpower != NULL) { free(fgm->gridpower); }
        if (fgm->data != NULL) { free(fgm->data); }
        free(fgm);
    } else {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
};

int freeFGM_MAP(FGM_MAP *fgm)
{
    if (fgm != NULL) {
        if (fgm->Ngrid != NULL) { free(fgm->Ngrid); }
        if (fgm->gridpower != NULL) { free(fgm->gridpower); }
        if (fgm->data != NULL) { free(fgm->data); }
        free(fgm);
    } else {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
};

int locateKeyWord(FILE *fid, char keyword[])
{
    char line[MAX_LINE_LENGTH];
    int keywordlength;
  
    keywordlength = strlen(keyword);
  
    rewind(fid);
    do {
        if (fgets(line,MAX_LINE_LENGTH,fid) == NULL) {
            return EXIT_FAILURE;
        }
    } while (strncmp(line,keyword,keywordlength) != 0);
  
    return EXIT_SUCCESS;
};

FGM *createFGM
(
    int Ncv, 
	int *Ngrid, 
	float *gridpower, 
	int Nvar, 
	char (*varname)[VAR_NAME_LENGTH], 
	float *data
)
{
	// Create new FGM instance
    FGM *newFGM = (FGM *)malloc(sizeof(FGM));
    if (newFGM == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for the new FGM\n");
        return NULL;
    }
	
	// Populate the new FGM struct
    newFGM->Ncv = Ncv;
    newFGM->Ngrid = (int *)malloc(Ncv * sizeof(int));
    newFGM->gridpower = (float *)malloc(Ncv * sizeof(float));
    newFGM->Nvar = Nvar;
    newFGM->varname = (char (*)[VAR_NAME_LENGTH])malloc(Nvar * VAR_NAME_LENGTH * sizeof(char));
	
	int dataSize = 1;
	for (int j = 0; j < Ncv; j++) {
		dataSize *= Ngrid[j]; // Calculate total points for all dimensions
	}
	dataSize *= Nvar;

    newFGM->data = (float *)malloc( dataSize * sizeof(float) );
	
    // Check for memory allocation failure for each array
    if (newFGM->Ngrid == NULL || newFGM->gridpower == NULL || newFGM->varname == NULL || newFGM->data == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for one or more of the new FGM variables\n");        
        freeFGM(newFGM);
        return NULL;
    }
	
    // Copy data into the new FGM struct
    memcpy(newFGM->Ngrid, Ngrid, Ncv * sizeof(int));
    memcpy(newFGM->gridpower, gridpower, Ncv * sizeof(float));
    for (int i = 0; i < Nvar; i++) {
        strcpy(newFGM->varname[i], varname[i]);
    }
    memcpy(newFGM->data, data, dataSize * sizeof(float));

    return newFGM;
}

FGM_MAP *createFGM_MAP
(
    int Ncv, 
	int *Ngrid, 
	float *gridpower, 
	float *mins, 
	float *maxs, 
	char (*varname)[VAR_NAME_LENGTH], 
	int *data
)
{
    // Create new FGM_MAP instance	
    FGM_MAP *newFGM_MAP = (FGM_MAP *)malloc(sizeof(FGM_MAP));
    if (newFGM_MAP == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for the new FGM\n");
        return NULL;
    }
	
	// Populate the new FGM_MAP struct
    newFGM_MAP->Ncv = Ncv;
    newFGM_MAP->Ngrid = (int *)malloc(Ncv * sizeof(int));
    newFGM_MAP->gridpower = (float *)malloc(Ncv * sizeof(float));
	newFGM_MAP->mins = (float *)malloc(Ncv * sizeof(float));
	newFGM_MAP->maxs = (float *)malloc(Ncv * sizeof(float));
    newFGM_MAP->Nvar = Ncv + 1;
    newFGM_MAP->varname = (char (*)[VAR_NAME_LENGTH])malloc(Ncv + 1 * VAR_NAME_LENGTH * sizeof(char));

	int dataSize = 1;
	for (int j = 0; j < Ncv; j++) {
		dataSize *= Ngrid[j]; // Calculate total points for all dimensions
	}
	dataSize *= Ncv + 1;
 
    newFGM_MAP->data = (int *)malloc( dataSize * sizeof(int) );

    // Check for memory allocation failure for each array
    if (newFGM_MAP->Ngrid == NULL || newFGM_MAP->gridpower == NULL || newFGM_MAP->varname == NULL || newFGM_MAP->data == NULL) {
        fprintf(stderr, "Error: Unable to allocate space for one or more of the new FGM variables\n");        
        freeFGM_MAP(newFGM_MAP);
        return NULL;
    }

    // Copy data into the new FGM struct
    memcpy(newFGM_MAP->Ngrid, Ngrid, Ncv * sizeof(int));
    memcpy(newFGM_MAP->gridpower, gridpower, Ncv * sizeof(float));

	memcpy(newFGM_MAP->mins, mins, Ncv * sizeof(float));
	memcpy(newFGM_MAP->maxs, maxs, Ncv * sizeof(float));

    for (int i = 0; i < (Ncv + 1); i++) {
        strcpy(newFGM_MAP->varname[i], varname[i]);
    }

    memcpy(newFGM_MAP->data, data, dataSize * sizeof(int));

    return newFGM_MAP;
}

int writeFGM
(
    const FGM *fgm, 
	const char *filename
)
{	
    printf("\n Writing FGM ascii file...\n\n");

    if (fgm == NULL || filename == NULL) {
        fprintf(stderr, "Error: fgm or filename is empty\n");        
        return EXIT_FAILURE;
    }

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Failed to open FGM database file\n");
        return EXIT_FAILURE;
    }

    // Write FGM data to the file
	fprintf(file, "[FGM]\n");
	fprintf(file, "NDRL\n");
	fprintf(file, "[DIMENSION]\n");
	fprintf(file, "%d\n", fgm->Ncv);
	fprintf(file, "[DATASIZE]\n");
	
	// Write the fgm->Ngrid values
	for (int i = 0; i < fgm->Ncv; i++) {
        fprintf(file, "%d ", fgm->Ngrid[i]);
    }
    fprintf(file, "\n"); // New line after the last Ngrid value
	
	fprintf(file, "[GRIDPOWER]\n");
	
	// Write the fgm->gridpower values
	for (int i = 0; i < fgm->Ncv; i++) {
		fprintf(file, "%f ", fgm->gridpower[i]); // Assuming gridpower is an array of floats
	}
	fprintf(file, "\n"); // New line after the last gridpower value
	
	fprintf(file, "[VARIABLES]\n");
	fprintf(file, "%d\n", fgm->Nvar);
	// Write the variable names, one per line
	for (int i = 0; i < fgm->Nvar; i++) {
		fprintf(file, "%s\n", fgm->varname[i]);
	}
	fprintf(file, "[END]\n");
	fprintf(file, "[DATA]\n");
	
	int totalGridPoints = 1;
	for (int i = 0; i < fgm->Ncv; i++) {
		totalGridPoints *= fgm->Ngrid[i];
	}
	
	// Write the fgm->data values
	for (int i = 0; i < totalGridPoints; i++) {
		for (int j = 0; j < fgm->Nvar; j++) {
			fprintf(file, "%16.10e ", fgm->data[i * fgm->Nvar + j]);
		}
		fprintf(file, "\n"); // New line after each set of Nvar variables
	}
	
	fprintf(file, "[END]\n");
	
    fclose(file);
	
	printf("\n FGM ascii file written\n\n");


    return EXIT_SUCCESS; // Success
}

int writeFGM_MAP
(
    const FGM_MAP *fgm, 
	const char *filename
)
{
	
    printf("\n Writing FGM mapping ascii file...\n\n");

    if (fgm == NULL || filename == NULL) {
        fprintf(stderr, "Error: fgm or filename is empty\n");        
        return EXIT_FAILURE;
    }

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Failed to open FGM database file\n");
        return EXIT_FAILURE;
    }

    // Write FGM data to the file
	fprintf(file, "[FGM]\n");
	fprintf(file, "NDRL\n");
	fprintf(file, "[DIMENSION]\n");
	fprintf(file, "%d\n", fgm->Ncv);
	fprintf(file, "[DATASIZE]\n");

	// Write the fgm->Ngrid values
	for (int i = 0; i < fgm->Ncv; i++) {
        fprintf(file, "%d ", fgm->Ngrid[i]);
    }
    fprintf(file, "\n"); // New line after the last Ngrid value

	fprintf(file, "[GRIDPOWER]\n");
	
	// Write the fgm->gridpower values
	for (int i = 0; i < fgm->Ncv; i++) {
		fprintf(file, "%f ", fgm->gridpower[i]); // Assuming gridpower is an array of floats
	}
	fprintf(file, "\n"); // New line after the last gridpower value

	fprintf(file, "[MINS]\n");
	
	// Write the fgm->mins values
	for (int i = 0; i < fgm->Ncv; i++) {
		fprintf(file, "%f ", fgm->mins[i]); // Assuming mins is an array of floats
	}
	fprintf(file, "\n"); // New line after the last mins value

	fprintf(file, "[MAXS]\n");
	
	// Write the fgm->mins values
	for (int i = 0; i < fgm->Ncv; i++) {
		fprintf(file, "%f ", fgm->maxs[i]); // Assuming maxs is an array of floats
	}
	fprintf(file, "\n"); // New line after the last maxs value

	
	fprintf(file, "[VARIABLES]\n");
	fprintf(file, "%d\n", fgm->Nvar);
	// Write the variable names, one per line
	for (int i = 0; i < fgm->Nvar; i++) {
		if (i < fgm->Ncv) {
			fprintf(file, "iCV%d\n",i+1);
		} else {
			fprintf(file, "i\n");
		}
	}
	fprintf(file, "[END]\n");
	fprintf(file, "[DATA]\n");

	int totalGridPoints = 1;
	for (int i = 0; i < fgm->Ncv; i++) {
		totalGridPoints *= fgm->Ngrid[i];
	}

	// Write the fgm->data values
	for (int i = 0; i < totalGridPoints; i++) {
		for (int j = 0; j < fgm->Nvar; j++) {
			fprintf(file, "%d ", fgm->data[i * fgm->Nvar + j]);
		}
		fprintf(file, "\n"); // New line after each set of Nvar variables
	}
	
	fprintf(file, "[END]\n");
	
    fclose(file);
	
	printf("\n FGM map ascii file written\n\n");


    return EXIT_SUCCESS; // Success
}

// *************************************************************************** //