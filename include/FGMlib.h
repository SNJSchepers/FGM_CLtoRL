/* --------------------------------------------------------------------------- *\

                             Create and read the FGM

\* --------------------------------------------------------------------------- */

#ifndef FGMLIB_H
#define FGMLIB_H

# define MAX_LINE_LENGTH 1024

# define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; })
# define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })
	  
#define VAR_NAME_LENGTH 50

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
#include <stdbool.h>

# include <unistd.h>
# include <sys/stat.h>

// Structure representing the FGM data
typedef struct {
    int Ncv; 			// Number of control variables
    int *Ngrid;			// Array storing the number of grid points for each control variable
    double *gridpower; 	// Array storing grid power values
    int Nvar;			// Number of variables
    char (*varname)[VAR_NAME_LENGTH]; // Array of variable names
    double *data;		// Pointer to the data array
    double *mins; 		// Array to store minimum values for each control variable
    double *maxs; 		// Array to store maximum values for each control variable
} FGM;

// Structure representing the FGM mapping data
typedef struct {
    int Ncv; 			// Number of control variables
    int *Ngrid;			// Array storing the number of grid points for each control variable
    double *gridpower; 	// Array storing grid power values
    int Nvar;			// Number of variables
    char (*varname)[VAR_NAME_LENGTH]; // Array of variable names
    int *data;		// Pointer to the data array
	double *mins; 		// Array to store minimum values for each control variable
    double *maxs; 		// Array to store maximum values for each control variable
} FGM_MAP;

// Read FGM data from a file and return a pointer to the FGM struct.
// @param filename: The name of the file containing FGM data.
// @return: A pointer to the FGM struct containing the data.
FGM * readFGM(const char filename[]);

// Free memory allocated for the FGM struct and its members.
// @param fgm: Pointer to the FGM struct to be freed.
// @return: 0 if successful (or if fgm is NULL), -1 if there was an error.
int freeFGM
(
    FGM *fgm
);

// Free memory allocated for the FGM_MAP struct and its members.
// @param fgm: Pointer to the FGM_MAP struct to be freed.
// @return: 0 if successful (or if fgm is NULL), -1 if there was an error.
int freeFGM_MAP
(
    FGM_MAP *fgm
);

// Search for a keyword in a file.
// @param fid: Pointer to the file to search.
// @param keyword: The keyword to search for.
// @return: 0 if the keyword is found, -1 if not found, or an error occurred.
int locateKeyWord
(
    FILE *fid, 
	char keyword[]
);

// Create a new FGM struct and populate its members.
// @param Ncv: Number of control variables.
// @param Ngrid: Array specifying the number of grid points in each dimension.
// @param gridpower: Array of grid power values for each control variable.
// @param Nvar: Number of variables.
// @param varname: Array of variable names.
// @param data: Pointer to the data array.
// @return: A pointer to the newly created FGM struct.
FGM *createFGM
(
    int Ncv, 
	int *Ngrid, 
	double *gridpower, 
	int Nvar, 
	char (*varname)[VAR_NAME_LENGTH], 
	double *data
);

// Create a new FGM_MAP struct and populate its members.
// @param Ncv: Number of control variables.
// @param Ngrid: Array specifying the number of grid points in each dimension.
// @param gridpower: Array of grid power values for each control variable.
// @param mins: Array of minimum values for each control variable.
// @param maxs: Array of maximum values for each control variable.
// @param varname: Array of variable names.
// @param data: Pointer to the map data array.
// @return: A pointer to the newly created FGM_MAP struct.
FGM_MAP *createFGM_MAP
(
    int Ncv, 
	int *Ngrid, 
	double *gridpower, 
	double *mins, 
	double *maxs, 
	char (*varname)[VAR_NAME_LENGTH], 
	int *data
);

// Write the contents of an FGM struct to a file.
// @param fgm: Pointer to the FGM struct to be written.
// @param filename: The name of the file where FGM data will be written.
// @return: 0 if successful, -1 if there was an error.
int writeFGM
(
    const FGM *fgm, 
	const char *filename
);

// Write the contents of an FGM_MAP struct to a file.
// @param fgm: Pointer to the FGM_MAP struct to be written.
// @param filename: The name of the file where FGM_MAP data will be written.
// @return: 0 if successful, -1 if there was an error.
int writeFGM_MAP
(
    const FGM_MAP *fgm, 
	const char *filename
);

#endif // FGM_H

// *************************************************************************** //