#include <stdio.h> /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <ctype.h> /* needed for isspace() */
#include <string.h> /* needed for memset() */
#include <glpk.h> /* the linear programming toolkit */
#include <stdbool.h> /* needed for boolean data type */

/* global variables */
int numRows; /* number of rows in the puzzle */
int length; /* length of each row in the puzzle */
int *lengths; /* array containing the lengths of each break (0 for no brick) */
int *values; /* array containing the values of each break (undef if no brick) */
int *solution; /* array containing the solution (1 for brick taken) */
int debug; /* flag for debug mode; 1 means debug mode, 0 means debug off */

/* prototypes of functions */

int readInput(char *filename);
/* reads puzzle from file */
/* readInput creates and fills the global variables as needed */
/* it returns 0 if all is okay and 1 otherwise */

void computeSolution(void); /* computes a solution and stores it in solution */

void countBricksAndConstraints(int *totalBricks, int *numConstraints); /* Iterate through input arrays counting bricks and constraints */

void setVariableBoundsAndCoefficients(glp_prob *lp, int totalBricks); /* For every brick - Set its variable, variable bounds, and variable coefficient */

void addConstraints(glp_prob *lp, int numConstraints, int totalBricks);
/* Add Constraints - For every brick that has brick(s) above it, you have to pick the brick(s) above before you can pick the current brick */

void solveLinearProgram(glp_prob *lp); /* Solve the linear program, and update the solutions array to indicate which bricks are picked */

int getCurrentBrickPosition(int rowIdx, int columnIdx); /* Compute position of the brick given its row and column index */

bool areAllBricksNegative(void); /* Checks if all bricks in input have negative values */

bool areAllBricksPositive(void); /* Checks if all bricks in input have positive values */

void takeAllBricks(void); /* Select all bricks - Update solution array with 1s to indicate all bricks are picked */

bool doesPositionContainBrick(int currentBrickPosition); /* Checks if a brick is present at a specified position */

int getBrickLeftBound(int rowIdx, int columnIdx); /* Get the left bound (start) of a brick */

int *createBrickPosToVarIdxMap(void); /* Create an array to map brick's position to its variable index */

/*
 * computeSolution - Function to actually solve the problem
				   - Models the problem as a linear programming problem, and solves using GNU Linear Programming Kit (GLPK)
				   - Uses the solutions obtained by LP to update solutions arrays - 1 if brick taken, 0 if not
 * @param void - Function doesn't take any parameter
 * @return void - Function doesn't return anything
*/
void computeSolution(void) {
	// Edge Case 1: if all bricks have negative values (and/or 0s), it is optimal to not take any bricks
	bool allBricksNegative = areAllBricksNegative();
	if (allBricksNegative) {
		// No bricks should be picked. Solution array by default is initialised with all zeros. Don't modify solutions array and exit function 
		return;
	}

	// Edge Case 2: if all bricks have positive values, it is optimal to take all of them
	bool allBricksPositive = areAllBricksPositive();
	if (allBricksPositive) {
		// All bricks should be picked. Update solution arrays with 1s to indicate all bricks are picked and exit function 
		takeAllBricks();
		return;
	}

	// Linear programming problem
	glp_prob *lp;

	// Create LP problem instance
	lp = glp_create_prob();
	// Set objective for LP - Maximise the weights in our closure
	glp_set_obj_dir(lp, GLP_MAX);
		
	int totalBricks = 0;	// Total number of bricks - gives the LP the total number of variables (columns)
	int numConstraints = 0; // Total number of constraints our LP has to satisfy (rows)

	// Iterate through input arrays - count bricks and constraints
	countBricksAndConstraints(&totalBricks, &numConstraints);

	// For every brick - Set its variable, variable bounds, and variable coefficient for the linear program
	setVariableBoundsAndCoefficients(lp, totalBricks);

	// For every brick that has brick(s) above it - add a constraint such that you have to pick the brick(s) above before you can pick the current brick
	addConstraints(lp, numConstraints, totalBricks);

	// Solve the LP - Update solutions array to indicate which bricks are picked
	solveLinearProgram(lp);

	// Release memory used for LP
	glp_delete_prob(lp);
}

/* countBricksAndConstraints - Iterate through input arrays counting bricks and constraints
 * @param *totalBricks - Pointer to total bricks variable which will be incremented for every new brick
 * @param *numConstraints - Pointer to total constraints variable which will be incremented for every new constraint
 * @return void - Function doesn't return anything
*/
void countBricksAndConstraints(int *totalBricks, int *numConstraints) {
	// For each row
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		// For each column position in that row
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {
			
			// Get potential brick position
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);

			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition); 

			if (positionContainsBrick) {
				// Brick Found
				(*totalBricks)++;
				if (debug) { fprintf(stdout, "DEBUG: Brick Found at Position = %d\n", currentBrickPosition); }
				
				// Bricks from the 2nd row onwards will have bricks above them
				if (rowIdx > 0) {

					// Check for constraints: Check which brick(s) are above the current one

					// Run over bricks before current brick and add lengths to compute beginning (left) and end (right)
					int left = getBrickLeftBound(rowIdx, columnIdx);
					int right = left + lengths[currentBrickPosition];

					// Check bricks in row above:
					int brickAbovePos = 0;

					// For each brick positions in the row above
					for (int aboveColumnIdx = 0; aboveColumnIdx < length; aboveColumnIdx++) {
						// Position of the (potential) brick above the current brick
						int aboveBrickPosition = getCurrentBrickPosition(rowIdx-1, aboveColumnIdx);
						int aboveBrickLength = lengths[aboveBrickPosition];

						// Check if a brick is found above current brick
						if (aboveBrickLength > 0) {
							
							// Check if bricks Overlap - 2 conditions
							// (1) To overlap, the position of where brick above ends cannot be before the left
							// (2) To overlap, the position of where brick above starts cannot be beyond the right
							if (
								(brickAbovePos + aboveBrickLength > left) && (brickAbovePos < right)
							) {
								// Current brick has a brick above it.
								// You cannot pick current brick unless you pick the brick above
								// This gives us a constraint. Increment the number of constraints
								(*numConstraints)++;
							}
							brickAbovePos += aboveBrickLength;
						}
					}
				}
			}
		}
	}

	if (debug) {
		fprintf(stdout, "\nDEBUG: Total number of bricks = %d\n", *totalBricks);
		fprintf(stdout, "DEBUG: Total number of constraints = %d\n\n", *numConstraints);
	}

}

/*
 * setVariableBoundsAndCoefficients - Set the variables, variable bounds, and variable coefficients for the linear program
 * @param *lp - Pointer to the linear program
 * @param totalBricks - Total number of bricks
 * @return void - Function doesn't return anything
*/
void setVariableBoundsAndCoefficients(glp_prob *lp, int totalBricks) {
	// Set number of variables (columns)
	// Create a variable for each brick
	// Variable is 1 if brick is taken, and 0 otherwise 
	glp_add_cols(lp, totalBricks);
	
	if (debug) {
		fprintf(stdout, "DEBUG: Setting total number of variables (columns) for the LP. A variable per each brick, so a total of %d variables \n", totalBricks);
		fprintf(stdout, "DEBUG: LP maximises the sum of the following terms: \n");
	}

	// Set bounds and objective coefficients for the brick variables
	int brickVarIdx = 1;
	// For each row
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		// For each column position in that row
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);

			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);

			if (positionContainsBrick) {
				// Brick Found

				// Set Bounds: Range of values for the variable brickVarIdx
				// Bounds: 0 <= variable <= 1
				glp_set_col_bnds(lp, brickVarIdx, GLP_DB, 0.0, 1.0);

				// Set objective coefficient: LP multiplies the coefficient to the brickVarIdx variable
				// Coefficient = value (cost) of the brick
				// LP maximises sum of (variables * their objective coefficients)
				glp_set_obj_coef(lp, brickVarIdx, values[currentBrickPosition]);

				if (debug) {
					fprintf(stdout, "%d * x%d \n", values[currentBrickPosition], brickVarIdx); 
				}

				brickVarIdx++;
			}
		}
	}

	if (debug) {
		fprintf(stdout, "DEBUG: where for every variable: 0 <= variable <=1 \n\n");
	}

}

/*
 * addConstraints - For every brick that has brick(s) above it: add a constraint such that you have to pick the brick(s) above before you can pick the current brick
 * @param *lp - Pointer to the linear program
 * @param numConstraints - Total number of constraints
 * @param totalBricks - Total number of bricks
 * @return void - Function doesn't return anything
*/
void addConstraints(glp_prob *lp, int numConstraints, int totalBricks) {
	
	// Check if at least 1 constraint exists
	bool atLeastOneConstraint = numConstraints > 0;
	
	// If there are constraint(s), add them:
	if (atLeastOneConstraint) {
		// LP needs a row per constraints
		glp_add_rows(lp, numConstraints);

		// Create Arrays for constraint matrix (GLPK uses 1-based indexing)
		// These are used to define the LHS for the constraints

		// Every constraint contains 2 variables (multiplied with their coefficient) - one for the brick above and another for the brick below
		// 1-indexed so start from 1
		int constraintMatrixArraySize = 1 + 2 * numConstraints;
		
		int *ia = malloc(sizeof(int) * (constraintMatrixArraySize));      // Row index - Constraint Number. Constraint enforces the brick above is picked first
		int *ja = malloc(sizeof(int) * (constraintMatrixArraySize));      // Column index - Variable Number. LP needs variable per brick to indicate if its picked
		double *ar = malloc(sizeof(double) * (constraintMatrixArraySize));// Coefficient value - Coefficient value can be +1 or -1. It's multiplied with variable

		// Array to map the position of brick to the variable number of the brick
		int *brickPosToVarIdx = createBrickPosToVarIdxMap();

		// Adding constraints:
		int currConstraint = 1;           // Keep track of constraint number (1-indexed so start from 1)
		int constraintMatrixIdx = 1;      // Matrix entry number - For every variable (and it's coefficient +1/-1) we need a new entry in our matrix

		// For each row, starting from the 2nd
		for (int rowIdx = 1; rowIdx < numRows; rowIdx++) {
			// For each column position in that row
			for (int columnIdx = 0; columnIdx < length; columnIdx++) {
				
				// Get potential brick position
				int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);

				// Check if brick present in current position
				bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);

				if (positionContainsBrick) {
					// Brick Found

					// Add constraints: Check which brick(s) are above the current one:

					// Run over bricks before current brick and add lengths to compute beginning (left) and end (right)
					int left = getBrickLeftBound(rowIdx, columnIdx);
					int right = left + lengths[currentBrickPosition];

					// Check each brick in row above
					int brickAbovePos = 0;

					// For each brick positions in the row above
					for (int aboveColumnIdx = 0; aboveColumnIdx < length; aboveColumnIdx++) {

						// Position of the (potential) brick above the current brick
						int aboveBrickPosition = getCurrentBrickPosition(rowIdx-1, aboveColumnIdx);
						int aboveBrickLength = lengths[aboveBrickPosition];

						// Check if a brick is found above the current brick
						if (aboveBrickLength > 0) {

							// Check if bricks Overlap - 2 conditions
							// (1) To overlap, the position of where brick above ends cannot be before the left
							// (2) To overlap, the position of where brick above starts cannot be beyond the right
							if (
								(brickAbovePos + aboveBrickLength > left) && (brickAbovePos < right)
							) {
								// Add constraint: Current brick has a brick above it

								// If you want to pick the current brick, you first have to ensure you picked the brick above it
								// So if the above brick is picked (1)
								// The brick below may be picked (1) or it may not be picked (0)
								
								// This means, currentBrick <= aboveBrick
								// And hence, currentBrick - aboveBrick <= 0

								// Set the RHS to <= 0
								glp_set_row_bnds(lp, currConstraint, GLP_UP, 0.0, 0.0); 

								// Set the LHS of the constraints using the arrays for the constraint matrix
								ia[constraintMatrixIdx] = currConstraint;                         // current constraint (equation) number
								ja[constraintMatrixIdx] = brickPosToVarIdx[currentBrickPosition]; // brick variable number
								ar[constraintMatrixIdx] = 1.0;                                    // current brick coefficient: +1
								constraintMatrixIdx++;                                            

								ia[constraintMatrixIdx] = currConstraint;                        // current constraint (equation) number
								ja[constraintMatrixIdx] = brickPosToVarIdx[aboveBrickPosition];  // brick variable number
								ar[constraintMatrixIdx] = -1.0;                                  // above brick coefficient: -1
								constraintMatrixIdx++;

								if (debug) {
									fprintf(
										stdout, "DEBUG: Constraint %d: x%d - x%d <= 0\n",
										currConstraint,
										brickPosToVarIdx[currentBrickPosition],
										brickPosToVarIdx[aboveBrickPosition]
									);
								}

								currConstraint++;
							}

							// Move onto the next brick in the row above
							brickAbovePos += aboveBrickLength;
						}
					}
				}
			}
		}

		// Every constraint contains 2 variables (multiplied with their coefficient) - one for the brick above and another for the brick below
		// Total number of elements in the matrix = Total number of constraints * Number of variables per constraint = numConstraints * 2
		int constraintMatrixElementCount = 2 * numConstraints;

		// Load the constraint matrix containing all the constraints (LHS)
		glp_load_matrix(lp, constraintMatrixElementCount, ia, ja, ar);

		// Now that constraint matrix is loaded, free the memory that was allocated by malloc for constraint matrix arrays
		free(ia);
		free(ja);
		free(ar);
		free(brickPosToVarIdx);
	}
}


/*
 * solveLinearProgram - Solve the linear program using GLPK's simplex method and update the solutions array to indicate which bricks are picked
 * @param *lp - Pointer to the linear program
 * @return void - Function doesn't return anything
*/
void solveLinearProgram(glp_prob *lp) {
	if (debug) {
		glp_term_out(1); // Enable debug output from GLPK
	} else {
		glp_term_out(0); // Disable debug output from GLPK
	}

	// solve the LP using the simplex method
	glp_simplex(lp, NULL);
	
	// Current brick variable index
	int brickVarIdx = 1;

	// For each row
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		// For each column position in that row
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {
			
			// Get current brick position
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);

			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);

			if (positionContainsBrick) {

				// Obtain the answer for the current brick variable
				double brickVarSolution = glp_get_col_prim(lp, brickVarIdx);

				if (debug) {
					fprintf(stdout, "DEBUG: Solution for brick %d = %f\n", brickVarIdx, brickVarSolution);
				}

				// The answer for the brick variable, brickVarSolution has an allowed range of [0,1]
				// For Maximum-Weight Closure Problem problems, the brickVarSolution will really only have 2 possible values: 1 or 0
				// 1 indicates the brick is selected, and 0 indicates the brick is not
				// To be extra cautious, the brickVarSolution can be rounded to the closest integer - 1 if it's greater than 0.5 and 0 otherwise
				if (brickVarSolution > 0.5) {
					solution[currentBrickPosition] = 1;
				} else {
					solution[currentBrickPosition] = 0;
				}

				brickVarIdx++;
			}
		}
	}

	if (debug) {
		fprintf(stdout, "DEBUG: Solution array updated with variable solutions - 1 indicates the brick is taken and 0 if not\n\n");
	}

}

/*
 * getCurrentBrickPosition - Compute position of the brick given its row and column index
                           - Using the brick's position, you can get the length and value of the (columnIdx+1)th brick in the (rowIdx+1)th row
                           - The brick's position is used to get its length from the lengths array, and value from the values array
 * @param rowIdx - Row index
 * @param columnIdx - Column index
 * @return - Integer representing the current brick's position
*/
int getCurrentBrickPosition(int rowIdx, int columnIdx) {
	int currentBrickPosition = rowIdx*length + columnIdx;
	return currentBrickPosition;
}

/*
* areAllBricksNegative - Iterates through all bricks and checks if all of them have a negative value
* @param void - Function doesn't take any parameter
* @return - Returns true if all bricks are negative (or 0) and false otherwise
*/
bool areAllBricksNegative(void) {
	// For each row
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		// For each column position in that row
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {
			
			// Position of a potential brick
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);
			
			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);
			
			if (positionContainsBrick) {

				// Get value of current brick
				int currentBrickValue = values[currentBrickPosition];

				// If the current brick has a positive value, it is impossible for all bricks to be negative
				if (currentBrickValue > 0) {
					return false;
				}
			}
		}
	}

	if (debug) {
		fprintf(stdout, "DEBUG: All bricks are negative. Optimal to not take any bricks\n\n");
	}

	// No positive bricks were found, this means all the bricks were negative (or 0)
	return true;
}

/*
* areAllBricksPositive - Iterates through all bricks and checks if all of them have a positive value
* @param void - Function doesn't take any parameter
* @return - Returns true if all bricks are positive and false otherwise
*/
bool areAllBricksPositive(void) {
	// For each row
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		// For each column position in that row
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {
			
			// Position of a potential brick
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);
			
			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);
			
			if (positionContainsBrick) {

				// Get value of current brick
				int currentBrickValue = values[currentBrickPosition];

				// If the current brick has a negative value, it is impossible for all bricks to be positive
				if (currentBrickValue < 0) {
					return false;
				}
			}
		}
	}

	if (debug) {
		fprintf(stdout, "DEBUG: All bricks are positive. Optimal to take all bricks\n");
	}

	// No negative bricks were found, this means all the bricks were positive
	return true;
}


/*
* takeAllBricks - Select all bricks: Update solution array with 1s to indicate all bricks are picked
* @param void - Function doesn't take any parameter
* @return void - Function doesn't return anything 
*/
void takeAllBricks(void) {
	if (debug) {
		fprintf(stdout, "DEBUG: Selecting all bricks ...\n");
	}

	// For each row
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		// For each column position in that row
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {
			
			// Get potential brick position
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);

			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);

			if (positionContainsBrick) {

				// Select the brick. Update its position in the solution array with 1 to indicate that the brick is selected
				solution[currentBrickPosition] = 1;

				if (debug) {
					fprintf(stdout, "DEBUG: Selecting the brick at position = %d\n", currentBrickPosition);
				}

			}
		}
	}

	if (debug) {
		fprintf(stdout, "DEBUG: All bricks selected - updated solution array with all 1s\n\n");
	}
}

/*
* doesPositionContainBrick - Checks if a brick is present at a specified position
* @param currentBrickPosition - Integer representing position of brick
* @return - Returns true if a brick is present at current position and false otherwise
*/
bool doesPositionContainBrick(int currentBrickPosition){

	// Brick exists at that position if the length array has a non-zero entry at that index
	bool positionContainsBrick = lengths[currentBrickPosition] > 0;

	return positionContainsBrick;
}

/*
* getBrickLeftBound - Get the left bound (start) of a brick
* @param rowIdx - Row index of current brick
* @param columnIdx - Column index of current brick
* @return - Integer specifying the starting position of the current brick
*/
int getBrickLeftBound(int rowIdx, int columnIdx) {
	int left = 0;
	
	for (int posIdx = 0; posIdx < columnIdx; posIdx++) {
		int brickPosition = getCurrentBrickPosition(rowIdx, posIdx);
		left += lengths[brickPosition];
	}

	return left;
}

/*
* createBrickPosToVarIdxMap - Create an array to map brick's position to its variable index
* @param void - Function doesn't take any parameter
* @return - Pointer to array of Integers mapping brick's position to its variable index
*/
int *createBrickPosToVarIdxMap(void) {

	// Variable to represent the brick. 1st brick is x1, 2nd brick is x2, 3rd brick is x3, etc.
	int brickVarIdx = 1;
	// Array to map the position of brick to the variable number of the brick
	int *brickPosToVarIdx = malloc(sizeof(int) * numRows * length);

	// Mapping brick positions in grid to variable numbers
	// Having this mapping is efficient to get the brick's variable number instead of computing it on the fly for every brick
	for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
		for (int columnIdx = 0; columnIdx < length; columnIdx++) {

			// Obtain position of current brick
			int currentBrickPosition = getCurrentBrickPosition(rowIdx, columnIdx);

			// Check if brick present in current position
			bool positionContainsBrick = doesPositionContainBrick(currentBrickPosition);

			if (positionContainsBrick) {
				// Brick Found, map currentBrickPosition to the brick's variable number
				brickPosToVarIdx[currentBrickPosition] = brickVarIdx;

				if (debug) { fprintf(stdout, "DEBUG: Brick at position %d is represented as x%d \n", currentBrickPosition, brickVarIdx);}

				brickVarIdx++;
			} else {
				// No brick at current position, map currentBrickPosition to 0
				brickPosToVarIdx[currentBrickPosition] = 0;
			}		
		}
	}

	return brickPosToVarIdx;
}



/* int canBeTaken(int row, int brick) checks if the brick at position brick */
/* in row can be taken by checking the bricks immediately above (only those) */
int canBeTaken(int row, int brick) {
	int i; /* loop variables to run over row */
	int left, right; /* beginning and end of brick */
	int pos; /* position that is currently checked; */

	if ( row==0 ) {
		return 1; /* every brick in row 0 can be taken */
	}

	/* run over bricks before brick and add lenghts to compute beginning and end */
	for ( i=0, left=0; i<brick; i++ ) {
		left += lengths[row*length+i];
	}
	right = left + lengths[row*length+brick];

	/* run over bricks above and check that all relevant bricks are taken */
	for ( i=0, pos=0; pos<length; i++ ) {
		if ( (pos>=left) && (pos<right) ) {
			/* brick i is above brick; check if taken */
			if ( solution[(row-1)*length+i]==0 ) {
				/* brick not taken => brick cannot be taken */
				return 0;
			}
		}
		/* move to next brick */
		pos += lengths[(row-1)*length+i];
	}
	return 1; /* passed all tests => can be taken */
}

void printSolution(void) {
	int i, j; /* loop variables to go over rows and bricks in each row */
	int value = 0; /* stores the value of the solution */
	int taken; /* flags if at least one brick was taken */

	/* go over each row to print which bricks are taken */
	for ( i=0; i<numRows; i++ ) {
		taken = 0; /* no brick taken yet */
		fprintf(stdout, "Bricks taken from row %d: ", i+1);
		/* go over all bricks to print the ones that are taken */
		for ( j=0; j<length; j++ ) {
			if ( solution[i*length+j]==1 ) {
				taken=1; /* remember that a brick was taken */
				if ( canBeTaken(i, j) && ( lengths[i*length+j]>0 ) ) {
					/* brick is taken and this is okay */
					value += values[i*length+j]; /* add up value of taken bricks */
					fprintf(stdout, "%d (value %d) ", j+1, values[i*length+j]);
				} else {
					/* brick could not be taken! */
					fprintf(stdout, "Brick %d (value %d) was taken illegaly!\n",
						j+1, values[i*length+j]);
					return; /* stop printing solution */
				}
			} else {
				if ( solution[i*length+j]!=0 ) { /* only 0 and 1 are allowed values */
					fprintf(stdout, "Brick %d (value %d) has an illegal solution value (%d).\n",
						j+1, values[i*length+j], solution[i*length+j]);
					return; /* stop printing solution */
				}
			}
		}
		if ( taken ) {
			fprintf(stdout, "\n");
		} else {
			fprintf(stdout, "(none)\n");
		}
	}
	fprintf(stdout, "Value of the solution: %d\n", value);
}

int main(int argc, char **argv) {
	int i; /* used to run over the command line parameters */

	if ( argc<2 ) { /* no command line parameter given */
		fprintf(stderr, "Usage: %s [file1] [file2] [file3] [...]\n"
      "Where each [file] is the name of a file with a puzzle.\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	if ( argv[1][0]=='-' && argv[1][1]=='d' && argv[1][2]==0 ) {
    /* If the first parameter is -d we activate debug mode. */
		debug=1; /* switch debug mode on */
		fprintf(stdout, "DEBUG: Debug mode activated\n"); /* be explicit about it */
	} else {
		debug=0; /* switch debug mode off */
	}

  for ( i=1+debug; i<argc; i++ ) { /* go over remaining command line parameters */
    if ( readInput(argv[i]) ) { /* try to read file */
      /* returned with error message */
      fprintf(stderr, "%s: Cannot read puzzle with filename %s. Skipping it.\n",
        argv[0], argv[i]);
    } else { /* input read successfully */
			fprintf(stdout, "Looking at puzzle from %s\n", argv[i]);
			computeSolution(); /* compute a solution */
			printSolution(); /* print the solution including its value */
      /* free memory for next input; avoid a memory leak */
      free(lengths);
			free(values);
      free(solution);
    }
  }
	return EXIT_SUCCESS;
}


/* readInput(*char filename) reads the input and stores it */
/* return value 1 indicates an error; otherwise 0 is returned */
int readInput(char *filename) {
  FILE *fh;
	int i, j; /* loop variables to go over rows and bricks */
	int checkLen; /* add length in each row to check */
	int countSpace; /* two count space */
	char line[4096]; /* one line of text */
	char *inLine; /* pointer to current position in line */

	/* open file for reading */
	fh = fopen("testFile123.txt", "wt");
	if ( fh==NULL ) {
		fprintf(stderr, "Cannot create file testFile123.txt for writing.\n");
	} else {
		fprintf(fh, "Test, 1, 2, 3.\n");
		fclose(fh);
	}
  if ( ( fh = fopen(filename, "rt") ) == NULL ) {
    return 1;
  }
	/* read line containing number of rows and length */
	if ( fgets(line, 4096, fh)==NULL ) {
		if ( debug ) {
			fprintf(stdout, "DEBUG: Unable to read first line of input file.\n");
		}
		fclose(fh);
		return 1;
	}
	/* extract number of rows and length */
	if ( sscanf(line, "%d %d", &numRows, &length)!=2 ) {
		if ( debug ) {
			fprintf(stdout, "DEBUG: Unable to read number of rows and length.\n");
		}
		fclose(fh);
		return 1;
	}
	/* create arrays to store input and solution */
	if ( (lengths = (int *)malloc(sizeof(int)*numRows*length))==NULL ) {
		if ( debug ) {
			fprintf(stdout, "DEBUG: Unable to allocate memory for lengths.\n");
		}
		fclose(fh);
		return 1;
	}
	if ( (values = (int *)malloc(sizeof(int)*numRows*length))==NULL ) {
		if ( debug ) {
			fprintf(stdout, "DEBUG: Unable to allocate memory for values.\n");
		}
		free(lengths);
		fclose(fh);
		return 1;
	}
	if ( (solution = (int *)malloc(sizeof(int)*numRows*length))==NULL ) {
		if ( debug ) {
			fprintf(stdout, "DEBUG: Unable to allocate memory for solution.\n");
		}
		free(lengths);
		free(values);
		fclose(fh);
		return 1;
	}
	if ( debug ) {
		fprintf(stdout, "DEBUG: rows=%d, length=%d\n", numRows, length);
	}
	/* initialise solution with all 0s */
	memset(solution, 0, sizeof(int)*numRows*length);
	/* read rows one after the other */
	for ( i=0; i<numRows; i++ ) {
		/* read next line */
		if ( fgets(line, 4096, fh)==NULL ) {
			if ( debug ) {
				fprintf(stdout, "DEBUG: Unable to read row %d of input file.\n", i+1);
			}
			free(lengths);
			free(values);
			free(solution);
			fclose(fh);
			return 1;
		}
		inLine=line;
		checkLen=0; /* read bricks until length is reached */
		j=0; /* start with brick 8 */
		while ( checkLen<length ) {
			if ( sscanf(inLine, "%d %d", &lengths[i*length+j], &values[i*length+j])!=2 ) {
				if ( debug ) {
					fprintf(stdout, "DEBUG: Unable to read next two values.\n");
				}
				free(lengths);
				free(values);
				free(solution);
				fclose(fh);
				return 1;
			}
			if ( debug ) {
				fprintf(stdout, "DEBUG: brick %d in row %d has length=%d, value=%d\n",
					i, j, lengths[i*length+j], values[i*length+j]);
			}
			if ( lengths[i*length+j]<=0 ) {
				if ( debug ) {
					fprintf(stdout, "DEBUG: length needs to be positive\n");
				}
				free(lengths);
				free(values);
				free(solution);
				fclose(fh);
				return 1;
			}
			checkLen += lengths[i*length+j];
			/* remove two integers from line and continue */
			for ( countSpace=0; countSpace<2; countSpace++ ) {
				/* skip all whitespace */
				while ( isspace(*inLine) ) {
					inLine++;
				}
				/* skip one integer */
				while ( !isspace(*inLine ) ) {
					inLine++;
				}
			}
			/* move to next brick */
			j++;
		}
		if ( checkLen>length ) {
			if ( debug ) {
				fprintf(stdout, "DEBUG: total length %d too long (%d)\n", checkLen, length);
			}
			free(lengths);
			free(values);
			free(solution);
			fclose(fh);
			return 1;
		}
	}
  fclose(fh); /* close file after reading the input */
  return 0; /* signal all went well */
}
