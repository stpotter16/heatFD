/*--------------------------------------------------------------------------------------
 * Program for solving the 1D and 2D heat equation with 2nd and 4th order finite difference iterative solvers.
 *
 * Sam Potter
 *
 *--------------------------------------------------------------------------------------
 *  main.c: Main driver program.
 *-------------------------------------------------------------------------------------*/

#include "utilities.h"

int main(int argc, char *argv[])
{
	/* Declare problem struct */
	Heat Problem;
		
	grvy_timer_begin("Main Program");

	/* Get input file */
	/* Check if there's an actual input */
	if (argc < 2){
		grvy_printf(GRVY_ERROR, "\nUsage error: heatFD [input-file]\n");
		exit(1);
	}

	/* Read input file (log level, timer settings, problem setup, output file */
	int input_return;
	input_return = input(&Problem, argv[1]);
	if(input_return){
		grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", input_return);
		/* Code to clean up the struct */
		free((Problem).method);
		free((Problem).output_file);
		exit(1);
	}

	/* Sanitize input. If returns non-zero, clean up the struct */
	int sanitize_return;
	sanitize_return = sanitize(&Problem);
	if(sanitize_return){
		grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", sanitize_return);
		/* Code to clean up struct */
		free((Problem).method);
		free((Problem).output_file);
		exit(1);
	}

	/* Call PETSC if defined */
	#ifdef INCLUDE_PETSC
		int petsc_return;
		petsc_return = petsc(&Problem);
		if(petsc_return){
			grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", petsc_return);
			cleanup(&Problem);
			exit(1);
		}
	#else

		/* Initialize linear system matrices and assembly them */
		int assemble_return;
		assemble_return = assemble(&Problem);
		if(assemble_return){
			grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", assemble_return);
			/* Code to clean up struct */
			cleanup(&Problem);
			exit(1);
		}

		/* Solve the linear system */
		int solver_return;
		solver_return = solver(&Problem);
		if(solver_return){
			grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", solver_return);
			/* Code to clean up struct */
			cleanup(&Problem);
			exit(1);
		}
	#endif

	/* Write results to file (verify if needed, write settings too) */
	int output_return;
	output_return = hdf5_output(&Problem);
	if(output_return){
		grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", output_return);
		/* Code to clean up strcut */
		cleanup(&Problem);
		exit(1);
	}

	/* Clean up variables */
	cleanup(&Problem);
	
	/* Timer output */
	grvy_timer_end("Main Program");
	grvy_timer_finalize();
	grvy_timer_summarize();
	
	/* Return 0 on clean exit */
	return 0;
}
