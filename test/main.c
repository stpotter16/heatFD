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
	/* Declare variables */
	char *method;
	int inter_max;
	double tol;
	int order;
	int dimension;
	double kappa;
	int n;
	double Txmin;
	double Txmax;
	double Tymin;
	double Tymax;
	double Ax;
	double Bx;
	int verification;
	int out_mode;
	int timer;
	char *output_file;
	
	grvy_timer_begin("Main Program");
	/* Get log level and timer setting*/
	int log_return;
	log_return = get_settings(&out_mode, &timer);
	if(log_return){
		printf("Something went wrong...search for return code %d\n", log_return);
		exit(1);
	}
	
	/* Set log level */
	if(out_mode==0){
		grvy_log_setlevel(GRVY_NOLOG);
	} else if(out_mode==1){
		grvy_log_setlevel(GRVY_INFO);
	} else if (out_mode==2){
		grvy_log_setlevel(GRVY_DEBUG);
	} else {
		grvy_log_setlevel(GRVY_DEBUG):
	}	
	
	grvy_printf(GRVY_INFO, "Solving the Steady State Heat Equation with Finite Difference \n \n");
	grvy_printf(GRVY_INFO, "Parsing Log/Timing Settings...DONE\n");
	
	/* Call input file parsing */
	int input_return;
	input_return = input(&method, &iter_max, &tol, &order, &dimension, &kappa, &n, &Txmin, &Txmax, &Tymin, &Tymax, &Ax, &Bx, &verification, &output_file);
	if(input_return){
		grvy_printf(GRVY_INFO,"Something went wrong...search for return code %d\n", input_return);
		free(method);
		free(output_file);
		exit(1);
	}
	
	if(argc > 1){
		int cl_return;
		cl_return = cl_parse(argc, argv, &order, &dimension, &verification, &out_mode, &timer, &output_file);
		if(cl_return){
			grvy_printf(GRVY_INFO, "Something went wrong...search for return code %d\n", cl_return);
			free(method);
			free(output_file);
			exit(1);
		}
	}

	/* Assemble system */
	/* Allocate memory for matrix A, vector T, and vector Q */
	double *A;
	double *T;
	double *Q;
	int size;
	/* Only allocate space for interior nodes which are the unknowns */
	if(dimension==1){
		size = (n-2);
	} else if(dimension==2){
		size = (n-2) * (n-2);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %i or %i\n", 1, 2);
		free(method);
		free(output_file);
		exit(1);
	}
	A=malloc(size * size * sizeof(double));
	T=malloc(size * sizeof(double));
	Q=malloc(size * sizeof(double));
	if(!(A && T && Q){
		grvy_printf(GRVY_INFO, "Something went wrong with malloc!\n");
		free(A);
		free(T);
		free(Q);
		exit(1)
	}
	int assemble_return;
	assemble_return = assemble(A, T, Q, dimension, order, n, size, kappa, Txmin, Txmax, Tymin, Tymax, Ax, Bx, verification);
	if(assmble_return){
		grvy_printf(GRVY_INFO,"Something went wrong...search for return code %d\n", assemble_return);
		cleanup(A, T, Q, &method, &output_file);
		exit(1);
	}

	/* Call solver */
	int solver_return;
	solver_return = solver(A, T, Q, size, tol, itermax, method);
	if(solver_return){
		grvy_printf(GRVY_INFO,"Something went wrong...search for return code %d\n", solver_return);
		cleanup(A, T, Q, &method, &output_file);
		exit(1);
	}

	/* Write output */
	int output_return;
	output_return = output(T, n, size, dimension, kappa, Ax, By, output_file, verification);
	if(output_return){
		grvy_printf(GRVY_INFO,"Something went wrong...search for return code %d\n", output_return);
		cleanup(A, T, Q, &method, &output_file);
		exit(1);
	}

	/* Free variables */
	cleanup(A, T, Q, &method, &output_file);
	
	/* Timing */
	grvy_timer_end("Main Program");
	grvy_timer_finalize();
	if(debug>=1 && timer){
		grvy_timer_summarize();
	}
	
	/* Return 0 on clean exit */
	return 0;
}
