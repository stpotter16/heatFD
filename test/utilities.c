/*--------------------------------------------------------------------------------------
 * Program for solving the 1D and 2D heat equation with 2nd and 4th order finite difference iterative solvers.
 *
 * Sam Potter
 *
 *--------------------------------------------------------------------------------------
 *  utilities.c: Contains functions used in main.c 
 *-------------------------------------------------------------------------------------*/

#include "utilities.h"

/*--------------------------------------------------------------------------------------
 *  get_settings: function for getting log level.
 *-------------------------------------------------------------------------------------*/

int get_settings(int *out_mode, int *timer)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Declare variables
	*File open error */
	int erropen;

	/* Open file */
	erropen = grvy_input_fopen("./input.dat");
	if(erropen){
		return 200;
	}
	if(!(grvy_input_fread_int("output/out_mode",&out_mode))){
		return 201;
	}
	if(!(grvy_input_fread_int("output/timer", &timer))){
		return 202;
	}
	/* Timing */
	grvy_timer_end(__func__);
	
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  input: Routine for parsing input file and returning variables for rest of the code.
 *-------------------------------------------------------------------------------------*/

int input(char *solver, int *iter_max, double *tol, int *order, int *dimension, double *kappa, int *n, double *Txmin, double *Tymin, double *Txmax, double *Tymax, double *Ax, double *By, int *verification, char *output_file)
{
	
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Parsing problem details...\n");
	
	/* Declare variable
	*File open error */
	int erropen;

	/* Open file */
	erropen = grvy_input_fopen("./input.dat");
	if(erropen){
		grvy_printf(GRVY_DEBUG,"Unable to open file %-10s\n","./input.dat");
		return 100;
	}

	/* Get solver options */
	if(!(grvy_input_fread_char("solver/method",&solver))){
		grvy_printf(GRVY_DEBUG,"Unable to read char value %-10s\n","solver/method"); 
		grvy_input_fclose();
		return 101;
	 }

	if(!(grvy_input_fread_int("solver/iter_max",&iter_max))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n", "solver/iter_max");
		grvy_input_fclose();
		return 102;
	}

	if(!(grvy_input_fread_double("solver/tol",&tol))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n","solver/tol");
		grvy_input_fclose();	
		return 103;
	}

	if(!(grvy_input_fread_int("solver/order", &order))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n", "solver/order");
		grvy_input_fclose();	
		return 104;
	}

	/* Get problem details */
	if(!(grvy_input_fread_int("problem/dimension",&dimension))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n","problem/dimension"); 
		grvy_input_fclose();
		return 105;
	}

	if(!(dimension == 1 || dimension == 2)){	
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		grvy_input_fclose();
		return 112;
	}
	
	if(!(order == 2 || order == 4)){
		grvy_printf(GRVY_INFO, "Incorrect finite difference order! Must be either %d or %d\n", 2, 4);
		grvy_input_fclose();
		return 118;
	}
	if(!(grvy_input_fread_double("problem/kappa", &kappa))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/kappa"); 
		grvy_input_fclose();
		return 106;
	}

	if(!(grvy_input_fread_int("problem/n",&n))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n","problem/n");
		grvy_input_fclose();
		return 107;
	}
	/*
	if(!(grvy_input_fread_float("problem/xmin", &xmin))){
		printf("Unable to read float value %-10s", "problem/xmin");
		return 8;
	}

	if(!(grvy_input_fread_float("problem/xmax", &xmax))){
		printf("Unable to read float value %-10s", "problem/xmax");
		return 9;
	}
	*/
	if(!(grvy_input_fread_double("problem/Txmin", &Txmin))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Txmin");
		grvy_input_fclose();
		return 108;
	}

	if(!(grvy_input_fread_double("problem/Txmax", &Txmax))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Txmax");
		grvy_input_fclose();
		return 109;
	}
	/*
	if(!(grvy_input_fread_int("problem/ny",&ny))){
		printf("Unable to read int value %-10s","problem/ny");
		return 12;
	}

	if(!(grvy_input_fread_float("problem/ymin", &ymin))){
		printf("Unable to read float value %-10s", "problem/ymin");
		return 13;
	}

	if(!(grvy_input_fread_float("problem/ymax", &ymax))){
		printf("Unable to read float value %-10S", "problem/ymax");
		return 14;
	}
	*/
	if(!(grvy_input_fread_double("problem/Tymin", &Tymin))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Tymin");
		grvy_input_fclose();
		return 110;
	}

	if(!(grvy_input_fread_double("problem/Tymax", &Tymax))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Tymax");
		grvy_input_fclose();
		return 111;
	}

	/* Get MASA options */
	if(!(grvy_input_fread_int("MASA/Ax", &Ax))){
		grvy_printf(GRVY_DEBUG, "Unable to read double value %-10s\n", "Ax");
		grvy_input_fclose();
		return 120;
	}

	if(!(grvy_input_fread_int("MASA/By", &By))){
		grvy_printf(GRVY_DEBUG, "Unable to read double value %-10s\n", "By");
		grvy_input_fclose();
		return 121;
	}
	/* Get output options */
	if(!(grvy_input_fread_int("output/verification",&verify))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n", "verification");
		grvy_input_fclose();
		return 113;
	}
	
	if(!(grvy_input_fread_char("output/out_file", &output_file))){
		grvy_printf(GRVY_DEBUG, "Unable to read char value %-10s\n", "output_file");
		grvy_input_fclose();
		return 114;
	}

	/* Close file */
	grvy_input_fclose();
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Parsing problem details...DONE\n");
	
	/* Writing settings to stdout */
	grvy_printf(GRVY_INFO, "\n Input file: %-10s\n", "input.dat");
	
	grvy_printf(GRVY_INFO, "\n Mesh Settings:\n");
	grvy_printf(GRVY_INFO,"--Dimension:%d\n", dimension);
	if(dimension==1){
		grvy_printf(GRVY_INFO,"--nx:%d\n", n);
		grvy_printf(GRVY_INFO,"--(xmin, xmax): (%f, %f)\n", 0.0, 1.0);
	} else if(dimension==2){
		grvy_printf(GRVY_INFO,"--nx:%d\n", n);
		grvy_printf(GRVY_INFO,"--(xmin, xmax): (%f, %f)\n", 0.0, 1.0);
		grvy_printf(GRVY_INFO,"--ny:%d\n", n);
		grvy_printf(GRVY_INFO,"--(ymin, ymax): (%f, %f)\n", 0.0, 1.0);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 115;
	}
	
	grvy_printf(GRVY_INFO, "\n Boundary Conditions:\n");
	if(dimension==1){
		grvy_printf(GRVY_INFO,"--(Txmin, Txmax): (%f, %f)\n", Txmin, Txmax);
	} else if(dimension==2){
		grvy_printf(GRVY_INFO,"--(Txmin, Txmax): (%f, %f)\n", Txmin, Txmax);
		grvy_printf(GRVY_INFO,"--(Tymin, Tymax): (%f, %f)\n", Tymin, Tymax);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 115;
	}

	grvy_printf(GRVY_INFO, "\n MASA Settings:\n");
	if(dimension==1){
		grvy_printf(GRVY_INFO,"--Ax: %f\n", Ax);
	} else if(dimension==2){
		grvy_printf(GRVY_INFO,"--Ax: %f\n", Ax);
		grvy_printf(GRVY_INFO,"--By: %f\n", By);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 115;
	}

	grvy_printf(GRVY_INFO, "\n Solver Settings:\n");
	grvy_printf(GRVY_INFO,"--Method:%-10s\n", method);
	grvy_printf(GRVY_INFO,"--Max iterations:%d\n", iter_max);
	grvy_printf(GRVY_INFO,"--Convergence Tolerance:%f\n", tol);
	grvy_printf(GRVY_INFO,"--Finite Difference Order:%d\n", order);
	
	grvy_printf(GRVY_INFO, "\n Output Settings:\n");
	grvy_printf(GRVY_INFO,"--Verification:%d (0=no, 1=yes)\n", verification);
	grvy_printf(GRVY_INFO,"--Output Mode:%d (0=silent, 1=standard, 2=debug)\n", out_mode);
	grvy_printf(GRVY_INFO,"--Timer:%d (0=don't write, 1=write)\n", timer);
	grvy_printf(GRVY_INFO,"--Output File:%-10s\n \n", output_file);

	/* Timing */
	grvy_timer_end(__func__);

	/* Return 0 if function exited successfully */
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  cl_parse: Parse command line input
 *-------------------------------------------------------------------------------------*/

int cl_parse(int argc, char *argv[], int *order, int *dimension, int *verification, int *out_mode, int *timer, char *output_file)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Parsing command line options...\n");
	if(argc > 3){
		/* Too many inputs. Expecting something of the form: setting value */
		return 800;
	}
	
	if(strcmp(argv[1], "order"){
		if(!(argv[2]==2 || argv[2]==4){
			/* Invalid input */
			return 801;
		}
		*order = argv[2];
	} else if(strcmp(argv[1], "dimension"){
		if(!(argv[2]==1 || argv[2]==2){
			/* Invalid input */
			return 802;			
		}
		*dimension = argv[2];
	} else if(strcmp(argv[1], "verification"){
		if(!(argv[2]==0 || argv[2]==1){
			/* Invalid input */
			return 803;
		}
		*verification = argv[2];
	} else if(strcmp(argv[1], "out_mode"){
		if(!(argv[2]==0 || argv[2]==1 || argv[2]==2){
			/* Invalid input */
			return 804;
		}
		*out_mode = argv[2];
	} else if(strcmp(argv[1], "timer"){
		if(!(argv[2]==0 || argv[2]==1){
			/* Invalid input */
			return 805;
		}
		*timer = argv[2];
	} else if(strcmp(argv[1], "output_file"){
		*output_file = argv[2];
	} else{
		/* Invalid setting requested */
		return 806;
	}
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Parsing command line options...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	/* Clean exit */
	return 0;
} 
/*--------------------------------------------------------------------------------------
 *  assemble: Routine for creating system matrix and vectors
 *-------------------------------------------------------------------------------------*/

int assemble(double *A, double *T, double *Q, int dimension, int order, int n, int size, double kappa, double Txmin, double Txmax, double Tymin, double Tymax, int verification, double Ax, double By)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Assembling system...\n");
	
	/* Mesh spacing calculations */
	double h;
	h = 1/(n-1) * (1.0 - 0.0);
	grvy_printf(GRVY_DEBUG, "Mesh spacing: %f\n", h);
	
	/* Iteration variables */
	int i;
	int j;
	int diag;
	
	/* Prefill A, Q, T with zeros */
	/* Suboptimal, but easy to get going */
	/* Thrashing by not writing sequentially? */
	for(i=0; i<size, i++){
		for(j=0; j<size; j++){
			A[i*size + j] = 0.0;
		}
		T[i] = 0.0;
		Q[i] = 0.0;	
	}
	
	if(dimension==1){
		/* MASA */
		if(verification){
			masa_init("verify", "heateq_1d_steady_const");
			masa_set_param("k_0", kappa);
			masa_set_param("A_x", Ax);
		}
		if(order==2){
			for(i=0; i<size; i++){
				/* Fill A */
				/* Get first interior node */
				if(i==0){
					A[i*size] = 2.0;
					A[i*size + 1] = -1.0;
				/* Get last interior node */
				} else if(i==size-1){
					A[i*size + size-1] = 2.0;
					A[i*size + size-2] = -1.0;
				/* Get the rest of the interior nodes */
				} else {
					diag = i;
					A[i*size + diag] = 2.0;
					A[i*size + diag + 1] = -1.0;
					A[i*size + diag - 1] = -1.0;
				}
				 
				/* Fill Q */
				if(verification){
					/* Bump the index by one to account for boundary node index */
					Q[i] += masa_eval_1d_source_t((i+1)*h);
				}
				/* Multiply by h^2/kappa */
				Q[i] = h*h/kappa * Q[i];
				
			}
			/* Fill Q with BC T values */
			Q[0] += Txmin;
			Q[size-1] += Txmax;
		} else if(order==4){
			for(i=0; i<size; i++){
				/* Fill A */
				/* Get first interior node */
				} else if(i==0){
					A[i*size + 0] = 2.0;
					A[i*size + 1] = -1.0;
				/* Get second interior node */
				} else if(i==1){
					A[i*size + 0] = -4.0/3.0;
					A[i*size + 1] = 5.0/2.0;
					A[i*size + 2] = -4.0/3.0;
					A[i*size + 3] = 1.0/12.0;
				/* Get second to last interior node */
				} else if(i==size-2){
					A[i*size + size-1] = -4.0/3.0;
					A[i*size + size-2] = 5.0/2.0;
					A[i*size + size-3] = -4.0/3.0;
					A[i*size + size-4] = 1.0/12.0;
				/* Get last interior node */
				} else if(i==size-1){
					A[i*size + size-1] = 2.0;
					A[i*size + size-2] = -1.0;
				/* Get the rest of the interior nodes */
				} else {
					diag = i;
					A[i*size + diag] = 5.0/2.0;
					A[i*size + diag + 1] = -4.0/3.0;
					A[i*size + diag - 1] = -4.0/3.0;
					A[i*size + diag + 2] = 1.0/12.0;
					A[i*size + diag - 2] = 1.0/12.0;
				}
				 
				/* Fill Q */
				if(verification){
					/* Bump the index by one to account for boundary node index */
					Q[i] += masa_eval_1d_source_t((i+1)*h);
				}
				/* Multiply by h^2/kappa */
				Q[i] = h*h/kappa * Q[i];
			}
			/* Fill Q with BC T values */
			Q[0] += Txmin;
			Q[1] -= 1/12 * Txmin;
			Q[size-2] -= 1/12 * Txmax;
			Q[size-1] += Txmax;			
		} else {
			grvy_printf(GRVY_INFO, "Incorrect finite difference order! Must be either %d or %d\n", 2, 4);
			return 601;
		}
		
	} else if(dimension==2){
		/* Number of interior nodes is n-2 */
		int interior;
		interior = n-2;
		/* MASA */
		if(verification){
			masa_init("verify", "heateq_2d_steady_const");
			masa_set_param("k_0", kappa);
			masa_set_param("A_x", Ax);
			masa_set_param("B_y", By);
		}
		if(order==2){
			/* Fill A */
			for(i=0; i<size; i=i+interior;){
				for(j=0; j<interior; j++){
					/* Special case: first node row*/
					if(i==0){
						if(j==0){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j + 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j - 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
						} else {
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j - 1] = -1.0;
							A[(i+j)*size + j + 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
						}
					/* Special case: last node row */
					} else if(i==size-interior){
						if(j==0){
							A[(i+j)*size + size - interior] = 4.0;
							A[(i+j)*size + size - interior + 1] = -1.0;
							A[(i+j)*size + size - 2 * interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + size - 1] = 4.0;
							A[(i+j)*size + size - 2] = -1.0;
							A[(i+j)*size + size - 1 - interior] = -1.0;
						} else {
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j - 1] = -1.0;
							A[(i+j)*size + i + j + 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
						}
					/* Interior nodes */
					} else {
						if(j==0){
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j + 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
							A[(i+j)*size + i + j + interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j - 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
							A[(i+j)*size + i + j + interior] = -1.0;
						} else {
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j - 1] = -1.0;
							A[(i+j)*size + i + j + 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
							A[(i+j)*size + i + j + interior] = -1.0;
						}
						
					}
				}
			}
			/* Fill Q */
			if(verification){
				for(i=0; i<interior; i++){
					for(j=0; j<interior; j++){
						/* Bump index because these are interior nodes */
						Q[i*interior + j] += masa_eval_2d_source_t((j+1)*h, (i+1)*h);
					}
				}
			}
			for(i=0; i<size; i++){
				/* Multiply by h^2/kappa */
				Q[i] = h*h/kappa * Q[i];	
			}
			/* Fill Q with BC T values */
			/* Bottom left and top right nodes */
			Q[0] += Txmin + Tymin;
			Q[size-1] += Txmax + Tymax;
			/* Top left and bottom right nodes */
			Q[0*interior + interior - 1] += Txmax + Tymin;
			Q[(interior-1) * interior] += Txmin + Tymax;			
			/* Fill left and right inner column */
			for(i=1;i<interior-1;i++){
				Q[i*interior + 0] += Txmin;
				Q[i*interior + interior - 1] += Txmax;
			}
			/* Fill top and bottom row */
			for(i=1; i<interior-1; i++){
				Q[0*interior + i] += Tymin;
				Q[(interior-1)*interior + i] += Tymax;
			}
		} else if(order==4){
			/* Fill A */
			for(i=0; i<size; i=i+interior;){
				for(j=0; j<interior; j++){
					/* Special case: first node row*/
					if(i==0){
						if(j==0){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j + 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j - 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
						} else {
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j - 1] = -1.0;
							A[(i+j)*size + j + 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
						}
					/* Special case: second node row*/
					if(i==interior){
						if(j==0){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j + 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
							A[(i+j)*size + j - interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j - 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
							A[(i+j)*size + j - interior] = -1.0;
						} else if(j==1){
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j + 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j + 2*interior] = -1.0/12.0;
						} else if(j==interior-2){
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j - 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j + 2*interior] = -1.0/12.0;
						} else {
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j - 2] = -1.0/12.0;
							A[(i+j)*size + j + 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j + 2*interior] = -1.0/12.0;
						}
					/* Special case: second to last node row */
					} else if(i==size-2*interior){
						if(j==0){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j + 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
							A[(i+j)*size + j - interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + j] = 4.0;
							A[(i+j)*size + j - 1] = -1.0;
							A[(i+j)*size + j + interior] = -1.0;
							A[(i+j)*size + j - interior] = -1.0;
						} else if(j==1){
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j + 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j - 2*interior] = -1.0/12.0;
						} else if(j==interior-2){
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j - 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j - 2*interior] = -1.0/12.0;
						} else {
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j - 2] = -1.0/12.0;
							A[(i+j)*size + j + 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j - 2*interior] = -1.0/12.0;
					/* Special case: last node row*/
					} else if(i==size-interior){
						if(j==0){
							A[(i+j)*size + size - interior] = 4.0;
							A[(i+j)*size + size - interior + 1] = -1.0;
							A[(i+j)*size + size - 2 * interior] = -1.0;
						} else if(j==interior-1){
							A[(i+j)*size + size - 1] = 4.0;
							A[(i+j)*size + size - 2] = -1.0;
							A[(i+j)*size + size - 1 - interior] = -1.0;
						} else {
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j - 1] = -1.0;
							A[(i+j)*size + i + j + 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
						}
					/* Interior nodes */
					} else {
						if(j==0){
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j + 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
							A[(i+j)*size + i + j + interior] = -1.0;
						} else if(j==1){
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j + 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j - 2*interior] = -1.0/12.0;
							A[(i+j)*size + j + 2*interior] = -1.0/12.0;
						} else if(j==interior-2){
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j - 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j - 2*interior] = -1.0/12.0;
							A[(i+j)*size + j + 2*interior] = -1.0/12.0;
						} else if(j==interior-1){
							A[(i+j)*size + i + j] = 4.0;
							A[(i+j)*size + i + j - 1] = -1.0;
							A[(i+j)*size + i + j - interior] = -1.0;
							A[(i+j)*size + i + j + interior] = -1.0;
						} else {
							A[(i+j)*size + j] = 5;
							A[(i+j)*size + j - 1] = -4.0/3.0;
							A[(i+j)*size + j + 1] = -4.0/3.0;
							A[(i+j)*size + j - 2] = -1.0/12.0;
							A[(i+j)*size + j + 2] = -1.0/12.0;
							A[(i+j)*size + j + interior] = -4.0/3.0;
							A[(i+j)*size + j - interior] = -4.0/3.0;
							A[(i+j)*size + j - 2*interior] = -1.0/12.0;
							A[(i+j)*size + j + 2*interior] = -1.0/12.0;
						}
						
					}
				}
			}
			/* Fill Q */
			if(verification){
				for(i=0; i<interior; i++){
					for(j=0; j<interior; j++){
						/* Bump index because these are interior nodes */
						Q[i*interior + j] += masa_eval_2d_source_((j+1)*h, (i+1)*h);
					}
				}
			}
			for(i=0; i<size; i++){
				/* Multiply by h^2/kappa */
				Q[i] = h*h/kappa * Q[i];	
			}
			/* Fill Q with BC T values */
			/* First the 2D order nodes
			/* Bottom left and top right nodes */
			Q[0] += Txmin + Tymin;
			Q[size-1] += Txmax + Tymax;
			/* Top left and bottom right nodes */
			Q[0*interior + interior - 1] += Txmax + Tymin;
			Q[(interior-1) * interior] += Txmin + Tymax;			
			/* Fill left and right inner column */
			for(i=1;i<interior-1;i++){
				Q[i*interior + 0] += Txmin;
				Q[i*interior + interior - 1] += Txmax;
			}
			/* Fill top and bottom row */
			for(i=1; i<interior-1; i++){
				Q[0*interior + i] += Tymin;
				Q[(interior-1)*interior + i] += Tymax;
			}
			/* Now the 4th order nodes */
			/* Bottom left and top right nodes */
			Q[interior + 1] += -1.0/12.0*(Txmin + Tymin);
			Q[size-1 - interior - 1] += -1.0/12.0 * (Txmax + Tymax);
			/* Top left and bottom right nodes */
			Q[1*interior + interior - 2] += -1.0/12.0 * (Txmax + Tymin);
			Q[(interior-2)*interior + 1] += -1.0/12.0 * (Txmin + Tymax);
			/* Fill left and right inner column */
			for(i=2; i<interior-2; i++){
				Q[i*interior + 1] += -1.0/12.0 * (Txmin);
				Q[i*interior + interior - 2] += -1.0/12.0 * (Txmax);
			}
			/* Fill top and bottom row */
			for(i=2; i<interior-2; i++){
				Q[1*interior + i] += -1.0/12.0 * (Tymin);
				Q[(interior-2)*interior + i] += -1.0/12.0 * (Tymax);
			}
		} else {
			grvy_printf(GRVY_INFO, "Incorrect finite difference order! Must be either %d or %d\n", 2, 4);
			return 601;
		}
		
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 600;
	}
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Assembling system...DONE\n");
	grvy_printf(GRVY_DEBUG, "A=\n");
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			grvy_printf(GRVY_DEBUG, "%.3f ", A[i*size + j]);
		}
		grvy_printf(GRVY_DEBUG, "\n");
	}
	grvy_printf(GRVY_DEBUG, "Q=\n");
	for(i=0; i<size; i++){
		grvy_printf(GRVY_DEBUG, "%.3f\n", Q[i]);
	}
	/* Timing */
	grvy_timer_end(__func__);
	
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  solver: Routine for solving linear system
 *-------------------------------------------------------------------------------------*/

int solver(double *A, double *T, double *Q, int size, int tol, int itermax, char method)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Solving system...\n");
	
	/* Variables needed for iterating */
	float converge;
	converge = 1.0;
	int k;
	k = 1;
	int i;
	int j;
	
	/* Variables needed for error and updating */
	double sigma;
	double l2;
	l2 = 0.0;
	double delta;
	if(strcmp(method, "Gauss-Seidal"){
		while(converge >= tol){
			for(i=0; i<size; i++){
				sigma = 0.0;
				for(j=0; j<size; j++){
					if(i != j){
						sigma += A[i*size + j]*T[j];
					}
				}
				/* Grab amount delta to ith entry */
				delta = 1/A[i*size + i] * (Q[i] - sigma);
				/* Add it to iterate solution */
				T[i] = delta;
				/* Compute the sum of change squared for L2 error */
				l2 += delta * delta;
			}
			/* Compute final l2 error */
			l2 = sqrt(1/size * l2);
			grvy_printf(GRVY_DEBUG, "Iteration: %d\t Iteration Norm: %f\n", k, l2);
			if (l2 <= tol){
				grvy_printf(GRVY_INFO, "Converged at iteration: %d", k);
				break;
			}
			/* Check for max iterations and increment */
			if (k >= iter_max){
				grvy_printf(GRVY_INFO, "Maximum number of iterations (%d) reached. Iteration Norm: %f", iter_max, l2);
				return 500;
			}
			/* Increment iteration counter */
			k++;
			/* Set L2 error back to zero */
			l2 = 0.0;
		}
	} else if(strcmp(method, "Jacobi"){
		/* Allocate array for holding k iteration values */
		double *Tk;
		Tk = malloc(size * sizeof(double));
		/* Fill with zeros */
		for(i=0; i<size; i++){
			Tk[i] = 0.0;
		}
		while(converge >= tol){
			for(i=0; i<size; i++){
				sigma = 0.0;
				for(j=0; j<size; j++){
					if(i != j){
						sigma += A[i*size + j]*Tk[j];
					}
				}
				/* Grab amount delta to ith entry */
				delta = 1/A[i*size + i] * (Q[i] - sigma);
				/* Add it to iterate solution */
				T[i] = delta;
				/* Compute the sum of change squared for L2 error */
				l2 += delta * delta;	
			}
			/* Compute final l2 error */
			l2 = sqrt(1/size * l2);
			grvy_printf(GRVY_DEBUG, "Iteration: %d\t Iteration Norm: %f\n", k, l2);
			if (l2 <= tol){
				grvy_printf(GRVY_INFO, "Converged at iteration: %d", k);
				free(Tk);
				break;
			}
			/* Check for max iterations and increment */
			if (k >= iter_max){
				grvy_printf(GRVY_INFO, "Maximum number of iterations (%d) reached. Iteration Norm: %f", iter_max, l2);
				free(Tk);
				return 500;
			}
			/* Increment iteration counter */
			k++;
			/* Set L2 error back to zero */
			l2 = 0.0;
			/* Pass values of Tk to T so it becomes the k-1 iteration */
			for(i=0; i<size; i++){
				Tk[i] = T[i];
			}
		}		
	} else {
		return 501;
	}	
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Solving system...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  errornorm: Routine for computing error between MASA and computed solution
 *-------------------------------------------------------------------------------------*/

double errornorm(double *T, int size, int n, int dimension, double kappa, double Ax, double By)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Computing L2 error norm against MASA soln...\n");
	
	/* Mesh spacing calculations */
	double h;
	h = 1/(n-1) * (1.0 - 0.0);
	
	/* Iterators */
	int i;
	int j;
	/* Error calcs */
	double l2;
	l2 = 0.0;
	double delta;
	
	if(dimension==1){
		masa_init("error", "heateq_1d_steady_const");
		masa_set_param("A_x", Ax);
		masa_set_param("k_0", kappa);
		for(i=0; i<size; i++){
			delta = T[i] - masa_eval_1d_exact_t((i+1)*h);
			l2 += delta * delta;
		}
	} else{
		masa_init("error", "heateq_2d_steady_const");
		masa_set_param("A_x", Ax);
		masa_set_param("B_y", By);
		masa_set_param("k_0", kappa);
		for(i=0; i<n-2; i++){
			for(j=0; j<n-2; j++){
				/* not declaring a variable interior=n-2 for this case */
				delta = T[i*(n-2) + j] - masa_eval_2d_exact_t((j+1)*h, (i+1)*h);
				l2 += delta * delta;
			}
		}
	}
	
	l2 = sqrt(1/size * l2);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Computing L2 error norm against MASA soln\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	return l2;
}

/*--------------------------------------------------------------------------------------
 *  output: Routine for writing solution to file
 *-------------------------------------------------------------------------------------*/

int output(double *T, int n, int size, int dimension, double kappa, double Ax, double By, char output_file, int verification)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Writing results to file...\n");
	
	FILE *fp;
	fp = fopen(output_file, "w");
	if(fp==NULL){
		return 400;
	}
	
	grvy_printf(GRVY_INFO, "Writing to %-10s", output_file);
	/* Write values */
	if(dimension==1){
		fprintf(fp, "Node:\t X Coordinate:\t Temp:\n");
		int i;
		float xcoord;
		for(i=0; i<n; i++){
			xcoord = i/(n - 1);
			fprintf(fp, "%d\t %f\t %f\n", i, xcoord, T[i]);
		}
	} else if(dimension==2){
		fprintf(fp, "Node:\t X Coordinate:\t Y Coordinate:\t Temp:\n");
		int i;
		int j;
		float xcoord;
		float ycoord;
		int node;
		for(i=0; i<n; i++){
			ycoord = i/(n-1);
			for(j=0; j<n; j++){
				node = i*n + j;
				xcoord = j/(n-1);
				fprintf(fp, "%d\t %f\t %f\t %f\n", node, xcoord, ycoord, T[node]);
			}
		}
	}
	
	/* Compute L2 error if necessary */
	if(verification){
		double errnorm;
		errnorm = errornorm(T, size, n, dimension, kappa, Ax, By);
		grvy_printf(GRVY_INFO, "L2 error norm: %f", errnorm);
		fprintf(fp, "\n L2 Error Norm:\n");
		fprintf(fp,"f",errnorm);
	}
	
	/* Close and free */
	fclose(fp);
	free(fp);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Writing results to file...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  cleanup: Function for cleaning up and freeing dynamically allocated memory
 *-------------------------------------------------------------------------------------*/
 
void cleanup(double *A, double *T, double *Q, char *method, char *output_file)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Freeing variables...\n");

	free(A);
	free(T);
	free(Q);
	free(method);
	free(output_file);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Freeing variables...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
}
