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
 *  input: Routine for parsing input file and returning variables for rest of the code.
 *-------------------------------------------------------------------------------------*/

int input(Heat* problem, char* inputfile)
{

	/* Timing */
	grvy_timer_begin(__func__);

	/* Declare variables for file open */
	int erropen;

	/* Open file */
	erropen = grvy_input_fopen(inputfile);
	if(!(erropen)){
		grvy_printf(GRVY_DEBUG,"Unable to open file %-10s\n",*inputfile);
		return 100;
	}

	/* Parsing logging and timing details */
	if(!(grvy_input_fread_int("output/out_mode", &(*problem).out_mode))){
		return 101;
	}
	
	/* Set log level */
	if((*problem).out_mode == 0){
		grvy_log_setlevel(GRVY_NOLOG);
	} else if((*problem).out_mode == 1){
		grvy_log_setlevel(GRVY_INFO);
	} else if((*problem).out_mode == 2){
		grvy_log_setlevel(GRVY_DEBUG);
	} else {
		grvy_log_setlevel(GRVY_DEBUG);
	}

	grvy_printf(GRVY_INFO, "Solving the Steady State Heat Equation with Finite Difference \n \n");
	grvy_printf(GRVY_INFO, "Parsing Log settings...\n");
	grvy_printf(GRVY_INFO, "Parsing Log settings...DONE\n");
	
	/* Read Problem details */
	grvy_printf(GRVY_INFO, "Parsing problem details...\n");

	if(!(grvy_input_fread_char("solver/method", &(*problem).method))){
		grvy_printf(GRVY_DEBUG,"Unable to read char value %-10s\n","solver/method"); 
		grvy_input_fclose();
		return 102;
	 }

	if(!(grvy_input_fread_int("solver/iter_max", &(*problem).iter_max))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n", "solver/iter_max");
		grvy_input_fclose();
		return 103;
	}

	if(!(grvy_input_fread_double("solver/tol", &(*problem).tol))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n","solver/tol");
		grvy_input_fclose();	
		return 104;
	}

	if(!(grvy_input_fread_int("solver/order", &(*problem).order))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n", "solver/order");
		grvy_input_fclose();	
		return 105;
	}

	if(!(grvy_input_fread_int("problem/dimension", &(*problem).dimension))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n","problem/dimension"); 
		grvy_input_fclose();
		return 106;
	}

	if(!(grvy_input_fread_double("problem/kappa", &(*problem).kappa))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/kappa"); 
		grvy_input_fclose();
		return 107;
	}

	if(!(grvy_input_fread_int("problem/n", &(*problem).n))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n","problem/n");
		grvy_input_fclose();
		return 108;
	}
	if(!(grvy_input_fread_double("problem/Txmin", &(*problem).Txmin))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Txmin");
		grvy_input_fclose();
		return 109;
	}

	if(!(grvy_input_fread_double("problem/Txmax", &(*problem).Txmax))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Txmax");
		grvy_input_fclose();
		return 110;
	}

	if(!(grvy_input_fread_double("problem/Tymin", &(*problem).Tymin))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Tymin");
		grvy_input_fclose();
		return 111;
	}

	if(!(grvy_input_fread_double("problem/Tymax", &(*problem).Tymax))){
		grvy_printf(GRVY_DEBUG,"Unable to read double value %-10s\n", "problem/Tymax");
		grvy_input_fclose();
		return 112;
	}
	
	/* MASA */
	if(!(grvy_input_fread_double("MASA/Ax", &(*problem).Ax))){
		grvy_printf(GRVY_DEBUG, "Unable to read double value %-10s\n", "Ax");
		grvy_input_fclose();
		return 113;
	}

	if(!(grvy_input_fread_double("MASA/By", &(*problem).By))){
		grvy_printf(GRVY_DEBUG, "Unable to read double value %-10s\n", "By");
		grvy_input_fclose();
		return 114;
	}

	/* Output settings */
	if(!(grvy_input_fread_int("output/verification", &(*problem).verification))){
		grvy_printf(GRVY_DEBUG,"Unable to read int value %-10s\n", "verification");
		grvy_input_fclose();
		return 115;
	}

	if(!(grvy_input_fread_char("output/out_file", &(*problem).output_file))){
		grvy_printf(GRVY_DEBUG, "Unable to read char value %-10s\n", "output_file");
		grvy_input_fclose();
		return 116;
	}

	grvy_input_fclose();

	/* Echo settings to stdout */
	grvy_printf(GRVY_INFO, "Parsing problem details...DONE\n");
	
	grvy_printf(GRVY_INFO, "\n Input file: %-10s\n", inputfile);
	
	grvy_printf(GRVY_INFO, "\n Mesh Settings:\n");
	grvy_printf(GRVY_INFO,"--Dimension: %d\n", (*problem).dimension);
	if((*problem).dimension==1){
		grvy_printf(GRVY_INFO,"--nx: %d\n", (*problem).n);
		grvy_printf(GRVY_INFO,"--(xmin, xmax): (%f, %f)\n", 0.0, 1.0);
	} else if((*problem).dimension==2){
		grvy_printf(GRVY_INFO,"--nx: %d\n", (*problem).n);
		grvy_printf(GRVY_INFO,"--(xmin, xmax): (%f, %f)\n", 0.0, 1.0);
		grvy_printf(GRVY_INFO,"--ny: %d\n", (*problem).n);
		grvy_printf(GRVY_INFO,"--(ymin, ymax): (%f, %f)\n", 0.0, 1.0);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 117;
	}
	
	grvy_printf(GRVY_INFO, "\n Boundary Conditions:\n");
	if((*problem).dimension==1){
		grvy_printf(GRVY_INFO,"--(Txmin, Txmax): (%f, %f)\n", (*problem).Txmin, (*problem).Txmax);
	} else if((*problem).dimension==2){
		grvy_printf(GRVY_INFO,"--(Txmin, Txmax): (%f, %f)\n", (*problem).Txmin, (*problem).Txmax);
		grvy_printf(GRVY_INFO,"--(Tymin, Tymax): (%f, %f)\n", (*problem).Tymin, (*problem).Tymax);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 117;
	}

	grvy_printf(GRVY_INFO, "\n MASA Settings:\n");
	if((*problem).dimension==1){
		grvy_printf(GRVY_INFO,"--Ax: %.15f\n", (*problem).Ax);
	} else if((*problem).dimension==2){
		grvy_printf(GRVY_INFO,"--Ax: %.15f\n", (*problem).Ax);
		grvy_printf(GRVY_INFO,"--By: %.15f\n", (*problem).By);
	} else {
		grvy_printf(GRVY_INFO,"Incorrect number of dimensions! Must be either %d or %d\n", 1, 2);
		return 117;
	}

	grvy_printf(GRVY_INFO, "\n Solver Settings:\n");
	grvy_printf(GRVY_INFO,"--Method: %-10s\n", (*problem).method);
	grvy_printf(GRVY_INFO,"--Max iterations: %d\n", (*problem).iter_max);
	grvy_printf(GRVY_INFO,"--Convergence Tolerance: %E\n", (*problem).tol);
	grvy_printf(GRVY_INFO,"--Finite Difference Order: %d\n", (*problem).order);
	
	grvy_printf(GRVY_INFO, "\n Output Settings:\n");
	grvy_printf(GRVY_INFO,"--Verification: %d (0=no, 1=yes)\n",(*problem).verification);
	grvy_printf(GRVY_INFO,"--Output Mode: %d (0=silent, 1=standard, 2=debug)\n", (*problem).out_mode);
	grvy_printf(GRVY_INFO,"--Output File: %-10s\n \n", (*problem).output_file);

	grvy_timer_end(__func__);

	/* Return zero on clean exit */
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  sanitize: Sanitize input (Sanity checks on dimension, order, etc.)
 *-------------------------------------------------------------------------------------*/

int sanitize(Heat* problem)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Sanity checking input...\n");

	/* Check order */
	if(!((*problem).order == 2 || (*problem).order == 4)){
		grvy_printf(GRVY_DEBUG, "Incorrect problem order! Must be either %d or %d\n", 2, 4);
		return 800;
	}

	/* Check dimension */
	if(!((*problem).dimension == 1 || (*problem).dimension == 2)){
		grvy_printf(GRVY_DEBUG, "Incorrect problem dimension! Must be either %d or %d\n", 1, 2);
		return 801;
	}

	/* Check that the number of nodes makes sense */
	if(!((*problem).n >= 3)){
		grvy_printf(GRVY_DEBUG, "Must have at least %d nodes\n", 3);
		return 801;
	}	

	/* Check that the file name ends with .h5 */
	char* endptr = (*problem).output_file + strlen((*problem).output_file) - 2;
	char ptr5[3];
	strcpy(ptr5, "h5");
	grvy_printf(GRVY_DEBUG, "endptr: %-10s\n", endptr);
	if(strcmp(endptr, ptr5)){
		grvy_printf(GRVY_DEBUG, "Output file must be HDF5 compatible\n");
		return 801;
	}	
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Sanity checking input...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	/* Clean exit */
	return 0;
} 
/*--------------------------------------------------------------------------------------
 *  assemble: Routine for creating system matrix and vectors
 *-------------------------------------------------------------------------------------*/

int assemble(Heat* problem)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Assembling system...\n");

	/* Allocate memory */
	if((*problem).dimension == 1){
		if((*problem).order == 2){
			(*problem).size = (*problem).n - 2;
		} else if((*problem).order == 4){
			(*problem).size = (*problem).n - 4;
		} else {
			return 901;
		}
	} else if((*problem).dimension == 2){
		if((*problem).order == 2){
			(*problem).size = ((*problem).n - 2) * ((*problem).n - 2);
		} else if((*problem).order == 4){
			(*problem).size = ((*problem).n - 4) * ((*problem).n - 4);
		} else {
			return 902;
		}
	} else {
		return 900;
	}

	(*problem).A = malloc((*problem).size * (*problem).size * sizeof(double));
	(*problem).T = malloc((*problem).size * sizeof(double));		
	(*problem).Q = malloc((*problem).size * sizeof(double));
	
	/* Mesh spacing calculations */
	double h;
	h = 1.0/((*problem).n-1) * (1.0 - 0.0);
	grvy_printf(GRVY_INFO, "Mesh spacing: %f\n", h);
	
	/* Iteration variables */
	int i;
	int j;
	int diag;
	
	/* Prefill A, Q, T with zeros */
	/* Suboptimal, but easy to get going */
	/* Thrashing by not writing sequentially? */
	for(i=0; i<(*problem).size; i++){
		for(j=0; j<(*problem).size; j++){
			(*problem).A[i*(*problem).size + j] = 0.0;
		}
		(*problem).T[i] = 0.0;
		(*problem).Q[i] = 0.0;	
	}
	if((*problem).dimension==1){
		/* MASA */
		if((*problem).verification){
			masa_init("verify", "heateq_1d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
			masa_display_param();
		}
		if((*problem).order==2){
			for(i=0; i<(*problem).size; i++){
				/* Fill A */
				/* Get first node */
				if(i==0){
					(*problem).A[i*(*problem).size] = 2.0;
					(*problem).A[i*(*problem).size + 1] = -1.0;
				/* Get last  node */
				} else if(i==(*problem).size-1){
					(*problem).A[(i+1)*(*problem).size - 1] = 2.0;	
					(*problem).A[(i+1)*(*problem).size - 2] = -1.0;
				/* Get the interior nodes */
				} else {
					diag = i;
					(*problem).A[i*(*problem).size + diag] = 2.0;
					(*problem).A[i*(*problem).size + diag + 1] = -1.0;
					(*problem).A[i*(*problem).size + diag - 1] = -1.0;
				}
				 
				/* Fill Q */
				if((*problem).verification){
					/* Bump the index by one to account for boundary node index */
					(*problem).Q[i] += masa_eval_1d_source_t((i+1)*h);
				}
				/* Multiply by h^2/kappa */
				(*problem).Q[i] = h*h/((*problem).kappa) * (*problem).Q[i];
				
			}
			/* Fill Q with BC T values */
			(*problem).Q[0] += (*problem).Txmin;
			(*problem).Q[(*problem).size-1] += (*problem).Txmax;
		} else if((*problem).order==4){
			for(i=0; i<(*problem).size; i++){
				/* Fill A */
				/* Get first interior node */
				/*
				if(i==0){
					(*problem).A[i*(*problem).size + 0] = 2.0;
					(*problem).A[i*(*problem).size + 1] = -1.0;
				*/
				if(i==0){
					(*problem).A[i*(*problem).size + 0] = 5.0/2.0;
					(*problem).A[i*(*problem).size + 1] = -4.0/3.0;
					(*problem).A[i*(*problem).size + 2] = 1.0/12.0;
				/* Get second interior node */
				} else if(i==1){
					(*problem).A[i*(*problem).size + 0] = -4.0/3.0;
					(*problem).A[i*(*problem).size + 1] = 5.0/2.0;
					(*problem).A[i*(*problem).size + 2] = -4.0/3.0;
					(*problem).A[i*(*problem).size + 3] = 1.0/12.0;
				/* Get second to last interior node */
				} else if(i==(*problem).size-2){
					(*problem).A[i*(*problem).size +(*problem).size-1] = -4.0/3.0;
					(*problem).A[i*(*problem).size +(*problem).size-2] = 5.0/2.0;
					(*problem).A[i*(*problem).size +(*problem).size-3] = -4.0/3.0;
					(*problem).A[i*(*problem).size +(*problem).size-4] = 1.0/12.0;
				/* Get last interior node */
				/*
				} else if(i==(*problem).size-1){
					(*problem).A[i*(*problem).size +(*problem).size-1] = 2.0;
					(*problem).A[i*(*problem).size +(*problem).size-2] = -1.0;
				*/
				} else if(i==(*problem).size-1){
					(*problem).A[i*(*problem).size +(*problem).size - 1] = 5.0/2.0;
					(*problem).A[i*(*problem).size +(*problem).size - 2] = -4.0/3.0;
					(*problem).A[i*(*problem).size +(*problem).size - 3] = 1.0/12.0;
				/* Get the rest of the interior nodes */
				} else {
					diag = i;
					(*problem).A[i*(*problem).size + diag] = 5.0/2.0;
					(*problem).A[i*(*problem).size + diag + 1] = -4.0/3.0;
					(*problem).A[i*(*problem).size + diag - 1] = -4.0/3.0;
					(*problem).A[i*(*problem).size + diag + 2] = 1.0/12.0;
					(*problem).A[i*(*problem).size + diag - 2] = 1.0/12.0;
				}
				 
				/* Fill Q */
				if((*problem).verification){
					/* Bump the index by one to account for boundary node index */
					(*problem).Q[i] += masa_eval_1d_source_t((i+2)*h);
				}
				/* Multiply by h^2/kappa */
				(*problem).Q[i] = h*h/((*problem).kappa) * (*problem).Q[i];
			}
			/* Fill Q with BC T values */
			/*
			(*problem).Q[0] += (*problem).Txmin;
			(*problem).Q[1] -= 1.0/12.0 * (*problem).Txmin;
			(*problem).Q[(*problem).size-2] -= 1.0/12.0 * (*problem).Txmax;
			(*problem).Q[(*problem).size-1] += (*problem).Txmax;	
			*/
			(*problem).Q[0] += -1.0/12.0 * masa_eval_1d_exact_t(0.0) + 4.0/3.0 * masa_eval_1d_exact_t(h);
			(*problem).Q[1] += -1.0/12.0 * masa_eval_1d_exact_t(0.0);
			(*problem).Q[(*problem).size - 2] += -1.0/12.0 * masa_eval_1d_exact_t(1.0);
			(*problem).Q[(*problem).size - 1] += -1.0/12.0 * masa_eval_1d_exact_t(1.0) + 4.0/3.0 * masa_eval_1d_exact_t(1.0 - h);
		}		
	} else if((*problem).dimension==2){
		/* Number of interior nodes is n-2 */
		int interior;
		interior = (*problem).n-2;
		/* MASA */
		if((*problem).verification){
			masa_init("verify", "heateq_2d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
			masa_set_param("B_y", (*problem).By);
		}
		if((*problem).order==2){
			/* Fill A */
			for(i=0; i<(*problem).size; i=i+interior){
				for(j=0; j<interior; j++){
					/* Special case: first node row*/
					if(i==0){
						if(j==0){
							(*problem).A[(i+j)*(*problem).size + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + interior] = -1.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + interior] = -1.0;
						}
					/* Special case: last node row */
					} else if(i==(*problem).size-interior){
						if(j==0){
							(*problem).A[(i+j)*(*problem).size +(*problem).size - interior] = 4.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - interior + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 2 * interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 1] = 4.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 2] = -1.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 1 - interior] = -1.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
						}
					/* Interior nodes */
					} else {
						if(j==0){
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -1.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -1.0;
						}
						
					}
				}
			}
			/* Fill Q */
			if((*problem).verification){
				for(i=0; i<interior; i++){
					for(j=0; j<interior; j++){
						/* Bump index because these are interior nodes */
						(*problem).Q[i*interior + j] += masa_eval_2d_source_t((j+1)*h, (i+1)*h);
					}
				}
			} 
			for(i=0; i<(*problem).size; i++){
				/* Multiply by h^2/kappa */
				(*problem).Q[i] = h*h/((*problem).kappa )* (*problem).Q[i];	
			}
			
			if((*problem).verification){
				/* Fill (*problem).Q with BC T values */
				/* Bottom left and top right nodes */
				(*problem).Q[0] += masa_eval_2d_exact_t(0.0, h) + masa_eval_2d_exact_t(h, 0.0);
				(*problem).Q[(*problem).size-1] += masa_eval_2d_exact_t(1.0, interior*h) + masa_eval_2d_exact_t(interior * h, 1.0);
				/* Top left and bottom right nodes */
				(*problem).Q[0*interior + interior - 1] += masa_eval_2d_exact_t(h, 1.0) + masa_eval_2d_exact_t(0.0, interior * h);
				(*problem).Q[(interior-1) * interior] += masa_eval_2d_exact_t(1.0, h) + masa_eval_2d_exact_t(interior * h, 0.0);			
				/* Fill left and right inner column */
				for(i=1;i<interior-1;i++){
					(*problem).Q[i*interior + 0] += masa_eval_2d_exact_t(0.0, i*h);
					(*problem).Q[i*interior + interior - 1] += masa_eval_2d_exact_t(1.0, i*h);
				}
				/* Fill top and bottom row */
				for(i=1; i<interior-1; i++){
					(*problem).Q[0*interior + i] += masa_eval_2d_exact_t(i*h, 0.0);
					(*problem).Q[(interior-1)*interior + i] += masa_eval_2d_exact_t(i*h, 1.0);
				}
			} else {
				/* Fill (*problem).Q with BC T values */
				/* Bottom left and top right nodes */
				(*problem).Q[0] += (*problem).Txmin + (*problem).Tymin;
				(*problem).Q[(*problem).size-1] += (*problem).Txmax + (*problem).Tymax;
				/* Top left and bottom right nodes */
				(*problem).Q[0*interior + interior - 1] += (*problem).Txmax + (*problem).Tymin;
				(*problem).Q[(interior-1) * interior] += (*problem).Txmin + (*problem).Tymax;			
				/* Fill left and right inner column */
				for(i=1;i<interior-1;i++){
					(*problem).Q[i*interior + 0] += (*problem).Txmin;
					(*problem).Q[i*interior + interior - 1] += (*problem).Txmax;
				}
				/* Fill top and bottom row */
				for(i=1; i<interior-1; i++){
					(*problem).Q[0*interior + i] += (*problem).Tymin;
					(*problem).Q[(interior-1)*interior + i] += (*problem).Tymax;
				}
			}
		} else if((*problem).order==4){
			/* Fill A */
			for(i=0; i<(*problem).size; i=i+interior){
				for(j=0; j<interior; j++){
					/* Special case: first node row*/
					if(i==0){
						if(j==0){
							(*problem).A[(i+j)*(*problem).size + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + interior] = -1.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + j + interior] = -1.0;
						}
					/* Special case: second node row*/
					} else if(i==interior){
						if(j==0){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 4.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 4.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -1.0;
						} else if(j==1){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 5.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 2*interior] = 1.0/12.0;
						} else if(j==interior-2){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 5.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 2*interior] = 1.0/12.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 5.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)- 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 2*interior] = 1.0/12.0;
						}
					/* Special case: second to last node row */
					} else if(i==(*problem).size-2*interior){
						if(j==0){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 4.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)+ interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 4.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -1.0;
						} else if(j==1){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 5.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)+ 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)+ 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)+ interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)- interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)- 2*interior] = 1.0/12.0;
						} else if(j==interior-2){
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 5.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 2*interior] = 1.0/12.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + (i+j)] = 5.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j) - 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)+ 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)+ interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)- interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + (i+j)- 2*interior] = 1.0/12.0;
						}
					/* Special case: last node row*/
					} else if(i==(*problem).size-interior){
						if(j==0){
							(*problem).A[(i+j)*(*problem).size +(*problem).size - interior] = 4.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - interior + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 2 * interior] = -1.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 1] = 4.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 2] = -1.0;
							(*problem).A[(i+j)*(*problem).size +(*problem).size - 1 - interior] = -1.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
						}
					/* Interior nodes */
					} else {
						if(j==0){
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -1.0;
						} else if(j==1){
							(*problem).A[(i+j)*(*problem).size + i + j] = 5.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 2*interior] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 2*interior] = 1.0/12.0;
						} else if(j==interior-2){
							(*problem).A[(i+j)*(*problem).size + i + j] = 5.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 2*interior] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 2*interior] = 1.0/12.0;
						} else if(j==interior-1){
							(*problem).A[(i+j)*(*problem).size + i + j] = 4.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -1.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -1.0;
						} else {
							(*problem).A[(i+j)*(*problem).size + i + j] = 5.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 1] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 2] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - interior] = -4.0/3.0;
							(*problem).A[(i+j)*(*problem).size + i + j - 2*interior] = 1.0/12.0;
							(*problem).A[(i+j)*(*problem).size + i + j + 2*interior] = 1.0/12.0;
						}
						
					}
				}
			}
			/* Fill Q */
			if((*problem).verification){
				for(i=0; i<interior; i++){
					for(j=0; j<interior; j++){
						/* Bump index because these are interior nodes */
						(*problem).Q[i*interior + j] += masa_eval_2d_source_t((j+1)*h, (i+1)*h);
					}
				}
			}
			for(i=0; i<(*problem).size; i++){
				/* Multiply by h^2/kappa */
				(*problem).Q[i] = h*h/((*problem).kappa) * (*problem).Q[i];	
			}
			/* Fill (*problem).Q with BC T values */
			/* First the 2D order nodes
			* Bottom left and top right nodes */
			(*problem).Q[0] += (*problem).Txmin + (*problem).Tymin;
			(*problem).Q[(*problem).size-1] += (*problem).Txmax + (*problem).Tymax;
			/* Top left and bottom right nodes */
			(*problem).Q[0*interior + interior - 1] += (*problem).Txmax + (*problem).Tymin;
			(*problem).Q[(interior-1) * interior] += (*problem).Txmin + (*problem).Tymax;			
			/* Fill left and right inner column */
			for(i=1;i<interior-1;i++){
				(*problem).Q[i*interior + 0] += (*problem).Txmin;
				(*problem).Q[i*interior + interior - 1] += (*problem).Txmax;
			}
			/* Fill top and bottom row */
			for(i=1; i<interior-1; i++){
				(*problem).Q[0*interior + i] += (*problem).Tymin;
				(*problem).Q[(interior-1)*interior + i] += (*problem).Tymax;
			}
			/* Now the 4th order nodes */
			/* Bottom left and top right nodes */
			(*problem).Q[interior + 1] += -1.0/12.0*((*problem).Txmin + (*problem).Tymin);
			(*problem).Q[(*problem).size-1 - interior - 1] += -1.0/12.0 * ((*problem).Txmax + (*problem).Tymax);
			/* Top left and bottom right nodes */
			(*problem).Q[1*interior + interior - 2] += -1.0/12.0 * ((*problem).Txmax + (*problem).Tymin);
			(*problem).Q[(interior-2)*interior + 1] += -1.0/12.0 * ((*problem).Txmin + (*problem).Tymax);
			/* Fill left and right inner column */
			for(i=2; i<interior-2; i++){
				(*problem).Q[i*interior + 1] += -1.0/12.0 * ((*problem).Txmin);
				(*problem).Q[i*interior + interior - 2] += -1.0/12.0 * ((*problem).Txmax);
			}
			/* Fill top and bottom row */
			for(i=2; i<interior-2; i++){
				(*problem).Q[1*interior + i] += -1.0/12.0 * ((*problem).Tymin);
				(*problem).Q[(interior-2)*interior + i] += -1.0/12.0 * ((*problem).Tymax);
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
	for(i=0; i<(*problem).size; i++){
		for(j=0; j<(*problem).size; j++){
			grvy_printf(GRVY_DEBUG, "%.3f ", (*problem).A[i*(*problem).size + j]);
		}
		grvy_printf(GRVY_DEBUG, "\n");
	}
	grvy_printf(GRVY_DEBUG, "Q=\n");
	for(i=0; i<(*problem).size; i++){
		grvy_printf(GRVY_DEBUG, "%.3f\n", (*problem).Q[i]);
	}
	/* Timing */
	grvy_timer_end(__func__);
	
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  petsc: Routine for solving linear system with petsc
 *-------------------------------------------------------------------------------------*/
#ifdef INCLUDE_PETSC
	int petsc(Heat *problem)
	{
		/* Timing */
		grvy_timer_begin(__func__);
		/* Logging */
		grvy_printf(GRVY_INFO, "Solving system with PETSC...\n");
		/* Declare petsc variables */
		PetscErrorCode petstatus;
		Mat A;
		Vec Rhs, Sol;
		KSP Solver;
		petstatus=MatCreate(PETSC_COMM_SERIAL, &A);
		petstatus=MatSetSizes(&A, PETSC_DECIDE, PETSC_DECIDE, (*problem).size, (*problem).size);
		
		/* Just showing prototype code for 1d 2nd for now until petsc compilation fixed */
		int i;
		int j;
		/* MASA */
		if((*problem).verification){
			masa_init("verify", "heateq_1d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
			masa_display_param();
		
		for(i=0; i<(*problem).size; i++){
			for(j=0; j<(*problem).size; j++){
				petstatus=MatSetSizes(&A, 1, &i, 1, &j, 0.0, INSERT_VALUES);
				/* Fill A */
				/* Get first node */
				if(i==0){
					petstatus=MatSetSizes(&A, 0, &i, 0, &j, 2.0, INSERT_VALUES);
					petstatus=MatSetSizes(&A, 0, &i, 1, &j, -1.0, INSERT_VALUES);	
				/* Get last  node */
				} else if(i==(*problem).size-1 && j==(problem).size - 1){
					petstatus=MatSetSizes(&A, i, &i, j, &j, 2.0, INSERT_VALUES);
					petstatus=MatSetSizes(&A, i, &i, j-1, &j, -1.0, INSERT_VALUES);
				/* Get the interior nodes */
				} else {
					petstatus=MatSetSizes(&A, i, &i, j, &j, 2.0, INSERT_VALUES);
					petstatus=MatSetSizes(&A, i, &i, j+1, &j, -1.0, INSERT_VALUES);
					petstatus=MatSetSizes(&A, i, &i, j-1, &j, -1.0, INSERT_VALUES);
				}
			}
		}
		petstatus=MatView(A,0);
		petstatus=MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
		petstatus=MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
		/* Solver setup */
		petstatus = KSPCreate(PETSC_COMM_SERIAL, &Solver);
		/* RHS and Sol setup */
		petstatus=VecCreate(PETS_COMM_SELF, &Rhs);
		petstatus=VecSetSize(RHS, PETSC_DECIDE, (*problem).size);
		pestatus=VecDuplicate(RHS, &Sol);
		for(i=0; i<(*problem).size; i++){
			/* Fill Q */ 
			if((*problem).verification){
				/* Bump the index by one to account for boundary node index */
				double val;
				val = h*h/((*problem).kappat) * masa_eval_1d_source_t((i+1)*h);
				petstatus=VecSetValues(Rhs, 1, &i, &val, INSERT_VALUES);
			}
		}
		/* Fill Q with BC T values */
		petstatus=VecSetValues(RHS, 1, 0, &(*problem).Txmin, ADD_VALUES);
		petstatus=VecSetValues(RHS, 1, (*problem).size - 1, &(*problem).Txmax, ADD_VALUES);
		petstatus=VecAssemblyBegin(RHS);
		petstatus=VecAssemblyEnd(RHS);
		petstatus=VecAssemblyBegin(Sol);
		petstatus=VecAssemblyEnd(Sol);

		/* Solve system */
		petstatus=KSPSolve(Solver, Rhs, Sol);
		
		/* Write solution to T member */
		
		
		/* Clean up */
		petstatus=VecDestroy(Rhs);
		petstatus=VecDestroy(Sol);
		petstatus=KSPDestroy(Solver);
		petstatus=MatDestroy(&A);
		
		
		/* Logging */
		grvy_printf(GRVY_INFO, "Solving system with PETSC...DONE\n");
		/* Timing */
		grvy_timer_end(__func__);
		
		return 0;
	}

#endif

/*--------------------------------------------------------------------------------------
 *  solver: Routine for solving linear system
 *-------------------------------------------------------------------------------------*/

int solver(Heat* problem)
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
	double new;
	double delta;
	if(!(strcmp((*problem).method, "Gauss-Seidel"))){
		grvy_printf(GRVY_DEBUG, "Using solver method Gauss-Seidel\n");
		while(converge >= (*problem).tol){
			for(i=0; i<(*problem).size; i++){
				sigma = 0.0;
				for(j=0; j<(*problem).size; j++){
					if(i != j){
						sigma += (*problem).A[i*(*problem).size + j]*(*problem).T[j];
					}
				}
				/* Grab amount delta to ith entry */
				new = 1.0/((*problem).A[i*(*problem).size + i]) * ((*problem).Q[i] - sigma);
				/* Compute difference between new value and old value */
				delta = new - (*problem).T[i];
				/* Set new  iterate solution value */
				(*problem).T[i] = new;
				/* Compute the sum of change squared for L2 error */
				l2 += delta * delta;
			}
			/* Compute final l2 error */
			l2 = sqrt(1.0/((*problem).size) * l2);
			grvy_printf(GRVY_DEBUG, "Iteration: %d\t Iteration Norm: %E\t Tolerance: %E\n", k, l2, (*problem).tol);
			if (l2 <= (*problem).tol){
				grvy_printf(GRVY_INFO, "Converged at iteration: %d\n", k);
				break;
			}
			/* Check for max iterations and increment */
			if (k >= (*problem).iter_max){
				grvy_printf(GRVY_INFO, "Maximum number of iterations (%d) reached. Iteration Norm: %E\n", (*problem).iter_max, l2);
				return 500;
			}
			/* Increment iteration counter */
			k++;
			/* Set L2 error back to zero */
			l2 = 0.0;
		}
	} else if(!(strcmp((*problem).method, "Jacobi"))){
		grvy_printf(GRVY_DEBUG, "Using solver method Jacobi\n");
		/* Allocate array for holding k iteration values */
		double *Tk;
		Tk = malloc((*problem).size * sizeof(double));
		/* Fill with zeros */
		for(i=0; i<(*problem).size; i++){
			Tk[i] = 0.0;
		}
		while(converge >= (*problem).tol){
			for(i=0; i<(*problem).size; i++){
				sigma = 0.0;
				for(j=0; j<(*problem).size; j++){
					if(i != j){
						sigma += (*problem).A[i*(*problem).size + j]*Tk[j];
					}
				}
				/* Grab amount delta to ith entry */
				delta = 1.0/((*problem).A[i*(*problem).size + i]) * ((*problem).Q[i] - sigma);
				/* Add it to iterate solution */
				(*problem).T[i] = delta;
				/* Compute the sum of change squared for L2 error */
				l2 += delta * delta;	
			}
			/* Compute final l2 error */
			l2 = sqrt(1.0/((*problem).size) * l2);
			grvy_printf(GRVY_DEBUG, "Iteration: %d\t Iteration Norm: %f\n", k, l2);
			if (l2 <= (*problem).tol){
				grvy_printf(GRVY_INFO, "Converged at iteration: %d", k);
				free(Tk);
				break;
			}
			/* Check for max iterations and increment */
			if (k >= (*problem).iter_max){
				grvy_printf(GRVY_INFO, "Maximum number of iterations (%d) reached. Iteration Norm: %f", (*problem).iter_max, l2);
				free(Tk);
				return 500;
			}
			/* Increment iteration counter */
			k++;
			/* Set L2 error back to zero */
			l2 = 0.0;
			/* Pass values of Tk to T so it becomes the k-1 iteration */
			for(i=0; i<(*problem).size; i++){
				Tk[i] = (*problem).T[i];
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

double errornorm(Heat* problem)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Computing L2 error norm against MASA soln...\n");
	
	/* Mesh spacing calculations */
	double h;
	h = 1.0/((*problem).n-1) * (1.0 - 0.0);
	
	/* Iterators */
	int i;
	int j;
	/* Error calcs */
	double l2;
	l2 = 0.0;
	double delta;
	
	if((*problem).dimension==1){
		masa_init("error", "heateq_1d_steady_const");
		masa_set_param("A_x", (*problem).Ax);
		masa_set_param("k_0", (*problem).kappa);
		for(i=0; i<(*problem).size; i++){
			if((*problem).order == 2){
				delta = (*problem).T[i] - masa_eval_1d_exact_t((i+1)*h);
			} else {
				delta = (*problem).T[i] - masa_eval_1d_exact_t((i+2)*h);
			}
			l2 += delta * delta;
		}
	} else{
		masa_init("error", "heateq_2d_steady_const");
		masa_set_param("A_x", (*problem).Ax);
		masa_set_param("B_y", (*problem).By);
		masa_set_param("k_0", (*problem).kappa);
		for(i=0; i<(*problem).n-2; i++){
			for(j=0; j<(*problem).n-2; j++){
				/* not declaring a variable interior=n-2 for this case */
				delta = (*problem).T[i*((*problem).n-2) + j] - masa_eval_2d_exact_t((j+1)*h, (i+1)*h);
				l2 += delta * delta;
			}
		}
	}
	
	l2 = sqrt(1.0/(*problem).size * l2);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Computing L2 error norm against MASA soln...\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	return l2;
}

/*--------------------------------------------------------------------------------------
 *  hdf5_output: Routine for writing solution to an hdf5 file
 *-------------------------------------------------------------------------------------*/

int hdf5_output(Heat* problem)
{

	/* Timing */
	grvy_timer_begin(__func__);
	
	/* Compute L2 error if necessary */
	double errnorm;
	if((*problem).verification){
		errnorm = errornorm(problem);
		grvy_printf(GRVY_INFO, "L2 error norm: %f\n", errnorm);
	}
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Writing results to file...\n");
	grvy_printf(GRVY_INFO, "Writing to...%-10s\n", (*problem).output_file);
	
	
	/* Setup hdf5 file identifier*/
	hid_t file;
	herr_t status;
	
	/* Write settings */
	file = H5Fcreate((*problem).output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hsize_t dimsints[1];
	dimsints[0] = 1;

	hid_t dsp_ints; /* dataspace definitions */
	hid_t dst_dimension, dst_order, dst_nodes; /* dataset definitions */

	dsp_ints = H5Screate_simple(1, dimsints, NULL);

	dst_dimension = H5Dcreate2(file, "dimension", H5T_NATIVE_INT, dsp_ints, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dst_dimension, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*problem).dimension);
	status = H5Dclose(dst_dimension);

	dst_order = H5Dcreate2(file, "order", H5T_NATIVE_INT, dsp_ints, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dst_order, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*problem).order);
	status = H5Dclose(dst_order);

	dst_nodes = H5Dcreate2(file, "nodes", H5T_NATIVE_INT, dsp_ints, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dst_nodes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*problem).n);
	status = H5Dclose(dst_nodes);


	status = H5Sclose(dsp_ints);

	/* Write data */
	hsize_t dimsdata[2];
	hid_t dsp_datax, dsp_datay, dsp_datat, dsp_datan;
	hid_t dst_datax, dst_datay, dst_datat, dst_datan;
	int rows;
	int columns;
	if((*problem).dimension == 1){
		rows = (*problem).n;
		columns = 1 ;
		dimsdata[0] = (*problem).n;
		dimsdata[1] = 1;
	} else if((*problem).dimension == 2){
		rows = (*problem).n * (*problem).n;
		columns = 1 ;
		dimsdata[0] = (*problem).n * (*problem).n;
		dimsdata[1] = 1;
	}
	
	double datax[rows][columns];
	double datay[rows][columns];
	double datat[rows][columns];
	int datan[rows][columns];

	if((*problem).dimension==1){
		dsp_datax = H5Screate_simple(2, dimsdata, NULL);
		dsp_datat = H5Screate_simple(2, dimsdata, NULL);
		dsp_datan = H5Screate_simple(2, dimsdata, NULL);
		dst_datax = H5Dcreate2(file, "X", H5T_NATIVE_DOUBLE, dsp_datax, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dst_datat = H5Dcreate2(file, "T", H5T_NATIVE_DOUBLE, dsp_datat, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dst_datan = H5Dcreate2(file, "Node", H5T_NATIVE_INT, dsp_datan, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		int i;
		double xcoord;
		double val;
		if((*problem).verification){
			masa_init("verify", "heateq_1d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
		}
		for(i=0; i<(*problem).n; i++){
			xcoord =(1.0 * i) / ((*problem).n - 1);
			if(i==0){ /* Left Node */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(0.0);
				} else {
					val = (*problem).Txmin;
				}
			} else if(i==(*problem).n-1){ /* Right node */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(1.0);
				} else {
					val = (*problem).Txmax;
				}
			} else if((*problem).order == 4 && i==1){ /* Second from left, 4th order */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(xcoord);
				} else {
					val = (*problem).Txmin;
				}
			} else if((*problem).order == 4 && i==(*problem).n-2){ /* Soecond from right, 4th order */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(xcoord);
				} else {
					val = (*problem).Txmax;
				}
			} else {
				if((*problem).order == 2){
					val = (*problem).T[i-1];
				} else if((*problem).order == 4){
					val = (*problem).T[i-2];
				} else {
					return 905;
				}
			}
			datan[i][0] = i;
			datax[i][0] = xcoord;
			datat[i][0] = val;
		}

		status = H5Dwrite(dst_datax, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datax);
		status = H5Dclose(dst_datax);
		status = H5Sclose(dsp_datax);

		status = H5Dwrite(dst_datat, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datat);
		status = H5Dclose(dst_datat);
		status = H5Sclose(dsp_datat);

		status = H5Dwrite(dst_datan, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, datan);
		status = H5Dclose(dst_datan);
		status = H5Sclose(dsp_datan);
	} else if((*problem).dimension==2){
		dsp_datax = H5Screate_simple(2, dimsdata, NULL);
		dsp_datay = H5Screate_simple(2, dimsdata, NULL);
		dsp_datat = H5Screate_simple(2, dimsdata, NULL);
		dsp_datan = H5Screate_simple(2, dimsdata, NULL);
		dst_datax = H5Dcreate2(file, "X", H5T_NATIVE_DOUBLE, dsp_datax, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dst_datay = H5Dcreate2(file, "Y", H5T_NATIVE_DOUBLE, dsp_datay, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dst_datat = H5Dcreate2(file, "T", H5T_NATIVE_DOUBLE, dsp_datat, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dst_datan = H5Dcreate2(file, "Node", H5T_NATIVE_INT, dsp_datan, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		int i;
		int j;
		double xcoord;
		double ycoord;
		double val;
		int index;
		index = 0;
		int node;
		if((*problem).verification){
			masa_init("verify", "heateq_2d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
			masa_set_param("B_y", (*problem).By);
		}
		for(i=0; i<(*problem).n; i++){
			ycoord = (1.0 * i)/((*problem).n - 1);
			for(j=0; j<(*problem).n; j++){
				node = i*(*problem).n + j;
				xcoord = (1.0 * j)/((*problem).n - 1);
				if(j==0 && i==0){ /* Bottom left */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(0.0, 0.0);
					} else {
						val = (*problem).Txmin + (*problem).Tymin;
					}
				} else if (j==(*problem).n-1 && i==(*problem).n-1){ /* Top right */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(1.0, 1.0);
					} else {
						val = (*problem).Txmax + (*problem).Tymax;
					}
				} else if(j==(*problem).n-1 && i==0){ /* Bottom right */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(1.0, 0.0);
					} else {
						val = (*problem).Txmax + (*problem).Tymin;
					}
				} else if(j==0 && i==(*problem).n-1){ /* Top left */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(0.0, 1.0);
					} else {
						val = (*problem).Txmin + (*problem).Tymax;
					}	
				} else if(i==0 && j!=0 && j!=(*problem).n-1){ /* Bottom row */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(xcoord, 0.0);
					} else {
						val = (*problem).Tymin;
					}
				} else if(i==(*problem).n-1 && j!=0 && j!=(*problem).n-1){ /* Top row */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(xcoord, 1.0);
					} else {
						val = (*problem).Tymax;
					}
				} else if(j==0 && i!=0 && i!=(*problem).n-1){ /* Left column */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(0.0, ycoord);
					} else {
						val = (*problem).Txmin;
					}
				} else if(j==(*problem).n-1 && i!=0 && i!=(*problem).n-1){ /* Right column */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(1.0, ycoord);
					} else {
						val = (*problem).Txmax;	
					}
				} else {
					index++;
					val = (*problem).T[index];
				}
				datan[i*(*problem).n + j][0] = node;
				datax[i*(*problem).n + j][1] = xcoord;
				datay[i*(*problem).n + j][2] = ycoord;
				datat[i*(*problem).n + j][3] = val;
			}
		}
		status = H5Dwrite(dst_datax, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datax);
		status = H5Dclose(dst_datax);
		status = H5Sclose(dsp_datax);

		status = H5Dwrite(dst_datay, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datay);
		status = H5Dclose(dst_datay);
		status = H5Sclose(dsp_datay);


		status = H5Dwrite(dst_datat, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datat);
		status = H5Dclose(dst_datat);
		status = H5Sclose(dsp_datat);

		status = H5Dwrite(dst_datan, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, datan);
		status = H5Dclose(dst_datan);
		status = H5Sclose(dsp_datan);
	}
	
	
	/* Write error */
	hsize_t dimserr[1];
	hid_t dsp_err;
	hid_t dst_err;
	dimserr[0] = 1;
	if((*problem).verification){
		dsp_err = H5Screate_simple(1, dimserr, NULL);
		dst_err = H5Dcreate(file, "L2 Error Norm", H5T_NATIVE_DOUBLE, dsp_err, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dst_err, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &errnorm);
		status = H5Dclose(dst_err);
		status = H5Sclose(dsp_err);
	}

	status = H5Fclose(file);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Writing results to file...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);

	return 0;

}

/*--------------------------------------------------------------------------------------
 *  output: Routine for writing solution to an ascii file
 *-------------------------------------------------------------------------------------*/

int output(Heat* problem, char* inputfile)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Writing results to file...\n");
	grvy_printf(GRVY_INFO, "Writing to...%-10s\n", (*problem).output_file);
	
	
	/*Dump settings */
	int erropen;
	erropen = grvy_input_fopen(inputfile);
	if(!(erropen)){
		grvy_printf(GRVY_DEBUG, "Unable to open file %-10s\n", *inputfile);
		return 401;
	} 

	grvy_input_fdump_file("%", (*problem).output_file);
	
	grvy_input_fclose();
	
	FILE *fp;
	fp = fopen((*problem).output_file, "a");
	if(fp==NULL){
		return 400;
	}
	

	/* Write values */
	if((*problem).dimension==1){
		fprintf(fp, "Node:\t X Coordinate:\t Temp:\n");
		int i;
		double xcoord;
		double val;
		if((*problem).verification){
			masa_init("verify", "heateq_1d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
		}
		for(i=0; i<(*problem).n; i++){
			xcoord =(1.0 * i) / ((*problem).n - 1);
			if(i==0){ /* Left Node */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(0.0);
				} else {
					val = (*problem).Txmin;
				}
			} else if(i==(*problem).n-1){ /* Right node */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(1.0);
				} else {
					val = (*problem).Txmax;
				}
			} else if((*problem).order == 4 && i==1){ /* Second from left, 4th order */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(xcoord);
				} else {
					val = (*problem).Txmin;
				}
			} else if((*problem).order == 4 && i==(*problem).n-2){ /* Soecond from right, 4th order */
				if((*problem).verification){
					val = masa_eval_1d_exact_t(xcoord);
				} else {
					val = (*problem).Txmax;
				}
			} else {
				if((*problem).order == 2){
					val = (*problem).T[i-1];
				} else if((*problem).order == 4){
					val = (*problem).T[i-2];
				} else {
					return 905;
				}
			}
			fprintf(fp, "%d\t %f\t %f\n", i, xcoord, val);
		}
	} else if((*problem).dimension==2){
		fprintf(fp, "Node:\t X Coordinate:\t Y Coordinate:\t Temp:\n");
		int i;
		int j;
		double xcoord;
		double ycoord;
		double val;
		int index;
		index = 0;
		int node;
		if((*problem).verification){
			masa_init("verify", "heateq_2d_steady_const");
			masa_set_param("k_0", (*problem).kappa);
			masa_set_param("A_x", (*problem).Ax);
			masa_set_param("B_y", (*problem).By);
		}
		for(i=0; i<(*problem).n; i++){
			ycoord = (1.0 * i)/((*problem).n - 1);
			for(j=0; j<(*problem).n; j++){
				node = i*(*problem).n + j;
				xcoord = (1.0 * j)/((*problem).n - 1);
				if(j==0 && i==0){ /* Bottom left */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(0.0, 0.0);
					} else {
						val = (*problem).Txmin + (*problem).Tymin;
					}
				} else if (j==(*problem).n-1 && i==(*problem).n-1){ /* Top right */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(1.0, 1.0);
					} else {
						val = (*problem).Txmax + (*problem).Tymax;
					}
				} else if(j==(*problem).n-1 && i==0){ /* Bottom right */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(1.0, 0.0);
					} else {
						val = (*problem).Txmax + (*problem).Tymin;
					}
				} else if(j==0 && i==(*problem).n-1){ /* Top left */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(0.0, 1.0);
					} else {
						val = (*problem).Txmin + (*problem).Tymax;
					}	
				} else if(i==0 && j!=0 && j!=(*problem).n-1){ /* Bottom row */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(xcoord, 0.0);
					} else {
						val = (*problem).Tymin;
					}
				} else if(i==(*problem).n-1 && j!=0 && j!=(*problem).n-1){ /* Top row */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(xcoord, 1.0);
					} else {
						val = (*problem).Tymax;
					}
				} else if(j==0 && i!=0 && i!=(*problem).n-1){ /* Left column */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(0.0, ycoord);
					} else {
						val = (*problem).Txmin;
					}
				} else if(j==(*problem).n-1 && i!=0 && i!=(*problem).n-1){ /* Right column */
					if((*problem).verification){
						val = masa_eval_2d_exact_t(1.0, ycoord);
					} else {
						val = (*problem).Txmax;	
					}
				} else {
					index++;
					val = (*problem).T[index];
				}
				fprintf(fp, "%d\t %f\t %f\t %f\n", node, xcoord, ycoord, val);
			}
		}
	}
	
	/* Compute L2 error if necessary */
	if((*problem).verification){
		double errnorm;
		errnorm = errornorm(problem);
		grvy_printf(GRVY_INFO, "L2 error norm: %f\n", errnorm);
		fprintf(fp, "\n L2 Error Norm: %E\n", errnorm);
	}
	
	/* Close and free */
	fclose(fp);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Writing results to file...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
	
	return 0;
}

/*--------------------------------------------------------------------------------------
 *  cleanup: Function for cleaning up and freeing dynamically allocated memory
 *-------------------------------------------------------------------------------------*/
 
void cleanup(Heat* problem)
{
	/* Timing */
	grvy_timer_begin(__func__);
	/* Logging */
	grvy_printf(GRVY_INFO, "Freeing variables...\n");

	free((*problem).A);
	free((*problem).T);
	free((*problem).Q);
	free((*problem).method);
	free((*problem).output_file);
	
	/* Logging */
	grvy_printf(GRVY_INFO, "Freeing variables...DONE\n");
	/* Timing */
	grvy_timer_end(__func__);
	return;
}
