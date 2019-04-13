/*--------------------------------------------------------------------------------------
 * Program for solving the 1D and 2D heat equation with 2nd and 4th order finite difference iterative solvers.
 *
 * Sam Potter
 *
 *--------------------------------------------------------------------------------------
 *  utilities.h: header file for functions defined in utilities.c
 *-------------------------------------------------------------------------------------*/
 
 #include <stdio.h>
 #include <stdlib.h>
 #include <grvy.h>
 #include <masa.h>
 #include <math.h>
 #include <string.h>
 #include <hdf5.h>
 #include <config.h>
 #ifdef INCUDE_PETSC
 	#include <petsc.h>
 #endif
 
 #ifndef UTILITIES_H
 #define UTILITIES_H

 /* Struct definition */

 typedef struct
 {
	char *method;
	int iter_max;
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
	double By;
	int verification;
	int out_mode;
	char *output_file;
	double *A;
	double *T;
	double *Q;
	int size;

 } Heat;
 
 /* Function prototypes */
 int input(Heat* problem, char* inputfile);
 int sanitize(Heat* problem);
 int petsc(Heat* problem);
 int assemble(Heat* problem);
 int solver(Heat* problem);
 double errornorm(Heat* problem);
 int hdf5_output(Heat* problem);
 int output(Heat* problem, char* inputfile);
 void cleanup(Heat* problem);
 #endif
