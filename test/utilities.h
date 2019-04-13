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
 
 #ifndef UTILITIES_H
 #define UTILITIES_H
 
 /* Function prototypes */
 int get_settings(int *out_mode, int *timer);
 int input(char *solver, int *iter_max, double *tol, int *order, int *dimension, double *kappa, int *n, double *Txmin, double *Tymin, double *Txmax, double *Tymax, double *Ax, double *By, int *verification, char *output_file);
 int cl_parse(int argc, char *argv[], int *order, int *dimension, int *verification, int *out_mode, int *timer, char *output_file);
 int assemble(double *A, double *T, double *Q, int dimension, int order, int n, int size, double kappa, double Txmin, double Txmax, double Tymin, double Tymax, int verification, double Ax, double By);
 int solver(double *A, double *T, double *Q, int size, int tol, int itermax, char method);
 double errornorm(double *T, int size, int n, int dimension, double kappa, double Ax, double By);
 int output(double *T, int n, int size, int dimension, double kappa, double Ax, double By, char output_file, int verification);
 void cleanup(double *A, double *T, double *Q, char *method, char *output_file);
 #endif