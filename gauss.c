/* Implements Gauss-Seidel iterative algorithm for
 * solving linear equations
 * Implements SOR iterative algorithm for solving
 * linear equations
 * Lab Assignment 06, PB09203226
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXREPT	409600
void lsolve_gauss(double *pA, double *Y, double *X,
		  int n, double epsilon, int *pnstep)
{
	int i, j, k;
	double sum, norm_inf;
	double (*A)[n] = (double (*)[n]) pA;
	double *old = malloc(n * sizeof(*old));

	for (i = 0; i < MAXREPT; i++) {
		for (j = 0; j < n; j++) {
			old[j] = X[j];
		}
		for (j = 0; j < n; j++) {
			sum = Y[j] / A[j][j];
			for (k = 0; k < n; k++) {
				if (k != j) {
					sum -= A[j][k] * X[k] / A[j][j];
				}
			}
			X[j] = sum;
		}
		norm_inf = 0;
		for (j = 0; j < n; j++) {
			if (fabs(X[j] - old[j]) > norm_inf)
				norm_inf = fabs(X[j] -old[j]);
		}
		if (norm_inf < epsilon)
			break;
	}

	*pnstep = i+1;
	free(old);
}

#define EPSILON	1e-13
int main(void)
{
	double A[][9] = {{31, -13, 0, 0, 0, -10, 0, 0, 0},
			  {-13, 35, -9, 0, -11, 0, 0, 0, 0},
			  {0, -9, 31, -10, 0, 0, 0, 0, 0},
			  {0, 0, -10, 79, -30, 0, 0, 0, -9},
			  {0, 0, 0, -30, 57, -7, 0, -5, 0},
			  {0, 0, 0, 0, -7, 47, -30, 0, 0},
			  {0, 0, 0, 0, 0, -30, 41, 0, 0},
			  {0, 0, 0, 0, -5, 0, 0, 27, -2},
			  {0, 0, 0, -9, 0, 0, 0, -2, 29}};
	double Y[9] = {-15, 27, -23, 0, -20, 12, -7, 7, 10};
	double X[9] = {1,1,1,1,1,1,1,1,1};
	int i, nstep;

	lsolve_gauss(A, Y, X, sizeof(X)/sizeof(*X), EPSILON, &nstep);
	printf("Roots = \n");
	for (i = 0; i < sizeof(X)/sizeof(*X); i++)
		printf("%.13f\n", X[i]);
	printf("Steps = %d\n", nstep);

	return 0;
}
