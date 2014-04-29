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
void lsolve_sor(double *pA, double *Y, double *X, double omega,
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
			X[j] = (1 - omega) * old[j] + omega * sum;
		}
		norm_inf = 0;
		for (j = 0; j < n; j++) {
			if (fabs(X[j] - old[j]) > norm_inf)
				norm_inf = fabs(X[j] - old[j]);
		}
		if (norm_inf < epsilon)
			break;
	}

	*pnstep = i+1;
	free(old);
}

void lsolve_gauss(double *pA, double *Y, double *X,
		  int n, double epsilon, int *pnstep)
{
	lsolve_sor(pA, Y, X, 1, n, epsilon, pnstep);
}

#define arr_len(x)	(sizeof(x)/sizeof(*(x)))
#define EPSILON	1e-13

#define test_sor(omega, ns)		do {				\
		double __X[9];						\
		int __i;						\
		for (__i = 0; __i < arr_len(__X); __i++) {		\
			__X[__i] = 0;					\
		}							\
		lsolve_sor((double *) A, Y, __X, omega,			\
			   arr_len(__X), EPSILON, &ns);			\
		printf("Omega = %.2f, Steps = %d\n", omega, ns);	\
	} while (0)

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
	double X[9] = {0};
	int i, nstep, min_nstep = MAXREPT;
	double min_omega;

	printf("Gauss-Seidel Iteration:\n");
	lsolve_gauss((double *) A, Y, X, arr_len(X), EPSILON, &nstep);
	printf("Roots = \n");
	for (i = 0; i < arr_len(X); i++)
		printf("%.13f\n", X[i]);
	printf("Steps = %d\n", nstep);

	printf("\nSOR Iteration:\n");
	for (i = 1; i <= 99; i++) {
		test_sor((double) i / 50.0, nstep);
		if (min_nstep > nstep) {
			min_nstep = nstep;
			min_omega = (double) i / 50.0;
		}
	}
	printf("Best omega = %.2f\n", min_omega);

	return 0;
}
