/* Implements column-major Gauss elimination for solving systems
 * of linear equations
 */

#undef NDEBUG
#include <assert.h>

#include <stdio.h>
#include <math.h>

/**
 * Solves X in AX = Y
 * @pA : (A|Y)
 * @X : X
 * Modifies X[]
 */
void lsolve_colmaj(double *pA, double *X, int n)
{
	assert(pA != NULL && X != NULL && n > 1);

	double (*A)[n+1] = (double (*)[n+1]) pA;
	double tmp, factor;
	int i, j, k, max;

	for (i = 0; i < n; i++) {
		max = i;	/* locate major row */
		for (j = i+1; j < n; j++) {
			if (fabs(A[j][i]) > fabs(A[max][i]))
				max = j;
		}
		if (i != max) { /* exchange i-th row and major row */
			for (j = i; j < n; j++) {
				tmp = A[i][j];
				A[i][j] = A[max][j];
				A[max][j] = tmp;
			}
		}
		for (j = i+1; j < n; j++) { /* elimination */
			factor = - A[j][i] / A[i][i];
			for (k = i; k < n+1; k++)
				A[j][k] += factor * A[i][k];
		}
	}
	for (i = n-1; i >= 0; i--) {
		X[i] = A[i][n] / A[i][i];
		for (j = i-1; j >= 0; j--)
			A[j][n] -= X[i] * A[j][i];
	}
}

int main(void)
{
	double A[][10] = {{31, -13, 0, 0, 0, -10, 0, 0, 0, -15},
			   {-13, 35, -9, 0, -11, 0, 0, 0, 0, 27},
			   {0, -9, 31, -10, 0, 0, 0, 0, 0, -23},
			   {0, 0, -10, 79, -30, 0, 0, 0, -9, 0},
			   {0, 0, 0, -30, 57, -7, 0, -5, 0, -20},
			   {0, 0, 0, 0, -7, 47, -30, 0, 0, 12},
			   {0, 0, 0, 0, 0, -30, 41, 0, 0, -7},
			   {0, 0, 0, 0, -5, 0, 0, 27, -2, 7},
			   {0, 0, 0, -9, 0, 0, 0, -2, 29, 10}};
	double X[9];
	lsolve_colmaj((double *)A, X, 9);

	printf("Roots = \n");
	int i;
	for (i = 0; i < 9; i++)
		printf("%.13f\n", X[i]);

	return 0;
}
