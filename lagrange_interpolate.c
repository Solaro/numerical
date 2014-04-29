/*
 * Lab Exercise 02, Sept. 10
 * Implements Lagrange Interpolate
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Do Lagrange interpolation */
long double lagrange_interpolate(
	long double x, int n, long double vx[], long double vy[])
{
	int i, j;
	long double res = 0.0;
	long double base;

	for (i = 0; i <= n; i++) {
		base = 1.0;
		for (j = 0; j <= n; j++) {
			if (j != i)
				base *= (x - vx[j]) / (vx[i] - vx[j]);
		}
		res += base * vy[i];
	}
	return res;
}

/* Generates grid for interpolation
 * fgen(i, n) gives x[i], fcn(x[i]) gives y[i]
 */
void generate_grid(int n, long double vx[], long double vy[],
		   long double (*fgen)(int, int),
		   long double (*fcn)(long double))
{
	int i;

	for (i = 0; i <= n; i++) {
		vx[i] = fgen(i, n);
		vy[i] = fcn(vx[i]);
	}
}

/* Computes and returns maximal error after
 * generating grid and doing interpolation
 * fgen(i, n) gives x[i], fcn(x[i]) gives y[i] (interpolating point)
 * fgeny(j) gives y[j] (sampling error)
 * On error, returns a negtive number
 */
long double get_max_error(int n, long double (*fgen)(int, int),
			  long double (*fcn)(long double),
			  int jmax, long double (*fgeny)(int))
{
	int j;
	long double max_err = 0.0;
	long double err = 0.0;
	long double y = 0.0;
	long double *vx = malloc((n+1) * sizeof(*vx));
	long double *vy = malloc((n+1) * sizeof(*vy));

	if (vx == NULL || vy == NULL)
		return -1.0;
	generate_grid(n, vx, vy, fgen, fcn);
	for (j = 0; j <= jmax; j++) {
		y = fgeny(j);
		err = fabs(fcn(y) - lagrange_interpolate(y, n, vx, vy));
		if (err > max_err)
			max_err = err;
	}
	free(vx);
	free(vy);
	return max_err;
}

long double fgen_1(int i, int n)
{
	return -5.0 + 10.0 * i / n;
}

#define PI	(atan(1.0) * 4.0)
long double fgen_2(int i, int n)
{
	return -5.0 * cos(PI * (2*i+1) / (2*n+2));
}

long double fcn(long double x)
{
	return 1.0 / (1.0 + x * x);
}

long double fgen_y(int i)
{
	return -5.0 + 0.1 * i;
}

int main(void)
{
	int i;

	for (i = 5; i <= 40; i *= 2) {
		printf("N = %d\nMax Error of grid (1): %.13Lf\n"
		       "Max Error of grid (2): %.13Lf\n", i,
		       get_max_error(i, fgen_1, fcn, 100, fgen_y),
		       get_max_error(i, fgen_2, fcn, 100, fgen_y));
	}

	return 0;
}
