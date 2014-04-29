/* Implements numerical composite trapezoidal and Simpson integration
 * in both 1D and 2D
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#undef NDEBUG
#include <assert.h>

/**
 * generate_sample:
 * @f: pointer to sampling function
 * @a: left boundary of sampling interval
 * @h: step
 * @n: number of partition intervals
 * @v: array of f[(a, cutting points, b)]
 *
 * Sample a function on [a, a + n * h]
 * at points x[i] = a + i * h, i = 0..n.
 *
 * Modifies v.	No return value.
 */
void generate_sample(double (*f)(double),
		     double a, double h, int n, double v[n+1])
{
	int i;

	assert(v != NULL);
	for (i = 0; i <= n; i++)
		v[i] = f(a + h * i);
}

/**
 * nintegrate_trapezodial:
 * @f : pointer to integrand
 * @a : left boundary of integrating interval
 * @b : right boundary of integrating interval
 * @n : number of partition intervals
 *
 * Performs numerical integration with
 * composite trapezodial algorithm.
 *
 * Returns: Result of integration
 */
double nintegrate_trapezodial(double (*f)(double),
			      double a, double b, int n)
{
	int i;
	double *v = malloc((n+1) * sizeof(*v));
	double sum = 0.0;
	double h = (b - a) / n;

	assert(f != NULL);
	assert(v != NULL);
	generate_sample(f, a, h, n, v);

	for (i = 1; i < n; i++)
		sum += v[i];
	sum += v[0] / 2.0;
	sum += v[n] / 2.0;

	free(v);

	return sum * h;
}

/**
 * nintegrate_simpson:
 * @f : pointer to integrand
 * @a : left boundary of integrating interval
 * @b : right boundary of integrating interval
 * @n : number of partiton intervals
 *
 * Performs numerical integration with
 * composite Simpson algorithm.
 *
 * Returns: Result of integration
 */
double nintegrate_simpson(double (*f)(double),
			  double a, double b, int n)
{
	int i;
	double *v = malloc((n+1) * sizeof(*v));
	double sum_odd = 0.0, sum_even = 0.0;
	double res;
	double h = (b - a) / n;

	assert(f != NULL);
	assert(v != NULL);
	generate_sample(f, a, h, n, v);

	for (i = 1; i < n; i += 2)
		sum_odd += v[i];
	for (i = 2; i < n; i += 2)
		sum_even += v[i];
	res = (v[0] + v[n] + 4.0*sum_odd + 2.0*sum_even) * h / 3.0;

	free(v);

	return res;
}

/**
 * generate_sample_2d:
 * @f: pointer to sampling function
 * @a: left x boundary of sampling interval
 * @h: x step
 * @n: number of partition intervals in x direction
 * @b: left y boundary of sampling interval
 * @k: y step
 * @m: number of partition intervals in y direction
 * @v: array of f[(a, cutting points, c), (b, cutting points, d)]
 *
 * Sample a function on [a, a + n * h] * [b, b + m * k]
 * at points (x_i, y_j) = (a + i * h, b + j * k), i = 0..n, j = 0..m
 * Store result in (n+1)*(m+1) array v with following layout:
 *	 y_0   y_j   y_m
 *  x_0	[		]
 *	[		]
 *	[		]
 *  x_i	[    v[i][j]	]
 *	[		]
 *	[		]
 *  x_n	[		]
 *
 * Modifies v.	No return value.
 */
void generate_sample_2d(double (*f)(double, double),
			double a, double h, int n,
			double b, double k, int m,
			double v[n+1][m+1])
{
	int i, j;

	assert(v != NULL);
	for (i = 0; i <= n; i++)
		for (j = 0; j <= m; j++)
			v[i][j] = f(a + h * i, b + k * j);
}

/**
 * nintegrate_trapezodial_2d:
 * @f : pointer to integrand
 * @a : left x boundary of integrating interval
 * @b : right x boundary of integrating interval
 * @n : number of partition intervals in x direction
 * @c : left y boundary of integrating interval
 * @d : right y boundary of integrating interval
 * @m : number of partition intervals in y direction
 *
 * Performs numerical integration with
 * composite trapezodial algorithm.
 *
 * Returns: Result of integration
 */
double nintegrate_trapezodial_2d(double (*f)(double, double),
				 double a, double b, int n,
				 double c, double d, int m)
{
	int i, j;
	double (*v)[m+1] = malloc((n+1) * (m+1) * sizeof(**v));
	double sum = 0.0, bsum = 0.0;
	double h = (b - a) / n, k = (d - c) / m;

	assert(f != NULL);
	assert(v != NULL);
	generate_sample_2d(f, a, h, n, c, k, m, v);

	/* interior */
	for (i = 1; i < n; i++)
		for (j = 1; j < m; j++)
			sum += v[i][j];
	/* boundary */
	for (i = 1; i < n; i++)
		bsum += v[i][0];
	for (i = 1; i < n; i++)
		bsum += v[i][m];
	for (j = 1; j < m; j++)
		bsum += v[0][j];
	for (j = 1; j < m; j++)
		bsum += v[n][j];
	sum += 0.5 * bsum;
	/* 4 vertexes */
	sum += (v[0][0] + v[n][0] + v[0][m] + v[n][m]) * 0.25;

	free(v);

	return sum * h * k;
}

/**
 * nintegrate_simpson_2d:
 * @n : number of partiton intervals
 * @f : pointer to integrand
 * @a : left boundary of integrating interval
 * @b : right boundary of integrating interval
 *
 * Performs numerical integration with
 * composite Simpson algorithm.
 *
 * Returns: Result of integration
 */
double nintegrate_simpson_2d(double (*f)(double, double),
			     double a, double b, int n,
			     double c, double d, int m)
{
	int i, j;
	double (*v)[m+1] = malloc((n+1) * (m+1) * sizeof(**v));
	double sum_odd_odd = 0.0, sum_odd_even = 0.0;
	double sum_even_odd = 0.0, sum_even_even = 0.0;
	double bsum_odd = 0.0, bsum_even = 0.0;
	double res;
	double h = (b - a) / n, k = (d - c) / m;

	assert(f != NULL);
	assert(v != NULL);
	generate_sample_2d(f, a, h, n, c, k, m, v);

	/* interior */
	for (i = 1; i < n; i += 2) {
		for (j = 1; j < m; j += 2)
			sum_odd_odd += v[i][j];
		for (j = 2; j < m; j += 2)
			sum_odd_even += v[i][j];
	}
	for (i = 2; i < n; i += 2) {
		for (j = 1; j < m; j += 2)
			sum_even_odd += v[i][j];
		for (j = 2; j < m; j += 2)
			sum_even_even += v[i][j];
	}

	/* boundary */
	for (i = 1; i < n; i += 2)
		bsum_odd += v[i][0] + v[i][m];
	for (i = 2; i < n; i += 2)
		bsum_even += v[i][0] + v[i][m];
	for (j = 1; j < m; j += 2)
		bsum_odd += v[0][j] + v[n][j];
	for (j = 2; j < m; j += 2)
		bsum_even += v[0][j] + v[n][j];

	res = (v[0][0] + v[0][m] + v[n][0] + v[n][m]
	       + 16.0 * sum_odd_odd + 8.0 * sum_odd_even
	       + 8.0 * sum_even_odd + 4.0 * sum_even_even
	       + 4.0 * bsum_odd + 2.0 * bsum_even) * h * k / 9.0;

	free(v);
	return res;
}

double f(double x, double y)
{
	return 1.0 / (x + y);
}

int main(void)
{
	int i;
	double true_value;
	double nres;

	true_value = cos(1.0) - cos(5.0);
	printf("Composite trapezoidal integration:\n");
	for (i = 1; i <= 1 << 12; i *= 2) {
		nres = nintegrate_trapezodial(sin, 1.0, 5.0, i);
		printf("n = %d, result = %.13f, error = %.13f\n",
		       i, nres, fabs(nres - true_value));
	}

	printf("\nComposite Simpson integration:\n");
	for (i = 2; i <= 1 << 12; i *= 2) {
		nres = nintegrate_simpson(sin, 1.0, 5.0, i);
		printf("n = %d, result = %.13f, error = %.13f\n",
		       i, nres, fabs(nres - true_value));
	}

	printf("\n2D Composite trapezodial integration:\n");
	true_value = 2.911031660324;
	for (i = 1; i <= 1 << 12; i *= 2) {
		nres = nintegrate_trapezodial_2d(f, 1, 5, i, 1, 5, i);
		printf("n = %d, result = %.13f, error = %.13f\n",
		       i, nres, fabs(nres - true_value));
	}

	printf("\n2D Composite Simpson integration:\n");
	for (i = 1; i <= 1 << 12; i *= 2) {
		nres = nintegrate_simpson_2d(f, 1, 5, i, 1, 5, i);
		printf("n = %d, result = %.13f, error = %.13f\n",
		       i, nres, fabs(nres - true_value));
	}

	return 0;
}
