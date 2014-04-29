/* Implements numerical composite trapezoidal and Simpson integration
 * Lab experiment 03, PB09203226
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * generate_sample:
 * @n: number of partition intervals
 * @v: array of f[(a, cutting points, b)]
 * @f: pointer to sampling function
 * @a: left boundary of sampling interval
 * @h: step
 *
 * Sample a function on [a, a + n * h]
 * at points x[i] = a + i * h, i = 0..n.
 *
 * Modifies v.	No return value.
 */
void generate_sample(int n, double v[], double (*f)(double),
		     double a, double h)
{
	int i;

	for (i = 0; i <= n; i++)
		v[i] = f(a + h * i);
}

/**
 * nintegrate_trapezodial:
 * @n : number of partition intervals
 * @f : pointer to integrand
 * @a : left boundary of integrating interval
 * @b : right boundary of integrating interval
 *
 * Performs numerical integration with
 * composite trapezodial algorithm.
 *
 * Returns: Result of integration
 */
double nintegrate_trapezodial(int n, double (*f)(double),
			      double a, double b)
{
	int i;
	double *v = malloc((n+1) * sizeof(*v));
	double sum = 0.0;
	double h = (b - a) / n;

	generate_sample(n, v, f, a, h);

	for (i = 1; i < n; i++)
		sum += v[i];
	sum += v[0] / 2.0;
	sum += v[n] / 2.0;

	free(v);

	return sum * h;
}

/**
 * nintegrate_simpson:
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
double nintegrate_simpson(int n, double (*f)(double),
			  double a, double b)
{
	int i;
	double *v = malloc((n+1) * sizeof(*v));
	double sum_odd = 0.0, sum_even = 0.0;
	double res;
	double h = (b - a) / n;

	generate_sample(n, v, f, a, h);

	for (i = 1; i < n; i += 2)
		sum_odd += v[i];
	for (i = 2; i < n; i += 2)
		sum_even += v[i];
	res = (v[0] + v[n] + 4.0*sum_odd + 2.0*sum_even) * h / 3.0;

	free(v);

	return res;
}

int main(void)
{
	int i;
	double true_value = cos(1.0) - cos(5.0);
	double nres;

	printf("Composite trapezoidal integration:\n");
	for (i = 1; i <= 1 << 12; i *= 2) {
		nres = nintegrate_trapezodial(i, sin, 1.0, 5.0);
		printf("n = %d, result = %.13f, error = %.13f\n",
		       i, nres, fabs(nres - true_value));
	}

	printf("\nComposite Simpson integration:\n");
	for (i = 2; i <= 1 << 12; i *= 2) {
		nres = nintegrate_simpson(i, sin, 1.0, 5.0);
		printf("n = %d, result = %.13f, error = %.13f\n",
		       i, nres, fabs(nres - true_value));
	}

	return 0;
}
