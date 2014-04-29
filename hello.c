/* Numerical Methods, Lab Exercise 01
 * PB09203226 */

#include <stdio.h>

/* calculates $\Psi(x) = \sum^\infty_{n=1} \frac{1}{n(n+x)}$
 * returns result
 * end-of-loop criteria: $$
 * \left| \left( \sum^\infty_{n=1} - \sum^N_{n=1} \right)
 *    \frac{1}{n(n+x)} \right| < \sum^\infty_{n=N+1} \frac{1}{n^2}
 *    < \sum^\infty_{n=N+1} \frac{1}{(n-1)n}
 *    = \frac{1}{N} < \epsilon
 * $$
 */
#define EPSILON		1e-6
long double get_psi(double x)
{
	int n;
	long double sum = 0.0;

	for (n = 1; n <= 1/EPSILON; n++)
		sum += 1.0 / (n * (n + x));

	return sum;
}

int main(void)
{
	float x;

	for (x = 0.0; x < 1.01; x += 0.1) /* x <= 1.0 won't work (fp issue) */
		printf("%.2f\t%.13Lf\n", x, get_psi(x));
	for (x = 10.0; x <= 300.0; x += 10.0)
		printf("%.2f\t%.13Lf\n", x, get_psi(x));

	return 0;
}
