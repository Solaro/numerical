/* Implements Newton's iterative algorithm for solving non-linear equations
 * Lab Assignment 04, PB09203226
 */

#include <stdio.h>
#include <math.h>

#define MAXREPT	1024

double nsolve_newton(double (*f)(double), double (*fprime)(double),
		     double initv, double epsilon, int *nr_iter)
{
	int i;
	double res, prev = initv;

	for (i = 0; i < MAXREPT; i++) {
		res = prev - f(prev) / fprime(prev);
		if (fabs(res - prev) < epsilon) {
			*nr_iter = i+1;
			return res;
		}
		prev = res;
	}
	*nr_iter = MAXREPT;
	return nan("no root");
}

double nsolve_secant(double (*f)(double), double initv1,
		     double initv2, double epsilon, int *nr_iter)
{
	int i;
	double res = initv2, prev = initv1, oldres;

	for (i = 0; i < MAXREPT; i++) {
		oldres = res;
		res = res - f(res) * (res - prev) / (f(res) - f(prev));
		if (fabs(oldres - res) < epsilon) {
			*nr_iter = i+1;
			return res;
		}
		prev = oldres;
	}
	*nr_iter = MAXREPT;
	return nan("no root");
}

double f(double x)
{
	return x * x * x / 3.0 - x;
}

double fprime(double x)
{
	return x * x - 1.0;
}

#define EPSILON	1e-13
#define test_nsolve_newton(x)	do {				\
	int _nr_iter;						\
	double _res;						\
	_res = nsolve_newton(f, fprime, x, EPSILON, &_nr_iter);	\
	printf("Initial value = %.1f, Steps of"			\
	" iteration = %d, root = %.13f\n",			\
	       x, _nr_iter, _res);				\
	} while (0)
#define test_nsolve_secant(x1, x2)	do {			\
	int _nr_iter;						\
	double _res;						\
	_res = nsolve_secant(f, x1, x2, EPSILON, &_nr_iter);	\
	printf("Initial value = (%.1f, %.1f), Steps of"		\
	" iteration = %d, root = %.13f\n",			\
	       x1, x2, _nr_iter, _res);				\
	} while (0)

int main(void)
{
	printf("Newton's Method:\n");
	test_nsolve_newton(0.1);
	test_nsolve_newton(0.2);
	test_nsolve_newton(0.9);
	test_nsolve_newton(9.0);
	printf("Secant Method:\n");
	test_nsolve_secant(0.0, 0.1);
	test_nsolve_secant(0.1, 0.2);
	test_nsolve_secant(0.2, 0.9);
	test_nsolve_secant(8.0, 9.0);

	return 0;
}
