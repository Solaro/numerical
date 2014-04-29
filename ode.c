/* Lab Assignment 07
 * PB09203226
 */

#include <stdio.h>
#include <math.h>

#define trace(x, d)	do {			\
		printf(#x " = %" #d "\n", x);	\
	} while (0)

double ndsolve_runge(double (*f)(double, double), double a,
		     double b, double h, double initv)
{
	int i, n = (b - a) / h;
	double k1, k2, k3, k4;
	double x, y = initv;

	trace(n, d);
	for (i = 0; i < n; i++) {
		trace(i, d);
		x = a + i * h;
		trace(x, f);
		trace(y, f);
		k1 = f(x, y);
		k2 = f(x + h / 2, y + h * k1 / 2);
		k3 = f(x + h / 2, y + h * k2 / 2);
		k4 = f(x + h, y + h * k3);
		y += (k1 + 2*k2 + 2*k3 + k4) * h / 6;
	}

	return y;
}

double ndsolve_adams(double (*f)(double, double), double a,
		     double b, double h, double initv)
{
	int i, n = (b - a) / h;
	double xnp1, xn, xnm1, xnm2;
	double ynp1, yn, ynm1, ynm2;

	ynm2 = initv;
	ynm1 = ndsolve_runge(f, a, a + h, h, initv);
	yn = ndsolve_runge(f, a, a + 2 * h, h, initv);
	for (i = 0; i < n-2; i++) {
		xnp1 = a + (i + 3) * h;
		xn = a + (i + 2) * h;
		xnm1 = a + (i + 1) * h;
		xnm2 = a + i * h;
		ynp1 = yn + (23 * f(xn, yn) - 16 * f(xnm1, ynm1)
			     + 5 * f(xnm2, ynm2)) * h / 12;
		ynp1 = yn + (5 * f(xnp1, ynp1) + 8 * f(xn, yn)
			     - f(xnm1, ynm1)) * h / 12;
		ynm2 = ynm1;
		ynm1 = yn;
		yn = ynp1;
	}
	return ynp1;
}

double f(double x, double y)
{
	return - x * x * y * y;
}

double y(double x)
{
	return 3.0 / (1.0 + x * x * x);
}

#define test_ndsolve(h, method)		do {				\
		double __res = ndsolve_##method(&f, 0, 1.5, (h), 3);	\
		printf("Step = %f, Result = %.13f, Error = %.13f\n",	\
		       (h), __res, fabs(__res - y(1.5)));		\
	} while (0)
#define test_runge(h)	test_ndsolve((h), runge)
#define test_adams(h)	test_ndsolve((h), adams)

int main(void)
{
	puts("Runge-Kutta Method:");
	test_runge(0.1);
	test_runge(0.05);
	test_runge(0.025);
	test_runge(0.0125);
	puts("\nAdams Method:");
	test_adams(0.1);
	test_adams(0.05);
	test_adams(0.025);
	test_adams(0.0125);

	return 0;
}
