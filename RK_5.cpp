#include <iostream>
#include "RK_5.h"

using namespace std;

/* 5th-Order Explicit Runge-Kutta Method
 * "Classical fifth-, sixth-, seventh-, and eighth-order
 * Runge-Kutta formulas with stepsize control"
 * By:  Erwin Fehlberg      1968
 */

/* Time-step coefficients */
static double A[] = {
	0.0,
	1.0 / 3.0,
	2.0 / 5.0,
	1.0,
	2.0 / 3.0,
	4.0 / 5.0
};

/* Solution weighting coefficients */
static double C[] = {
	23.0 / 192.0,
	0.0,
	125.0 / 192.0,
	0.0,
	-27.0 / 64.0,
	125.0 / 192.0
};

/* Simulation weighting coefficients */
static double B0[] = {0.0};
static double B1[] = {
	1.0 / 3.0,
	0.0, 0.0, 0.0, 0.0
};
static double B2[] = {
	4.0 / 25.0,
	6.0 / 25.0,
	0.0, 0.0, 0.0
};
static double B3[] = {
	1.0 / 4.0,
	-3.0,
	15.0 / 4.0,
	0.0, 0.0
};
static double B4[] = {
	2.0 / 27.0,
	10.0 / 9.0,
	-50.0 / 81.0,
	8.0 / 81.0,
	0.0
};
static double B5[] = {
	2.0 / 25.0,
	12.0 / 25.0,
	2.0 / 15.0,
	8.0 / 75.0,
	0.0
};
static double *B[] = {B0, B1, B2, B3, B4, B5};

static const int nStage = 6;

/* Actual integration step happens here */
void rk5step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim) {
double dt = tUpp - tLow;
	printf("\n\ndt = %4.4f   tLow = %4.4f    tUpp = %4.4f \n", dt, tLow, tUpp);

	/// Allocate memory:
	double t[nStage];
	double** z = new double*[nStage];
	for (int i = 0; i < nStage; i++) {
		z[i] = new double[nDim];
	}
	double** f = new double*[nStage];
	for (int i = 0; i < nStage; i++) {
		f[i] = new double[nDim];
	}

	/// Populate time grid:
	for (int iStage = 0; iStage < nStage; iStage++) {
		t[iStage] = tLow + dt * A[iStage];
	}

	/// Initial State:
	int iStage = 0;
	for (int iDim = 0; iDim < nDim; iDim++) {
		z[iStage][iDim] = zLow[iDim];
	}

	/// Dynamics at initial point:
	dynFun(t[0], z[0], f[0]);

	/// March through each stage:
	double sum;
	for (int iStage = 1; iStage < nStage; iStage++) {
		printf("t = %4.4f \nB[%d]:", t[iStage], iStage);
		for (int iDim = 0; iDim < nDim; iDim++) {
			sum = 0.0;
			for (int j = 0; j < iStage; j++) {
				sum = sum + B[iStage][j] * f[iStage - 1][iDim];
				if (iDim == 0) printf("  %4.4f", B[iStage][j]);
			}
			if (iDim == 0) printf("\n");
			z[iStage][iDim] = zLow[iDim] + dt * sum;
		}
		dynFun(t[iStage], z[iStage], f[iStage]);
	}

	/// Compute the final estimate:
	for (int iDim = 0; iDim < nDim; iDim++) {
		sum = 0.0;
		for (int iStage = 0; iStage < nStage; iStage++) {
			sum = sum + C[iStage] * f[iStage][iDim];
		}
		zUpp[iDim] = zLow[iDim] + dt * sum;
	}

	/// Debugging:
	printf("C: ");
	for (int iStage = 0; iStage < nStage; iStage++) {
		printf("%4.4f    ", C[iStage]);
	}	
	for (int iDim = 0; iDim < nDim; iDim++) {
		printf("\nz[%d]:  ", iDim);
		for (int iStage = 0; iStage < nStage; iStage++) {
			printf("%4.4f    ", z[iStage][iDim]);
		}
	}
	for (int iDim = 0; iDim < nDim; iDim++) {
		printf("\nf[%d]:  ", iDim);
		for (int iStage = 0; iStage < nStage; iStage++) {
			printf("%4.4f    ", f[iStage][iDim]);
		}
	}
	printf("\nzLow: ");
	for (int iDim = 0; iDim < nDim; iDim++) {
		printf("%4.4f    ", zLow[iDim]);
	}
	printf("\nzUpp: ");
	for (int iDim = 0; iDim < nDim; iDim++) {
		printf("%4.4f    ", zUpp[iDim]);
	}
	printf("\n");

	/// Release memory:
	for (int i = 0; i < nStage; i++) {
		delete [] z[i];
	}
	delete [] z;
	for (int i = 0; i < nStage; i++) {
		delete [] f[i];
	}
	delete [] f;
}