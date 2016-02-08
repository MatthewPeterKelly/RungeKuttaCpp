#include <iostream>
#include "RK_4A.h"
#include "integrator.h"

using namespace std;

/* 4th-Order Explicit "Classical" Runge-Kutta Method */

/* Time-step coefficients */
static double A[] = {
	0.0,
	0.5,
	0.5,
	1.0
};

/* Solution weighting coefficients */
static double C[] = {
	1.0 / 6.0,
	1.0 / 3.0,
	1.0 / 3.0,
	1.0 / 6.0
};

/* Simulation weighting coefficients */
static double B[] = {
	0.5,
	0.0, 0.5,
	0.0, 0.0, 1.0
};

static const int nStage = 4;

/* Actual integration step happens here */
void rk4Astep(DynFun dynFun, 
	double tLow, double tUpp, double zLow[], double zUpp[], int nDim) 
{
	
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
	double dt = tUpp - tLow;
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
	int idx;
	for (int iStage = 1; iStage < nStage; iStage++) {
		for (int iDim = 0; iDim < nDim; iDim++) {
			sum = 0.0;
			for (int j = 0; j < iStage; j++) {
				idx = iStage*(iStage-1)/2 + j;   // Triangle numbers   
				sum = sum + B[idx]*f[iStage - 1][iDim];
			}
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