#include <iostream>
#include "RK_4A.h"
#include "integrator.h"

using namespace std;

/* 4th-Order Explicit "Classical" Runge-Kutta Method */

/* Time-step coefficients */
static double RK4A__A[] = {
	0.0,
	0.5,
	0.5,
	1.0
};

/* Solution weighting coefficients */
static double RK4A__C[] = {
	1.0 / 6.0,
	1.0 / 3.0,
	1.0 / 3.0,
	1.0 / 6.0
};

/* Simulation weighting coefficients */
static double RK4A__B[] = {
	0.5,
	0.0, 0.5,
	0.0, 0.0, 1.0
};

static const int RK4A__nStage = 4;

/* Actual integration step happens here */
void rk4Astep(DynFun dynFun,
              double tLow, double tUpp, double zLow[], double zUpp[], int nDim)
{
	RK_STEP( dynFun,
	         tLow,  tUpp,  zLow,  zUpp,  nDim,
	         RK4A__A,  RK4A__B,  RK4A__C,  RK4A__nStage);
}