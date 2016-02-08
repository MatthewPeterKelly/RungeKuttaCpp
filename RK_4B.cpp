#include <iostream>
#include "RK_4B.h"
#include "integrator.h"

using namespace std;

/* 4th-Order Explicit "3/8 Rule" Runge-Kutta Method */

/* Time-step coefficients */
static double RK4B__A[] = {
	0.0,
	1.0/3.0,
	2.0/3.0,
	1.0
};

/* Solution weighting coefficients */
static double RK4B__C[] = {
	1.0 / 8.0,
	3.0 / 8.0,
	3.0 / 8.0,
	1.0 / 8.0
};

/* Simulation weighting coefficients */
static double RK4B__B[] = {
	1.0/3.0,
	-1.0/3.0, 	1.0,
	1.0 		-1.0, 	1.0
};

static const int RK4B__nStage = 4;

/* Actual integration step happens here */
void rk4Bstep(DynFun dynFun,
              double tLow, double tUpp, double zLow[], double zUpp[], int nDim)
{
	RK_STEP( dynFun,
	         tLow,  tUpp,  zLow,  zUpp,  nDim,
	         RK4B__A,  RK4B__B,  RK4B__C,  RK4B__nStage);
}