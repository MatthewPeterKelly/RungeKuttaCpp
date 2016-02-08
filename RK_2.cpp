#include <iostream>
#include "RK_2.h"

#include "integrator.h"

using namespace std;

/* 4th-Order Explicit "Classical" Runge-Kutta Method */

/* Time-step coefficients */
static double RK2__A[] = {0.0, 0.5};

/* Solution weighting coefficients */
static double RK2__C[] = {0.0, 1.0};

/* Simulation weighting coefficients */
static double RK2__B[] = {0.5};

static const int RK2__nStage = 2;

/* Actual integration step happens here */
void rk2step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim){
		RK_STEP( dynFun,
	         tLow,  tUpp,  zLow,  zUpp,  nDim,
	         RK2__A,  RK2__B,  RK2__C,  RK2__nStage);
}
