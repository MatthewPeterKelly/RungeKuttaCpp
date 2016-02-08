#include <iostream>
#include "RK_5.h"

#include "integrator.h"

using namespace std;

/* 5th-Order Explicit Runge-Kutta Method
 * "Classical fifth-, sixth-, seventh-, and eighth-order
 * Runge-Kutta formulas with stepsize control"
 * By:  Erwin Fehlberg      1968
 */

/* Time-step coefficients */
static double RK5__A[] = {
	0.0,
	1.0/3.0,
	2.0/5.0,
	1.0,
	2.0/3.0,
	4.0/5.0
};

/* Solution weighting coefficients */
static double RK5__C[] = {
	23.0/192.0,
	0.0,
	125.0/192.0,
	0.0,
	-27.0/64.0,
	125.0/192.0
};

/* Simulation weighting coefficients */
static double RK5__B[] = {
	1.0/3.0,
	4.0/25.0, 	6.0/25.0,
	1.0/4.0,		-3.0,			15.0/4.0,
	2.0/27.0, 	10.0/9.0,  	-50.0/81.0,   8.0/81.0,
	2.0/25.0, 	12.0/25.0,  	2.0/15.0,  	8.0/75.0,  0.0
};

static const int RK5__nStage = 6;

/* Runge-Kutta 5th-order method
 * tLow = time at beginning of the step
 * tUpp = time at the end of the step
 * zLow = state at the beginning of the step
 * zUpp = state at the end of the step (unknown  --  Computed by this function)
 * nDim = dimension of the state space*/
 void rk5step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim) {
	RK_STEP( dynFun,
	         tLow,  tUpp,  zLow,  zUpp,  nDim,
	         RK5__A,  RK5__B,  RK5__C,  RK5__nStage);
}