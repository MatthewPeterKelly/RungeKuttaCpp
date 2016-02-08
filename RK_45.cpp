#include <iostream>
#include "RK_45.h"

#include "integrator.h"

using namespace std;

/* Runge-Kutta-Fehlberg 
 * https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
 */
static double RK45__A[] = {
	0.0,
	1.0 / 4.0,
	3.0 / 8.0,
	12.0 / 13.0,
	1.0,
	1.0 / 2.0
};

static double RK45__C[] = {
	16.0 / 135.0,
	0.0,
	6656.0 / 12825.0,
	28561.0 / 56430.0,
	-9.0 / 50.0,
	2.0 / 55.0
};

static double RK45__B[] = {
	1.0 / 4.0,
	3.0 / 32.0,			9.0 / 32.0,
	1932.0 / 2197.0,	-7200.0 / 2197.0,	7296.0 / 2197.0,
	439.0 / 216.0,		-8.0,				3680.0 / 513.0,		-845.0 / 4104.0,
	-8.0 / 27.0,		2.0,				-3544.0 / 2565.0,	1859.0 / 4104.0,	-11.0 / 40.0
};

static const int RK45__nStage = 6;

/* Runge-Kutta-Fehlberg Integration method
 * tLow = time at beginning of the step
 * tUpp = time at the end of the step
 * zLow = state at the beginning of the step
 * zUpp = state at the end of the step (unknown  --  Computed by this function)
 * nDim = dimension of the state space*/
void rk45step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim) {
	RK_STEP( dynFun,
	         tLow,  tUpp,  zLow,  zUpp,  nDim,
	         RK45__A,  RK45__B,  RK45__C,  RK45__nStage);
}