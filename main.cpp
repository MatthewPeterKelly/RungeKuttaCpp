#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

#include "integrator.h"


/* Test dynamics function --  simple pendulum*/
void simplePendulum(double t, double z[], double dz[]) {
	double x = z[0];
	double v = z[1];
	double dx = v;
	double dv = -0.1 * v - sin(x);
	dz[0] = dx;
	dz[1] = dv;
}


/* Test dynamics function -- driven damped pendulum*/
void drivenDampedPendulum(double t, double z[], double dz[]) {
	double x = z[0];
	double v = z[1];
	double u = cos(t);
	double dx = v;
	double dv = u - 0.1 * v - sin(x);
	dz[0] = dx;
	dz[1] = dv;
}

int main()
{
	double t0 = 0.0;
	double t1 = 10 * (2.0 * M_PI);
	double z0[2];
	double z1[2];
	int nDim = 2;
	double dt = 0.2;
	int nStep = ceil((t1 - t0) / dt);

	DynFun dynFun = simplePendulum;
	// DynFun dynFun = drivenDampedPendulum;

	z0[0] = 1.9;
	z0[1] = -4.5;

	// IntegrationMethod method = Euler;      	// 1st-order method
	// IntegrationMethod method = MidPoint;  	// 2nd-order method
	// IntegrationMethod method = RungeKutta;  // Classical 4th-order Runge-Kutta
	// IntegrationMethod method = RK_2;  		// 2nd-order Runge-Kutta (MidPoint Method)
	// IntegrationMethod method = RK_4A;  		// Classical 4th-order Runge-Kutta
	// IntegrationMethod method = RK_4B;  		// Classical 4th-order Runge-Kutta
	IntegrationMethod method = RK_45; 			// Runge-Kutta-Fehlberg
	// IntegrationMethod method = RK_5;  		// 5th-order Runge-Kutta
	// IntegrationMethod method = RK_10;  		// 10th-order Runge-Kutta

	simulate(dynFun, t0, t1, z0, z1, nDim, nStep, method);

}

