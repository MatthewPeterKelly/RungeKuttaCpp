#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

typedef void (*DynFun)(double, double[], double[]);

enum IntegrationMethod {
	Euler,
	MidPoint,
	RungeKutta,
	RK_2,
	RK_4A,
	RK_4B,
	RK_45,
	RK_5,
	RK_10
};

void RK_STEP(DynFun dynFun,
             double tLow, double tUpp, double zLow[], double zUpp[], int nDim,
             double A[], double B[], double C[], int nStage);

void simulate(DynFun dynFun,double t0, double t1, 
	double z0[], double z1[], int nDim, int nStep, 
	IntegrationMethod method);

#endif