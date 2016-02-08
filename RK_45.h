#ifndef __RK45_H__
#define __RK45_H__

#include "integrator.h"

void rk45step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim);

#endif