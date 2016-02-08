#ifndef __RK10_H__
#define __RK10_H__

#include "integrator.h"

void rk10step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim);

#endif