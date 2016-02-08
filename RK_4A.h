#ifndef __RK4A_H__
#define __RK4A_H__

#include "integrator.h"

void rk4Astep(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim);

#endif