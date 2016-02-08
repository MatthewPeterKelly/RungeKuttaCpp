#ifndef __RK4B_H__
#define __RK4B_H__

#include "integrator.h"

void rk4Bstep(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim);

#endif