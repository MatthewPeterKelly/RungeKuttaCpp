#ifndef __RK2_H__
#define __RK2_H__

#include "integrator.h"

void rk2step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim);

#endif