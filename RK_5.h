#ifndef __RK5_H__
#define __RK5_H__

#include "integrator.h"

void rk5step(DynFun dynFun, double tLow, double tUpp, double zLow[], double zUpp[], int nDim);

#endif