#ifndef PERIODIC_BOUNDARY_CUH
#define PERIODIC_BOUNDARY_CUH


#ifdef __cplusplus
extern "C" {
#endif

#include "general_config.h"

    extern BoundaryCondition h_bc;
    void set_condition(int cond);
    void periodic_boundary(double *dx, double *dy);

#ifdef __cplusplus
}
#endif

#endif // PERIODIC_BOUNDARY_CUH
