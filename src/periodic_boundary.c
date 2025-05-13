#include "periodic_boundary.h"

#include <general_config.h>
#include <tgmath.h>

BoundaryCondition h_bc;

void set_condition(int cond) {
    h_bc = (BoundaryCondition)cond;
}

void periodic_boundary(double *dx, double *dy) {

    if ((h_bc & PERIODIC_X) != 0) {
        *dx = fmod(*dx + 0.5 * Lx, Lx);
        if (*dx < 0) *dx += Lx;
        *dx -= 0.5 * Lx;
    }
    if ((h_bc & PERIODIC_Y) != 0) {
        *dy = fmod(*dy + 0.5 * Ly, Ly);
        if (*dy < 0) *dy += Ly;
        *dy -= 0.5 * Ly;
    }
}
