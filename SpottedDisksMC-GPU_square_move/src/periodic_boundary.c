#include "periodic_boundary.h"
#include "square_config.h"

BoundaryCondition boundary_condition = NO_BOUNDARY;

#include <math.h>

void periodic_boundary(double *dx, double *dy) {
    if (boundary_condition == PERIODIC_X || boundary_condition == PERIODIC_BOTH) {
        *dx = fmod(*dx + 0.5 * Lx, Lx);
        if (*dx < 0)
            *dx += Lx;
        *dx -= 0.5 * Lx;
    }
    if (boundary_condition == PERIODIC_Y || boundary_condition == PERIODIC_BOTH) {
        *dy = fmod(*dy + 0.5 * Ly, Ly);
        if (*dy < 0)
            *dy += Ly;
        *dy -= 0.5 * Ly;
    }
}

void set_condition(BoundaryCondition condition) {
     boundary_condition = condition;
}