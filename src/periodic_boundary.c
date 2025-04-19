#include "periodic_boundary.h"
#include "config.h"

BoundaryCondition boundary_condition = NO_BOUNDARY;

#include <math.h>

void periodic_boundary(double *dx, double *dy) {
    if (boundary_condition == PERIODIC_X || boundary_condition == PERIODIC_BOTH) {
        *dx = fmod(*dx + 0.5 * GL_CFG->Lx, GL_CFG->Lx);
        if (*dx < 0)
            *dx += GL_CFG->Lx;
        *dx -= 0.5 * GL_CFG->Lx;
    }
    if (boundary_condition == PERIODIC_Y || boundary_condition == PERIODIC_BOTH) {
        *dy = fmod(*dy + 0.5 * GL_CFG->Ly, GL_CFG->Ly);
        if (*dy < 0)
            *dy += GL_CFG->Ly;
        *dy -= 0.5 * GL_CFG->Ly;
    }
}

void set_condition(BoundaryCondition condition) {
    boundary_condition = condition;
}


double apply_boundary(double coord, double L, int periodic) {
    if (periodic) {
        if (coord < 0)
            coord += L;
        else if (coord > L)
            coord -= L;
    } else {
        if (coord < 0)
            coord = 0;
        else if (coord > L)
            coord = L;
    }
    return coord;
}
