#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H

#include "generate_config.h"
typedef enum {
    NO_BOUNDARY = 0,
    PERIODIC_X = 1,
    PERIODIC_Y = 2,
    PERIODIC_BOTH = 3
} BoundaryCondition;

extern BoundaryCondition boundary_condition;

void set_condition(BoundaryCondition condition);
Particle adjust_circle_for_periodic(Particle ref, Particle sp);
void periodic_boundary(double *dx, double *dy);
int is_overlapping_pbc(const Particle new_p, const Particle particles[], int count);
#endif
