#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef enum {
    NO_BOUNDARY = 0,
    PERIODIC_X = 1,
    PERIODIC_Y = 2,
    PERIODIC_BOTH = 3
} BoundaryCondition;

extern BoundaryCondition boundary_condition;

void periodic_boundary(double *dx, double *dy);
void set_condition(BoundaryCondition condition);
#endif //PERIODIC_BOUNDARY_H
