#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H

typedef enum {
    NO_BOUNDARY = 0,
    PERIODIC_X = 1,
    PERIODIC_Y = 2,
    PERIODIC_BOTH = 3
} BoundaryCondition;

extern BoundaryCondition boundary_condition;

void periodic_boundary(double *, double *);

void set_condition(BoundaryCondition);

double apply_boundary(double, double, int);
#endif // PERIODIC_BOUNDARY_H
