#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H

#include "generate_config.h"

double distance(Particle p1, Particle p2);
int is_near_boundary(Particle p);
int check_overlap(Particle copy, Particle particles[], int n);
void apply_periodic_boundary(Particle particles[], int *n, int mode);
#endif