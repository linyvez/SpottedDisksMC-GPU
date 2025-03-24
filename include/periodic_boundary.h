#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H

#include "generate_config.h"
#include "cell_linked_algorithm.h"

double distance(Particle p1, Particle p2);
int is_near_boundary(Particle p);
int check_overlap(Particle copy, int n);
void apply_periodic_boundary(int *n, int mode);
void generate_random_with_pbc(int mode);
#endif