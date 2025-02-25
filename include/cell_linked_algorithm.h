//
// Created by Андрій Білий on 25.02.2025.
//

#ifndef CELL_LINKED_ALGORITHM_H
#define CELL_LINKED_ALGORITHM_H
#include "generate_config.h"

#define CELL_SIZE SIGMA
#define Mx ((int)(Lx / CELL_SIZE))
#define My ((int)(Ly / CELL_SIZE))
#define NUM_CELLS (Mx * My)

Particle particles[N];
int head[NUM_CELLS];
int list[N];


int get_cell_index(double x, double y);
int is_overlapping_algo(Particle particle, int num_particles);
void insert_particle_in_cell(int p_index, Particle particle);

#endif //CELL_LINKED_ALGORITHM_H
