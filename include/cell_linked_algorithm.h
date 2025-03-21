
#ifndef CELL_LINKED_ALGORITHM_H
#define CELL_LINKED_ALGORITHM_H
#include "generate_config.h"

#define CELL_SIZE SIGMA
// number of cells along the x and y axes respectively
#define Mx ((int)(Lx / CELL_SIZE))
#define My ((int)(Ly / CELL_SIZE))
#define NUM_CELLS (Mx * My)


// maximum number of particles in a cell


extern int particles_idx[];
extern Particle particles[];
extern int parts_in_cells[][4];


int get_cell_index(double x, double y);
int is_overlapping_algo(Particle particle);
int insert_particle_in_cell(int p_index, Particle particle);
void move_particle(int p_index, double newX, double newY);
void generate_random_algo();
void initialize_cells();
#endif //CELL_LINKED_ALGORITHM_H
