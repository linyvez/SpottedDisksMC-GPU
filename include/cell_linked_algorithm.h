
#ifndef CELL_LINKED_ALGORITHM_H
#define CELL_LINKED_ALGORITHM_H
#include "generate_config.h"

#define CELL_SIZE SIGMA
// number of cells along the x and y axes respectively
#define Mx ((int)(Lx / CELL_SIZE))
#define My ((int)(Ly / CELL_SIZE))
#define NUM_CELLS (Mx * My)

extern int head[];  // Declare, but do not define
extern Particle particles[];
extern int list[];


int get_cell_index(double x, double y);
int is_overlapping_algo(Particle particle);
void insert_particle_in_cell(int p_index, Particle particle);
void generate_random_algo();
#endif //CELL_LINKED_ALGORITHM_H
