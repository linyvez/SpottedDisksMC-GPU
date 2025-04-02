#ifndef CIRCLE_CONFIG_H
#define CIRCLE_CONFIG_H

#include "general_config.h"

#define SIGMA   1
#define PATCH_NUMBER    4

#define Mx_C ((int)(Lx / SIGMA))
#define My_C ((int)(Ly / SIGMA))
#define NUM_CELLS_C (Mx_C * My_C)

#define MAX_DISPLACEMENT SIGMA / 2

#define PATCH_RADIUS_C SIGMA / 4
#define PATCH_DELTA_C SIGMA / 4

typedef struct {
    double x, y;
    double q[4];

    Patch patches[PATCH_NUMBER];
} CircleParticle;

extern int particles_idx[];
extern CircleParticle particles[];
extern int parts_in_cells[][4];

extern CircleParticle circles[];



int get_cell_index(double x, double y);
int is_overlapping_circle(CircleParticle particle);
int insert_particle_in_cell(int p_index, CircleParticle *particle);
void move_particle(int p_index, double newX, double newY);
int generate_random_circles();
void compute_circle_patch_global_position(const CircleParticle sp, const Patch patch, double *global_x, double *global_y);
void initialize_cells();

#endif;