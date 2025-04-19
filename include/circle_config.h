#ifndef CIRCLE_CONFIG_H
#define CIRCLE_CONFIG_H

#include "general_config.h"


typedef struct {
    double x, y;
    double q[4];

    Patch patches[MAX_PATCHES];
} CircleParticle;

extern int *particles_idx;
extern CircleParticle *particles;
extern int (*parts_in_cells)[MAX_NEIGHBOURS];



int get_cell_index(double x, double y);
int is_overlapping_circle(CircleParticle particle);
void insert_particle_in_cell(int p_index, const CircleParticle *particle);
void move_particle(int p_index, double newX, double newY);
int generate_random_circles();
void compute_circle_patch_global_position(const CircleParticle sp, const Patch patch, double *global_x, double *global_y);
void initialize_cells();

#endif // CIRCLE_CONFIG_H