#ifndef SQUARE_CONFIG_H
#define SQUARE_CONFIG_H

#include "general_config.h"


typedef struct {
    double x, y;
    double q[2];
    double side;

    int cell_min_x, cell_min_y;
    int cell_max_x, cell_max_y;

    double patch_radius;

    Patch patches[MAX_PATCHES];
} SquareParticle;


extern SquareParticle squares[N];

int is_overlapping_square(SquareParticle, int);

void update_square_AABB(SquareParticle *);

void insert_square_in_cells(int, SquareParticle);

int generate_random_squares(int, const double (*)[2]);

void compute_patch_global_position(SquareParticle, Patch, double *, double *);

#endif // SQUARE_CONFIG_H
