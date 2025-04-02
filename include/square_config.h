#ifndef SQUARE_CONFIG_H
#define SQUARE_CONFIG_H

#include "general_config.h"

#define SQUARE_SIDE 2.0

#define CELL_SIZE SQUARE_SIDE * 1.41421356
#define Mx ((int)((Lx + 1) / CELL_SIZE))
#define My ((int)((Ly + 1) / CELL_SIZE))
#define NUM_CELLS Mx *My

#define PATCH_RADIUS SQUARE_SIDE / 4
#define PATCH_DELTA SQUARE_SIDE / 2

static const double CORNERS_MODEL[4][2] = {
    {PATCH_DELTA, -PATCH_DELTA},
    {-PATCH_DELTA, PATCH_DELTA},
    {-PATCH_DELTA, -PATCH_DELTA},
    {PATCH_DELTA, PATCH_DELTA}
};

static const double SIDES_MODEL[4][2] = {
    {0, -PATCH_DELTA},
    {0, PATCH_DELTA},
    {-PATCH_DELTA, 0},
    {PATCH_DELTA, 0}
};

static const double TSIDES_MODEL[2][2] = {
    {PATCH_DELTA, 0},
    {0, -PATCH_DELTA},
};

static const double TNPSIDES_MODEL[2][2] = {
    {PATCH_DELTA, 0},
    {-PATCH_DELTA, 0},
};







typedef struct
{
    double x, y, z;
    double q[4];
    double shape[3];

    int cell_min_x, cell_min_y;
    int cell_max_x, cell_max_y;

    Patch patches[4];
} SquareParticle;

typedef struct
{
    int squareIndex;
    int next;
} Node;

#define MAX_NODES (N * 16)
extern Node nodePool[MAX_NODES];
extern long nodeCount;
extern int freeList;

extern int npatches;

extern int head[NUM_CELLS];
extern SquareParticle squares[N];

int is_overlapping_square(SquareParticle, int);
void update_square_AABB(SquareParticle *);
void insert_square_in_cells(int, SquareParticle);
int generate_random_squares(int);

void compute_patch_global_position(const SquareParticle, const Patch, double *, double *);

#endif // SQUARE_CONFIG_H
