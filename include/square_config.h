#ifndef SQUARE_CONFIG_H
#define SQUARE_CONFIG_H

#include <math.h>

#define N 140
#define Lx 48.0
#define Ly 48.0
#define SQUARE_SIDE 2.0
#define MAX_ANGLE (2 * M_PI)


#define CELL_SIZE SQUARE_SIDE * 1.41421356
#define Mx ((int)((Lx + 1) / CELL_SIZE))
#define My ((int)((Ly + 1) / CELL_SIZE))
#define NUM_CELLS Mx * My

#define PATCH_STRENGTH 0
#define PATCH_RADIUS SQUARE_SIDE / 4


typedef struct {
    double rel_x, rel_y;
    double shape[3];
    double strength;
} Patch;

typedef struct {
    double x, y, z;
    double q[4];
    double shape[3];

    int cell_min_x, cell_min_y;
    int cell_max_x, cell_max_y;

    Patch patches[4];
} SquareParticle;


typedef struct {
    int squareIndex;
    int next;
} Node;

#define MAX_NODES (N * 16)
Node nodePool[MAX_NODES];
extern long nodeCount;
extern int freeList;


int head[NUM_CELLS];

SquareParticle squares[N];


int is_overlapping_square(SquareParticle, int);
void update_square_AABB(SquareParticle *);
void insert_square_in_cells(int, SquareParticle);
int generate_random_squares();

void compute_patch_global_position(const SquareParticle, const Patch, double *, double *);

#endif //SQUARE_CONFIG_H
