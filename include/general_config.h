#ifndef GENERAL_CONFIG_H
#define GENERAL_CONFIG_H

#include <math.h>


#define Z_AXIS 0.0
#define OVITO_MUL 0.5

#define PARTICLE_SIZE 2.0

#define N 3000
#define Lx 96.0
#define Ly 96.0

#define MAX_ATTEMPTS 100
#define ANIMATION_STEPS 2000

#define CELL_SIZE (PARTICLE_SIZE * 1.41421356)
#define Mx ((int)(Lx / CELL_SIZE + 0.99))
#define My ((int)(Ly / CELL_SIZE + 0.99))
#define NUM_CELLS Mx * My

#define MAX_ANGLE (2 * M_PI)
#define RROTMAX (M_PI / 10)
#define LRROTMAX (M_PI / 2)

#define RDISPMAX PARTICLE_SIZE
#define LRDISPMAX Lx / 2

#define RMOVE 0.5
#define LRMOVE 0.5

#define LFLAG 1

#define KT 0.1
#define PATCH_STRENGTH -4.2 * KT
#define PATCH_RADIUS (PARTICLE_SIZE / 4)
#define PATCH_DELTA (PARTICLE_SIZE / 2)

#define MAX_PATCHES 4
#define MAX_NODES (N * 16)


typedef struct {
    int squareIndex;
    int next;
} Node;


extern Node nodePool[MAX_NODES];
extern int nodeCount;
extern int freeList;

extern int head[NUM_CELLS];



typedef struct
{
    double rel_x, rel_y;
    double strength;
} Patch;


static const double PATCH_FOUR_CORNER_M[4][2] = {
    {PATCH_DELTA, -PATCH_DELTA},
    {-PATCH_DELTA, PATCH_DELTA},
    {-PATCH_DELTA, -PATCH_DELTA},
    {PATCH_DELTA, PATCH_DELTA}
};

static const double PATCH_FOUR_CENTER_M[4][2] = {
    {0, -PATCH_DELTA},
    {0, PATCH_DELTA},
    {-PATCH_DELTA, 0},
    {PATCH_DELTA, 0}
};

static const double PATCH_TWO_CENTER_M[2][2] = {
    {PATCH_DELTA, 0},
    {-PATCH_DELTA, 0},
};





#endif