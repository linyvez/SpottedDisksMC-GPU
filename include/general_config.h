#ifndef GENERAL_CONFIG_H
#define GENERAL_CONFIG_H

#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#include <corecrt_math_defines.h>
#endif

#define Z_AXIS 0.0
#define OVITO_MUL 0.5

#define PARTICLE_SIZE 2.0

#define N_PAR 9000
#define Lx 50.0
#define Ly 50.0

#define MAX_ATTEMPTS 1000
#define ANIMATION_STEPS 3000

#define CELL_SIZE (PARTICLE_SIZE * 2)
#define Mx ((int)(Lx / CELL_SIZE))
#define My ((int)(Ly / CELL_SIZE))
#define NUM_CELLS Mx * My

#define MAX_ANGLE (2 * M_PI)
#define RROTMAX (M_PI / 10)
#define LRROTMAX (M_PI / 2)

#define RDISPMAX PARTICLE_SIZE / 2
#define LRDISPMAX Lx / 2

#define RMOVE 0.5
#define LRMOVE 0.5

#define KT 1
#define PATCH_STRENGTH -4.2 * KT
#define PATCH_RADIUS (PARTICLE_SIZE / 4)
#define PATCH_DELTA (PARTICLE_SIZE / 2)

#define MAX_PATCHES 4

#define MAX_NEIGH 15


typedef struct {
    int (*neighbors)[MAX_NEIGH];
    int *count;
} LinkedCell;

typedef struct {
    double rel_x, rel_y;
} Patch;

typedef enum {
    NO_BOUNDARY = 0,
    PERIODIC_X = 1,
    PERIODIC_Y = 2,
    PERIODIC_BOTH = 3
} BoundaryCondition;




#endif
