#ifndef GENERAL_CONFIG_H
#define GENERAL_CONFIG_H

#include <math.h>

#define MAX_ATTEMPTS 10000
#define MAX_ANGLE (2 * M_PI)

#define MAX_NEIGHBOURS 4

#define RMOVE 0.5
#define LRMOVE 0.5

#define MAX_PATCHES 4

#define Z_AXIS 0.0

extern int lflag;

extern int *visited;

extern int current_epoch;

typedef enum {
    PARTICLE_DISK = 0,
    PARTICLE_SQUARE = 1
} ParticleType;


typedef struct
{
    double rel_x, rel_y;
} Patch;


void set_lflag(int);

#endif