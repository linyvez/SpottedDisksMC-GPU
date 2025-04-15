#ifndef GENERAL_CONFIG_H
#define GENERAL_CONFIG_H

#include <math.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "periodic_boundary.h"

#define N 300
#define Lx 25.0
#define Ly 25.0

#define MAX_ATTEMPTS 100
#define ANIMATION_STEPS 200

#define MAX_ANGLE (2 * M_PI)
#define MAX_ROTATION_ANGLE (M_PI / 10)
#define LMAX_ROTATION_ANGLE (M_PI / 2)
#define RMOVE 0.5
#define RLMOVE 0.5

#define LFLAG 0

#define KT 0.1

#define PATCH_STRENGTH -4.2 * KT


typedef struct {
    double x;
    double y;
} Coordinates;

typedef struct
{
    double rel_x, rel_y;
    double shape[3];
    double strength;
} Patch;

#endif