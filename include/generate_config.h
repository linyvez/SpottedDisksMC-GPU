#ifndef GENERATE_CONFIG_H
#define GENERATE_CONFIG_H

#define N 1000
#define Lx 50
#define Ly 50
#define SIGMA 1.0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// add patches here
typedef struct Particle
{
    double x, y;
    // int cell;
} Particle;

int is_overlaping(Particle particle, Particle particles[], int n);
void generate_random();
void generate_square_grid();
void generate_hexagonal_grid();

#endif