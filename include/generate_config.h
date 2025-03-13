#ifndef GENERATE_CONFIG_H
#define GENERATE_CONFIG_H

#define N 50
#define Lx 10
#define Ly 10
#define SIGMA 1.0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct Particle
{
    double x, y;
} Particle;

int is_overlaping(Particle particle, Particle particles[], int n);
void generate_random();
void generate_square_grid();
void generate_hexagonal_grid();

#endif