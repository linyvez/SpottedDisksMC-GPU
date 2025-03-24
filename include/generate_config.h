#define N 300
#define Lx 24.0
#define Ly 24.0
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
