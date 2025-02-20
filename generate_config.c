#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 350
#define Lx 24.0
#define Ly 24.0
#define SIGMA 1.0

typedef struct Particle
{
    double x, y;
} Particle;

int is_overlaping(Particle particle, Particle particles[], int n) {
    for (int i = 0; i < n; i++) {
        double x_diff = particle.x - particles[i].x;
        double y_diff = particle.y - particles[i].y;
        double distance = sqrt(x_diff*x_diff + y_diff*y_diff);
        if (distance < SIGMA) {
            return 1;
        }
    }
    return 0;
}

int main() {
    Particle particles[N];
    int count = 0;

    srand((unsigned)time(NULL));

    while (count < N) {
        Particle particle;
        particle.x = ((double)rand() / RAND_MAX) * Lx;
        particle.y = ((double)rand() / RAND_MAX) * Ly;

        if (!is_overlaping(particle, particles, count)) {
            particles[count] = particle;
            count++;
        }
    }

    FILE *config = fopen("configuration.pdb", "w");
    if (config == NULL) {
        printf("Error while opening file\n");
    }

    for (int i = 0; i < N; i++) {
        fprintf(config, "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
            i+1,
            "ATOM",
            "RES",
            "A",      
            1,         
            particles[i].x,
            particles[i].y,
            0.0
            );
    }
    
    fclose(config);
    return 0;
}
