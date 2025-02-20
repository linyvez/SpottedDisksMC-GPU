#define N 350
#define Lx 24.0
#define Ly 24.0
#define SIGMA 1.0


typedef struct Particle
{
    double x, y;
} Particle;

int is_overlaping(Particle particle, Particle particles[], int n);

