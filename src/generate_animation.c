#include "generate_animation.h"


double calculate_energy(double x, double y) {
    return x*x + y*y;
}

double random_uniform() {
    return ((double)rand()/RAND_MAX);
}

// creates an initial configuration using Linked Cell Algorithm
void initial_generate() {
    int count = 0;

    for (int i = 0; i < NUM_CELLS; i++) {
        head[i] = -1;
    }

    srand((unsigned)time(NULL));

    while (count < N) {
        Particle particle;
        particle.x = ((double)rand() / RAND_MAX) * Lx;
        particle.y = ((double)rand() / RAND_MAX) * Ly;

        if (!is_overlapping_algo(particle, count)) {
            particles[count] = particle;
            insert_particle_in_cell(count, particle);
            count++;
        }
    }

}

void write_to_file(int frame_number) {
    char filename[20];
    snprintf(filename, sizeof(filename), "particles_%04d.pdb", frame_number);
    FILE *config = fopen(filename, "w");
    if (config == NULL) {
        printf("Error while opening file\n");
        return;
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
}


// move all particles using Monte Carlo method.
// no optimization, only use for small N!
void no_optimization_move() {
    for (int i = 0; i < N; i++) {
        double x0 = particles[i].x;
        double y0 = particles[i].y;

        double x1 = fabs(x0 + random_uniform() - 0.5);
        double y1 = fabs(y0 + random_uniform() - 0.5);

        double w = exp(-(calculate_energy(x1,y1)-calculate_energy(x0, y0))/kB*T);
        if (w > 1) w = 1;

        if ((random_uniform() < w)) {
            // accept
            particles[i].x = x1;
            particles[i].y = y1;
        }
    }
}


int main() {
    int frame_count = 50;
    initial_generate();
    write_to_file(0);
    for (int i = 1; i < frame_count; i++) {
        no_optimization_move();
        write_to_file(i);
        
    }
}