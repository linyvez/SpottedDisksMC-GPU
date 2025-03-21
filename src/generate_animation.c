#include "generate_animation.h"


double calculate_energy(double x, double y) {
    return x*x + y*y;
}

double random_uniform() {
    return ((double)rand()/RAND_MAX);
}

int is_overlaping2(int i) {
    for (int j = 0; j < N; j++) {
        if (j == i) continue;
        double x_diff = particles[i].x - particles[j].x;
        double y_diff = particles[i].y - particles[j].y;
        double distance = sqrt(x_diff*x_diff + y_diff*y_diff);
        if (distance < SIGMA) {
            return 1;
        }
    }
    return 0;
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

        // radius: vector + distance
        double x1 = x0 + MAX_DISPLACEMENT * (random_uniform() - 0.5);
        double y1 = y0 + MAX_DISPLACEMENT * (random_uniform() - 0.5);

        x1 = (x1 > Lx) ? Lx : x1;
        y1 = (y1 > Ly) ? Ly : y1;
        // double w = exp(-(calculate_energy(x1,y1)-calculate_energy(x0, y0))/kB*T);
        // if (w > 1) w = 1;
        double w = 1;

        // if ((random_uniform() < w)) {
            // accept
        particles[i].x = x1;
        particles[i].y = y1;
        if (is_overlaping2(i)) {
            particles[i].x = x0;
            particles[i].y = y0;
        }
            
        // }
    }
}

void randomize_coordinates(Particle particle) {
    double x1 = fabs(particle.x + random_uniform() - 0.5);
    double y1 = fabs(particle.y + random_uniform() - 0.5);

    // make sure the particles do not leave the field
    x1 = (x1 > Lx) ? x1 - Lx : x1;
    y1 = (y1 > Ly) ? y1 - Ly : y1;

    particle.x = x1;
    particle.y = y1;
}


// optimized using linked cell algorithm
void linked_move() {
    for (int i = 0; i < N; i++) {
        double x0 = particles[i].x;
        double y0 = particles[i].y;

        // the coordinates do not change here... .... ....
        while(is_overlapping_algo(particles[i])) {
            printf("Particle %d: x = %f, y = %f\n", i, particles[i].x, particles[i].y);
            randomize_coordinates(particles[i]);
        }

        double w = exp(-(calculate_energy(particles[i].x,particles[i].y)-calculate_energy(x0, y0))/kB*T);
        if (w > 1) w = 1;

        if ((random_uniform() > w)) {
            // reject
            particles[i].x = x0;
            particles[i].y = y0;
        }
    }
}


int main() {
    int frame_count = 50;
    generate_random_algo();
    write_to_file(0);
    for (int i = 1; i < frame_count; i++) {
        no_optimization_move();
        write_to_file(i);
    }
}