#include "periodic_boundary.h"

double distance(Particle p1, Particle p2) {
    double dx = fabs(p1.x - p2.x);
    double dy = fabs(p1.y - p2.y);

    if (dx > Lx / 2) dx = Lx - dx;
    if (dy > Ly / 2) dy = Ly - dy;

    return sqrt(dx * dx + dy * dy);
}

int is_near_boundary(Particle p) {
    double left_edge = Lx / 3;
    double right_edge = 2 * (Lx / 3);
    double bottom_edge = Ly / 3;
    double top_edge = 2 * (Ly / 3);

    return (fabs(p.x - left_edge) <= SIGMA / 2 || fabs(p.x - right_edge) <= SIGMA / 2 ||
            fabs(p.y - bottom_edge) <= SIGMA / 2 || fabs(p.y - top_edge) <= SIGMA / 2);
}

int check_overlap(Particle copy, int n) {
    for (int j = 0; j < n; j++) {
        Particle *neighbor = &particles[j];

        if ((fabs(neighbor->x - copy.x) <= SIGMA) || (fabs(neighbor->y - copy.y) <= SIGMA)) {
            if (distance(copy, *neighbor) < SIGMA) {
                return 1;
            }
        }
    }
    return 0;
}

void apply_periodic_boundary(int *n, int mode) {
    int original_count = *n;

    for (int i = 0; i < original_count; i++) {
        Particle *p = &particles[i];

        if (mode == 2) {
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (dx == 0 && dy == 0) continue;

                    Particle copy = *p;
                    copy.x += dx * (Lx / 3);
                    copy.y += dy * (Ly / 3);

                    if (*n >= N) {
                        printf("Reached max number of particles, can't add more.\n");
                        return;
                    }

                    if (is_near_boundary(*p) && is_overlapping_algo(copy)) {
                        printf("Skip: (%f, %f)\n", copy.x, copy.y);
                        continue;
                    }

                    particles[(*n)++] = copy;
                    insert_particle_in_cell(*n-1, &copy);
                }
            }
        }
    }
}



void generate_random_with_pbc(int mode) {
    initialize_cells();
    // 1. Generate the center cell
    int count = 0;
    // int attempts = 0;
    // int max_attempts = 3000;

    srand((unsigned)time(NULL));

    double region_width = Lx / 3.0 - SIGMA;
    double region_height = Ly / 3.0 - SIGMA;
    double offset_x = (Lx - region_width) / 2.0;
    double offset_y = (Ly - region_height) / 2.0;


    while (count < ceil(N / 9.0)) {
        Particle particle;
        particle.x = ((double)rand() / RAND_MAX) * region_width + offset_x;
        particle.y = ((double)rand() / RAND_MAX) * region_height + offset_y;


        if (!is_overlapping_algo(particle)) {
            particles[count] = particle;
            insert_particle_in_cell(count, &particle);
            count++;}
            // attempts = 0;
        // } else {
        //     attempts++;
        // }

        // if (attempts >= max_attempts) {
        //     printf("Max attempts reached, stopping particle generation.\n");
        //     break;
        // }
    }

    printf("Generated %d particles in the central region\n", count);

    apply_periodic_boundary(&count, mode);

}