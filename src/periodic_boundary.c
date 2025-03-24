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

int check_overlap(Particle copy, Particle particles[], int n) {
    for (int j = 0; j < n; j++) {
        Particle neighbor = particles[j];

        if ((fabs(neighbor.x - copy.x) <= SIGMA) || (fabs(neighbor.y - copy.y) <= SIGMA)) {
            if (distance(copy, neighbor) < SIGMA) {
                return 1;
            }
        }
    }
    return 0;
}

void apply_periodic_boundary(Particle particles[], int *n, int mode) {
    int original_count = *n;

    for (int i = 0; i < original_count; i++) {
        Particle p = particles[i];

        if (mode == 2) {
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (dx == 0 && dy == 0) continue;

                    Particle copy = p;
                    copy.x += dx * (Lx / 3);
                    copy.y += dy * (Ly / 3);

                    if (*n >= N) {
                        printf("Reached max number of particles, can't add more.\n");
                        return;
                    }

                    if (is_near_boundary(p) && check_overlap(copy, particles, *n)) {
                        continue;
                    }

                    particles[(*n)++] = copy;
                }
            }
        }
    }
}
