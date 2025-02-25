#include "cell_linked_algorithm.h"


int get_cell_index(double x, double y) {
    int ix = (int)(x / CELL_SIZE);
    int iy = (int)(y / CELL_SIZE);
    if (ix >= Mx) ix = Mx - 1;
    if (iy >= My) iy = My - 1;
    return ix + Mx * iy;
}

int is_overlapping_algo(Particle particle, int num_particles) {
    int cell = get_cell_index(particle.x, particle.y);

    int ix = cell % Mx;
    int iy = cell / Mx;

    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int nx = ix + dx;
            int ny = iy + dy;
            if (nx < 0 || nx >= Mx || ny < 0 || ny >= My)
                continue;
            int ncell = nx + Mx * ny;
            for (int j = head[ncell]; j != -1; j = list[j]) {
                double x_diff = particle.x - particles[j].x;
                double y_diff = particle.y - particles[j].y;
                double distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                if (distance < SIGMA) {
                    return 1;
                }
            }
        }
    }
    return 0;
}


void insert_particle_in_cell(int p_index, Particle particle) {
    int cell = get_cell_index(particle.x, particle.y);
    list[p_index] = head[cell];
    head[cell] = p_index;
}