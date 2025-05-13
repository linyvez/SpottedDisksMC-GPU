#include "circle_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <patch.h>
#include <tgmath.h>

#include "periodic_boundary.h"
#include "shared_utilities.h"


inline int cell_index(double coord) {
    int raw = (int) (coord / CELL_SIZE);

    if (raw >= Mx) {
        raw -= 1;
    }
    return raw;
}


static CircleParticle adjust_circle_for_periodic(const CircleParticle ref, const CircleParticle cp) {
    CircleParticle sp_adj = cp;
    double dx = cp.x - ref.x;
    double dy = cp.y - ref.y;

    periodic_boundary(&dx, &dy);

    sp_adj.x = ref.x + dx;
    sp_adj.y = ref.y + dy;

    return sp_adj;
}


static int check_circles_overlap(CircleParticle a, CircleParticle b) {
    CircleParticle b_adj = adjust_circle_for_periodic(a, b);

    double x_diff = a.x - b_adj.x;
    double y_diff = a.y - b_adj.y;

    double distance = x_diff * x_diff + y_diff * y_diff;

    if (distance < PARTICLE_SIZE * PARTICLE_SIZE) {
        return 1;
    }
    return 0;
}

static int is_overlapping_circle(const CircleParticle cp, const CircleParticle * circles, LinkedCell cell_struct) {
    int cell_ix = cell_index(cp.x);
    int cell_iy = cell_index(cp.y);

    for (int ox = -1; ox <= 1; ox++) {
        for (int oy = -1; oy <= 1; oy++) {
            int ghost_ix = (cell_ix + ox + Mx) % Mx;
            int ghost_iy = (cell_iy + oy + My) % My;
            int ghost_cell = ghost_ix + ghost_iy * Mx;

            int cnt = cell_struct.count[ghost_cell];
            for (int i = 0; i < cnt; i++) {
                int nidx = cell_struct.neighbors[ghost_cell][i];
                if (check_circles_overlap(circles[nidx], cp)) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

static void insert_circle_in_cells(int circleIndex, const CircleParticle cp, LinkedCell cell_struct) {
    int cell = cell_index(cp.x) + cell_index(cp.y) * Mx;
    int cnt = cell_struct.count[cell];
    cell_struct.neighbors[cell][cnt] = circleIndex;
    cell_struct.count[cell]++;
}



int generate_random_circles(int patch, CircleParticle * circles, LinkedCell cell_struct) {
    int totalCircles = 0;
    int attempts = 0;

    h_patch = assign_patch_type(patch);

    while (totalCircles < N_PAR) {
        if (attempts++ > MAX_ATTEMPTS) {
            fprintf(stderr, "Error: Too many attempts to place circle.\n");
            break;
        }
        CircleParticle cp;
        cp.x = (double) rand() / RAND_MAX * Lx;
        cp.y = (double) rand() / RAND_MAX * Ly;

        double theta = (double) rand() / RAND_MAX * MAX_ANGLE;

        cp.q[0] = sin(theta / 2);
        cp.q[1] = cos(theta / 2);

        if (!is_overlapping_circle(cp, circles, cell_struct)) {
            circles[totalCircles] = cp;
            insert_circle_in_cells(totalCircles, cp, cell_struct);
            totalCircles++;
            attempts = 0;
        }
    }

    FILE *f = fopen("data/configuration_circle.xyz", "w");
    if (!f) {
        printf("Error opening configuration file\n");
        return 0;
    }
    write_file_c(f, circles, totalCircles);
    fclose(f);
    return totalCircles;
}
