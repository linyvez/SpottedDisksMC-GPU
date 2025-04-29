#include "square_move.h"
#include "square_config.h"
#include "periodic_boundary.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double compute_patch_interaction(const SquareParticle sp1, const SquareParticle sp2, int num_patches) {
    double energy = 0.0;
    for (int i = 0; i < num_patches; i++) {
        double gx1, gy1;
        compute_patch_global_position(sp1, sp1.patches[i], &gx1, &gy1);
        for (int j = 0; j < num_patches; j++) {
            double gx2, gy2;
            compute_patch_global_position(sp2, sp2.patches[j], &gx2, &gy2);
            double dx = gx1 - gx2;
            double dy = gy1 - gy2;
            double distance = sqrt(dx * dx + dy * dy);
            if (distance < 2 * sp1.patch_radius) {
                energy += sp1.patches[i].strength;
            }
        }
    }
    return energy;
}

double compute_patch_energy(const SquareParticle sp, int num_patches) {
    double sum = 0.0;

    int uniqueNeighbors[64];
    int uniqueCount = 0;

    for (int ix = sp.cell_min_x; ix <= sp.cell_max_x; ix++) {
        for (int iy = sp.cell_min_y; iy <= sp.cell_max_y; iy++) {
            for (int ox = -1; ox <= 1; ox++) {
                for (int oy = -1; oy <= 1; oy++) {
                    int ghost_ix = (ix + ox + Mx) % Mx;
                    int ghost_iy = (iy + oy + My) % My;
                    int ghost_cell = ghost_ix + ghost_iy * Mx;
                    if ((ghost_ix != ix || ghost_iy != iy) && (abs(ghost_ix - ix) <= 1 && abs(ghost_iy - iy) <= 1))
                        continue;
                    for (int nodeIndex = head[ghost_cell]; nodeIndex != -1; nodeIndex = nodePool[nodeIndex].next) {
                        int neighbor_idx = nodePool[nodeIndex].squareIndex;


                        int alreadyIncluded = 0;
                        for (int k = 0; k < uniqueCount; k++) {
                            if (uniqueNeighbors[k] == neighbor_idx) {
                                alreadyIncluded = 1;
                                break;
                            }
                        }
                        if (!alreadyIncluded) {
                            uniqueNeighbors[uniqueCount++] = neighbor_idx;
                        }
                    }
                }
            }
        }
    }

    for (int k = 0; k < uniqueCount; k++) {
        int neighbor_idx = uniqueNeighbors[k];

        sum += compute_patch_interaction(squares[neighbor_idx], sp, num_patches);
    }
    return sum;
}

void remove_square_from_cells(int squareIndex, const SquareParticle sp) {
    for (int ix = sp.cell_min_x; ix <= sp.cell_max_x; ix++) {
        for (int iy = sp.cell_min_y; iy <= sp.cell_max_y; iy++) {
            int cell = ix + Mx * iy;
            int prev = -1;
            int curr = head[cell];
            while (curr != -1) {
                if (nodePool[curr].squareIndex == squareIndex) {
                    if (prev == -1) {
                        head[cell] = nodePool[curr].next;
                    } else {
                        nodePool[prev].next = nodePool[curr].next;
                    }
                    nodePool[curr].next = freeList;
                    freeList = curr;
                    break;
                }
                prev = curr;
                curr = nodePool[curr].next;
            }
        }
    }
}

void rotate_square(SquareParticle *sp) {
    double current_angle = 2.0 * atan2(sp->q[0], sp->q[1]);
    double range = ((LFLAG != 0) && (rand() / (double) RAND_MAX < LRMOVE)) ? LRROTMAX : RROTMAX;

    double angle = rand() / (double) RAND_MAX * 2 * range - range;


    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0)
        current_angle += 2 * M_PI;

    sp->q[0] = sin(current_angle / 2);
    sp->q[1] = cos(current_angle / 2);
}

void move_square(SquareParticle *sp) {
    double range = ((LFLAG != 0) && (rand() / (double) RAND_MAX < LRMOVE)) ? LRDISPMAX : RDISPMAX;

    double dx = (rand() / (double) RAND_MAX * 2.0 - 1.0) * range;
    double dy = (rand() / (double) RAND_MAX * 2.0 - 1.0) * range;

    sp->x += dx;
    sp->y += dy;

    int periodic_x = (boundary_condition == PERIODIC_BOTH || boundary_condition == PERIODIC_X);
    int periodic_y = (boundary_condition == PERIODIC_BOTH || boundary_condition == PERIODIC_Y);

    sp->x = apply_boundary(sp->x, Lx, periodic_x);
    sp->y = apply_boundary(sp->y, Ly, periodic_y);
}


void animate_movement(int steps, int totalSquares, int num_patches) {
    FILE *f = fopen("data/square_animation.xyz", "w");
    if (!f) {
        printf("Error opening animation file\n");
        return;
    }

    FILE *energyFile = fopen("data/energy.dat", "w");

    if (!energyFile) {
        printf("Error opening energy file\n");
        exit(1);
    }

    double potential_energy = 0;

    for (int i = 0; i < totalSquares; i++) {
        for (int j = i + 1; j < totalSquares; j++) {
            potential_energy += compute_patch_interaction(squares[i], squares[j], num_patches);
        }
    }
    fprintf(energyFile, "%d %lf\n", 0, potential_energy);

    for (int step = 1; step <= steps; step++) {
        for (int attempt = 0; attempt < N; attempt++) {
            int idx = rand() % totalSquares;
            double r = (double) rand() / RAND_MAX;

            remove_square_from_cells(idx, squares[idx]);

            SquareParticle candidate = squares[idx];
            if (r < RMOVE)
                move_square(&candidate);
            else
                rotate_square(&candidate);

            update_square_AABB(&candidate);

            if (!is_overlapping_square(candidate, idx)) {
                double deltaE = compute_patch_energy(candidate, num_patches) - compute_patch_energy(
                                    squares[idx], num_patches);


                if (deltaE <= 0 || exp(-deltaE / KT) > (double) rand() / RAND_MAX) {
                    potential_energy += deltaE / 2;
                    squares[idx] = candidate;
                }
            }

            insert_square_in_cells(idx, squares[idx]);
        }

        fprintf(f, "%d\n", totalSquares + totalSquares * num_patches);
        fprintf(
            f,
            "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n",
            Lx, Ly);
        for (int i = 0; i < totalSquares; i++) {
            fprintf(f, "B %.6f %.6f %.6f  %.6f %.6f %.6f %.6f  %.6f %.6f %.6f\n",
                    squares[i].x,
                    squares[i].y,
                    Z_AXIS,
                    Z_AXIS,
                    Z_AXIS,
                    squares[i].q[0],
                    squares[i].q[1],
                    squares->side * OVITO_MUL,
                    squares->side * OVITO_MUL,
                    squares->side * pow(OVITO_MUL, 4));
        }

        for (int i = 0; i < totalSquares; i++) {
            for (int j = 0; j < num_patches; j++) {
                double global_x, global_y;
                compute_patch_global_position(squares[i], squares[i].patches[j], &global_x, &global_y);
                fprintf(f, "P %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                        global_x,
                        global_y,
                        Z_AXIS,
                        Z_AXIS,
                        Z_AXIS,
                        squares[i].q[0],
                        squares[i].q[1],
                        squares->patch_radius * OVITO_MUL,
                        squares->patch_radius * OVITO_MUL,
                        squares->patch_radius * pow(OVITO_MUL, 3));
            }
        }
        fprintf(energyFile, "%d %lf\n", step, potential_energy);
    }
    fclose(energyFile);
    fclose(f);
}
