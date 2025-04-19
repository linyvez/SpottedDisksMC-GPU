#include "square_move.h"
#include "square_config.h"
#include "periodic_boundary.h"

#include <math.h>
#include <print_error.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "general_config.h"
#include "patch.h"


double compute_patch_interaction(const SquareParticle sp1, const SquareParticle sp2) {
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
            if (distance < 2 * PH_CFG->radius) {
                energy += PH_CFG->interaction;
            }
        }
    }
    return energy;
}

double compute_patch_energy(const SquareParticle sp) {
    double sum = 0.0;

    int epoch = current_epoch++;

    for (int ix = sp.cell_min_x; ix <= sp.cell_max_x; ix++) {
        for (int iy = sp.cell_min_y; iy <= sp.cell_max_y; iy++) {
            for (int ox = -1; ox <= 1; ox++) {
                for (int oy = -1; oy <= 1; oy++) {
                    int ghost_ix = (ix + ox + CL_CFG->Mx) % CL_CFG->Mx;
                    int ghost_iy = (iy + oy + CL_CFG->My) % CL_CFG->My;
                    int ghost_cell = ghost_ix + ghost_iy * CL_CFG->Mx;

                    if ((ghost_ix != ix || ghost_iy != iy) && (abs(ghost_ix - ix) <= 1 && abs(ghost_iy - iy) <= 1))
                        continue;

                    for (int nodeIndex = head[ghost_cell]; nodeIndex != -1; nodeIndex = nodePool[nodeIndex].next) {
                        int neighbor_idx = nodePool[nodeIndex].squareIndex;

                        if (visited[neighbor_idx] != epoch) {
                            sum += compute_patch_interaction(squares[neighbor_idx], sp);
                            visited[neighbor_idx] = epoch;
                        }
                    }
                }
            }
        }
    }

    return sum;
}


void remove_square_from_cells(int squareIndex, const SquareParticle sp) {
    for (int ix = sp.cell_min_x; ix <= sp.cell_max_x; ix++) {
        for (int iy = sp.cell_min_y; iy <= sp.cell_max_y; iy++) {
            int cell = ix + CL_CFG->Mx * iy;
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
    double current_angle = 2.0 * atan2(sp->q[2], sp->q[3]);
    double range = ((lflag != 0) && (drand48() < LRMOVE)) ? MV_CFG->lrrotmax : MV_CFG->rrotmax;

    double angle = drand48() * 2 * range - range;


    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0)
        current_angle += 2 * M_PI;

    sp->q[2] = sin(current_angle / 2);
    sp->q[3] = cos(current_angle / 2);
}

void move_square(SquareParticle *sp) {
    double range = ((lflag != 0) && (drand48() < LRMOVE)) ? MV_CFG->lrdispmax : MV_CFG->rdispmax;

    double dx = ((drand48() * 2.0) - 1.0) * range;
    double dy = ((drand48() * 2.0) - 1.0) * range;

    sp->x += dx;
    sp->y += dy;

    int periodic_x = (boundary_condition == PERIODIC_BOTH || boundary_condition == PERIODIC_X);
    int periodic_y = (boundary_condition == PERIODIC_BOTH || boundary_condition == PERIODIC_Y);

    sp->x = apply_boundary(sp->x, GL_CFG->Lx, periodic_x);
    sp->y = apply_boundary(sp->y, GL_CFG->Ly, periodic_y);
}


void animate_movement(int steps, int totalSquares) {
    FILE *f = fopen("data/square_animation.xyz", "w");
    if (!f) {
        print_error(true, "Error opening animation file\n");
        return;
    }

    FILE *energyFile = fopen("data/energy.dat", "w");

    if (!energyFile) {
        print_error(true, "Error opening energy file\n");
        exit(1);
    }

    double potential_energy = 0;

    for (int i = 0; i < totalSquares; i++) {
        for (int j = i + 1; j < totalSquares; j++) {
            potential_energy += compute_patch_interaction(squares[i], squares[j]);
        }
    }
    fprintf(energyFile, "%d %lf\n", 0, potential_energy);

    for (int step = 1; step <= steps; step++) {
        for (int attempt = 0; attempt < GL_CFG->num_particles; attempt++) {
            int idx = (int) (drand48() * totalSquares);
            double r = drand48();

            remove_square_from_cells(idx, squares[idx]);

            SquareParticle candidate = squares[idx];
            if (r < RMOVE)
                move_square(&candidate);
            else
                rotate_square(&candidate);

            update_square_AABB(&candidate);

            if (!is_overlapping_square(candidate, idx)) {
                double deltaE = compute_patch_energy(candidate) - compute_patch_energy(squares[idx]);


                if (deltaE <= 0 || exp(-deltaE / GL_CFG->temperature) > drand48()) {
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
            GL_CFG->Lx, GL_CFG->Ly);
        for (int i = 0; i < totalSquares; i++) {
            fprintf(f, "B %.6f %.6f %.6f  %.6f %.6f %.6f %.6f  %.6f %.6f %.6f\n",
                    squares[i].x,
                    squares[i].y,
                    Z_AXIS,
                    squares[i].q[0],
                    squares[i].q[1],
                    squares[i].q[2],
                    squares[i].q[3],
                    GL_CFG->particle_size * 0.5,
                    GL_CFG->particle_size * 0.5,
                    GL_CFG->particle_size * 0.0625);
        }

        for (int i = 0; i < totalSquares; i++) {
            for (int j = 0; j < num_patches; j++) {
                double global_x, global_y;
                compute_patch_global_position(squares[i], squares[i].patches[j], &global_x, &global_y);
                fprintf(f, "P %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                        global_x,
                        global_y,
                        Z_AXIS,
                        squares[i].q[0],
                        squares[i].q[1],
                        squares[i].q[2],
                        squares[i].q[3],
                        PH_CFG->radius * 0.5,
                        PH_CFG->radius * 0.5,
                        PH_CFG->radius * 0.0625);
            }
        }
        fprintf(energyFile, "%d %lf\n", step, potential_energy);
    }
    fclose(energyFile);
    fclose(f);

    free(head);
    free(nodePool);
    free(squares);
}
