#include "circle_move.h"

#include <config.h>
#include <patch.h>
#include <periodic_boundary.h>
#include <print_error.h>
#include <stdio.h>
#include <_stdlib.h>


void rotate_circle(CircleParticle *c) {
    double current_angle = 2.0 * atan2(c->q[2], c->q[3]);
    double range = ((lflag != 0) && (drand48() < LRMOVE)) ? MV_CFG->lrrotmax : MV_CFG->rrotmax;

    double angle = drand48() * 2 * range - range;


    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0)
        current_angle += 2 * M_PI;

    c->q[2] = sin(current_angle / 2);
    c->q[3] = cos(current_angle / 2);
}

Coordinates calculate_new_coordinates(const CircleParticle *c) {
    Coordinates coordinates;
    double range = ((lflag != 0) && (drand48() < LRMOVE)) ? MV_CFG->lrdispmax : MV_CFG->rdispmax;

    double dx = ((drand48() * 2.0) - 1.0) * range;
    double dy = ((drand48() * 2.0) - 1.0) * range;

    double x = c->x + dx;
    double y = c->y + dy;

    int periodic_x = (boundary_condition == PERIODIC_BOTH || boundary_condition == PERIODIC_X);
    int periodic_y = (boundary_condition == PERIODIC_BOTH || boundary_condition == PERIODIC_Y);

    x = apply_boundary(x, GL_CFG->Lx, periodic_x);
    y = apply_boundary(y, GL_CFG->Ly, periodic_y);
    coordinates.x = x;
    coordinates.y = y;
    return coordinates;
}

double compute_circle_patch_interaction(const CircleParticle sp1, const CircleParticle sp2) {
    double energy = 0.0;
    // replace 4 by number of patches
    for (int i = 0; i < num_patches; i++) {
        double gx1, gy1;
        compute_circle_patch_global_position(sp1, sp1.patches[i], &gx1, &gy1);
        for (int j = 0; j < num_patches; j++) {
            double gx2, gy2;
            compute_circle_patch_global_position(sp2, sp2.patches[j], &gx2, &gy2);
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

double compute_circle_patch_energy(const CircleParticle c, int idx) {
    double sum = 0.0;

    int epoch = current_epoch++;


    int cell_idx = get_cell_index(c.x, c.y);
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int neighbour_cell = cell_idx + dx * CL_CFG->Mx + dy;
            // account for boundary conditions
            neighbour_cell = (neighbour_cell + CL_CFG->num_cells) % CL_CFG->num_cells;

            for (int i = 0; i < 4; i++) {
                if (parts_in_cells[neighbour_cell][i] == -1) break;

                int n_idx = parts_in_cells[neighbour_cell][i];
                if (n_idx == idx) continue;

                if (visited[n_idx] != epoch) {
                    sum += compute_circle_patch_interaction(particles[n_idx], c);
                    visited[n_idx] = epoch;
                }
            }
        }
    }

    return sum;
}


void animate_circle_movement(int steps, int totalCircles) {
    FILE *f = fopen("data/circle_animation.xyz", "w");
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
    for (int i = 0; i < totalCircles; i++) {
        for (int j = i + 1; j < totalCircles; j++) {
            potential_energy += compute_circle_patch_interaction(particles[i], particles[j]);
        }
    }
    fprintf(energyFile, "%d %lf\n", 0, potential_energy);

    for (int step = 1; step <= steps; step++) {
        for (int attempt = 0; attempt < GL_CFG->num_particles; attempt++) {
            int idx = (int) (drand48() * totalCircles);
            double r = drand48();

            CircleParticle candidate = particles[idx];
            Coordinates coordinates = calculate_new_coordinates(&candidate);

            double oldCx = candidate.x;
            double oldCy = candidate.y;
            double oldCq2 = candidate.q[2];
            double oldCq3 = candidate.q[3];
            double oldEnergy = compute_circle_patch_energy(particles[idx], idx);

            if (r < RMOVE) {
                move_particle(idx, coordinates.x, coordinates.y);
            } else
                rotate_circle(&particles[idx]);

            double deltaE = compute_circle_patch_energy(particles[idx], idx) - oldEnergy;
            if (deltaE != 0) {
                if (deltaE <= 0 || exp(-deltaE / GL_CFG->temperature) > drand48()) {
                    potential_energy += deltaE / 2;
                } else {
                    move_particle(idx, oldCx, oldCy);
                    particles[idx].q[2] = oldCq2;
                    particles[idx].q[3] = oldCq3;
                }
            }
        }

        fprintf(f, "%d\n", totalCircles + totalCircles * num_patches);
        fprintf(
            f,
            "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n",
            GL_CFG->Lx, GL_CFG->Ly);
        for (int i = 0; i < totalCircles; i++) {
            fprintf(f, "B %.6f %.6f %.6f  %.6f %.6f %.6f %.6f  %.6f %.6f %.6f\n",
                    particles[i].x,
                    particles[i].y,
                    Z_AXIS,
                    particles[i].q[0],
                    particles[i].q[1],
                    particles[i].q[2],
                    particles[i].q[3],
                    GL_CFG->particle_size * 0.5,
                    GL_CFG->particle_size * 0.5,
                    GL_CFG->particle_size * 0.125);
        }

        for (int i = 0; i < totalCircles; i++) {
            for (int j = 0; j < num_patches; j++) {
                double global_x, global_y;
                compute_circle_patch_global_position(particles[i], particles[i].patches[j], &global_x, &global_y);
                fprintf(f, "P %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                        global_x,
                        global_y,
                        Z_AXIS,
                        particles[i].q[0],
                        particles[i].q[1],
                        particles[i].q[2],
                        particles[i].q[3],
                        PH_CFG->radius * 0.5,
                        PH_CFG->radius * 0.5,
                        PH_CFG->radius * 0.125);
            }
        }
        fprintf(energyFile, "%d %lf\n", step, potential_energy);
    }
    fclose(energyFile);
    fclose(f);

    free(particles);
    free(particles_idx);
    free(parts_in_cells);
}
