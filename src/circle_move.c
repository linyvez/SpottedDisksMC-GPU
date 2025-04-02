#include "circle_move.h"

void rotate_circle(CircleParticle *c) {
    double current_angle = 2.0 * atan2(c->q[2], c->q[3]);
    double angle = ((rand() / (double)RAND_MAX) * 2 * MAX_ROTATION_ANGLE) - MAX_ROTATION_ANGLE;
    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0) {
        current_angle += 2 * M_PI;
    }
    c->q[2] = sin(current_angle / 2);
    c->q[3] = cos(current_angle / 2);
}

Coordinates calculate_new_coordinates(CircleParticle *c) {
    double dx = (rand() / (double)RAND_MAX) * 2 * SIGMA - SIGMA;
    double dy = (rand() / (double)RAND_MAX) * 2 * SIGMA - SIGMA;
    Coordinates coordinates;
    double x = c->x + dx;
    double y = c-> y + dy;
    if (boundary_condition == NO_BOUNDARY)
    {
        if (x < 0)
            x = 0;
        if (x > Lx)
            x = Lx;
        if (y < 0)
            y = 0;
        if (y > Ly)
            y = Ly;
    }
    else if (boundary_condition == PERIODIC_BOTH)
    {
        if (x < 0)
            x += Lx;
        if (x > Lx)
            x -= Lx;
        if (y < 0)
            y += Ly;
        if (y > Ly)
            y -= Ly;
    }
    else if (boundary_condition == PERIODIC_X)
    {
        if (x < 0)
            x += Lx;
        if (x > Lx)
            x -= Lx;
        if (y < 0)
            y = 0;
        if (y > Ly)
            y = Ly;
    }
    else if (boundary_condition == PERIODIC_Y)
    {
        if (x < 0)
            x = 0;
        if (x > Lx)
            x = Lx;
        if (y < 0)
            y += Ly;
        if (y > Ly)
            y -= Ly;
    }
    coordinates.x = x;
    coordinates.y = y;
    return coordinates;
}

double compute_circle_patch_interaction(const CircleParticle sp1, const CircleParticle sp2)
{
    double energy = 0.0;
    // replace 4 by number of patches
    for (int i = 0; i < PATCH_NUMBER; i++)
    {
        double gx1, gy1;
        compute_circle_patch_global_position(sp1, sp1.patches[i], &gx1, &gy1);
        for (int j = 0; j < PATCH_NUMBER; j++)
        {
            double gx2, gy2;
            compute_circle_patch_global_position(sp2, sp2.patches[j], &gx2, &gy2);
            double dx = gx1 - gx2;
            double dy = gy1 - gy2;
            double distance = sqrt(dx * dx + dy * dy);
            if (distance < 4 * PATCH_RADIUS_C)
            {
                energy += sp1.patches[i].strength;
            }
        }
    }
    return energy;
}

double compute_circle_patch_energy(const CircleParticle c, int idx) {
    double sum = 0.0;

    int uniqueNeighbors[64];
    int uniqueCount = 0;

    int cell_idx = get_cell_index(c.x, c.y);
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int neighbour_cell = cell_idx + dx * Mx_C + dy;
            // account for boundary conditions
            neighbour_cell = (neighbour_cell + NUM_CELLS_C) % NUM_CELLS_C;

            for (int i = 0; i < 4; i++) {
                if (parts_in_cells[neighbour_cell][i] == -1) break;

                int n_idx = parts_in_cells[neighbour_cell][i];
                if (n_idx == idx) continue;

                sum += compute_circle_patch_interaction(particles[n_idx], c);
            }
        }
    }
    // printf("%lf\n", sum);
    // for (int k = 0; k < uniqueCount; k++)
    // {
    //     int neighbor_idx = uniqueNeighbors[k];

    //     sum += compute_circle_patch_interaction(particles[neighbor_idx], c);
    // }
    
    return sum;
}


void animate_circle_movement(int steps, int totalCircles) {
    FILE *f = fopen("data/circle_animation.xyz", "w");
    if (!f)
    {
        printf("Error opening animation file\n");
        return;
    }

    FILE *energyFile = fopen("data/energy.dat", "w");

    if (!energyFile) {
        printf("Error opening energy file\n");
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
        for (int attempt = 0; attempt < N; attempt++) {
            int idx = rand() % totalCircles;
            double r = (double)rand() / RAND_MAX;

            CircleParticle candidate = particles[idx];
            Coordinates coordinates= calculate_new_coordinates(&candidate);
            
            double oldCx = candidate.x;
            double oldCy = candidate.y;
            double oldCq2 = candidate.q[2];
            double oldCq3 = candidate.q[3];
            double oldEnergy = compute_circle_patch_energy(particles[idx], idx);

            if (r < RMOVE) {
                move_particle(idx, coordinates.x, coordinates.y);
            }
            else
                rotate_circle(&particles[idx]);
            
            double deltaE = compute_circle_patch_energy(particles[idx], idx) - oldEnergy;
            if (deltaE != 0)
            if (deltaE <= 0 || exp(-deltaE / KT) > (double)rand() / RAND_MAX)
            {
                potential_energy += deltaE / 2;
            }
            else {
                move_particle(idx, oldCx, oldCy);
                particles[idx].q[2] = oldCq2;
                particles[idx].q[3] = oldCq3;
            }

            
        }

        fprintf(f, "%d\n", totalCircles + totalCircles * PATCH_NUMBER);
            fprintf(f, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n", Lx, Ly);
            for (int i = 0; i < totalCircles; i++)
            {
                fprintf(f, "B %.6f %.6f %.6f  %.6f %.6f %.6f %.6f  %.6f %.6f %.6f\n",
                        particles[i].x,
                        particles[i].y,
                        0.0,
                        particles[i].q[0],
                        particles[i].q[1],
                        particles[i].q[2],
                        particles[i].q[3],
                        SIGMA / 2.0,
                        SIGMA / 2.0,
                        SIGMA / 2.0);
            }

            for (int i = 0; i < totalCircles; i++)
            {
                for (int j = 0; j < PATCH_NUMBER; j++)
                {
                    double global_x, global_y;
                    compute_circle_patch_global_position(particles[i], particles[i].patches[j], &global_x, &global_y);
                    fprintf(f, "P %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                            global_x,
                            global_y,
                            0.0,
                            particles[i].q[0],
                            particles[i].q[1],
                            particles[i].q[2],
                            particles[i].q[3],
                            particles[i].patches[j].shape[0] / 2.0,
                            particles[i].patches[j].shape[1] / 2.0,
                            particles[i].patches[j].shape[2] / 2.0);
                }
            }
            fprintf(energyFile, "%d %lf\n", step, potential_energy);
        }
    fclose(energyFile);
    fclose(f);
}