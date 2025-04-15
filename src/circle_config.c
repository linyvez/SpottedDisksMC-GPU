#include "circle_config.h"


CircleParticle particles[N]; // all particles
int particles_idx[N]; // for each particle, store cell index
int parts_in_cells[NUM_CELLS_C][4]; // for each cell, store array of particle indices


void initialize_cells() {
    for (int i = 0; i < NUM_CELLS_C; i++) {
        for (int j = 0; j < 4; j++) {
            parts_in_cells[i][j] = -1;
        }
    }
}



int get_cell_index(double x, double y) {
    int ix = (int)(x / SIGMA);
    int iy = (int)(y / SIGMA);
    if (ix >= Mx_C) ix = Mx_C - 1;
    if (iy >= My_C) iy = My_C - 1;
    return ix + Mx_C * iy;
}



int is_overlapping_circle(CircleParticle particle) {
    int cell = get_cell_index(particle.x, particle.y);

    int ix = cell % Mx_C;
    int iy = cell / Mx_C;
    
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int nx = ix + dx;
            int ny = iy + dy;
            if (nx < 0 || nx >= Mx_C || ny < 0 || ny >= My_C)
                continue;
            int ncell = nx + Mx_C * ny;
            
            // for each particle in the cell, check overlap.
            for (int j = 0; j < 4; j++) {
                int p_i = parts_in_cells[ncell][j];
                if (p_i == -1) {
                    continue;
                };

                double x_diff = particle.x - particles[p_i].x;
                double y_diff = particle.y - particles[p_i].y;
                double distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                if (distance < SIGMA) {
                    return 1;
                }
            }
        }
    }

    return 0;
}


void insert_particle_in_cell(int p_index, CircleParticle *particle) {
    int cell = get_cell_index(particle->x, particle->y);

    // find free place in the cell
    for (int i = 0; i < 4; i++) {
        if (parts_in_cells[cell][i] == -1) {
            parts_in_cells[cell][i] = p_index;
            particles_idx[p_index] = cell;
            particles[p_index] = *particle;
            // printf("Inserted particle #%d (%f, %f) in cell #%d\n", p_index, particles[p_index].x,particles[p_index].y,particles_idx[p_index]);
            break;
        }
        // if (i == 3) printf("Failed to insert particle (%f, %f)\n", particle->x, particle->y);
    }
}


// remove -> check -> accept -> insert into new position
//                    reject -> insert into old position
void move_particle(int p_index, double newX, double newY) {
    CircleParticle *p = &particles[p_index];
    int oldCell = particles_idx[p_index];

    // remove
    for (int i = 0; i < 4; i++) {
        if (parts_in_cells[oldCell][i] == p_index) {
            parts_in_cells[oldCell][i] = -1;
            particles_idx[p_index] = -1;
            break;
        }
    }

    double oldX = p->x;
    double oldY = p->y;


    p->x = newX;
    p->y = newY;
    
    if (is_overlapping_circle(*p)) {
        p->x = oldX;
        p->y = oldY;
    }
    insert_particle_in_cell(p_index, p);
}

void compute_circle_patch_global_position(const CircleParticle sp, const Patch patch, double *global_x, double *global_y)
{

    double angle = 2.0 * atan2(sp.q[2], sp.q[3]);
    double cosA = cos(angle);
    double sinA = sin(angle);

    *global_x = sp.x + patch.rel_x * cosA - patch.rel_y * sinA;
    *global_y = sp.y + patch.rel_x * sinA + patch.rel_y * cosA;
}


int generate_random_circles() {
    initialize_cells();
    int count = 0;

    srand((unsigned)time(NULL));

    int attempts = 0;
    while (count < N) {
        if (attempts ++ > MAX_ATTEMPTS) {
            fprintf(stderr, "Error: Too many attempts to place a circle.\n");
            break;
        }
        CircleParticle particle;
        particle.x = ((double)rand() / RAND_MAX) * Lx;
        particle.y = ((double)rand() / RAND_MAX) * Ly;

        double theta = ((double)rand() / RAND_MAX) * MAX_ANGLE;
        particle.q[0] = 0.0;
        particle.q[1] = 0.0;
        particle.q[2] = sin(theta / 2);
        particle.q[3] = cos(theta / 2);

        if (!is_overlapping_circle(particle)) {
            const double radius = SIGMA / 2.0;
            const double angles[4] = {
                0,
                M_PI / 2,
                M_PI,
                3 * M_PI / 2
            };
            
            for (int i = 0; i < PATCH_NUMBER; i++) {
                particle.patches[i].rel_x = radius * cos(angles[i]);
                particle.patches[i].rel_y = radius * sin(angles[i]);
                particle.patches[i].shape[0] = PATCH_RADIUS_C;  
                particle.patches[i].shape[1] = PATCH_RADIUS_C;  
                particle.patches[i].shape[2] = PATCH_RADIUS_C;  
                particle.patches[i].strength = PATCH_STRENGTH;
            }


            particles[count] = particle;
            insert_particle_in_cell(count, &particle);
            count++;
            attempts = 0;
        }
    }

    FILE *f = fopen("data/configuration_circle.xyz", "w");
    if (!f) {
        printf("Error while opening configuration file\n");
        return 1;
    }
    printf("Number of circles laid: %d\n", count);

    fprintf(f, "%d\n", count + PATCH_NUMBER * count);
    fprintf(f, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n", Lx, Ly);

    for (int i = 0; i < count; i++)
    {
        fprintf(f, "B %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
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

    for (int i = 0; i < count; i++)
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
                    PATCH_RADIUS_C,
                    PATCH_RADIUS_C,
                    PATCH_RADIUS_C);
        }
    }

    fclose(f);
    return count;
}