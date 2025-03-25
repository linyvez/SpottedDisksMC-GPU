#include "cell_linked_algorithm.h"

Particle particles[N]; // all particles
int particles_idx[N]; // for each particle, store cell index
int parts_in_cells[NUM_CELLS][4]; // for each cell, store array of particle indices


void initialize_cells() {
    for (int i = 0; i < NUM_CELLS; i++) {
        for (int j = 0; j < 4; j++) {
            parts_in_cells[i][j] = -1;
        }
    }
}



int get_cell_index(double x, double y) {
    int ix = (int)(x / CELL_SIZE);
    int iy = (int)(y / CELL_SIZE);
    if (ix >= Mx) ix = Mx - 1;
    if (iy >= My) iy = My - 1;
    return ix + Mx * iy;
}



int is_overlapping_algo(Particle particle) {
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


int insert_particle_in_cell(int p_index, Particle *particle) {
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
    Particle *p = &particles[p_index];
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
    // printf("move: (%f, %f) --> (%f, %f)\n", oldX, oldY, newX, newY);
    // check
    if (is_overlapping_algo(*p)) {
        // reject
        p->x = oldX;
        p->y = oldY;
        // printf("rejected.\n");
    }
    else {
        // printf("accepted.\n");
    }
    insert_particle_in_cell(p_index, p);
}




void generate_random_algo() {
    initialize_cells();
    int count = 0;

    // for (int i = 0; i < NUM_CELLS; i++) {
    //     parts_in_cells[i] = -1;
    // }

    srand((unsigned)time(NULL));

    while (count < N) {
        Particle particle;
        particle.x = ((double)rand() / RAND_MAX) * Lx;
        particle.y = ((double)rand() / RAND_MAX) * Ly;

        if (!is_overlapping_algo(particle)) {
            particles[count] = particle;
            insert_particle_in_cell(count, &particle);
            count++;
        }
    }

    FILE *config = fopen("configuration_random_fast.pdb", "w");
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