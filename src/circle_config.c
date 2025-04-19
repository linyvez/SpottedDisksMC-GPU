#include "circle_config.h"

#include <periodic_boundary.h>
#include <stdio.h>

#include "patch.h"
#include "print_error.h"
#include <_stdlib.h>

#include "config.h"





CircleParticle *particles = NULL;
int *particles_idx = NULL;
int (*parts_in_cells)[MAX_NEIGHBOURS] = NULL;


static CircleParticle adjust_circles_for_periodic(const CircleParticle ref, const CircleParticle sp) {
    CircleParticle sp_adj = sp;
    double dx = sp.x - ref.x;
    double dy = sp.y - ref.y;

    periodic_boundary(&dx, &dy);

    sp_adj.x = ref.x + dx;
    sp_adj.y = ref.y + dy;

    return sp_adj;
}

void initialize_cells() {
    for (int i = 0; i < CL_CFG->num_cells; i++) {
        for (int j = 0; j < 4; j++) {
            parts_in_cells[i][j] = -1;
        }
    }
}



int get_cell_index(double x, double y) {
    int ix = (int)(x / GL_CFG->particle_size);
    int iy = (int)(y / GL_CFG->particle_size);
    if (ix >= CL_CFG->Mx) ix = CL_CFG->Mx - 1;
    if (iy >= CL_CFG->My) iy = CL_CFG->My - 1;
    return ix + CL_CFG->Mx * iy;
}


int is_overlapping_circle(CircleParticle particle) {
    int cell = get_cell_index(particle.x, particle.y);

    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int neighbour_cell = cell + dx * CL_CFG->Mx + dy;
            neighbour_cell = (neighbour_cell + CL_CFG->num_cells) % CL_CFG->num_cells;

            for (int j = 0; j < 4; j++) {
                int p_i = parts_in_cells[neighbour_cell][j];
                if (p_i == -1) {
                    continue;
                };
                CircleParticle adj_p_i = adjust_circles_for_periodic(particle, particles[p_i]);



                double x_diff = particle.x - adj_p_i.x;
                double y_diff = particle.y - adj_p_i.y;
                double distance = sqrt(x_diff * x_diff + y_diff * y_diff);
                if (distance < GL_CFG->particle_size) {
                    return 1;
                }
            }
        }
    }

    return 0;
}


void insert_particle_in_cell(int p_index, const CircleParticle *particle) {
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
    particles_idx = malloc( GL_CFG->num_particles * sizeof(int));
    if (!particles_idx) {
        print_error(true, "Failed to allocate memory for particles_idx");
        exit(EXIT_FAILURE);
    }

    parts_in_cells = malloc(CL_CFG->num_cells * sizeof(int[MAX_NEIGHBOURS]));
    if (!parts_in_cells) {
        print_error(true, "Failed to allocate memory for part_in_cells");
        exit(EXIT_FAILURE);
    }

    particles = malloc(GL_CFG->num_particles * sizeof(CircleParticle));
    if (!particles) {
        print_error(true, "Failed to allocate memory for particles");
        exit(EXIT_FAILURE);
    }

    visited = malloc(GL_CFG->num_particles * sizeof(int));
    if (!visited) {
        print_error(true, "Failed to allocate memory for visited");
        exit(EXIT_FAILURE);
    }


    initialize_cells();
    int count = 0;
    const double (*local)[2] = assign_patch_type();

    int attempts = 0;
    while (count < GL_CFG->num_particles) {
        if (attempts ++ > MAX_ATTEMPTS) {
            print_error(false, "Error: Too many attempts to place a circle.\n");
            break;
        }
        CircleParticle particle;
        particle.x = drand48() * GL_CFG->Lx;
        particle.y = drand48() * GL_CFG->Ly;




        double theta = drand48() * MAX_ANGLE;
        particle.q[0] = Z_AXIS;
        particle.q[1] = Z_AXIS;
        particle.q[2] = sin(theta / 2);
        particle.q[3] = cos(theta / 2);


        if (!is_overlapping_circle(particle)) {

            for (int i = 0; i < num_patches && local != NULL; i++) {

                particle.patches[i].rel_x = local[i][0];
                particle.patches[i].rel_y = local[i][1];
            }


            particles[count] = particle;
            insert_particle_in_cell(count, &particle);
            count++;
            attempts = 0;
        }
    }

    if (local != NULL) {
        free((void *)local);
    }

    FILE *f = fopen("data/configuration_circle.xyz", "w");
    if (!f) {
        print_error(true, "Error while opening configuration file\n");
        return 1;
    }
    printf("Number of circles laid: %d\n", count);

    fprintf(f, "%d\n", count + num_patches * count);
    fprintf(f, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n", GL_CFG->Lx, GL_CFG->Ly);

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
                GL_CFG->particle_size * 0.5,
                GL_CFG->particle_size * 0.5,
                GL_CFG->particle_size * 0.125);
    }

    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < num_patches; j++)
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
                    PH_CFG->radius * 0.5,
                    PH_CFG->radius * 0.5,
                    PH_CFG->radius * 0.125);
        }
    }

    fclose(f);
    return count;
}