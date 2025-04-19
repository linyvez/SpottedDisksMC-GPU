#include "square_config.h"
#include "periodic_boundary.h"
#include "general_config.h"

#include "config.h"
#include "patch.h"

#include <math.h>
#include <print_error.h>
#include <stdio.h>
#include <stdlib.h>


int nodeCount = 0;
int freeList = -1;

int *head = NULL;
Node *nodePool = NULL;
SquareParticle *squares = NULL;


static SquareParticle adjust_square_for_periodic(const SquareParticle ref, const SquareParticle sp) {
    SquareParticle sp_adj = sp;
    double dx = sp.x - ref.x;
    double dy = sp.y - ref.y;

    periodic_boundary(&dx, &dy);

    sp_adj.x = ref.x + dx;
    sp_adj.y = ref.y + dy;

    return sp_adj;
}

static void normalize(double *ax, double *ay) {
    const double len = sqrt((*ax) * (*ax) + (*ay) * (*ay));
    if (len > 1e-10) {
        *ax /= len;
        *ay /= len;
    }
}

static void project_polygon(const double corners[][2],
                            double ax, double ay,
                            double *minProj, double *maxProj) {
    const int num = 4;
    double dot = corners[0][0] * ax + corners[0][1] * ay;
    *minProj = *maxProj = dot;
    for (int i = 1; i < num; i++) {
        dot = corners[i][0] * ax + corners[i][1] * ay;
        if (dot < *minProj)
            *minProj = dot;
        if (dot > *maxProj)
            *maxProj = dot;
    }
}

static int intervals_overlap(double Amin, double Amax,
                             double Bmin, double Bmax) {
    if (Amax < Bmin)
        return 0;
    if (Bmax < Amin)
        return 0;
    return 1;
}

void compute_square_corners(const SquareParticle sp, double corners[4][2]) {
    double half = GL_CFG->particle_size * 0.5;
    double angle = 2.0 * atan2(sp.q[2], sp.q[3]);
    double cosA = cos(angle);
    double sinA = sin(angle);

    double local[4][2] = {
        {-half, half},
        {-half, -half},
        {half, -half},
        {half, half}
    };
    for (int i = 0; i < 4; i++) {
        corners[i][0] = sp.x + local[i][0] * cosA - local[i][1] * sinA;
        corners[i][1] = sp.y + local[i][0] * sinA + local[i][1] * cosA;
    }
}

static void get_square_axes(const double corners[4][2], double axes[2][2]) {
    double e0x = corners[1][0] - corners[0][0];
    double e0y = corners[1][1] - corners[0][1];
    axes[0][0] = -e0y;
    axes[0][1] = e0x;
    normalize(&axes[0][0], &axes[0][1]);

    double e1x = corners[2][0] - corners[1][0];
    double e1y = corners[2][1] - corners[1][1];
    axes[1][0] = -e1y;
    axes[1][1] = e1x;
    normalize(&axes[1][0], &axes[1][1]);
}

static int check_squares_overlap(const SquareParticle a, const SquareParticle b) {
    SquareParticle b_adj = adjust_square_for_periodic(a, b);
    double A[4][2], B[4][2];
    compute_square_corners(a, A);
    compute_square_corners(b_adj, B);

    double axesA[2][2];
    double axesB[2][2];
    get_square_axes(A, axesA);
    get_square_axes(B, axesB);

#define TEST_AXIS(ax, ay)                                \
    do                                                   \
    {                                                    \
        double Amin, Amax, Bmin, Bmax;                   \
        project_polygon(A, (ax), (ay), &Amin, &Amax); \
        project_polygon(B, (ax), (ay), &Bmin, &Bmax); \
        if (!intervals_overlap(Amin, Amax, Bmin, Bmax))  \
        {                                                \
            return 0;                                    \
        }                                                \
    } while (0)

    TEST_AXIS(axesA[0][0], axesA[0][1]);
    TEST_AXIS(axesA[1][0], axesA[1][1]);
    TEST_AXIS(axesB[0][0], axesB[0][1]);
    TEST_AXIS(axesB[1][0], axesB[1][1]);

    return 1;
}

void update_square_AABB(SquareParticle *sp) {
    double max_x, max_y;
    double corners[4][2];
    compute_square_corners(*sp, corners);
    double min_x = max_x = corners[0][0];
    double min_y = max_y = corners[0][1];
    for (int i = 1; i < 4; i++) {
        if (corners[i][0] < min_x)
            min_x = corners[i][0];
        if (corners[i][0] > max_x)
            max_x = corners[i][0];
        if (corners[i][1] < min_y)
            min_y = corners[i][1];
        if (corners[i][1] > max_y)
            max_y = corners[i][1];
    }

    sp->cell_min_x = (int) (min_x / CL_CFG->size);
    sp->cell_max_x = (int) (max_x / CL_CFG->size);
    sp->cell_min_y = (int) (min_y / CL_CFG->size);
    sp->cell_max_y = (int) (max_y / CL_CFG->size);

    if (sp->cell_min_x < 0)
        sp->cell_min_x = 0;
    if (sp->cell_min_y < 0)
        sp->cell_min_y = 0;
    if (sp->cell_max_x >= CL_CFG->Mx)
        sp->cell_max_x = CL_CFG->Mx - 1;
    if (sp->cell_max_y >= CL_CFG->My)
        sp->cell_max_y = CL_CFG->My - 1;
}


int is_overlapping_square(const SquareParticle sp, int candidate_index) {
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
                        int sqIdx = nodePool[nodeIndex].squareIndex;

                        if (sqIdx == candidate_index)
                            continue;
                        if (check_squares_overlap(squares[sqIdx], sp)) {
                            return 1;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

void insert_square_in_cells(int squareIndex, const SquareParticle sp) {
    for (int ix = sp.cell_min_x; ix <= sp.cell_max_x; ix++) {
        for (int iy = sp.cell_min_y; iy <= sp.cell_max_y; iy++) {
            int cell = ix + CL_CFG->Mx * iy;
            int n;
            if (freeList != -1) {
                n = freeList;
                freeList = nodePool[freeList].next;
            } else {
                n = nodeCount++;
                if (n >= CL_CFG->max_nodes) {
                    fprintf(stderr, "Error: nodePool overflow!\n");
                    exit(1);
                }
            }
            nodePool[n].squareIndex = squareIndex;
            nodePool[n].next = head[cell];
            head[cell] = n;
        }
    }
}

void compute_patch_global_position(const SquareParticle sp, const Patch patch, double *global_x, double *global_y) {
    double angle = 2.0 * atan2(sp.q[2], sp.q[3]);
    double cosA = cos(angle);
    double sinA = sin(angle);

    *global_x = sp.x + patch.rel_x * cosA - patch.rel_y * sinA;
    *global_y = sp.y + patch.rel_x * sinA + patch.rel_y * cosA;
}

int generate_random_squares() {
    head = malloc(CL_CFG->num_cells * sizeof(int));
    if (!head) {
        print_error(true, "Failed to allocate memory for head");
        exit(EXIT_FAILURE);
    }

    nodePool = malloc(CL_CFG->max_nodes * sizeof(Node));
    if (!nodePool) {
        print_error(true, "Failed to allocate memory for nodePool");
        exit(EXIT_FAILURE);
    }

    squares = malloc(GL_CFG->num_particles * sizeof(SquareParticle));
    if (!squares) {
        print_error(true, "Failed to allocate memory for squares");
        exit(EXIT_FAILURE);
    }

    visited = malloc(GL_CFG->num_particles * sizeof(int));
    if (!visited) {
        print_error(true, "Failed to allocate memory for visited");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < CL_CFG->num_cells; i++) {
        head[i] = -1;
    }
    nodeCount = 0;
    int totalSquares = 0;

    int attempts = 0;


    const double (*local)[2] = assign_patch_type();


    while (totalSquares < GL_CFG->num_particles) {
        if (attempts++ > MAX_ATTEMPTS) {
            fprintf(stderr, "Error: Too many attempts to place square.\n");
            break;
        }
        SquareParticle sp;
        sp.x = drand48() * GL_CFG->Lx;
        sp.y = drand48() * GL_CFG->Ly;

        double theta = drand48() * MAX_ANGLE;
        sp.q[0] = Z_AXIS;
        sp.q[1] = Z_AXIS;
        sp.q[2] = sin(theta / 2);
        sp.q[3] = cos(theta / 2);


        update_square_AABB(&sp);

        if (!is_overlapping_square(sp, -1)) {
            for (int i = 0; i < num_patches && local != NULL; i++) {
                sp.patches[i].rel_x = local[i][0];
                sp.patches[i].rel_y = local[i][1];
            }

            squares[totalSquares] = sp;
            insert_square_in_cells(totalSquares, sp);
            totalSquares++;
            attempts = 0;
        }
    }
    if (local != NULL) {
        free((void *) local);
    }

    FILE *f = fopen("data/configuration_square.xyz", "w");
    if (!f) {
        print_error(true, "Error opening configuration file\n");
        return 0;
    }
    fprintf(f, "%d\n", totalSquares + num_patches * totalSquares);
    fprintf(
        f,
        "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n",
        GL_CFG->Lx, GL_CFG->Ly);

    for (int i = 0; i < totalSquares; i++) {
        fprintf(f, "B %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                squares[i].x,
                squares[i].y,
                Z_AXIS,
                squares[i].q[0],
                squares[i].q[1],
                squares[i].q[2],
                squares[i].q[3],
                GL_CFG->particle_size * 0.5,
                GL_CFG->particle_size * 0.5,
                GL_CFG->particle_size * 0.125
        );
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
                    PH_CFG->radius * 0.125
            );
        }
    }
    fclose(f);
    return totalSquares;
}
