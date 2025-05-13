#include "square_config.h"

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

static void compute_square_corners(const SquareParticle sp, double corners[4][2]) {
    double half = PARTICLE_SIZE / 2.0;
    double angle = 2.0 * atan2(sp.q[0], sp.q[1]);
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


static int is_overlapping_square(const SquareParticle sp, const SquareParticle * squares, LinkedCell cell_struct) {
    int cell_ix = cell_index(sp.x);
    int cell_iy = cell_index(sp.y);

    for (int ox = -1; ox <= 1; ox++) {
        for (int oy = -1; oy <= 1; oy++) {
            int ghost_ix = (cell_ix + ox + Mx) % Mx;
            int ghost_iy = (cell_iy + oy + My) % My;
            int ghost_cell = ghost_ix + ghost_iy * Mx;

            int cnt = cell_struct.count[ghost_cell];
            for (int i = 0; i < cnt; i++) {
                int nidx = cell_struct.neighbors[ghost_cell][i];
                if (check_squares_overlap(squares[nidx], sp)) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

static void insert_square_in_cells(int squareIndex, const SquareParticle sp, LinkedCell cell_struct) {
    int cell = cell_index(sp.x) + cell_index(sp.y) * Mx;
    int cnt = cell_struct.count[cell];
    cell_struct.neighbors[cell][cnt] = squareIndex;
    cell_struct.count[cell]++;
}



int generate_random_squares(int patch, SquareParticle * squares, LinkedCell cell_struct) {
    int totalSquares = 0;
    int attempts = 0;

    h_patch = assign_patch_type(patch);

    while (totalSquares < N_PAR) {
        if (attempts++ > MAX_ATTEMPTS) {
            fprintf(stderr, "Error: Too many attempts to place square.\n");
            break;
        }
        SquareParticle sp;
        sp.x = (double) rand() / RAND_MAX * Lx;
        sp.y = (double) rand() / RAND_MAX * Ly;

        double theta = (double) rand() / RAND_MAX * MAX_ANGLE;

        sp.q[0] = sin(theta / 2);
        sp.q[1] = cos(theta / 2);

        if (!is_overlapping_square(sp, squares, cell_struct)) {
            squares[totalSquares] = sp;
            insert_square_in_cells(totalSquares, sp, cell_struct);
            totalSquares++;
            attempts = 0;
        }
    }

    FILE *f = fopen("data/configuration_square.xyz", "w");
    if (!f) {
        printf("Error opening configuration file\n");
        return 0;
    }
    write_file(f, squares, totalSquares);
    fclose(f);
    return totalSquares;
}
