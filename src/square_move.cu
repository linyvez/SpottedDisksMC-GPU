#include "square_move.cuh"
#include "curand_kernel.h"
#include <cuda_runtime_api.h>

#include <cmath>
#include <cstdio>
#include <patch.h>
#include <square_config.h>

#include "periodic_boundary.h"
#include "shared_utilities.h"
#include <random>
#include <algorithm>

__device__
inline int cell_index(double coord) {
    int raw = static_cast<int>(coord / CELL_SIZE);

    if (raw >= Mx) {
        raw -= 1;
    }
    return raw;
}

__device__ static
void d_periodic_boundary(double *dx, double *dy) {
    if ((d_bc & PERIODIC_X) != 0) {
        *dx = fmod(*dx + 0.5 * Lx, Lx);
        if (*dx < 0) *dx += Lx;
        *dx -= 0.5 * Lx;
    }
    if ((d_bc & PERIODIC_Y) != 0) {
        *dy = fmod(*dy + 0.5 * Ly, Ly);
        if (*dy < 0) *dy += Ly;
        *dy -= 0.5 * Ly;
    }
}

__device__ static SquareParticle d_adjust_square_for_periodic(const SquareParticle &ref, const SquareParticle &sp) {
    SquareParticle sp_adj = sp;
    double dx = sp.x - ref.x;
    double dy = sp.y - ref.y;

    d_periodic_boundary(&dx, &dy);

    sp_adj.x = ref.x + dx;
    sp_adj.y = ref.y + dy;

    return sp_adj;
}

__device__ static void d_normalize(double *ax, double *ay) {
    const double len = sqrt((*ax) * (*ax) + (*ay) * (*ay));
    if (len > 1e-10) {
        *ax /= len;
        *ay /= len;
    }
}

__device__ static void d_project_polygon(const double corners[][2],
                                         double ax, double ay,
                                         double *minProj, double *maxProj) {
    constexpr int num = 4;
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

__device__ static int d_intervals_overlap(double Amin, double Amax,
                                          double Bmin, double Bmax) {
    if (Amax < Bmin)
        return 0;
    if (Bmax < Amin)
        return 0;
    return 1;
}

__device__ static void d_compute_square_corners(const SquareParticle &sp, double corners[4][2]) {
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

__device__ static void d_get_square_axes(const double corners[4][2], double axes[2][2]) {
    double e0x = corners[1][0] - corners[0][0];
    double e0y = corners[1][1] - corners[0][1];
    axes[0][0] = -e0y;
    axes[0][1] = e0x;
    d_normalize(&axes[0][0], &axes[0][1]);

    double e1x = corners[2][0] - corners[1][0];
    double e1y = corners[2][1] - corners[1][1];
    axes[1][0] = -e1y;
    axes[1][1] = e1x;
    d_normalize(&axes[1][0], &axes[1][1]);
}

__device__ static int check_squares_overlap(const SquareParticle &a, const SquareParticle &b) {
    SquareParticle b_adj = d_adjust_square_for_periodic(a, b);
    double A[4][2], B[4][2];
    d_compute_square_corners(a, A);
    d_compute_square_corners(b_adj, B);

    double axesA[2][2];
    double axesB[2][2];
    d_get_square_axes(A, axesA);
    d_get_square_axes(B, axesB);

#define TEST_AXIS(ax, ay)                                \
     do                                                   \
     {                                                    \
         double Amin, Amax, Bmin, Bmax;                   \
         d_project_polygon(A, (ax), (ay), &Amin, &Amax); \
         d_project_polygon(B, (ax), (ay), &Bmin, &Bmax); \
         if (!d_intervals_overlap(Amin, Amax, Bmin, Bmax))  \
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


__device__ int static d_is_overlapping_square(const SquareParticle &sp, const SquareParticle *squares,
                                              LinkedCell cell_struct, int squareIdx) {
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

                if (nidx == squareIdx) {
                    continue;
                }

                if (check_squares_overlap(squares[nidx], sp)) {
                    return 1;
                }
            }
        }
    }
    return 0;
}


__device__ double compute_patch_interaction(
    const SquareParticle &sp1,
    const SquareParticle &sp2) {
    SquareParticle sp2_adj = d_adjust_square_for_periodic(sp1, sp2);

    double angle1 = 2.0 * atan2(sp1.q[0], sp1.q[1]);
    double cos1 = cos(angle1), sin1 = sin(angle1);
    double angle2 = 2.0 * atan2(sp2_adj.q[0], sp2_adj.q[1]);
    double cos2 = cos(angle2), sin2 = sin(angle2);

    int n = d_num_patches;
    double gx1[MAX_PATCHES], gy1[MAX_PATCHES];
    double gx2[MAX_PATCHES], gy2[MAX_PATCHES];

    for (int i = 0; i < n; ++i) {
        double rx = d_patch[i][0], ry = d_patch[i][1];
        gx1[i] = sp1.x + rx * cos1 - ry * sin1;
        gy1[i] = sp1.y + rx * sin1 + ry * cos1;
        gx2[i] = sp2_adj.x + rx * cos2 - ry * sin2;
        gy2[i] = sp2_adj.y + rx * sin2 + ry * cos2;
    }

    double energy = 0.0;
    constexpr double R2 = (2.0 * PATCH_RADIUS) * (2.0 * PATCH_RADIUS);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dx = gx1[i] - gx2[j];
            double dy = gy1[i] - gy2[j];
            d_periodic_boundary(&dx, &dy);
            double dist2 = dx * dx + dy * dy;
            if (dist2 < R2) {
                energy += PATCH_STRENGTH;
            }
        }
    }

    return energy;
}


__device__ static double compute_patch_energy(const SquareParticle &sp, const SquareParticle *squares,
                                              LinkedCell cell_struct, int squareIdx) {
    double sum = 0.0;

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

                if (nidx == squareIdx) {
                    continue;
                }

                sum += compute_patch_interaction(squares[nidx], sp);
            }
        }
    }
    return sum;
}




__device__ static void insert_square_in_cells(int squareIndex,
                                              const SquareParticle &sp,
                                              LinkedCell cell_struct) {
    int ix = cell_index(sp.x);
    int iy = cell_index(sp.y);
    int cell = ix + iy * Mx;

    int old_cnt = atomicAdd(&cell_struct.count[cell], 1);

    cell_struct.neighbors[cell][old_cnt] = squareIndex;
}


__device__ static void remove_square_from_cells(int squareIndex,
                                                const SquareParticle &sp,
                                                LinkedCell cell_struct) {
    int ix = cell_index(sp.x);
    int iy = cell_index(sp.y);
    int cell = ix + iy * Mx;

    int cnt = cell_struct.count[cell];
    int pos = -1;
    for (int s = 0; s < cnt; ++s) {
        if (cell_struct.neighbors[cell][s] == squareIndex) {
            pos = s;
            break;
        }
    }
    if (pos < 0) return;

    int old_cnt = atomicSub(&cell_struct.count[cell], 1);
    int last = old_cnt - 1;

    if (pos != last) {
        atomicExch(&cell_struct.neighbors[cell][pos],
                   cell_struct.neighbors[cell][last]);
    }
}


__device__ double apply_boundary(double coord, double L) {
    if (d_bc) {
        if (coord < 0)
            coord += L;
        else if (coord > L)
            coord -= L;
    } else {
        if (coord < 0)
            coord = 0;
        else if (coord > L)
            coord = L;
    }
    return coord;
}

__device__ static void rotate_square(SquareParticle *sp, curandState *rng) {
    double current_angle = 2.0 * atan2(sp->q[0], sp->q[1]);
    float u = curand_uniform(rng);

    double range = d_lflag != 0 && u < LRMOVE ? LRROTMAX : RROTMAX;

    float v = curand_uniform(rng);

    double angle = v * 2 * range - range;


    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0)
        current_angle += 2 * M_PI;

    sp->q[0] = sin(current_angle / 2);
    sp->q[1] = cos(current_angle / 2);
}

__device__ static void move_square(SquareParticle *sp, curandState *rng) {
    double range = d_lflag != 0 && (curand_uniform(rng) < LRMOVE) ? LRDISPMAX : RDISPMAX;

    double dx = (curand_uniform(rng) * 2.0 - 1.0) * range;
    double dy = (curand_uniform(rng) * 2.0 - 1.0) * range;

    sp->x += dx;
    sp->y += dy;

    sp->x = apply_boundary(sp->x, Lx);
    sp->y = apply_boundary(sp->y, Ly);
}


static void shuffle(int *array) {
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::shuffle(array, array + 4, gen);
}


__global__ void initRNG(curandState *states, unsigned long seed) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;

    curand_init(
        seed,
        tid,
        0,
        &states[tid]
    );
}



__global__ void mc_step_kernel(int color,
                               curandState *rngStates,
                               SquareParticle *squares,
                               LinkedCell cell_struct) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;

    unsigned cellX = tid % Mx;
    unsigned cellY = (tid / Mx) % My;
    unsigned myColor = (cellX & 1) << 1 | (cellY & 1);


    if (tid >= NUM_CELLS || myColor != color) {
        return;
    }

    curandState localState = rngStates[tid];


    int cnt = cell_struct.count[tid];
    if (cnt == 0) return;

    int r = static_cast<int>(curand_uniform_double(&localState) * cnt);
    if (r == cnt) r = cnt - 1;

    int squareIdx = cell_struct.neighbors[tid][r];


    SquareParticle candidate = squares[squareIdx];

    float rm = curand_uniform(&localState);

    if (rm < RMOVE) {
        move_square(&candidate, &localState);
    } else {
        rotate_square(&candidate, &localState);
    }

    if (!d_is_overlapping_square(candidate, squares, cell_struct, squareIdx)) {
        double myDE = compute_patch_energy(candidate, squares, cell_struct, squareIdx) - compute_patch_energy(
                          squares[squareIdx], squares, cell_struct, squareIdx);


        if (myDE <= 0 || exp(-myDE / KT) > curand_uniform(&localState)) {
            int new_cell = cell_index(candidate.x) + cell_index(candidate.y) * Mx;
            if (new_cell != tid) {
                remove_square_from_cells(squareIdx, squares[squareIdx], cell_struct);
                insert_square_in_cells(squareIdx, candidate, cell_struct);
            }

            squares[squareIdx] = candidate;
        }
    }
    rngStates[tid] = localState;
}

__constant__ int d_num_patches;
__constant__ BoundaryCondition d_bc;
__constant__ int d_lflag;
__constant__ double d_patch[MAX_PATCHES][2];


void animate_movement(const int totalSquares, const SquareParticle *squares, LinkedCell cell_struct, int lflag) {
    cudaMemcpyToSymbol(d_bc, &h_bc, sizeof(h_bc));
    cudaMemcpyToSymbol(d_num_patches, &h_num_patches, sizeof(h_num_patches));
    cudaMemcpyToSymbol(d_lflag, &lflag, sizeof(lflag));
    cudaMemcpyToSymbol(d_patch, h_patch, h_num_patches * 2 * sizeof(double));


    SquareParticle *d_squares;

    SquareParticle *d_snapshots = nullptr;
    size_t snapshotCount = static_cast<size_t>(ANIMATION_STEPS) * totalSquares;
    cudaMalloc(&d_snapshots, snapshotCount * sizeof(SquareParticle));

    cudaMalloc(&d_squares, totalSquares * sizeof(SquareParticle));
    cudaMemcpy(d_squares, squares, totalSquares * sizeof(SquareParticle), cudaMemcpyHostToDevice);

    int *d_count_mem;
    int *d_neigh_mem;
    LinkedCell d_cell_struct;


    cudaMalloc(&d_count_mem, sizeof(int) * NUM_CELLS);
    cudaMalloc(&d_neigh_mem, sizeof(int) * NUM_CELLS * MAX_NEIGH);

    cudaMemcpy(d_count_mem,
               cell_struct.count,
               sizeof(int) * NUM_CELLS,
               cudaMemcpyHostToDevice);

    cudaMemcpy(
        d_neigh_mem,
        &cell_struct.neighbors[0][0],
        sizeof(int) * NUM_CELLS * MAX_NEIGH,
        cudaMemcpyHostToDevice
    );

    d_cell_struct.count = d_count_mem;
    d_cell_struct.neighbors = reinterpret_cast<int (*)[MAX_NEIGH]>(d_neigh_mem);


    int threads = 32;
    int blocks = (NUM_CELLS + threads - 1) / threads;

    curandState *d_rngStates;
    size_t nThreads = blocks * threads;
    cudaMalloc(&d_rngStates, nThreads * sizeof(curandState));


    int colors[4] = {0, 1, 2, 3};


    initRNG<<<blocks, threads>>>(d_rngStates, 1234UL);
    cudaDeviceSynchronize();

    const clock_t start = clock();

    for (int step = 0; step < ANIMATION_STEPS; step++) {
        shuffle(colors);

        for (int color: colors) {
            mc_step_kernel<<<blocks,threads>>>(
                color,
                d_rngStates,
                d_squares,
                d_cell_struct
            );
        }
        SquareParticle *dst = d_snapshots + static_cast<size_t>(step) * totalSquares;

        cudaMemcpy(dst,
                   d_squares,
                   totalSquares * sizeof(SquareParticle),
                   cudaMemcpyDeviceToDevice);
    }

    const clock_t end = clock();

    double elapsed_ms = (double)(end - start) * 1000.0 / (double)CLOCKS_PER_SEC;
    printf("Elapsed time: %.3f ms\n", elapsed_ms);

    SquareParticle *h_snapshots = nullptr;
    cudaMallocHost(&h_snapshots, snapshotCount * sizeof(SquareParticle));

    cudaMemcpy(h_snapshots,
               d_snapshots,
               snapshotCount * sizeof(SquareParticle),
               cudaMemcpyDeviceToHost);


    FILE *f = fopen("data/square_animation.xyz", "w");
    if (!f) {
        printf("Error opening animation file\n");
        return;
    }

    for (int step = 0; step < ANIMATION_STEPS; ++step) {
        SquareParticle *frame = h_snapshots + static_cast<size_t>(step) * totalSquares;
        write_file(f, frame, totalSquares);
    }


    cudaFree(d_squares);
    cudaFree(d_count_mem);
    cudaFree(d_neigh_mem);
    cudaFree(d_rngStates);
    cudaFreeHost(h_snapshots);
    cudaFree(d_snapshots);

    fclose(f);
}
