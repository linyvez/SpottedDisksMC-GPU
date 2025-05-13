
#include "circle_move.cuh"
#include "curand_kernel.h"
#include <cuda_runtime_api.h>

#include <cmath>
#include <cstdio>
#include <patch.h>
#include <circle_config.h>

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
    if ((dc_bc & PERIODIC_X) != 0) {
        *dx = fmod(*dx + 0.5 * Lx, Lx);
        if (*dx < 0) *dx += Lx;
        *dx -= 0.5 * Lx;
    }
    if ((dc_bc & PERIODIC_Y) != 0) {
        *dy = fmod(*dy + 0.5 * Ly, Ly);
        if (*dy < 0) *dy += Ly;
        *dy -= 0.5 * Ly;
    }
}

__device__ static CircleParticle adjust_circle_for_periodic(const CircleParticle ref, const CircleParticle cp) {
    CircleParticle sp_adj = cp;
    double dx = cp.x - ref.x;
    double dy = cp.y - ref.y;

    d_periodic_boundary(&dx, &dy);

    sp_adj.x = ref.x + dx;
    sp_adj.y = ref.y + dy;

    return sp_adj;
}


__device__ static int check_circles_overlap(CircleParticle a, CircleParticle b) {
    CircleParticle b_adj = adjust_circle_for_periodic(a, b);

    double x_diff = a.x - b_adj.x;
    double y_diff = a.y - b_adj.y;

    double distance = x_diff * x_diff + y_diff * y_diff;

    if (distance < PARTICLE_SIZE * PARTICLE_SIZE) {
        return 1;
    }
    return 0;
}


__device__ static int d_is_overlapping_circle(const CircleParticle cp, CircleParticle * circles, LinkedCell cell_struct, int Idx) {
    int cell_ix = cell_index(cp.x);
    int cell_iy = cell_index(cp.y);

    for (int ox = -1; ox <= 1; ox++) {
        for (int oy = -1; oy <= 1; oy++) {
            int ghost_ix = (cell_ix + ox + Mx) % Mx;
            int ghost_iy = (cell_iy + oy + My) % My;
            int ghost_cell = ghost_ix + ghost_iy * Mx;

            int cnt = cell_struct.count[ghost_cell];
            for (int i = 0; i < cnt; i++) {
                int nidx = cell_struct.neighbors[ghost_cell][i];

                if (Idx == nidx) {
                    continue;
                }
                if (check_circles_overlap(circles[nidx], cp)) {
                    return 1;
                }
            }
        }
    }
    return 0;
}


__device__ static double compute_patch_interaction(
    const CircleParticle &a,
    const CircleParticle &b)
{
    CircleParticle b_adj = adjust_circle_for_periodic(a, b);

    double thetaA = 2.0 * atan2(a.q[0], a.q[1]);
    double thetaB = 2.0 * atan2(b_adj.q[0], b_adj.q[1]);
    double cA = cos(thetaA), sA = sin(thetaA);
    double cB = cos(thetaB), sB = sin(thetaB);

    double energy = 0.0;
    constexpr double cutoff2 =  PARTICLE_SIZE * PARTICLE_SIZE;

    for (int i = 0; i < dc_num_patches; ++i) {
        double rx1 = dc_patch[i][0], ry1 = dc_patch[i][1];
        double xi  = a.x + (rx1 * cA - ry1 * sA);
        double yi  = a.y + (rx1 * sA + ry1 * cA);

        for (int j = 0; j < dc_num_patches; ++j) {
            double rx2 = dc_patch[j][0], ry2 = dc_patch[j][1];
            double xj  = b_adj.x + (rx2 * cB - ry2 * sB);
            double yj  = b_adj.y + (rx2 * sB + ry2 * cB);

            double dx = xi - xj;
            double dy = yi - yj;
            d_periodic_boundary(&dx, &dy);

            if (dx*dx + dy*dy < cutoff2) {
                energy += PATCH_STRENGTH;
            }
        }
    }

    return energy;
}



__device__ static double compute_patch_energy(
    const CircleParticle &cp,
    const CircleParticle *circles,
    LinkedCell       cell_struct,
    int              circleIdx)
{
    double sum = 0.0;

    int ix = cell_index(cp.x);
    int iy = cell_index(cp.y);

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            int nix = (ix + dx + Mx) % Mx;
            int niy = (iy + dy + My) % My;
            int cell = nix + niy * Mx;

            int cnt = cell_struct.count[cell];
            for (int k = 0; k < cnt; ++k) {
                int j = cell_struct.neighbors[cell][k];
                if (j == circleIdx) continue;
                sum += compute_patch_interaction(cp, circles[j]);
            }
        }
    }

    return sum;
}





__device__ static void insert_circle_in_cells(int circleIndex,
                                              const CircleParticle &cp,
                                              LinkedCell cell_struct) {
    int ix = cell_index(cp.x);
    int iy = cell_index(cp.y);
    int cell = ix + iy * Mx;

    int old_cnt = atomicAdd(&cell_struct.count[cell], 1);

    cell_struct.neighbors[cell][old_cnt] = circleIndex;
}


__device__ static void remove_circle_from_cells(int circleIndex,
                                                const CircleParticle &cp,
                                                LinkedCell cell_struct) {
    int ix = cell_index(cp.x);
    int iy = cell_index(cp.y);
    int cell = ix + iy * Mx;

    int cnt = cell_struct.count[cell];
    int pos = -1;
    for (int s = 0; s < cnt; ++s) {
        if (cell_struct.neighbors[cell][s] == circleIndex) {
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


__device__ static double apply_boundary(double coord, double L) {
    if (dc_bc) {
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

__device__ static void rotate_circle(CircleParticle *sp, curandState *rng) {
    double current_angle = 2.0 * atan2(sp->q[0], sp->q[1]);
    float u = curand_uniform(rng);

    double range = dc_lflag != 0 && u < LRMOVE ? LRROTMAX : RROTMAX;

    float v = curand_uniform(rng);

    double angle = v * 2 * range - range;


    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0)
        current_angle += 2 * M_PI;

    sp->q[0] = sin(current_angle / 2);
    sp->q[1] = cos(current_angle / 2);
}



__device__ static void move_circle(CircleParticle *cp, curandState *rng) {
    double range = dc_lflag != 0 && (curand_uniform(rng) < LRMOVE) ? LRDISPMAX : RDISPMAX;

    double dx = (curand_uniform(rng) * 2.0 - 1.0) * range;
    double dy = (curand_uniform(rng) * 2.0 - 1.0) * range;

    cp->x += dx;
    cp->y += dy;

    cp->x = apply_boundary(cp->x, Lx);
    cp->y = apply_boundary(cp->y, Ly);
}


static void shuffle_c(int *array) {
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::shuffle(array, array + 4, gen);
}


__global__ static void initRNG_c(curandState *states, unsigned long seed) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;

    curand_init(
        seed,
        tid,
        0,
        &states[tid]
    );
}



__global__ static void mc_step_kernel_c(int color,
                               curandState *rngStates,
                               CircleParticle *circles,
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

    int circleIdx = cell_struct.neighbors[tid][r];


    CircleParticle candidate = circles[circleIdx];

    float rm = curand_uniform(&localState);

    if (rm < RMOVE) {
        move_circle(&candidate, &localState);
    } else {
        rotate_circle(&candidate, &localState);
    }

    if (!d_is_overlapping_circle(candidate, circles, cell_struct, circleIdx)) {
        double myDE = compute_patch_energy(candidate, circles, cell_struct, circleIdx) - compute_patch_energy(
                          circles[circleIdx], circles, cell_struct, circleIdx);


        if (myDE <= 0 || exp(-myDE / KT) > curand_uniform(&localState)) {
            int new_cell = cell_index(candidate.x) + cell_index(candidate.y) * Mx;
            if (new_cell != tid) {
                remove_circle_from_cells(circleIdx, circles[circleIdx], cell_struct);
                insert_circle_in_cells(circleIdx, candidate, cell_struct);
            }

            circles[circleIdx] = candidate;
        }
    }
    rngStates[tid] = localState;
}

__constant__ int dc_num_patches;
__constant__ BoundaryCondition dc_bc;
__constant__ int dc_lflag;
__constant__ double dc_patch[MAX_PATCHES][2];


void animate_movement_c(int totalCircles, const CircleParticle *circles, LinkedCell cell_struct, int lflag) {
    cudaMemcpyToSymbol(dc_bc, &h_bc, sizeof(h_bc));
    cudaMemcpyToSymbol(dc_num_patches, &h_num_patches, sizeof(h_num_patches));
    cudaMemcpyToSymbol(dc_lflag, &lflag, sizeof(lflag));
    cudaMemcpyToSymbol(dc_patch, h_patch, h_num_patches * 2 * sizeof(double));


    CircleParticle *d_circles;

    CircleParticle *d_snapshots = nullptr;
    size_t snapshotCount = static_cast<size_t>(ANIMATION_STEPS) * totalCircles;
    cudaMalloc(&d_snapshots, snapshotCount * sizeof(CircleParticle));

    cudaMalloc(&d_circles, totalCircles * sizeof(CircleParticle));
    cudaMemcpy(d_circles, circles, totalCircles * sizeof(CircleParticle), cudaMemcpyHostToDevice);

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


    initRNG_c<<<blocks, threads>>>(d_rngStates, 1234UL);
    cudaDeviceSynchronize();

    const clock_t start = clock();

    for (int step = 0; step < ANIMATION_STEPS; step++) {
        shuffle_c(colors);

        for (int color: colors) {
            mc_step_kernel_c<<<blocks,threads>>>(
                color,
                d_rngStates,
                d_circles,
                d_cell_struct
            );
        }
        CircleParticle *dst = d_snapshots + static_cast<size_t>(step) * totalCircles;

        cudaMemcpy(dst,
                   d_circles,
                   totalCircles * sizeof(CircleParticle),
                   cudaMemcpyDeviceToDevice);
    }

    const clock_t end = clock();

    double elapsed_ms = (double)(end - start) * 1000.0 / (double)CLOCKS_PER_SEC;
    printf("Elapsed time: %.3f ms\n", elapsed_ms);

    CircleParticle *h_snapshots = nullptr;
    cudaMallocHost(&h_snapshots, snapshotCount * sizeof(CircleParticle));

    cudaMemcpy(h_snapshots,
               d_snapshots,
               snapshotCount * sizeof(CircleParticle),
               cudaMemcpyDeviceToHost);


    FILE *f = fopen("data/circle_animation.xyz", "w");
    if (!f) {
        printf("Error opening animation file\n");
        return;
    }

    for (int step = 0; step < ANIMATION_STEPS; ++step) {
        CircleParticle *frame = h_snapshots + static_cast<size_t>(step) * totalCircles;
        write_file_c(f, frame, totalCircles);
    }


    cudaFree(d_circles);
    cudaFree(d_count_mem);
    cudaFree(d_neigh_mem);
    cudaFree(d_rngStates);
    cudaFreeHost(h_snapshots);
    cudaFree(d_snapshots);

    fclose(f);
}
