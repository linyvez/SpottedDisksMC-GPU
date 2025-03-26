#include "periodic_boundary.h"

BoundaryCondition boundary_condition = NO_BOUNDARY;

void set_condition(BoundaryCondition condition) {
    boundary_condition = condition;
}

Particle adjust_circle_for_periodic(const Particle ref, const Particle sp) {
    Particle sp_adj = sp;
    double dx = sp.x - ref.x;
    double dy = sp.y - ref.y;

    periodic_boundary(&dx, &dy);

    sp_adj.x = ref.x + dx;
    sp_adj.y = ref.y + dy;

    return sp_adj;
}

void periodic_boundary(double *dx, double *dy) {
    if (boundary_condition == PERIODIC_X || boundary_condition == PERIODIC_BOTH) {
        *dx = fmod(*dx + 0.5 * Lx, Lx);
        if (*dx < 0)
            *dx += Lx;
        *dx -= 0.5 * Lx;
    }
    if (boundary_condition == PERIODIC_Y || boundary_condition == PERIODIC_BOTH) {
        *dy = fmod(*dy + 0.5 * Ly, Ly);
        if (*dy < 0)
            *dy += Ly;
        *dy -= 0.5 * Ly;
    }
}

int is_overlapping_pbc(const Particle new_p, const Particle particles[], int count) {
    for (int i = 0; i < count; i++) {
        Particle neighbor = particles[i];

        Particle adjusted_neighbor = adjust_circle_for_periodic(new_p, neighbor);

        double dx = new_p.x - adjusted_neighbor.x;
        double dy = new_p.y - adjusted_neighbor.y;
        double dist = sqrt(dx * dx + dy * dy);

        if (dist < SIGMA) {
            return 1;
        }
    }
    return 0;
}
