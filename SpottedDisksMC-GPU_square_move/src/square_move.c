#include "square_move.h"
#include "square_config.h"
#include "periodic_boundary.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void remove_square_from_cells(int squareIndex, const SquareParticle sp)
{
    for (int ix = sp.cell_min_x; ix <= sp.cell_max_x; ix++)
    {
        for (int iy = sp.cell_min_y; iy <= sp.cell_max_y; iy++)
        {
            int cell = ix + Mx * iy;
            int prev = -1;
            int curr = head[cell];
            while (curr != -1)
            {
                if (nodePool[curr].squareIndex == squareIndex)
                {
                    if (prev == -1)
                    {
                        head[cell] = nodePool[curr].next;
                    }
                    else
                    {
                        nodePool[prev].next = nodePool[curr].next;
                    }
                    nodePool[curr].next = freeList;
                    freeList = curr;
                    break;
                }
                prev = curr;
                curr = nodePool[curr].next;
            }
        }
    }
}

void rotate_square(SquareParticle *sp)
{

    double current_angle = 2.0 * atan2(sp->q[2], sp->q[3]);
    double angle = ((rand() / (double)RAND_MAX) * 2 * MAX_ROTATION_ANGLE) - MAX_ROTATION_ANGLE;
    current_angle += angle;
    current_angle = fmod(current_angle, 2 * M_PI);
    if (current_angle < 0)
        current_angle += 2 * M_PI;

    sp->q[2] = sin(current_angle / 2);
    sp->q[3] = cos(current_angle / 2);
}

void move_square(SquareParticle *sp)
{

    double dx = (rand() / (double)RAND_MAX) * 2 * SQUARE_SIDE - SQUARE_SIDE;
    double dy = (rand() / (double)RAND_MAX) * 2 * SQUARE_SIDE - SQUARE_SIDE;
    sp->x += dx;
    sp->y += dy;
    if (boundary_condition == NO_BOUNDARY)
    {
        if (sp->x < 0)
            sp->x = 0;
        if (sp->x > Lx)
            sp->x = Lx;
        if (sp->y < 0)
            sp->y = 0;
        if (sp->y > Ly)
            sp->y = Ly;
    }
    else if (boundary_condition == PERIODIC_BOTH)
    {
        if (sp->x < 0)
            sp->x += Lx;
        if (sp->x > Lx)
            sp->x -= Lx;
        if (sp->y < 0)
            sp->y += Ly;
        if (sp->y > Ly)
            sp->y -= Ly;
    }
    else if (boundary_condition == PERIODIC_X)
    {
        if (sp->x < 0)
            sp->x += Lx;
        if (sp->x > Lx)
            sp->x -= Lx;
        if (sp->y < 0)
            sp->y = 0;
        if (sp->y > Ly)
            sp->y = Ly;
    }
    else if (boundary_condition == PERIODIC_Y)
    {
        if (sp->x < 0)
            sp->x = 0;
        if (sp->x > Lx)
            sp->x = Lx;
        if (sp->y < 0)
            sp->y += Ly;
        if (sp->y > Ly)
            sp->y -= Ly;
    }
}

void animate_movement(int steps, int totalSquares)
{
    FILE *f = fopen("square_animation.xyz", "w");
    if (!f)
    {
        printf("Error opening animation file\n");
        return;
    }

    for (int step = 0; step < steps; step++)
    {
        for (int attempt = 0; attempt < N; attempt++)
        {
            int idx = rand() % totalSquares;
            double r = (double)rand() / RAND_MAX;

            remove_square_from_cells(idx, squares[idx]);

            SquareParticle candidate = squares[idx];
            if (r < RMOVE)
                move_square(&candidate);
            else
                rotate_square(&candidate);

            update_square_AABB(&candidate);

            if (!is_overlapping_square(candidate, idx))
            {
                squares[idx] = candidate;
            }

            insert_square_in_cells(idx, squares[idx]);
        }

        fprintf(f, "%d\n", totalSquares + totalSquares * 4);
        fprintf(f, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n", Lx, Ly);
        for (int i = 0; i < totalSquares; i++)
        {

            fprintf(f, "B %.6f %.6f %.6f  %.6f %.6f %.6f %.6f  %.6f %.6f %.6f\n",
                    squares[i].x,
                    squares[i].y,
                    squares[i].z,
                    squares[i].q[0],
                    squares[i].q[1],
                    squares[i].q[2],
                    squares[i].q[3],
                    squares[i].shape[0] / 2.0,
                    squares[i].shape[1] / 2.0,
                    squares[i].shape[2] / 2.0);
        }

        for (int i = 0; i < totalSquares; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                double global_x, global_y;
                compute_patch_global_position(squares[i], squares[i].patches[j], &global_x, &global_y);
                fprintf(f, "P %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                        global_x,
                        global_y,
                        squares[i].z,
                        squares[i].q[0],
                        squares[i].q[1],
                        squares[i].q[2],
                        squares[i].q[3],
                        squares[i].patches[j].shape[0] / 2.0,
                        squares[i].patches[j].shape[1] / 2.0,
                        squares[i].patches[j].shape[2] / 2.0);
            }
        }
    }
    fclose(f);
}