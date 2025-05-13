#include "shared_utilities.h"
#include "circle_config.h"
#include <patch.h>
#include <tgmath.h>


void write_file(FILE *f, SquareParticle *squares, int totalSquares) {
    fprintf(f, "%d\n", totalSquares + h_num_patches * totalSquares);
    fprintf(
        f,
        "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3 Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n",
        Lx, Ly);

    for (int i = 0; i < totalSquares; i++) {
        double angle = 2.0 * atan2(squares[i].q[0], squares[i].q[1]);
        double cosA = cos(angle);
        double sinA = sin(angle);
        fprintf(f, "B %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                squares[i].x,
                squares[i].y,
                Z_AXIS,
                Z_AXIS,
                Z_AXIS,
                squares[i].q[0],
                squares[i].q[1],
                PARTICLE_SIZE * OVITO_MUL,
                PARTICLE_SIZE * OVITO_MUL,
                PARTICLE_SIZE * pow(OVITO_MUL, 4));

        double x0 = squares[i].x, y0 = squares[i].y;
        for (int j = 0; j < h_num_patches; j++) {
            double rx = h_patch[j][0], ry = h_patch[j][1];
            double gx = x0 + rx * cosA - ry * sinA;
            double gy = y0 + rx * sinA + ry * cosA;
            fprintf(f,
                    "P %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n",
                    gx, gy,
                    Z_AXIS, Z_AXIS, Z_AXIS,
                    squares[i].q[0], squares[i].q[1],
                    PATCH_RADIUS * OVITO_MUL,
                    PATCH_RADIUS * OVITO_MUL,
                    PATCH_RADIUS * pow(OVITO_MUL, 3));
        }
    }
}

void write_file_c(FILE *f, CircleParticle *circles, int totalCircles) {
    fprintf(f, "%d\n", totalCircles + h_num_patches * totalCircles);
    fprintf(f,
        "Properties=species:S:1:pos:R:3:"
        "orientation:R:4:aspherical_shape:R:3 "
        "Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0001\"\n",
        Lx, Ly);

    for (int i = 0; i < totalCircles; ++i) {

        double qz = circles[i].q[0];
        double qw = circles[i].q[1];
        double norm = sqrt(qz*qz + qw*qw);
        if (norm > 0.0) { qz /= norm;  qw /= norm; }

        double theta = 2.0 * atan2(qz, qw);
        double c = cos(theta), s = sin(theta);

        fprintf(f,
            "B %8.3f %8.3f %8.3f  "
            "%8.3f %8.3f %8.3f %8.3f  "
            "%8.3f %8.3f %8.3f\n",
            circles[i].x,
            circles[i].y,
            Z_AXIS,
            0.0, 0.0,
            qz, qw,
            PARTICLE_SIZE * OVITO_MUL,
            PARTICLE_SIZE * OVITO_MUL,
            PARTICLE_SIZE * OVITO_MUL
        );


        for (int j = 0; j < h_num_patches; ++j) {
            double rx = h_patch[j][0];
            double ry = h_patch[j][1];
            double gx = circles[i].x + rx * c - ry * s;
            double gy = circles[i].y + rx * s + ry * c;

            fprintf(f,
                "P %8.3f %8.3f %8.3f  "
                "%8.3f %8.3f %8.3f %8.3f  "
                "%8.3f %8.3f %8.3f\n",
                gx, gy,
                Z_AXIS,
                0.0, 0.0,
                qz, qw,
                PATCH_RADIUS * OVITO_MUL,
                PATCH_RADIUS * OVITO_MUL,
                PATCH_RADIUS * OVITO_MUL
            );
        }
    }
}



