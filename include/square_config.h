#ifndef SQUARE_CONFIG_H
#define SQUARE_CONFIG_H

#include "general_config.h"


typedef struct {
    double x, y;
    double q[2];

} SquareParticle;


int generate_random_squares(int patch, SquareParticle * squares, LinkedCell cell_struct);

#endif // SQUARE_CONFIG_H
