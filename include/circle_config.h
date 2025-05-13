#ifndef CIRCLE_CONFIG_H
#define CIRCLE_CONFIG_H

#include "general_config.h"


typedef struct {
    double x, y;

    double q[2];

} CircleParticle;


int generate_random_circles(int patch, CircleParticle * squares, LinkedCell cell_struct);

#endif //CIRCLE_CONFIG_H
