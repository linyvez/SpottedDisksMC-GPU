#ifndef CIRCLE_MOVE_H
#define CIRCLE_MOVE_H
#include "circle_config.h"

void rotate_circle(CircleParticle *c);
Coordinates calculate_new_coordinates(CircleParticle *c);
void animate_circle_movement(int steps, int totalCircles);

#endif // CIRCLE_MOVE_H