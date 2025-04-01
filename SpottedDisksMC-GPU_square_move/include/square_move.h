
#ifndef SQUARE_MOVE_H
#define SQUARE_MOVE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "square_config.h"

#define ANIMATION_STEPS 1000
#define MAX_ROTATION_ANGLE (M_PI / 10)
#define RMOVE 0.5

void rotate_square(SquareParticle *sp);
void move_square(SquareParticle *sp);
void animate_movement(int steps, int totalSquares);

#endif //SQUARE_MOVE_H
