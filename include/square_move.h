#ifndef SQUARE_MOVE_H
#define SQUARE_MOVE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "square_config.h"



#define MAX_TRAVERSE SQUARE_SIDE
#define LMAX_TRAVERSE Lx / 2



void rotate_square(SquareParticle *sp);
void move_square(SquareParticle *sp);
void animate_movement(int, int);

#endif // SQUARE_MOVE_H
