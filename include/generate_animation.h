#ifndef GENERATE_ANIMATION_H
#define GENERATE_ANIMATION_H

#define kB  1.0
#define T   0.1

#define MAX_DISPLACEMENT    1


#include "cell_linked_algorithm.h"

double calculate_energy(double x, double y);
double random_uniform();
void write_to_file(int frame_number);
void no_optimization_move();



#endif