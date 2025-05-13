#ifndef SHARED_UTILITIES_H
#define SHARED_UTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif
#include <square_config.h>
#include <stdio.h>

#include "circle_config.h"


    void write_file(FILE *f, SquareParticle *squares, int totalSquares);
    void write_file_c(FILE *f, CircleParticle *circles, int totalCircles);

#ifdef __cplusplus
}
#endif


#endif //SHARED_UTILITIES_H
