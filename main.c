#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <periodic_boundary.h>

#include "square_move.cuh"
#include <square_config.h>
#include <string.h>

#include "circle_config.h"
#include "circle_move.cuh"


int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Error, incorrect number of arguments.\n");
        return 1;
    }

    srand(time(NULL));
    int condition = atoi(argv[1]);

    set_condition(condition);

    int particle = atoi(argv[2]);


    int patch = atoi(argv[3]);

    int flag = atoi(argv[4]);


    SquareParticle squares[N_PAR];
    CircleParticle circles[N_PAR];
    int neighborCount[NUM_CELLS];
    int neighbors[NUM_CELLS][MAX_NEIGH];

    memset(neighborCount, 0, sizeof(neighborCount));

    LinkedCell cell_struct;
    cell_struct.count = neighborCount;
    cell_struct.neighbors = neighbors;

    if (particle == 1) {
        int totalSquares = generate_random_squares(patch, squares, cell_struct);
        animate_movement(totalSquares, squares, cell_struct, flag);

    } else if (particle == 0) {
        int totalCircles = generate_random_circles(patch, circles, cell_struct);
        animate_movement_c(totalCircles, circles, cell_struct, flag);
    }



    return 0;
}
