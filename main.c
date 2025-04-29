#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "square_config.h"
#include "square_move.h"
#include "periodic_boundary.h"
#include "circle_config.h"
#include "circle_move.h"
#include "patch.h"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Error, incorrect number of arguments.\n");
        return 1;
    }

    srand(time(NULL));
    const int condition = atoi(argv[1]);
    set_condition(condition);

    // const int particle_type = atoi(argv[2]);
    const int patch_type = atoi(argv[3]);

    set_patch_type(patch_type);


    int num_patches;


    const double (*patch_delta)[2] = assign_patch_type(&num_patches);


    const clock_t start = clock();
    int totalSquares = generate_random_squares(num_patches, patch_delta);
    animate_movement(ANIMATION_STEPS, totalSquares, num_patches);

    const clock_t end = clock();

    double cpu_time_used = (double) (end - start) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);

    return 0;
}
