#include "square_config.h"
#include "square_move.h"
#include "periodic_boundary.h"

int main(int argc, char *argv[])
{   
    // condition, particle type, patchy type
    if (argc != 4)
    {
        fprintf(stderr, "Error, incorrect number of arguments.\n");
        return 1;
    }

    srand((unsigned)time(NULL));
    const int condition = atoi(argv[1]);
    set_condition(condition);

    const int particle_type = atoi(argv[2]);
    const int patchy_type = atoi(argv[3]);

    const clock_t start = clock();
    if (particle_type == 0) {
        int totalCircles = generate_random_circles();
        animate_circle_movement(ANIMATION_STEPS, totalCircles);
    }
    else {
        int totalSquares = generate_random_squares(patchy_type);
        animate_movement(ANIMATION_STEPS, totalSquares);
    }
    
    const clock_t end = clock();

    double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);

    return 0;
}