#include "square_config.h"
#include "square_move.h"
#include "periodic_boundary.h"

int main(int argc, char *argv[]) {
    srand((unsigned)time(NULL));
    const int condition = atoi(argv[1]);
    set_condition(condition);

    const clock_t start = clock();
    int total_squares = generate_random_squares();
    animate_movement(ANIMATION_STEPS, total_squares);
    const clock_t end = clock();

    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", cpu_time_used);

    return 0;
}