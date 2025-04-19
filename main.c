#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "circle_move.h"
#include "square_move.h"
#include "periodic_boundary.h"
#include "errno.h"
#include "patch.h"
#include "config.h"
#include "print_error.h"



static int parse_int_arg(const char *arg, const char *name) {
    char *endptr;
    errno = 0;
    long tmp = strtol(arg, &endptr, 10);

    if (errno != 0 || *endptr != '\0' || tmp < INT_MIN || tmp > INT_MAX) {
        print_error(true, "Invalid integer for %s: '%s'\n", name, arg);
        exit(EXIT_FAILURE);
    }

    return (int) tmp;
}

static int parse_lflag(const char *arg) {
    int val = parse_int_arg(arg, "lflag");
    if (val != 0 && val != 1) {
        print_error(false, "Invalid value for lflag (must be 0 or 1): '%s'\n", arg);
        exit(EXIT_FAILURE);
    }
    return val;
}


int main(int argc, char *argv[]) {
    if (argc != 6) {
        print_error(false, "Incorrect number of arguments.");
        return EXIT_FAILURE;
    }

    srand48(time(NULL));

    const int condition = parse_int_arg(argv[1], "condition");
    set_condition(condition);

    const ParticleType particle_type = parse_int_arg(argv[2], "particle_type");

    const PatchType patch = parse_int_arg(argv[3], "patchy_type");
    set_patch_type(patch);

    const int lflag_value = parse_lflag(argv[4]);
    set_lflag(lflag_value);

    const char *config_file = argv[5];

    initialize_global_configs(config_file);


    int count = 0;
    clock_t start = clock();;
    switch (particle_type) {
        case PARTICLE_DISK:
            count = generate_random_circles();
            animate_circle_movement(GL_CFG->animation_steps, count);
            break;

        case PARTICLE_SQUARE:
            count = generate_random_squares();
            animate_movement(GL_CFG->animation_steps, count);
            break;

        default:
            print_error(false, "Unknown particle type: '%d'", particle_type);
            exit(EXIT_FAILURE);
    }
    clock_t end = clock();

    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);

    free((void *) GL_CFG);
    free((void *) MV_CFG);
    free((void *) CL_CFG);
    free((void *) PH_CFG);


    free(visited);

    return EXIT_SUCCESS;
}
