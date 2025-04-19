#include "patch.h"
#include "config.h"
#include "print_error.h"
#include <stdlib.h>

PatchType patch_type = PATCH_NONE;
int num_patches = 0;

void set_patch_type(PatchType patch) {
    patch_type = patch;
}

double (*assign_patch_type())[2] {
    double (*local)[2] = NULL;

    switch (patch_type) {
        case PATCH_NONE:
            num_patches = 0;
            return NULL;

        case PATCH_TWO_CENTER:
            num_patches = 2;
            local = malloc(sizeof(double[2]) * num_patches);
            if (!local) {
                print_error(true, "Memory allocation failed for PATCH_TWO_CENTER");
                exit(EXIT_FAILURE);
            }
            local[0][0] = -PH_CFG->delta; local[0][1] = 0.0;
            local[1][0] = PH_CFG->delta; local[1][1] = 0.0;
            break;

        case PATCH_FOUR_CENTER:
            num_patches = 4;
            local = malloc(sizeof(double[2]) * num_patches);
            if (!local) {
                print_error(true, "Memory allocation failed for PATCH_FOUR_CENTER");
                exit(EXIT_FAILURE);
            }
            local[0][0] = -PH_CFG->delta; local[0][1] = 0.0;
            local[1][0] = PH_CFG->delta;  local[1][1] = 0.0;
            local[2][0] = 0.0; local[2][1] = PH_CFG->delta;
            local[3][0] = 0.0; local[3][1] = -PH_CFG->delta;
            break;

        case PATCH_FOUR_CORNER:
            num_patches = 4;
            local = malloc(sizeof(double[2]) * num_patches);
            if (!local) {
                print_error(true, "Memory allocation failed for PATCH_FOUR_CORNER");
                exit(EXIT_FAILURE);
            }
            local[0][0] = -PH_CFG->delta; local[0][1] = -PH_CFG->delta;
            local[1][0] = PH_CFG->delta;  local[1][1] = -PH_CFG->delta;
            local[2][0] = PH_CFG->delta;  local[2][1] = PH_CFG->delta;
            local[3][0] = -PH_CFG->delta; local[3][1] = PH_CFG->delta;
            break;

        default:
            num_patches = 0;
            return NULL;
    }

    return local;
}
