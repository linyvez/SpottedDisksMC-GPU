#include "patch.h"
#include <stddef.h>
#include "general_config.h"

PatchType patch_type = PATCH_NONE;

void set_patch_type(PatchType patch) {
    patch_type = patch;
}

const double (*assign_patch_type(int *num_patches))[2] {
    switch (patch_type) {
        case PATCH_TWO_CENTER:
            *num_patches = 2;
        return PATCH_TWO_CENTER_M;

        case PATCH_FOUR_CORNER:
            *num_patches = 4;
        return PATCH_FOUR_CORNER_M;

        case PATCH_FOUR_CENTER:
            *num_patches = 4;
        return PATCH_FOUR_CENTER_M;

        default:
            *num_patches = 0;
        return NULL;
    }
}