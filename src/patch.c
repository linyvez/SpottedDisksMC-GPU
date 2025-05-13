#include "patch.h"

int h_num_patches;

const double (*h_patch)[2];


const double (*assign_patch_type(int patch_type))[2] {
    switch (patch_type) {
        case PATCH_TWO_CENTER:
            h_num_patches = 2;
            return PATCH_TWO_CENTER_M;

        case PATCH_FOUR_CORNER:
            h_num_patches = 4;
            return PATCH_FOUR_CORNER_M;

        case PATCH_FOUR_CENTER:
            h_num_patches = 4;
            return PATCH_FOUR_CENTER_M;

        default:
            h_num_patches = 0;
            return NULL;
    }
}
