#ifndef PATCH_H
#define PATCH_H

#include "general_config.h"

typedef enum {
    PATCH_NONE = 0,
    PATCH_TWO_CENTER = 1,
    PATCH_FOUR_CORNER = 2,
    PATCH_FOUR_CENTER = 3
} PatchType;

static const double PATCH_FOUR_CORNER_M[4][2] = {
    {PATCH_DELTA, -PATCH_DELTA},
    {-PATCH_DELTA, PATCH_DELTA},
    {-PATCH_DELTA, -PATCH_DELTA},
    {PATCH_DELTA, PATCH_DELTA}
};

static const double PATCH_FOUR_CENTER_M[4][2] = {
    {0, -PATCH_DELTA},
    {0, PATCH_DELTA},
    {-PATCH_DELTA, 0},
    {PATCH_DELTA, 0}
};

static const double PATCH_TWO_CENTER_M[2][2] = {
    {PATCH_DELTA, 0},
    {-PATCH_DELTA, 0},
};

#ifdef __cplusplus
extern "C" {
#endif

#include "general_config.h"

    extern int h_num_patches;
    extern const double (*h_patch)[2];

#ifdef __cplusplus
}
#endif



const double (*assign_patch_type(int))[2];

#endif //PATCH_H
