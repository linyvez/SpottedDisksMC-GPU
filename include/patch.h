#ifndef PATCH_H
#define PATCH_H

typedef enum {
    PATCH_NONE = 0,
    PATCH_TWO_CENTER = 1,
    PATCH_FOUR_CORNER = 2,
    PATCH_FOUR_CENTER = 3
} PatchType;

extern PatchType patch_type;

void set_patch_type(PatchType);

const double (*assign_patch_type(int *))[2];

#endif //PATCH_H
