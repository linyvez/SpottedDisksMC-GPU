#ifndef SQUARE_MOVE_CUH
#define SQUARE_MOVE_CUH
#include "general_config.h"
#include "square_config.h"

#ifdef __CUDACC__        // only NVCC defines this
# define CUDA_CONST __constant__
#else
# define CUDA_CONST       // plain C/C++ sees nothing
#endif


extern CUDA_CONST int d_lflag;
extern CUDA_CONST BoundaryCondition d_bc;
extern CUDA_CONST int d_num_patches;
extern CUDA_CONST double d_patch[MAX_PATCHES][2];

#ifdef __cplusplus
extern "C" {
#endif


void animate_movement(int, const SquareParticle *squares, LinkedCell cell_struct, int);

#ifdef __cplusplus
}
#endif


#endif // SQUARE_MOVE_CUH
