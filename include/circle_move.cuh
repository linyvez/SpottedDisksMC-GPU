#ifndef CIRCLE_MOVE_CUH
#define CIRCLE_MOVE_CUH

#include "general_config.h"
#include "circle_config.h"

#ifdef __CUDACC__        // only NVCC defines this
# define CUDA_CONST __constant__
#else
# define CUDA_CONST       // plain C/C++ sees nothing
#endif


extern CUDA_CONST int dc_lflag;
extern CUDA_CONST BoundaryCondition dc_bc;
extern CUDA_CONST int dc_num_patches;
extern CUDA_CONST double dc_patch[MAX_PATCHES][2];

#ifdef __cplusplus
extern "C" {
#endif


    void animate_movement_c(int, const CircleParticle *circles, LinkedCell cell_struct, int);

#ifdef __cplusplus
}
#endif



#endif //CIRCLE_MOVE_CUH
