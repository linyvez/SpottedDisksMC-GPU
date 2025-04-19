#ifndef CONFIG_H
#define CONFIG_H

typedef struct {
    int num_particles;
    int animation_steps;
    double Lx;
    double Ly;
    double particle_size;
    int patch_radius_coeff;
    double temperature;
    double particle_interact_coeff;
} Config;


typedef struct {
    double radius;
    double delta;
    double interaction;

} PatchConfig;

typedef struct {
    double rdispmax;
    double lrdispmax;
    double rrotmax;
    double lrrotmax;
} MoveConfig;


typedef struct {
    double size;
    int Mx;
    int My;
    int num_cells;
    int max_nodes;
} CellConfig;


void initialize_global_configs(const char *);

extern const Config *GL_CFG;
extern const MoveConfig *MV_CFG;
extern const CellConfig *CL_CFG;
extern const PatchConfig *PH_CFG;



#endif //CONFIG_H
