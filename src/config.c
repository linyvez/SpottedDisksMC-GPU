#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include "print_error.h"
#include "config.h"

const Config *GL_CFG = NULL;
const MoveConfig *MV_CFG = NULL;
const CellConfig *CL_CFG = NULL;
const PatchConfig *PH_CFG = NULL;

static char *trim(char *str) {
    while (isspace(*str)) str++;

    if (*str == '\0') return str;

    char *end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) {
        *end = '\0';
        end--;
    }
    return str;
}

static int parse_int(const char *value, const char *key, int line_number) {
    char *endptr;
    errno = 0;
    long val = strtol(value, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
        print_error(true, "Error parsing integer \"%s\" for key \"%s\" at line %d\n", value, key, line_number);
        exit(EXIT_FAILURE);
    }
    return (int) val;
}

static double parse_double(const char *value, const char *key, int line_number) {
    char *endptr;
    errno = 0;
    double val = strtod(value, &endptr);
    if (errno != 0 || *endptr != '\0') {
        print_error(true, "Error parsing double \"%s\" for key \"%s\" at line %d\n", value, key, line_number);
        exit(EXIT_FAILURE);
    }
    return val;
}


Config config_reader(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        print_error(true, "Error opening configuration file");
        exit(EXIT_FAILURE);
    }
    Config config = {0};

    char line[256];
    int line_number = 0;

    while (fgets(line, sizeof(line), fp)) {
        line_number++;

        line[strcspn(line, "\n")] = '\0';

        char *comment = strchr(line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }

        char *trimmed = trim(line);
        if (strlen(trimmed) == 0) {
            continue;
        }

        char *equal_sign = strchr(trimmed, '=');
        if (!equal_sign) {
            print_error(true, "Error in configuration file at line %d: missing '='.\n", line_number);
            exit(EXIT_FAILURE);
        }

        *equal_sign = '\0';
        char *key = trim(trimmed);
        char *value = trim(equal_sign + 1);

        if (strlen(key) == 0 || strlen(value) == 0) {
            print_error(true, "Error: Empty key or value at line %d.\n", line_number);
            exit(EXIT_FAILURE);
        }

        if (value[0] == '\"' && value[strlen(value) - 1] == '\"') {
            value[strlen(value) - 1] = '\0';
            value++;
        }

        if (strcmp(key, "num_particles") == 0) {
            config.num_particles = parse_int(value, key, line_number);
        } else if (strcmp(key, "animation_steps") == 0) {
            config.animation_steps = parse_int(value, key, line_number);
        } else if (strcmp(key, "patch_radius_coeff") == 0) {
            config.patch_radius_coeff = parse_int(value, key, line_number);
        } else if (strcmp(key, "particle_size") == 0) {
            config.particle_size = parse_double(value, key, line_number);
        } else if (strcmp(key, "temperature") == 0) {
            config.temperature = parse_double(value, key, line_number);
        } else if (strcmp(key, "Lx") == 0) {
            config.Lx = parse_double(value, key, line_number);
        } else if (strcmp(key, "Ly") == 0) {
            config.Ly = parse_double(value, key, line_number);
        } else if (strcmp(key, "particle_interact_coeff") == 0) {
            config.particle_interact_coeff = parse_double(value, key, line_number);
        } else {
            fprintf(stderr, "Warning: Unknown configuration key \"%s\" at line %d\n", key, line_number);
        }
    }

    fclose(fp);

#define CHECK_PARAM(x) \
        if (config.x == -1 || config.x == -1.0) { \
        print_error(true, "Missing required config parameter: %s\n", #x); \
        exit(EXIT_FAILURE); \
        }

    CHECK_PARAM(num_particles);
    CHECK_PARAM(animation_steps);
    CHECK_PARAM(particle_size);
    CHECK_PARAM(patch_radius_coeff);
    CHECK_PARAM(temperature);
    CHECK_PARAM(Lx);
    CHECK_PARAM(Ly);


    return config;
}

void initialize_global_configs(const char *filename) {
    Config *config = malloc(sizeof(Config));
    if (!config) {
        print_error(true, "Memory allocation failed for Config");
        exit(EXIT_FAILURE);
    }
    *config = config_reader(filename);
    GL_CFG = config;

    MoveConfig *mv_cfg = malloc(sizeof(MoveConfig));
    if (!mv_cfg) {
        print_error(true, "Memory allocation failed for MoveConfig");
        exit(EXIT_FAILURE);
    }
    mv_cfg->rdispmax = config->particle_size;
    mv_cfg->lrdispmax = config->Lx * 0.5;
    mv_cfg->rrotmax = M_PI / 10;
    mv_cfg->lrrotmax = M_PI / 2;
    MV_CFG = mv_cfg;

    CellConfig *cl_cfg = malloc(sizeof(CellConfig));
    if (!cl_cfg) {
        print_error(true, "Memory allocation failed for CellConfig");
        exit(EXIT_FAILURE);
    }
    cl_cfg->size = sqrt(2);
    cl_cfg->Mx = ceil(config->Lx / cl_cfg->size);
    cl_cfg->My = ceil(config->Ly / cl_cfg->size);
    cl_cfg->num_cells = cl_cfg->Mx * cl_cfg->My;
    cl_cfg->max_nodes = config->num_particles * 8;
    CL_CFG = cl_cfg;

    PatchConfig *ph_cfg = malloc(sizeof(PatchConfig));
    if (!ph_cfg) {
        print_error(true, "Memory allocation failed for PatchConfig");
        exit(EXIT_FAILURE);
    }
    ph_cfg->radius = config->particle_size / config->patch_radius_coeff;
    ph_cfg->delta = config->particle_size / 2.0;
    ph_cfg->interaction = config->temperature * config->particle_interact_coeff;
    PH_CFG = ph_cfg;
}
