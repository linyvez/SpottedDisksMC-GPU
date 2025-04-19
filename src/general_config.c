#include "general_config.h"

#include <stddef.h>


int lflag = 0;
int *visited = NULL;

int current_epoch = 1;



void set_lflag(const int disp_condition) {
    lflag = disp_condition;
}

