#include "periodic_boundary.h"
#include <time.h>



int is_overlaping(Particle particle, Particle particles[], int n) {
    for (int i = 0; i < n; i++) {
        double x_diff = particle.x - particles[i].x;
        double y_diff = particle.y - particles[i].y;
        double distance = sqrt(x_diff*x_diff + y_diff*y_diff);
        if (distance < SIGMA) {
            return 1;
        }
    }
    return 0;
}

void generate_random() {
    Particle particles[N];
    int count = 0;

    srand((unsigned)time(NULL));

    while (count < N) {
        Particle particle;
        particle.x = ((double)rand() / RAND_MAX) * Lx;
        particle.y = ((double)rand() / RAND_MAX) * Ly;

        if (!is_overlaping(particle, particles, count)) {
            particles[count] = particle;
            count++;
        }
    }

    FILE *config = fopen("configuration_random.pdb", "w");
    if (config == NULL) {
        printf("Error while opening file\n");
    }

    for (int i = 0; i < N; i++) {
        fprintf(config, "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
            i+1,
            "ATOM",
            "RES",
            "A",      
            1,         
            particles[i].x,
            particles[i].y,
            0.0
            );
    }
    
    fclose(config);
}


void generate_random_with_pbc(int mode) {
    Particle particles[N];
    int count = 0;
    int max_attempts = 50;
    int attempts = 0;

    srand((unsigned)time(NULL));

    double region_width = Lx / 3.0;
    double region_height = Ly / 3.0;
    double offset_x = (Lx - region_width) / 2.0;
    double offset_y = (Ly - region_height) / 2.0;


    while (count < N / 9.0) {
        Particle particle;
        particle.x = ((double)rand() / RAND_MAX) * region_width + offset_x;
        particle.y = ((double)rand() / RAND_MAX) * region_height + offset_y;

        if (!is_overlaping(particle, particles, count)) {
            particles[count++] = particle;
            attempts = 0;
        } else {
            attempts++;
        }

        if (attempts >= max_attempts) {
            printf("Max attempts reached, stopping particle generation.\n");
            break;
        }
    }

    printf("Generated %d particles in the central region\n", count);

    apply_periodic_boundary(particles, &count, mode);

    FILE *config = fopen("configuration_pbc.pdb", "w");
    if (config == NULL) {
        printf("Error while opening file\n");
        return;
    }

    for (int i = 0; i < count; i++) {
        fprintf(config, "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                i + 1, "ATOM", "RES", "A", 1, particles[i].x, particles[i].y, 0.0);
    }

    fclose(config);
}



int main(int argc, char *argv[]) {
    clock_t start, end;
    double cpu_time_used;
    if (argc != 2) {
        return 1;
    }

    int mode = atoi(argv[1]);
    if (mode < 0 || mode > 2) {
        printf("Invalid mode. Choose 0, 1, or 2.\n");
        return 1;
    }
    start = clock();
    generate_random_with_pbc(mode);
    // generate_random();
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken: %f s\n", cpu_time_used);
    return 0;
}
