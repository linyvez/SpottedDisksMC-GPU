
#include "generate_config.h"



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

// Task 1 : generate particles randomly 
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

// Task 2 : generate particles in a grid
void generate_square_grid() {

    Particle particles[N];
    int count = 0;

    // square grid

    int k = Lx/SIGMA; // number of particles in one row

    while (count < N) {
        Particle particle;
        particle.x = SIGMA/2 + count%k * SIGMA;
        particle.y = SIGMA/2 + SIGMA * (count/k);

        particles[count] = particle;
        count++;
    }

    FILE *config = fopen("configuration_squaregrid.pdb", "w");
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

void generate_hexagonal_grid() {
    Particle particles[N];
    int count = 0;
    double spacing = SIGMA;
    double vertical_spacing = spacing * sqrt(3.0) / 2.0;

    int num_cols = (int)(Lx / spacing);

    int row = 0;
    while (count < N) {
        for (int col = 0; col < num_cols && count < N; col++) {
            Particle particle;
            particle.x = SIGMA / 2.0 + col * spacing;
            if (row % 2 == 1) {
                particle.x += spacing / 2.0;
            }
            particle.y = vertical_spacing / 2.0 + row * vertical_spacing;
            particles[count++] = particle;
        }
        row++;
    }

    FILE *config = fopen("configuration_hexgrid.pdb", "w");
    if (config == NULL) {
        printf("Error while opening file\n");
        return;
    }

    for (int i = 0; i < N; i++) {
        fprintf(config, "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
            i + 1,
            "ATOM",
            "RES",
            "A",
            1,
            particles[i].x,
            particles[i].y,
            0.0);
    }

    fclose(config);
}


int main() {
    generate_random();
    generate_grid();
    return 0;
}
