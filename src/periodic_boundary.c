// #include "periodic_boundary.h"


//  void apply_periodic_boundary(double *dx, double *dy) {
//      if (*dx >  0.5 * X_SIZE) *dx -= X_SIZE;
//      if (*dx <= -0.5 * X_SIZE) *dx += X_SIZE;
//      if (*dy >  0.5 * Y_SIZE) *dy -= Y_SIZE;
//      if (*dy <= -0.5 * Y_SIZE) *dy += Y_SIZE;
//  }

//  int check_overlap_circle(double x, double y, int count) {
//      for (int i = 0; i < count; i++) {
//          double dx = x - particles[i].x;
//          double dy = y - particles[i].y;
//          apply_periodic_boundary(&dx, &dy);
//          double distance = sqrt(dx * dx + dy * dy);
//          if (distance < SIGMA) return 1;
//      }
//      return 0;
//  }

//  void generate_particles_circle() {
//      srand(time(NULL));
//      int count = 0;
//      double cell_width = X_SIZE / GRID_SIZE;
//      double cell_height = Y_SIZE / GRID_SIZE;

//      for (int gx = 0; gx < GRID_SIZE; gx++) {
//          for (int gy = 0; gy < GRID_SIZE; gy++) {
//              if (count >= N) return;

//              double x = (gx + 0.5) * cell_width + ((double)rand() / RAND_MAX - 0.5) * cell_width;
//              double y = (gy + 0.5) * cell_height + ((double)rand() / RAND_MAX - 0.5) * cell_height;

//              for (int dx = -1; dx <= 1; dx++) {
//                  for (int dy = -1; dy <= 1; dy++) {
//                      double new_x = x + dx * cell_width;
//                      double new_y = y + dy * cell_height;
//                      apply_periodic_boundary(&new_x, &new_y);
//                      if (!check_overlap_circle(new_x, new_y, count)) {
//                          particles[count].x = new_x;
//                          particles[count].y = new_y;
//                          count++;
//                      }
//                  }
//              }
//          }
//      }
//  }


// int check_overlap_sq(double x, double y, int count) {
//      for (int i = 0; i < count; i++) {
//          double dx = x - particles[i].x;
//          double dy = y - particles[i].y;
//          apply_periodic_boundary(&dx, &dy);

//          double distance = sqrt(dx * dx + dy * dy);
//          if (distance < 2 * RADIUS) {
//              return 1;
//          }
//      }
//      return 0;
//  }

// void generate_particles_sq() {
//      srand(time(NULL));
//      int count = 0;
//      double cell_width = X_SIZE / GRID_SIZE;
//      double cell_height = Y_SIZE / GRID_SIZE;

//      for (int gx = 0; gx < GRID_SIZE; gx++) {
//          for (int gy = 0; gy < GRID_SIZE; gy++) {
//              if (count >= N) return;

//              double x = (gx + 0.5) * cell_width + ((double)rand() / RAND_MAX - 0.5) * cell_width;
//              double y = (gy + 0.5) * cell_height + ((double)rand() / RAND_MAX - 0.5) * cell_height;

//              for (int dx = -1; dx <= 1; dx++) {
//                  for (int dy = -1; dy <= 1; dy++) {
//                      double new_x = x + dx * cell_width;
//                      double new_y = y + dy * cell_height;
//                      apply_periodic_boundary(&new_x, &new_y);
//                      if (!check_overlap_sq(new_x, new_y, count)) {
//                          particles[count].x = new_x;
//                          particles[count].y = new_y;
//                          count++;
//                      }
//                  }
//              }
//          }
//      }
//  }

//  void save_circle_to_file(const char *filename) {
//      FILE *file = fopen(filename, "w");
//      if (!file) {
//          perror("Failed to open file for writing");
//          return;
//      }
//      fprintf(file, "%d\n", N);
//      fprintf(file, "Properties=pos:R:2\n");
//      for (int i = 0; i < N; i++) {
//          fprintf(file, "%.6f %.6f\n", particles[i].x, particles[i].y);
//      }
//      fclose(file);
//  }

//  // int main() {
//  //     generate_particles_circle();
//  //     save_circle_to_file("particles_circle.xyz");
//  //
//  //     generate_particles_sq();
//  //     save_circle_to_file("particles_sqare.xyz");
//  //     return 0;
//  // }