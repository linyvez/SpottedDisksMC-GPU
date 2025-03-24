#ifndef SQUARE_CONFIG_H
#define SQUARE_CONFIG_H
#define N 80
#define Lx 24.0
#define Ly 24.0
#define SQUARE_SIDE 2.0
#define MAX_ANGLE (2 * M_PI)

#define CIRCUM_RADIUS(side) ((side) * sqrt(2.0) / 2.0)

#define CELL_SIZE SQUARE_SIDE * 1.41421356
#define Mx ((int)(Lx / CELL_SIZE))
#define My ((int)(Ly / CELL_SIZE))
#define NUM_CELLS Mx * My

typedef struct {
    double x, y, z;
    double q[4];
    double shape[3];
} SquareParticle;


typedef struct {
    int squareIndex;
    int next;
} Node;

#define MAX_NODES (N * 16)
Node nodePool[MAX_NODES];
extern int nodeCount;


int head[NUM_CELLS];

SquareParticle squares[N];

int get_cell_index(double x, double y);
int is_overlapping_square(SquareParticle sp);
void insert_square_in_cell(int index, SquareParticle sp);
void generate_random_squares();
static double get_angle(const SquareParticle sp);

#endif //SQUARE_CONFIG_H

