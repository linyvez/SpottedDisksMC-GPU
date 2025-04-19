 #ifndef CIRCLE_MOVE_H
 #define CIRCLE_MOVE_H
 #include "circle_config.h"

 typedef struct {
     double x;
     double y;
 } Coordinates;

 void rotate_circle(CircleParticle *c);
 Coordinates calculate_new_coordinates(const CircleParticle *c);
 void animate_circle_movement(int steps, int totalCircles);

 #endif // CIRCLE_MOVE_H