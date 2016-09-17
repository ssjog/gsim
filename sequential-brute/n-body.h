#ifndef _NBODY_HDR_
#define _NBODY_HDR_

#include <stdlib.h>

#define PI 3.141592653589793238462643383279502884197
#define FRAND(x) (rand()/((double)RAND_MAX + 1) * x)

#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 768

double EPSILON;
double G;
int RUN;

/* velocity, force and direction are vectors, 
 * not points. this is just to make the code
 * a bit easier to read. */
typedef struct {
  double x;
	double y;
} point_t, vector_t;


/* A body has a position in the plane, velocity, force and mass */
typedef struct {
  point_t p;
  vector_t v;
	vector_t f;
  double mass;
} body_t;

body_t* bodies;
size_t num_bodies;
size_t num_steps;

/* initialize body properties */
void body_init(body_t* b, int id);

void compute_forces(); /* calculate all the forces in the system */
void compute_positions(); /* move all the bodies in the system */

#ifdef GRAPHICAL_SIMULATION
#include <SDL/SDL.h>
#include <GL/gl.h>

void gui_init();
void update_gfx();

#endif /*GRAPHICAL_SIMULATION*/

#endif /*_NBODY_HDR_*/
