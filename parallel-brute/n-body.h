/* N-Body problem: Sequential Brute Force method 
   blabla
*/
#ifndef _NBODY_HDR_
#define _NBODY_HDR_
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 768

#define PI 3.141592653589793238462643383279502884197

#define FRAND(x) (rand()/((double)RAND_MAX + 1) * x)

int RUN;
double EPSILON;
double G;

/********************************************************
                  DATA STRUCTURES
*********************************************************/

/* A simple struct, representing a point in the xy-plane */
typedef struct {
  double x, y;
} point_t;

/* velocity, force and direction are vectors, 
 * not points. this is just to make the code
 * easier to read. */
typedef point_t vector_t; 

/* a body has a position in the plane, velocity, force and mass */
typedef struct {
  point_t p;
  vector_t v;
//	vector_t f;
  double mass;
} body_t; 

/* all the bodies in the simulation. 
   the bodies are allocated on startup. */
body_t* bodies;

vector_t** forces;

/* the number of bodies */
int num_bodies;

/* the number of steps to do the simulation */
int num_steps;


int num_threads;
int num_arrived;

pthread_mutex_t time_mtx; // for protecting the shared simulation_time variable
pthread_mutex_t barrier_lock;
pthread_cond_t calculations_done;
pthread_t* threads;



double simulation_time; // total simulation time

/* initialize body properties */
void body_init(body_t* b, int id);

/* calculate forces for all bodies */
void compute_forces();

/* update velocities and move bodies */
void compute_positions();

void* work(void* data);

void barrier(int id, int mode);

#ifdef GRAPHICAL_SIMULATION
#include <SDL/SDL.h>
#include <GL/gl.h>

void gui_init();
void update_gfx();

#endif/*GRAPHICAL_SIMULATION*/

#endif /*_NBODY_HDR_*/
