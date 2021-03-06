#ifndef _NBODY_HDR_
#define _NBODY_HDR_

#ifdef GFX
#include <SDL/SDL.h>
#include <GL/gl.h>
#endif

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>

#define FRAND(x) (rand()/((double)RAND_MAX + 1) * x)
#define DEBUG 0

#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 768

#define PI 3.141592653589793238462643383279502884197

double THETA;
double EPSILON;
double G;
int RUN;

/** DATA STRUCTURES **/

typedef struct {
  double x;
	double y;
} point_t, vector_t;

typedef struct {
	double x1;
	double y1;
	double x2;
	double y2;
} rect_t;

typedef struct {
	int id;
  double mass;
  point_t p;
  vector_t v;
  vector_t f;
} body_t;

struct node {
	int id;
  body_t* b; 
  double m;
  point_t com;
  rect_t q;
  struct node** c;
};

typedef struct node node_t;

int num_steps;
size_t num_threads;
size_t num_bodies;

int num_arrived;

double simulation_time;

pthread_mutex_t time_mtx;
pthread_mutex_t barrier_lock;
pthread_mutex_t tree_lock;

pthread_cond_t calculations_done;
pthread_cond_t tree_done;

pthread_t* threads;

body_t* bodies;
node_t* tree;


/* initialize body properties */
void body_init(body_t* b, int id);

/* allocates a new node on heap memory, initializes it's
 * quadrant with given rectangle */
node_t* new_node(rect_t q);

/* given a the child vector of a node (preferably tree root)
 * recursively traverse the tree and deallocate the nodes */
void delete_tree(node_t** children);


#ifdef GFX
/* setup a gui window with SDL/OpenGL */
void gui_init();

/* draw something in OpenGL */
void draw_body(GLfloat x, GLfloat y);

/* refreshes GL buffers and also listens for SDL events
 * that exits the program */
void update_gfx();
#endif

/* given a rectangle and a point and assuming 
 * the following "quadrant indexes":
 *
 * north west quadrant has index 0,
 * north east quadrant has index 1,
 * south west quadrant has index 2,
 * south east quadrant has index 3,
 * 
 * return the index for which quadrant the given point lies in.
 */
int in_which_quadrant(const rect_t parent, point_t p);

/* given a "quadrant index" and a node, create and return 
 * the quadrant in that node that matches given quadrant index. */
rect_t get_appropriate_quadrant(size_t quad_index, const node_t* parent_node);

/* given a set of bodies, compute the smallest square covering those bodies.
 * 
 * (this is done because if one uses a statically sized root rectangle when (re)building
 * the barnes-hut tree, some bodies might have travelled outside it's bounds and would
 * thus be excluded from the simulation!) */
rect_t minimal_square_cover(const body_t* bodies, size_t len);

/* given a set of bodies, compute a barnes-hut tree */
void compute_tree(const body_t* b, size_t num_bodies, node_t** tree);

/* computes force applied by a body or a group of bodies to given body
 * according to barnes-hut algorithm */
void compute_force(body_t* b, const node_t* n);

/* given an interval of the set of bodies (which is globally accessed)
 * compute the forces applied to them by other bodies/groups of bodies */
void compute_forces(size_t first, size_t last);

/* based on newly computed forces, compute a bodies next position.
 * reset the force when done */
void move_body(body_t* b);

/* given an interval of the set of bodies (which is globally accessed)
 * compute their positions */
void compute_positions(size_t first, size_t last);

/* recursively computes center of mass for all nodes in
 * a barnes-hut tree */
void compute_com(node_t* n);

/* inserts a body in the barnes-hut tree according to description given at
 * http://www.cs.princeton.edu/courses/archive/fall03/cs126/assignments/barnes-hut.html */
void insert(const body_t* b, node_t* tree);

/* barrier used to force all threads to wait until thread 0 has (re)built the
 * barnes-hut tree */
void barrier_build_tree(int id);

/* general barrier used to syncronize threads */
void barrier(int id, int mode);
/* this is basically the "main"-loop, adapted to be run
 * in parallel by N number of threads. the work is
 * sliced up in equally sized segments */
void* work(void* data);


#endif /*_NBODY_HDR_*/
