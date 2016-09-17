/* physics.c: Physics component of the simulation.
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"
#include <math.h>
#include <omp.h>

void compute_forces(int id)
{
  double distance, magnitude;
  int i, j, c = 0, alt_increment = 2 * (num_threads - id) - 1;
  
  for(i = id; i < num_bodies - 1 ; i += alt_increment) {

    /* "reverse stripes" thing */
    alt_increment = c++ % 2 ? 2 * id + 1 : 2 * (num_threads - id) - 1;
    
    double dx, dy;
    for(j = i + 1; j < num_bodies; j++) { 
      /* forces follow Newtons third law of motion, so
			 * there is no need to calculate f_ji after f_ij */
      
      /* calculate the distance between two bodies (and the direction of movement) */
      dx = bodies[i].p.x - bodies[j].p.x;
      dy = bodies[i].p.y - bodies[j].p.y;
			
			/* epsilon is to avoid the very strong accelerations at close distances */
      distance = sqrt(dx * dx + dy * dy) + EPSILON; 
      
			/* Newtons law of universal gravitation */
      magnitude = (G * bodies[i].mass * bodies[j].mass) / (distance * distance); 
      
			forces[id][i].x -= magnitude * dx / distance;
			forces[id][j].x += magnitude * dx / distance;
			forces[id][i].y -= magnitude * dy / distance;
			forces[id][j].y += magnitude * dy / distance;
    }
  }
}

void compute_positions(int id)
{
  size_t i, k;
	int c = 0, alt_increment = 2 * (num_threads - id) - 1;
  
  for(i = id; i < num_bodies; i += alt_increment) {
    /* "reverse stripes" thing */
    alt_increment = c++ % 2 ? 2 * id + 1: 2 * (num_threads - id) - 1;


		// sum the forces on body i and reset force column
		vector_t f = {0., 0.};    
		for (k = 0; k < num_threads; ++k) {
			f.x += forces[k][i].x;
			f.y += forces[k][i].y;
			
			forces[k][i].x = 0;
			forces[k][i].y = 0;
		}

		// dv = f/m
//    vector_t deltav = {bodies[i].f.x / bodies[i].mass, bodies[i].f.y / bodies[i].mass}; 
    vector_t deltav = {f.x / bodies[i].mass, f.y / bodies[i].mass}; 

    
		// dp = (v + dv/2)
		point_t deltapos = {(bodies[i].v.x + deltav.x / 2), (bodies[i].v.y + deltav.y/2)};
    
		/* update velocity */
    bodies[i].v.x = bodies[i].v.x + deltav.x;
    bodies[i].v.y = bodies[i].v.y + deltav.y;

    /* update position */
    bodies[i].p.x = bodies[i].p.x + deltapos.x;
    bodies[i].p.y = bodies[i].p.y + deltapos.y;

    /* and finally reset the force vector */
//    bodies[i].f.x = 0.0;
//    bodies[i].f.y = 0.0;
  }
}

void barrier(int id, int mode) {
  pthread_mutex_lock(&barrier_lock);
  ++num_arrived;
  if (num_arrived == num_threads) {
    num_arrived = 0;
    pthread_cond_broadcast(&calculations_done);
  } else {
    pthread_cond_wait(&calculations_done, &barrier_lock);
  }

	#ifdef GRAPHICAL_SIMULATION
	/* thread 0 updates graphics */
	if (id == 0 && mode == 1) update_gfx();

	if (RUN == 0 && mode == 2) {
		pthread_mutex_unlock(&barrier_lock);
		pthread_exit(NULL);
	}
	#endif

  pthread_mutex_unlock(&barrier_lock);
}

void* work(void* data)
{
  int id = (int)data;
  double local_simulation_time = 0.;

	/* we let thread 0 be responsible for graphics */
  #ifdef GRAPHICAL_SIMULATION
  if (id == 0)
    gui_init();
  #endif

  size_t i;
  for(i = 0; i < num_steps; ++i) {
    double s1 = omp_get_wtime(); 
    compute_forces(id);
    double e1 = omp_get_wtime(); 
    barrier(id, 0);
    
    double s2 = omp_get_wtime(); 
    compute_positions(id);
    double e2 = omp_get_wtime(); 
    barrier(id, 1);
    
    local_simulation_time += (e1 - s1) + (e2 - s2);

    #ifdef GRAPHICAL_SIMULATION
    /* update_gfx polls SDL events. if got "quit program"-event,
     * global variable RUN is set to zero but if we just exit thread
     * here, other threads will deadlock. a barrier is needed. */
    barrier(id, 2);
    #endif
  }

  pthread_mutex_lock(&time_mtx);
	if (simulation_time < local_simulation_time)
	  simulation_time = local_simulation_time;
  pthread_mutex_unlock(&time_mtx);
  
  pthread_exit(NULL);
}

/* initialize body properties */
void body_init(body_t* b, int id) {
	double max_radius = 100;

	double alpha = 0.4;
	double radius =  (1 - alpha) * max_radius + alpha * FRAND(max_radius);

	double v = FRAND(2 * PI);
	double x = radius * cos(v);
	double y = radius * sin(v);	
	point_t p = {x + SCREEN_WIDTH / 2, y + SCREEN_HEIGHT / 2};

	b->mass = 200 / radius;
	b->p = p;

//	b->f.x = 0;
//	b->f.y = 0;

	double dx = p.x - SCREEN_WIDTH / 2;
	double dy = p.y - SCREEN_HEIGHT / 2;
	
	double angle = atan(dy / dx);
	double cos_beta = cos(angle + PI / 4);
	double sin_beta = sin(angle + PI / 4);

	double rdiff = (max_radius - radius) / max_radius * 2;
	rdiff *= rdiff;

	if (dx > 0) {
    b->v.x = rdiff * cos_beta;
    b->v.y = rdiff * sin_beta;
	} else {
		b->v.x = rdiff * (-cos_beta);
		b->v.y = rdiff * (-sin_beta);
	}
}

