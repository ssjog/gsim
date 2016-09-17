/* physics.c: Physics component of the simulation.
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"
#include <math.h>

/* create a new body */
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

	b->f.x = 0;
	b->f.y = 0;

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

void compute_forces()
{
  double distance, magnitude;
  size_t i, j;
	
	/* calculate the forces between all the bodies */
  for(i = 0; i < num_bodies - 1 ; i++) 
		/* forces follow Newtons third law of motion, so
		 * there is no need to calculate f_ji after f_ij */
    for(j = i+1; j < num_bodies; j++) { 
      /* calculate the distance between two bodies */
			double dx = bodies[j].p.x - bodies[i].p.x;
			double dy = bodies[j].p.y - bodies[i].p.y;

			/* epsilon is to avoid the very strong accelerations at close distances */
      distance = sqrt(dx * dx + dy * dy) + EPSILON; 

			/* Newtons law of universal gravitation */
      magnitude = G * bodies[i].mass * bodies[j].mass / (distance * distance); 
      bodies[i].f.x += magnitude * dx / distance; 
      bodies[j].f.x -= magnitude * dx / distance; 
      bodies[i].f.y += magnitude * dy / distance; 
      bodies[j].f.y -= magnitude * dy / distance;
    }
}

void compute_positions()
{
  int i;
  for(i = 0; i < num_bodies; i++){
    vector_t deltav = {bodies[i].f.x / bodies[i].mass, bodies[i].f.y / bodies[i].mass};
		
		// dp = (v + dv/2)
    point_t deltapos = {(bodies[i].v.x + deltav.x/2), (bodies[i].v.y + deltav.y/2)}; 

    /* update v */
    bodies[i].v.x += deltav.x;
    bodies[i].v.y += deltav.y;

    /* update position */
    bodies[i].p.x = bodies[i].p.x + deltapos.x;
    bodies[i].p.y = bodies[i].p.y + deltapos.y;

    /* and finally reset the force vector */
    bodies[i].f.x = 0.0;
    bodies[i].f.y = 0.0;
  }
}

