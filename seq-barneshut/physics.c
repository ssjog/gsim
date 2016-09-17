/* physics.c: Physics component of the simulation.
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>

inline int in_which_quadrant(rect_t parent, point_t p) {
	double vertical_center = parent.y1 + (parent.y2 - parent.y1) / 2;
	double horizontal_center = parent.x1 + (parent.x2 - parent.x1) / 2;
	
	if (p.x >= parent.x1 && p.x < horizontal_center) {
		if (p.y >= parent.y1 && p.y < vertical_center)
			return 2;
		if (p.y >= vertical_center && p.y < parent.y2)
			return 0;
	} else if (p.x >= horizontal_center && p.x < parent.x2) {
		if (p.y >= parent.y1 && p.y < vertical_center)
			return 3;
		if (p.y >= vertical_center && p.y < parent.y2)
			return 1;
	}
	
	return -1;
}


int node_count = 0, body_count = 0;

void body_init(body_t* b, int id) {
	double max_radius = 100;

	double alpha = 0.4;
	double radius =  (1 - alpha) * max_radius + alpha * FRAND(max_radius);

	double v = FRAND(2 * PI);
	double x = radius * cos(v);
	double y = radius * sin(v);	
	point_t p = {x + SCREEN_WIDTH / 2, y + SCREEN_HEIGHT / 2};

	b->id = id;
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

node_t* new_node(rect_t q) {
	node_t* n = (node_t*)calloc(1, sizeof(node_t));
	n->q = q;
	n->id = node_count++;
	return n;
}


rect_t get_appropriate_quadrant(size_t quad_index, const node_t* parent_node) {
	rect_t parent = parent_node->q;

	double vertical_center = parent.y1 + (parent.y2 - parent.y1) / 2;
	double horizontal_center = parent.x1 + (parent.x2 - parent.x1) / 2;

	if (quad_index == 0) {
		rect_t q0 = { parent.x1, vertical_center, horizontal_center, parent.y2 };
		return q0;
	} else if (quad_index == 1) {
		rect_t q1 = { horizontal_center, vertical_center, parent.x2, parent.y2 };
		return q1;
	} else if (quad_index == 2) {
		rect_t q2 = { parent.x1, parent.y1, horizontal_center, vertical_center };
		return q2;
	} else if (quad_index == 3) {
		rect_t q3 = { horizontal_center, parent.y1, parent.x2, vertical_center };
		return q3;
	}

	fprintf(stderr, "weird trouble\n");
	exit(1);
}

void insert(const body_t* b, node_t* n) {
	// node has children => node is internal
	if (n->c) {
		if (DEBUG) printf("%d node is internal\n", b->id);

		// determine in what quadrant body should be placed.
		int quad_index = in_which_quadrant(n->q, b->p);
		if (quad_index != -1) {
			if (DEBUG) printf("%d calling insert on quadrant %d\n", b->id, quad_index);
			
			/* we only allocate child nodes when it's neccesary, in order to avoid
			 * excessive memory use */
			if (n->c[quad_index] == NULL) {
				n->c[quad_index] = new_node(get_appropriate_quadrant(quad_index, n));
			}
			
			/* now call insert again but with the child node corresponding to
			 * the quadrant matching the body's position */
			insert(b, n->c[quad_index]);
		} else {
			fprintf(stderr, "insert: point_t not inside parent quadrant\n");
		}

	// node has no children => node is external/leaf
	} else {

		// if there's no a body here already, just place given body here and be done
		if (n->b == NULL) {
			if (DEBUG) printf("%d empty leaf, inserting body\n", b->id);
			n->b = (body_t*)b;
		
		// otherwise, we have to split this leaf node
		} else {
			if (DEBUG) printf("%d not empty leaf, splitting\n", b->id);
		
			n->c = (node_t**)calloc(4, sizeof(node_t*));

			/* there's a body here already, this has to be moved down the tree now, since
			 * internal nodes can't contain bodies */

			int existing_body_quadrant_index = in_which_quadrant(n->q, n->b->p);			
			int body_quadrant_index = in_which_quadrant(n->q, b->p);
	
			if (DEBUG) printf("%d existing body %d moves into quadrant %d\n", b->id, n->b->id, existing_body_quadrant_index);
			body_t* existing_body = n->b;
			n->b = NULL;

			/* again, we only allocate child nodes when they're needed.
			 * this is the child node where we put the body that was previously _here_ */			
			if (n->c[existing_body_quadrant_index] == NULL) {
				n->c[existing_body_quadrant_index] = new_node(get_appropriate_quadrant(existing_body_quadrant_index, n));
			}

			/* insert again but for the child node for the body that was _here_ */
			insert(existing_body, n->c[existing_body_quadrant_index]);

			if (DEBUG) printf("%d body moves into quadrant %d\n", b->id, body_quadrant_index);
			
			/* ..and the same thing for new body */
			if (n->c[body_quadrant_index] == NULL) {
				n->c[body_quadrant_index] = new_node(get_appropriate_quadrant(body_quadrant_index, n));
			}

			insert(b, n->c[body_quadrant_index]);			
		}
	}
}

void compute_com(node_t* n) {
	// if node has children	
	if (n->c) {
		double com_x = 0;
		double com_y = 0;
		double mass_sum = 0;

		size_t i;

		for (i = 0; i < 4; ++i)	{
			if (n->c[i]) {
				// call compute_com for those children recursively
				compute_com(n->c[i]);

				// at each branch in three compute the center of mass + mass sum
				com_x += n->c[i]->com.x * n->c[i]->m;
				com_y += n->c[i]->com.y * n->c[i]->m;
				mass_sum += n->c[i]->m;
			}
		}

		// set mass
		n->m = mass_sum;
	
		// set center of mass
		n->com.x = com_x / mass_sum;
		n->com.y = com_y / mass_sum;
	} else {
	/* node is a leaf, center of mass is body center	
	 * the barnes-hut algo guarantees that there's a body
   * so NULL-check isn't needed */
		n->m = n->b->mass;
		n->com = n->b->p;
	}
}

void delete_tree(node_t** children) {
	// children is a nodes children but we don't know how many
	if (children) {
		size_t i;
		// call delete_tree recursively on the children
		for (i = 0; i < 4; ++i) {
			if (children[i])
				delete_tree(children[i]->c);
	
			// at each branch in three, free nodes
			if (children[i]) {
				free(children[i]);
				children[i] = NULL;
			}			
		}
		free(children);
		children = NULL;
	}
}

rect_t minimal_square_cover(const body_t* bodies, size_t len) {
	double max_x = bodies[0].p.x;
	double max_y = bodies[0].p.y;
	double min_x = max_x;
	double min_y = max_y;

	size_t i;
	for (i = 1; i < len; ++i) {
		if (bodies[i].p.x > max_x)
			max_x = bodies[i].p.x;		

		if (bodies[i].p.x < min_x)
			min_x = bodies[i].p.x;

		if (bodies[i].p.y > max_y)
			max_y = bodies[i].p.y;
		
		if (bodies[i].p.y < min_y)
			min_y = bodies[i].p.y;
	}

	/* modify the covering rectangle so that it becomes
 	 * a square with the side equal to the largest of side of the rectangle */
	double diff = ((max_x - min_x) - (max_y - min_y)) / 2.;
	rect_t cover;
	if (diff > 0.) {
		rect_t q = { min_x - 1., min_y - diff - 1., max_x + 1., max_y + diff + 1.};
		cover = q;
	} else {
		rect_t q = { min_x + diff - 1., min_y - 1., max_x - diff + 1., max_y + 1. };
		cover = q;
	}
		
	return cover;
}

void compute_force(body_t* b, const node_t* n) {	
	// node has children => node is internal
	if (n->c) {
		// find the direction of movement (and distance)
		double dx = n->com.x - b->p.x;
		double dy = n->com.y - b->p.y;

		/* pythagoras + fix to prevent possible division by zero 
		 * and extreme accelerations for close bodies */
		double distance_to_com = sqrt(dx * dx + dy * dy) + EPSILON;

		// x2-x1 = width since node quadrants always are squares
		double quadrant_width = n->q.x2 - n->q.x1;
		double theta = quadrant_width / distance_to_com;

		/* if the quadrant width - distance ratio is below THETA
 		 * approximate the force excerted using center of mass for given node */
		if (theta < THETA) {
			// Newtons law of universal gravitation
      double magnitude = (G * b->mass * n->m) / (distance_to_com * distance_to_com); 
      b->f.x += magnitude * dx / distance_to_com; 
      b->f.y += magnitude * dy / distance_to_com; 			

		/* otherwise traverse down the nodes children and do the samething there */
		} else {
			size_t i;
			for (i = 0; i < 4; ++i) {
				if (n->c[i])
					compute_force(b, n->c[i]);
			}
		}

	/* node is external/leaf which means we compute the force the normal way */
	} else {
    // calculate the distance between two bodies
		double dx = n->b->p.x - b->p.x;
		double dy = n->b->p.y - b->p.y;

		// epsilon is to avoid the very strong accelerations at close distances
    double distance = sqrt(dx * dx + dy * dy) + EPSILON; 

		// Newtons law of universal gravitation
    double magnitude = (G * b->mass * n->b->mass) / (distance * distance); 
    b->f.x += magnitude * dx / distance; 
    b->f.y += magnitude * dy / distance; 
	}
}

void move_body(body_t* b) {
	// dv = f/m*DT
  vector_t deltav = {b->f.x / b->mass, b->f.y / b->mass}; 

	// dp = (v + dv/2)*DT
  point_t deltapos = {(b->v.x + deltav.x / 2), (b->v.y + deltav.y / 2)}; 

  /* update velocity */
  b->v.x = b->v.x + deltav.x;
  b->v.y = b->v.y + deltav.y;

  /* update position */
  b->p.x = b->p.x + deltapos.x;
  b->p.y = b->p.y + deltapos.y;

  /* and finally reset the force vector */
  b->f.x = 0.0;
  b->f.y = 0.0;
}

void compute_tree(const body_t* bodies, size_t num_bodies, node_t** tree) {
	rect_t space = minimal_square_cover(bodies, num_bodies);
	*tree = new_node(space);

	size_t i;
	for (i = 0; i < num_bodies; ++i) {
		insert(&bodies[i], *tree);
	}

	compute_com(*tree);
}

void compute_forces(body_t* bodies, size_t num_bodies, const node_t* tree) {
	size_t i;
	for (i = 0; i < num_bodies; ++i) {
		compute_force(&bodies[i], tree);
	}
}

void compute_positions(body_t* bodies, size_t num_bodies) {
	size_t i;
	for (i = 0; i < num_bodies; ++i) {
		move_body(&bodies[i]);
	}
}
