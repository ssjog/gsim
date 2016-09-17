/* Sequential Barnes-Hut N-Body Simulation
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

int main(int argc, char* argv[])
{
  /* Set up command line arguments */
  if(argc > 3) {
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
		THETA = atof(argv[3]);
		if (THETA > 1.0) THETA = 1.0;
		if (THETA < 0.0) THETA = 0.0;
  } else {
    printf("Error, invalid arguments.\nUsage: n-body num_bodies num_steps theta\n");
    return 1;
  }

  body_t* bodies = (body_t*)calloc(num_bodies, sizeof(body_t));
  
  EPSILON = 20;
  G = 0.0667;
  RUN = 1;
	
  srand(time(NULL));
  
  size_t i;
  for (i = 0; i < num_bodies; ++i)
    body_init(&bodies[i], i);

  printf("Sequential (Barnes-Hut) N-Body Problem\nby Tobias Eriksson <tobier@kth.se> and Sebastian Sjögren <ssjog@kth.se>\n\n");
  printf("The number of bodies are %zu, and the simulation will run for %zu steps. Theta is %g.\n\n", num_bodies, num_steps, THETA);
	
  #ifdef GFX
  gui_init();
  #endif

  node_t* tree;

  double simulation_time = 0.0;

  while (--num_steps) {
    double start_time = omp_get_wtime(); 
    compute_tree(bodies, num_bodies, &tree);
    compute_forces(bodies, num_bodies, tree);
    compute_positions(bodies, num_bodies);
    
    delete_tree(tree->c);
    free(tree);
    
    double end_time = omp_get_wtime(); 

    #ifdef GFX
    update_gfx(bodies, num_bodies);
    if (RUN == 0)
      break;
    #endif
    
    simulation_time += (end_time - start_time); 
  }
  
  free(bodies);

  printf("The simulation time was %g seconds.\n", simulation_time);
  return 0;
}
