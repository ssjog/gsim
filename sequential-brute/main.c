/* Sequential Brute Force N-Body Simulation
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

int main(int argc, char* argv[]) 
{
  /* Set up command line arguments */
  if(argc > 2) {
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
  } else {
    printf("Error, invalid arguments.\nUsage: n-body num_bodies num_steps\n");
    return 1;
  }
  
  bodies = (body_t*)malloc(sizeof(body_t) * num_bodies);
  
  EPSILON = 20;
  G = 0.0667;
  RUN = 1;

  srand(time(NULL));
  
  /* initialize bodies */
  size_t i;  
  for (i = 0; i < num_bodies; i++)
    body_init(&bodies[i], i);
  
  printf("Sequential N-Body Problem\nby Tobias Eriksson <tobier@kth.se> and Sebastian Sjögren <ssjog@kth.se>\n\n");
  printf("The number of bodies are %zu, and the simulation will run for %zu steps.\n\n", num_bodies, num_steps);
  
  double simulation_time = 0.;
  
  #ifdef GRAPHICAL_SIMULATION
  gui_init();
  #endif

  /* Run the simulation for numSteps times */
  for (i = 0; i < num_steps; ++i) {
    double start = omp_get_wtime();
    compute_forces();
    compute_positions();
    double end = omp_get_wtime();
    
    simulation_time += (end - start);
    
    #ifdef GRAPHICAL_SIMULATION
    update_gfx();
    if (RUN == 0)
      break;
    #endif
  }
 
  printf("Execution time was %g seconds.\n", simulation_time);
  return 0;
}
