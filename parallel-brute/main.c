/* Parallel Brute Force N-Body Simulation
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"

int main(int argc, char** argv) 
{
  /* Set up command line arguments */
  if(argc > 3) {
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    num_threads = atoi(argv[3]);
  } else {
    printf("Error, invalid arguments.\nUsage: n-body num_bodies num_steps numWorkers\n");
    return 1;
  }

  /* Allocate the threads */
  threads = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);

  /* And allocate the bodies */
  bodies = (body_t*)malloc(sizeof(body_t) * num_bodies);

	/* Allocate force matrix */
	forces = (vector_t**)malloc(sizeof(vector_t*) * num_threads);

  size_t i;  
	for (i = 0; i < num_threads; ++i)
		forces[i] = (vector_t*)malloc(sizeof(vector_t) * num_bodies);

  
  srand(time(NULL));
  
  RUN = 1;
  EPSILON = 20;
  G = 0.0667;
  
  /* initialize bodies */
  for (i = 0; i < num_bodies; i++)
    body_init(&bodies[i], i);
  
  printf("Parallel N-Body Problem\nby Tobias Eriksson <tobier@kth.se> and Sebastian Sjögren <ssjog@kth.se>\n\n");
  printf("The number of bodies are %zu, the simulation will run for %zu steps using %zu threads.\n\n", num_bodies, num_steps, num_threads);
  
  simulation_time = 0.0;
  
  pthread_mutex_init(&time_mtx, NULL);
  pthread_mutex_init(&barrier_lock, NULL);
  pthread_cond_init(&calculations_done, NULL);
  
  for(i = 0; i < num_threads; i++)
    pthread_create(&threads[i], NULL, work, (void*)i);
  
    for (i = 0; i < num_threads; ++i) {
    pthread_join(threads[i], NULL);
  }
  
  printf("The total simulation time was %g seconds.\n", simulation_time);

	for (i = 0; i < num_threads; ++i)
		free(forces[i]);

  free(bodies);
  free(threads);
  
  return 0;
}
