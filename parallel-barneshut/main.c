/* Parallel Barnes-Hut N-Body Simulation
   Authors: Tobias Eriksson <tobier@kth.se>, Sebastian Sjögren <ssjog@kth.se>
*/
#include "n-body.h"

int main(int argc, char* argv[])
{
  /* Set up command line arguments */
  if(argc > 4) {
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    num_threads = atoi(argv[3]);
		THETA = atof(argv[4]);
		if (THETA > 1.0) THETA = 1.0;
		if (THETA < 0.0) THETA = 0.0;
  } else {
    printf("Error, invalid arguments.\nUsage: n-body num_bodies num_steps numWorkers theta\n");
    return 1;
  }
  
  threads = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
  bodies = (body_t*)calloc(num_bodies, sizeof(body_t));
  
  srand(time(NULL));
  
  RUN = 1;
  EPSILON = 20;
  G = 0.0667;
  
  printf("Parallel (Barnes-Hut) N-Body Problem\nby Tobias Eriksson <tobier@kth.se> and Sebastian Sjögren <ssjog@kth.se>\n\n");
  printf("The number of bodies are %zu, the simulation will run for %zu steps using %zu threads. Thets is %g.\n\n", num_bodies, num_steps, num_threads, THETA);

  size_t i;
  for (i = 0; i < num_bodies; ++i) {
    body_init(&bodies[i], i);
  }
  
  pthread_mutex_init(&time_mtx, NULL);
  pthread_mutex_init(&barrier_lock, NULL);
  pthread_mutex_init(&tree_lock, NULL);
  pthread_cond_init(&calculations_done, NULL);
  pthread_cond_init(&tree_done, NULL);

  simulation_time = 0.0;
  
  for(i = 0; i < num_threads; i++)
    pthread_create(&threads[i], NULL, work, (void*)i);
  
  for (i = 0; i < num_threads; ++i)
    pthread_join(threads[i], NULL);
  
  free(bodies);
  free(threads);

  printf("The simulation time was %g seconds.\n", simulation_time);
  return 0;
}
