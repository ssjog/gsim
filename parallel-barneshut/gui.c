/* gui.c: Graphical representation of the N-Body Simulation. 
   Author: Tobias Eriksson <tobier@kth.se>
   
   This is an optional graphical component for the N-Body problem
   that simply draws the bodies on the screen in each step of the simulation.

   This component depends on having the SDL and OpenGL libraries installed
   on the system. The conditional inclusion of this component is done with
   the preprocessor in the main-program code.
*/
#include "n-body.h"

void gui_init()
{
  // Initialize SDL
  const SDL_VideoInfo* info = NULL;
  SDL_Init(SDL_INIT_VIDEO);
  info = SDL_GetVideoInfo();
  int vidFlags = SDL_OPENGL | SDL_GL_DOUBLEBUFFER;
  if (info->hw_available) {vidFlags |= SDL_HWSURFACE;}
  else {vidFlags |= SDL_SWSURFACE;}
  SDL_SetVideoMode(SCREEN_WIDTH, SCREEN_HEIGHT, info->vfmt->BitsPerPixel, vidFlags);

  // Initialize OpenGL for 2d
  glViewport( 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
  glMatrixMode( GL_PROJECTION );
  glOrtho( 0, SCREEN_WIDTH, SCREEN_HEIGHT, 0, -1, 1 );
  glMatrixMode( GL_MODELVIEW );
  glDisable(GL_DEPTH_TEST);
}

void draw_body(GLfloat x, GLfloat y)
{
  glPushMatrix();
  glTranslatef(x, y, 0);
  glBegin(GL_TRIANGLE_FAN);
  glVertex2i(0, 0);
  glVertex2i(1, 0);
  glVertex2i(1, 1);
  glVertex2i(0, 1);
  glVertex2i(0, 0);
  glEnd();
  glPopMatrix();
}

void update_gfx()
{
	SDL_Event event;
  while(SDL_PollEvent(&event)){      
	  switch (event.type) {
      case SDL_QUIT:
			case SDL_KEYDOWN:
				RUN = 0;
				SDL_Quit();
				return;
    }
  }

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

  size_t i;
  for(i = 0; i < num_bodies; i++)
    draw_body(bodies[i].p.x, bodies[i].p.y);

  SDL_GL_SwapBuffers();
}
