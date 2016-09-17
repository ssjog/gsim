#
# Makefile for smaller C/C++ projects
# by Tobias Eriksson <tobier@kth.se>
#

#
# Define compiler and any compiler flags here.
#
# gcc for C projects, g++ for c++ projects.
#
CC = gcc
CFLAGS = -Wall -fopenmp
LDFLAGS= -lpthread -lm
#
# Define source-files and output here.
#
SRC= main.c physics.c
PROG = nbody

########################################
# Don't edit anything beyond this line.
########################################

OBJ= $(addsuffix .o, $(basename $(SRC)))

default: $(PROG)

$(PROG): $(SRC)
	$(CC) $(CFLAGS) -o $(PROG) $(SRC) $(LDFLAGS)

gfx: $(SRC)
	$(CC) $(CFLAGS) -o $(PROG) $(SRC) gui.c $(LDFLAGS) -lSDL -lGL -D GRAPHICAL_SIMULATION=1

clean:
	rm -f *.o $(PROG) *~
