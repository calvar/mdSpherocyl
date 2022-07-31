
PROGS	=	DSCmd PDSCmd graphic

HEADS	=	Matrix.hpp classes.hpp functions.hpp jacobi.h
CCF	=	classes.cpp functions.cpp jacobi.c

OPTS	=	-Wall -g  # -m64 -O3 -fno-inline 
LIBS	=	-lgsl -lgslcblas -lgomp
LIBSGL	=	-lGL -lGLU -lSDL

CC	=	g++


DSCmd: main.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ main.cpp $(CCF) $(LIBS)

PDSCmd: main.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ main.cpp $(CCF) $(LIBS) -fopenmp

TimeCorrel: TimeCorrel.cpp Matrix.hpp
	$(CC) $(OPTS) -o $@ TimeCorrel.cpp

confUpdate: confUpdate.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ confUpdate.cpp $(CCF) $(LIBS)

graphic: graphic.cpp Matrix.hpp gl2ps.h
	$(CC) $(OPTS) -w -o $@ graphic.cpp gl2ps.c $(LIBSGL)


clean: 
	rm $(PROGS) Error* Etest* TabBC* track* core* trayect* dist*
