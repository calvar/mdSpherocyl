
PROGS	=	DSCmd PDSCmd graphic

HEADS	=	Matrix.hpp classes_mc.hpp functions_mc.hpp jacobi.h
CCF	=	classes_mc.cpp functions_mc.cpp jacobi.c

OPTS	=	-Wall -O3  # -m64 -g -fno-inline 
LIBS	=	-lgsl -lgslcblas -lgomp
LIBSGL	=	-lGL -lGLU -lSDL

CC	=	g++


DSCmd: main_mc.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ main_mc.cpp $(CCF) $(LIBS)

PDSCmd: main_mc.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ main_mc.cpp $(CCF) $(LIBS) -fopenmp

TimeCorrel: TimeCorrel.cpp Matrix.hpp
	$(CC) $(OPTS) -o $@ TimeCorrel.cpp

confUpdate: confUpdate.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ confUpdate.cpp $(CCF) $(LIBS)

graphic: graphic.cpp Matrix.hpp gl2ps.h
	$(CC) $(OPTS) -w -o $@ graphic.cpp gl2ps.c $(LIBSGL)


clean: 
	rm $(PROGS) Error* Etest* TabBC* track* core* trayect* dist*
