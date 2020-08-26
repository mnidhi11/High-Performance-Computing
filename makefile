SOURCES = wada.c rk4.c rk4.h pendulum.c pendulum.h          #Specified Source files
OBJECTS = wada.o rk4.o pendulum.o                           #Specified Object files
LIBS = -lm                                                  #Math library
CC = mpicc                                                    #Specified the Compiler
OPTIONS = -O -c                                             	#-c for compile only
CFLAGS = -Wall                              #Wall flag O3 flag for vectorization.

wada: $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LIBS) -o wada
#Specified the dependent header files for object files
wada.o: rk4.h                                              
wada.o: pendulum.h
rk4.o:rk4.h
pendulum.o:pendulum.h
pendulum.o:rk4.h

clean:
	@rm -rf $(OBJECTS) a.out core wada
