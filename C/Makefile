# 				M A K E  F I L E

CC = gcc

LIBS = -lm
FLAGS = -g -O0

SOURCES = Rat.c Vertex.c Polynf.c bitlist.c fitness.c population.c evolution.c
OBJECTS = $(SOURCES:.c=.o)

all: poly_gen gen_poly

poly_gen: $(OBJECTS)
	$(CC) -shared -o poly_gen.so $(OBJECTS)

gen_poly: gen_poly.o  $(OBJECTS)  Global.h
	$(CC) $(FLAGS) -o  gen_poly.x  gen_poly.o  $(OBJECTS) $(LIBS)

bitlist.o:      Global.h 
fitness.o:      Global.h 
population.o:   Global.h 
evolution.o:    Global.h 
Rat.o:          Global.h Rat.h
Vertex.o:       Global.h Rat.h
Polynf.o:       Global.h Rat.h
gen_poly.o:     Global.h 
