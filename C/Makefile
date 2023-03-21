# 				M A K E  F I L E

SOURCES= Rat.c Vertex.c Polynf.c Fibration.c bitlist.c fitness.c population.c evolution.c
OBJECTS= $(SOURCES:.c=.o)

LIBS = -lm

FLAGS = -g -O0

CC=cc 

gen_poly: gen_poly.o  $(OBJECTS)  Global.h
	$(CC) $(FLAGS) -o  gen_poly.x  gen_poly.o  $(OBJECTS) $(LIBS)

bitlist.o:      Global.h 
fitness.o:      Global.h 
population.o:   Global.h 
evolution.o:    Global.h 
Rat.o:          Global.h Rat.h
Vertex.o:       Global.h Rat.h
Polynf.o:       Global.h Rat.h
Fibration.o:    Global.h
gen_poly.o:     Global.h 
