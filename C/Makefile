# 				M A K E  F I L E

SOURCES= Rat.c Vertex.c bitlist.c fitness.c population.c evolution.c
OBJECTS= $(SOURCES:.c=.o)

CC=cc 

gen_poly: gen_poly.o  $(OBJECTS)  Global.h
	$(CC) -o  gen_poly.x  gen_poly.o  $(OBJECTS)

bitlist.o:      Global.h 
fitness.o:      Global.h 
population.o:   Global.h 
evolution.o:    Global.h 
Rat.o:          Global.h Rat.h
Vertex.o:       Global.h Rat.h
gen_poly.o:     Global.h 
