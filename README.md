# Reflexive Polytope Generation

This package is designed to generate reflexive polytopes of a certain type using a genetic algorithm. 

Authors: Per Berglund, Yang-Hui He, Elli Heyes, Edward Hirst, Vishnu Jejjala, Andre Lukas


# Installation:
To install this package and use it on your own machine follow these simple steps:
- Step 0: If you do not have a C compiller and/or CMake installed on your machine start by downloading these.
- Step 1: Download the zipped file for this package from GitHub.
- Step 2: Unzip the file and move it somewhere convenient in your computer.
- Step 3: From the command line move to the file directory. 
- Step 4: Type 'make' and hit Enter. This should compile all the necessary files in the package.
- Step 5: To run the main function within the gen_poly.c file simply type './gen_poly.x' and hit Enter. 


# General Notes:
- Currently the main function in gen_poly.c is set to randomly initialise a population, evolve it, 
  extract the terminal states (i.e. reflexive polytopes), reduce the list using the normal form,
  and repeat this several times. The reduced list of output reflexive polytopes are written to a file 
  and the number of terminal states generated after each evolution, plus the total number of reduced
  terminal states are also written to a file.
  
- To specify the type of reflexive polytopes you would like to generate, such as the polytope dimension,
  number of points, etc, one can change the global variables in the Global.h file accordingly. 
  Here one can also change the global genetic algorithm variables like the number of generations, 
  population size, mutation rate etc. In fact, in order to find the best performance one usually needs to
  play around with these variables quite a bit.  

- The content in the files Polynf.c, Rat.c, Rat.h and Vertex.c have been copied from the source code of 
  PALP software. More information on this package can be found on the website: http://hep.itp.tuwien.ac.at/~kreuzer/CY/.
  
- The content of the bitlist.c, evolution.c and population.c have been copied from ...

# bitlist.c

The functions defined in this file deal with bitlists, including converting a polytope defined by 
integer vertex coordinates into a bitlist, generating a random bitlist state, determining whether 
two bitlists define the same polytope, printing a bitlist to a file, etc.


# population.c

The functions defined in this file deal with populations of bitlists, for instance generating a random
population, randomly the bitlists in a population mutating, ranking individuals in a population by their 
fitness, and updating a population to the next generation.


# fitness.c

This file contains the fitness function of the genetic algorithm, which determines how close a bitlist 
is from defining a reflexive polytope. 

There are several components of the fitness function that can be turned
on or off by the weight parameters. The two main components which should always remain on are the distance of 
the facets from the origin and the IP property. Their associated weights are defined in Global.h as DIST_WEIGHT
and IP_WEIGHT respectively. If DIST_WEIGHT > 0 then a penalty is added fitness for the average difference of the 
distance of the facets from the origin and 1. If IP_WEIGHT > 0 then a penalty is added to the fitness if the polytope
does not satisfy the interior point property.

On top of these two main components there are additional components one can include in the fitness to direct the 
search for reflexive polytopes with specific properties. To turn on these additional components, amend the global
parameters Global.h accordingly.

  x NVERTS_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope
    does not have NVERTS number of vertices. This causes the genetic algorithm to generate 
    polytopes of a given number of vertices.
    
  x NPOINTS_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope
    does not have NPTS number of points. This causes the genetic algorithm to generate polytopes 
    of a given number of points.
    
  x H11_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope does
    not have a h11 hodge number equal to H11. This causes the genetic algorithm to generate polytopes 
    of a given hodge number. Similarly for H12_WEIGHT, H13_WEIGHT, H22_WEIGHT, EULER_WEIGHT.


# evolution.c


# Vertex.c


# Rat.c & Rat.h


# Polynf.c


# Global.h


# gen_poly.c
    

