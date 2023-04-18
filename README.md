# Polytope-Generation

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
  
  

- The fitness weights, which are defined in the Global.h file, define what is meant
  by a terminal state.
    
  x DIST_WEIGHT: if this is > 0 then a penalty is added to the fitness for the average 
    distance of the facets from the origin.
    
  x IP_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope
    does not satisfy the interior point property.

  x NVERTS_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope
    does not have NVERTS number of vertices. NVERTS is also defined in 'Global.h'. This 
    causes the genetic algorithm to generate polytopes of a given number of vertices.
    
  x NPOINTS_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope
    does not have NPTS number of points. NPTS is also defined in 'Global.h'. This 
    causes the genetic algorithm to generate polytopes of a given number of points.
    
  x H11_WEIGHT: if this is > 0 then a penalty is added to the fitness if the polytope does
    not have a h11 hodge number equal to H11 (defined in 'Global.h'). This causes the 
    genetic algorithm to generate polytopes of a given hodge number. Similarly for H12_WEIGHT,
    H13_WEIGHT, H22_WEIGHT, EULER_WEIGHT.
    

