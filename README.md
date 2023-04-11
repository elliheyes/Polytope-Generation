# Polytope-Generation

This package is designed to generate reflexive polytopes of a certain type using a genetic algorithm. 

Authors: Per Berglund, Yang-Hui He, Elli Heyes, Edward Hirst, Vishnu Jejjala, Andre Lukas

Notes:
- Currently the main function is set to randomly initialise a population, evolve it, 
  extract the terminal states (i.e. reflexive polytopes), reduce the list using the normal form,
  and then write the reduced list of terminal states to a file.
  
- You can change the global variables in the 'Global.h' file. This includes the polytope
  variables like the polytope dimension, the number of points, etc and also the 
  genetic algorithm variables like the number of generations, population size, etc.

- The fitness weights, which are defined in the 'Global.h' file, define what is meant
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
    

