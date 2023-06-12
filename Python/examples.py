from poly_genetic import *
import numpy as np
                        
# generate a random state and look at it features
random_bitlist = randomstate()
print("Random State: ")
print("Fitness: ",random_bitlist.fitness)
print("Terminal: ",random_bitlist.terminal)
print("Bits:")
print(random_bitlist.bits[:])

# convert bitlist into point list
random_pointlist = bts2pts(random_bitlist)
print("Number of Points: ",random_pointlist.len)
print("Poiints:")
M = []
for i in range(POLYDIM):
    M.append(random_pointlist.points[i][:])
print(np.transpose(M))

# convert point list back into bitlist
random_bitlist_2 = pts2bts(random_pointlist)

# check whether the original bitlist and the converted bitlist are equivalent
print("Equivalent: ",bitlistequiv(random_bitlist,random_bitlist_2))

# generate a random population and look at it features
random_pop = randompop()
print("Random Population: ")
print("Maximum Fitness: ",random_pop.maxfitness)
print("Average Fitness: ",random_pop.avfitness)
print("Number of Terminal States: ",random_pop.nterm)

# evolve the random population
evol = evolvepop(random_pop)

# get reduced list of terminal states
bl1 = termstatesred(evol)

# run a search for 10 populations
bl2 = searchenv(10, 1)

