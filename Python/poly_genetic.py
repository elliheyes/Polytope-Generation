import ctypes as ct
from ctypes import *
import os
import subprocess

# define path to directory
path='.../Polytope-Generation-main/'
os.chdir(path+'C/')

#%% Define the global parameters

# polytope parameters
MIN = -2                            # minimum vertex coefficient value
MAXNVRTS = 5                        # maximum number of polytope vertices 
POLYDIM = 4                         # polytope dimension
BINLEN = 2                          # set length of binary number

# genetic algorithm parameters
POPSIZE = 400                       # size of population 
NUMGEN = 100                        # number of generations in each evolution
NUMCUTS = 1                         # number of cuts made to bitlist during crossover
METHOD = 1                          # selection method: 0 = ranking 1 = roulette
MUTRATE = 0.005                     # mutation rate
ALPHA = 3.0                         # value of alpha in selection probability
KEEPFITEST = 1                      # cary fitest individual into next population 
MONITOR = 1                         # print progress 

# fitenss parameters
NVERTS_WEIGHT = 0                   # search for polytopes with specific number of vertices
NVERTS = 5                          # number of vertices looking for
NPTS_WEIGHT = 0                     # search for polytopes with specific number of points
NPTS = 6                            # number of points looking for
H11_WEIGHT = 0                      # search for polytopes with specific h11
H11 = 1                             # h11 looking for
H12_WEIGHT = 0                      # search for polytopes with specific h12
H12 = 1                             # h12 looking for
H13_WEIGHT = 0                      # search for polytopes with specific h13
H13 = 1                             # h13 looking for
H22_WEIGHT = 0                      # search for polytopes with specific h22
H22 = 1                             # h22 looking for
EULER_WEIGHT = 0                    # search for polytopes with specific euler number
EULER = 1                           # euler number looking for

#%% Make and call the C library 

# read the Global.h file 
f = open("Global.h", "r")
lines = f.readlines()
f.close()

# edit the Global.h file 
f = open("Global.h", "w")
for i in range(len(lines)):
    line = lines[i]
    if "#define MIN" == line[:11]:
        line = "#define MIN " + str(MIN) + "\n"
    if "#define MAXNVRTS" == line[:16]:
        line = "#define MAXNVRTS " + str(MAXNVRTS) + "\n"
    if "#define POLYDIM" == line[:15]:
        line = "#define POLYDIM " + str(POLYDIM) + "\n"
    if "#define BINLEN" == line[:14]:
        line = "#define BINLEN " + str(BINLEN) + "\n"
    if "#define POPSIZE" == line[:15]:
        line = "#define POPSIZE " + str(POPSIZE) + "\n"
    if "#define NUMGEN" == line[:14]:
        line = "#define NUMGEN " + str(NUMGEN) + "\n"
    if "#define NUMCUTS" == line[:15]:
        line = "#define NUMCUTS " + str(NUMCUTS) + "\n"
    if "#define METHOD" == line[:14]:
        line = "#define METHOD " + str(METHOD) + "\n"
    if "#define MUTRATE" == line[:15]:
        line = "#define MUTRATE " + str(MUTRATE) + "\n"
    if "#define ALPHA" == line[:13]:
        line = "#define ALPHA " + str(ALPHA) + "\n"
    if "#define KEEPFITEST" == line[:18]:
        line = "#define KEEPFITEST " + str(KEEPFITEST) + "\n"
    if "#define MONITOR" == line[:15]:
        line = "#define MONITOR " + str(MONITOR) + "\n"
    if "#define NVERTS_WEIGHT" == line[:21]:
        line = "#define NVERTS_WEIGHT " + str(NVERTS_WEIGHT) + "\n"
    if "#define NVERTS " == line[:15]:
        line = "#define NVERTS " + str(NVERTS) + "\n"
    if "#define NPTS_WEIGHT" == line[:19]:
        line = "#define NPTS_WEIGHT " + str(NPTS_WEIGHT) + "\n"
    if "#define NPTS " == line[:13]:
        line = "#define NPTS " + str(NPTS) + "\n"
    if "#define H11_WEIGHT" == line[:18]:
        line = "#define H11_WEIGHT " + str(H11_WEIGHT) + "\n"
    if "#define H11 " == line[:12]:
        line = "#define H11 " + str(H11) + "\n"
    if "#define H12_WEIGHT" == line[:18]:
        line = "#define H12_WEIGHT " + str(H12_WEIGHT) + "\n"
    if "#define H12 " == line[:12]:
        line = "#define H12 " + str(H12) + "\n"
    if "#define H13_WEIGHT" == line[:18]:
        line = "#define H13_WEIGHT " + str(H13_WEIGHT) + "\n"
    if "#define H13 " == line[:12]:
        line = "#define H13 " + str(H13) + "\n"
    if "#define H22_WEIGHT" == line[:18]:
        line = "#define H22_WEIGHT " + str(H22_WEIGHT) + "\n"
    if "#define H22 " == line[:12]:
        line = "#define H22 " + str(H22) + "\n"
    if "#define EULER_WEIGHT" == line[:20]:
        line = "#define EULER_WEIGHT " + str(EULER_WEIGHT) + "\n"
    if "#define EULER " == line[:14]:
        line = "#define EULER " + str(EULER) + "\n"
    lines[i]=line
f.truncate(0)
f.writelines(lines)
f.close()

# make the c library
subprocess.run('make') 

# call the c library
so_file = 'poly_gen.so'
poly_gen = CDLL(so_file)

#%% Define data structures and functions

class pointlist(Structure):
    _fields_ = [("len", c_int),                                    # number of points
                ("points", c_int * MAXNVRTS * POLYDIM)]            # 2d array of points

class bitlist(Structure):
    _fields_ = [("len", c_int),                                    # length of bitlist
                ("bits", c_int * (MAXNVRTS * POLYDIM * BINLEN)),   # 1d array of bits
                ("fitness", c_float),                              # fitness of bitlist
                ("terminal", c_int)]                               # terminal (1) or not (0)
    
class population(Structure):
    _fields_ = [("size", c_int),                                   # population size
                ("maxfitness", c_float),                           # maximum fitness in population
                ("avfitness", c_float),                            # average fitness of population
                ("nterm", c_int),                                  # number of terminal states in population
                ("bl", bitlist * POPSIZE)]                         # list of bitlists 
    
class normal_form(Structure):
    _fields_ = [("Nv", c_int),                                     # number of vertices
                ("M", c_longlong * POLYDIM * 128)]                 # 2d array of normal form
    
# define the data types for function outputs and arguments
poly_gen.randomstate.restype = bitlist

poly_gen.pts2bts.restype = bitlist
poly_gen.pts2bts.argtypes = [pointlist]

poly_gen.bts2pts.restype = pointlist
poly_gen.bts2pts.argtypes = [bitlist]

poly_gen.normalform.restype = normal_form
poly_gen.normalform.argtypes = [bitlist]

poly_gen.bitlistsequiv.restype = c_int
poly_gen.bitlistsequiv.argtypes = [bitlist, bitlist]

poly_gen.randompop.restype = population
poly_gen.randompop.argtypes = [c_int]

poly_gen.nextpop.restype = c_void_p
poly_gen.nextpop.argtypes = [population, POINTER(population), c_int, c_int, c_int, c_float, c_float]

poly_gen.evolvepop.restype = POINTER(population * NUMGEN)
poly_gen.evolvepop.argtypes = [population, c_int, c_int, c_int, c_int, c_float, c_float, c_int]

poly_gen.termstates.restype = POINTER(bitlist)
poly_gen.termstates.argtypes = [POINTER(population * NUMGEN), c_int, POINTER(c_int)]

poly_gen.termstatesred.restype = POINTER(bitlist)
poly_gen.termstatesred.argtypes = [POINTER(population * NUMGEN), c_int, POINTER(c_int)]

poly_gen.removeredundancy.restype = c_void_p
poly_gen.removeredundancy.argtypes = [POINTER(bitlist), POINTER(c_int)]

poly_gen.searchenv.restype = POINTER(bitlist)
poly_gen.searchenv.argtypes = [c_int, c_int, c_int, c_int, c_int, c_int, c_float, c_float, c_int, POINTER(c_int)]

# define pthe python functions

# function that returns the bitlist of a randomly generated polytope  
def randomstate():
    return poly_gen.randomstate()

# function that, given a pointlist, returns the corresponding bitlist
def pts2bts(pl):
    return poly_gen.pts2bts(pl)

# function that, given a bitlist, returns the corresponding pointlist
def bts2pts(bl):
    return poly_gen.bts2pts(bl)

# function that, given a bitlist, returns the corresponding normal form 
def normalform(bl):
    return poly_gen.normalform(bl)

# function that determines whether two bitlists describe equivalent polytopes
def bitlistequiv(bl1,bl2):
    return poly_gen.bitlistsequiv(bl1,bl2)

# function that returns a randomly generated population
def randompop():
    return poly_gen.randompop(POPSIZE)

# function that, given a population, returns the next mutated population
def nextpop(pop):
    nextpop = randompop()
    nexpopptr = ct.pointer(nextpop)
    poly_gen.nextpop(pop, nexpopptr, METHOD, NUMCUTS, KEEPFITEST, MUTRATE, ALPHA)
    return nextpop

# function that, evolves an initial population over NUMGEN generations and returns the list of populations
def evolvepop(initialpop):
    return poly_gen.evolvepop(initialpop,NUMGEN,METHOD,NUMCUTS,KEEPFITEST,MUTRATE,ALPHA,MONITOR)

# function that returns a list of terminal states extracted from a list of populations
def termstates(pops):
    n = ct.c_int(0)
    nptr = ct.pointer(n)
    return poly_gen.termstates(pops,NUMGEN,nptr),nptr.contents.value

# function that returns a reduced list of terminal states extracted from a list of populations
def termstatesred(pops):
    n = ct.c_int(0)
    nptr = ct.pointer(n)
    return poly_gen.termstatesred(pops,NUMGEN,nptr),nptr.contents.value

# function that removes redundancy and returns a reduced list of bitlists
def removeredundancy(bl,n):
    return poly_gen.removeredundancy(bl,n)

# function that repeatedly evolves random initial populations, extracting terminal states, removing redundancy and writing to file
def searchenv(numrun):
    n = ct.c_int(0)
    nptr = ct.pointer(n)
    return poly_gen.searchenv(numrun, NUMGEN, POPSIZE, METHOD, NUMCUTS, KEEPFITEST, MUTRATE, ALPHA, MONITOR, nptr);
