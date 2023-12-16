# author: elli heyes (elli.heyes@city.ac.uk)
# date last edited: 15/12/2023
import copy
import random
import numpy as np
from cytools import Polytope
from fitness import fitness


class state:
    """"Polytope state object
         .points is the points matrix of the polytope
         .dim is the dimension of the polytope
         .max_coeff is the absolute maximum coefficient value
         .fitness is the fitness of the triangulation
         .terminal is 1 if the state is good and 0 otherwise"""
    def __init__(self, points, dim, max_coeff):
        self.points = points
        self.dim = dim
        self.max_coeff = max_coeff
        self.max_verts = len(points)
        
        self.fitness = 0.
        self.terminal = 0
    
    def compute_fitness(self, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
        self.fitness = fitness(self.points, h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
        if self.fitness == 0:
            self.terminal = 1
        else:
            self.terminal = 0


class population:
    """"Population:
        .state_list is the list of states
        .pop_size is the population size
        .max_fitness is the maximum fitness score
        .av_fitness is the average fitness score
        .num_term is the number of terminal states"""
    def __init__(self, state_list):        
        self.state_list = state_list
        self.pop_size = len(state_list)
        
        self.max_fitness = 0.
        self.av_fitness = 0.
        self.num_term = 0.
        
        self.compute_max_fitness()
        self.compute_av_fitness()
        self.compute_num_term()
        
    def compute_max_fitness(self):
        max_fit = self.state_list[0].fitness
        for i in range(1,self.pop_size):
            if self.state_list[i].fitness > max_fit:
                max_fit = self.state_list[i].fitness
        self.max_fitness = max_fit 
    
    def compute_av_fitness(self):
        sum_fit = 0
        for i in range(self.pop_size):
            sum_fit += self.state_list[i].fitness
        self.av_fitness = sum_fit / self.pop_size
        
    def compute_num_term(self):
        total_term = 0
        for i in range(self.pop_size):
            if self.state_list[i].terminal:
                total_term += 1
        self.num_term = total_term
        
    def mutate_pop(self, mut_rate, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
        # determine the number of mutations to perform
        num_points = self.pop_size*self.state_list[0].max_verts*self.state_list[0].dim
        num_mut = round(num_points*mut_rate)
    
        for i in range(num_mut):
            # determine a random state in the population
            state_pos = random.choice(range(self.pop_size))
            
            # determine a random point in the polytope
            point_pos = random.choice(range(len(self.state_list[state_pos].points)))
            
            # determine a random coordinate in the point
            coord_pos = random.choice(range(1,len(self.state_list[state_pos].points[point_pos])))
            
            # mutate the point
            self.state_list[state_pos].points[point_pos][coord_pos] = random.choice(range(-self.state_list[state_pos].max_coeff,self.state_list[state_pos].max_coeff))
             
            # update fitness
            self.state_list[state_pos].compute_fitness(h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
        
        # update max fitness, average fitness and number of terminal states
        self.compute_max_fitness()
        self.compute_av_fitness()
        self.compute_num_term()
        

def crossover(state1, state2, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Cross two states."""
    # determine a random point in the polytopes
    point_pos = random.choice(range(state1.max_verts))
            
    # determine a random coordinate in the point
    coord_pos = random.choice(range(state1.dim))
    
    # swap the relevant parts of the bit lists
    new_points1 = copy.deepcopy(state1.points)
    new_points2 = copy.deepcopy(state2.points)
    for i in range(point_pos,len(new_points1)):
        if i == point_pos:
            start = coord_pos
        else:
            start = 0
        for j in range(start,len(new_points1[i])):
            new_points1[i][j] = state2.points[i][j]
            new_points2[i][j] = state1.points[i][j]
    
    # define the new states
    new_state1 = state(points=new_points1, dim=state1.dim, max_coeff=state1.max_coeff)
    new_state2 = state(points=new_points2, dim=state2.dim, max_coeff=state2.max_coeff)
    
    # update the fitness
    new_state1.compute_fitness(h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
    new_state2.compute_fitness(h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)

    return new_state1, new_state2
  

def random_state(dim, max_verts, max_coeff, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Generate a random polytope state."""
    points = []
    for i in range(max_verts):
        point = []
        for j in range(dim):
            point.append(random.choice(range(-max_coeff,max_coeff)))
        points.append(point)
        
    S = state(points=points, dim=dim, max_coeff=max_coeff)
    S.compute_fitness(h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
    
    return S


def random_pop(pop_size, dim, max_verts, max_coeff, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Generate a random population."""
    state_list = []
    for i in range(pop_size):
        state_list.append(random_state(dim, max_verts, max_coeff, h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav))
    
    pop = population(state_list)
    pop.compute_max_fitness()
    pop.compute_av_fitness()
    pop.compute_num_term()
    
    return pop


def sort_pop(pop):
    """"Sort a population based on the fitness scores."""
    fitness_scores = [pop.state_list[i].fitness for i in range(pop.size)]
    sorted_indices = np.argsort(np.array(fitness_scores))
    sorted_state_list = [pop.state_list[i] for i in sorted_indices]
    return population(sorted_state_list)


def select_and_cross(pop, p, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Select two pairs of individuals from the population and cross them."""
    indices = random.choices(range(pop.pop_size), p, k=2)
    state1, state2 = crossover(pop.state_list[indices[0]], pop.state_list[indices[1]], h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
    return [state1, state2]
    

def next_pop(pop, mut_rate=0.01, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Update a population by performing selection, crossover and mutation."""
    df = pop.max_fitness - pop.av_fitness
    if df <= 0:
        p = [1/pop.pop_size for i in range(pop.pop_size)]
    else:
        p = [((3-1)*(pop.state_list[i].fitness-pop.av_fitness)+df)/df/pop.pop_size for i in range(pop.pop_size)]

    state_list = []
    for i in range(int(pop.pop_size/2)):
        state_list = state_list + select_and_cross(pop, p, h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
    new_pop = population(state_list)
    
    new_pop.mutate_pop(mut_rate=mut_rate, h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)

    return new_pop
            
    
def term_states(pop):
    """"Select terminal states from a population."""
    states = [] 
    for i in range(pop.pop_size):
        if pop.state_list[i].terminal == 1:   
            P = Polytope(pop.state_list[i].points)
            states.append(P)
    return states


def remove_redundancy(term_states):
    """"Remove redundancy in a list of terminal states."""
    reduced_states = []
    for i in range(len(term_states)):
        equiv = 0
        for j in range(len(reduced_states)):
            if term_states[i].is_linearly_equivalent(reduced_states[j]):
                equiv = 1
                break
        if not equiv:
            reduced_states.append(term_states[i])
    return reduced_states
    

def evol_pop(pop, num_gen, mut_rate, monitor=True, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Evolve a population over generations and extract terminal states."""    

    term_states_list = remove_redundancy(term_states(pop))

    if monitor:
        print("Total # Terminal States    Average Fitness    Maximum Fitness")
    for i in range(num_gen):
        pop = next_pop(pop, mut_rate=mut_rate, h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
        
        term_states_list = term_states_list + term_states(pop)
        
        if monitor:
            print("    "+str(len(term_states_list))+"                        "+str(round(pop.av_fitness,2))+
                  "                  "+str(round(pop.max_fitness,2)))
    
    term_states_list = remove_redundancy(term_states_list)
    print("Total # of reduced terminal states: "+str(len(term_states_list)))

    return term_states_list


def search(num_run, num_gen, pop_size, mut_rate, dim, max_coeff, max_verts, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    """"Evolve several random populations over generations and extract terminal states."""    
    
    terminal_states = []
    print("Total # Terminal States")
    for i in range(num_run):
        initial_pop = random_pop(pop_size, polydim, maxcoeff, maxvert)
        
        terminal_states = terminal_states + evol_pop(initial_pop, num_gen, mut_rate, monitor=False, h11=h11, h12=h12, h13=h13, h22=h22, chi=chi, fav=fav)
        terminal_states = remove_redundancy(terminal_states)
        
        print("     "+str(len(terminal_states)))

    return terminal_states


