# author: elli heyes (elli.heyes@city.ac.uk)
# date last edited: 15/12/2023
import numpy as np
from cytools import Polytope

def IP(poly):
    pts = [list(item) for item in poly.points()]
    origin = [0 for i in range(poly.dim())]
    return len(poly.interior_points()) == 1 and origin in pts

def lattice_dists(poly):
    return poly.inequalities()[:,-1]

def fitness(points, h11=None, h12=None, h13=None, h22=None, chi=None, fav=False):
    w1 = w2 = w3 = w4 = w5 = w6 = w7 = w8 = 1

    # delete duplicate points 
    points_red = []
    for point in points:
        if not point in points_red:
            points_red.append(point)
    
    poly = Polytope(points_red)
    
    term1 = w1 * (IP(poly) - 1)
    term2 = - w2 * np.sum(np.abs(d:=lattice_dists(poly) - 1)) / len(d)
    
    term3 = term4 = term5 = term6 = term7 = term8 = 0
    if h11 != None:
        if poly.is_reflexive():
            term3 = - w3 * abs(h11 - poly.h11(lattice="N"))
        else:
            term3 = -100
    if h12 != None:
        if poly.is_reflexive():
            term4 = - w4 * abs(h12 - poly.h12(lattice="N"))
        else:
            term4 = -100
    if h13 != None:
        if poly.is_reflexive():
            term5 = - w5 * abs(h13 - poly.h13(lattice="N"))
        else:
            term5 = -100
    if h22 != None:
        if poly.is_reflexive():
            term6 = - w6 * abs(h22 - poly.h22(lattice="N"))
        else:
            term6 = -100
    if chi != None:
        if poly.is_reflexive():
            term7 = - w7 * abs(chi - poly.chi(lattice="N"))
        else:
            term7 = -100
    if fav == True:
        if poly.is_reflexive():
            if not poly.is_favorable(lattice="N"):
                term8 = - w8 * 1
        else:
            term8 = -100
    
    return term1 + term2 + term3 + term4 + term5 + term6 + term7

