from stn import STN, loadSTNfromJSONfile, invcdf_norm, invcdf_uniform
from util import STNtoDCSTN, PriorityQueue
from dc_stn import DC_STN
import empirical
import random
import json
from algorithm import *
from math import floor, ceil


# def sample_maxgain(, size):
#     for i in range(size)





def maxgain(inputstn,
         debug=False,
         returnAlpha=True,
         lb=0.0,
         ub=0.999):
    stncopy = inputstn.copy()
    # dictionary of alphas for binary search

    alphas = {i: i / 1000.0 for i in range(1001)}

    # bounds for binary search
    lower = ceil(lb * 1000) - 1
    upper = floor(ub * 1000) + 1

    result = None

    # First run binary search on alpha
    # list of untightened edges
    tedges = stncopy.contingentEdges
    print(tedges)
    result = None
    while len(tedges) != 0:
        #lowest alpha
        workingAlpha = []
        
        while upper - lower > 1:
            alpha = alphas[(upper + lower) // 2]

            if debug:
                print('trying alpha = {}'.format(alpha))
            stncopy = alphaUpdate(stncopy, tedges, alpha)
            if alpha == .999:
                return stncopy
            #
            # 
            # bounds is a dictionary of contigent edges you can change
            # contigent
            #weight sum of the cycle 
            dc, conflicts, bounds, weight = DC_Checker(stncopy)
            if dc:
                upper = (upper + lower) // 2
                result = (alpha, stncopy)
                print("xixi")
            else:
                lower = (upper + lower) // 2
            # finished our search, load the smallest alpha decoupling
            if upper - lower <= 1:
                if result is not None:
                    loweststn = alphaUpdate(stncopy, tedges, lower)
                    dc, conflicts, bounds, weight = DC_Checker(loweststn)
                    if dc:
                        if debug:
                            print('The system is dynamically controllable')
                        return loweststn
                    else:
                        # find the tightest contingent edges
                        tightest = bounds['contingent']
                        for i,j in list(tightest.keys()):
                            edge, bound = tightest[i,j]
                            tedges.remove(edge)
                        lower = ceil(lb * 1000) - 1
    
                        
    # skip the rest if there was no decoupling at all
    if result is None:
        if debug:
            print('could not produce dynamically controllable STNU.')
        return None

    # Fail here
    assert(False)

def alphaUpdate(inputstn, tedges, alpha):
    stncopy = inputstn.copy()
    # update the edges based on the list of target edges
    if alpha <= 0:
        return stncopy
    else:
        for (i, j), edge in list(tedges.items()):
            # if edge.distribution() == "gaussian":
            #     p_ij = invcdf_norm(1.0 - alpha * 0.5, edge.mu, edge.sigma)
            #     p_ji = -invcdf_norm(alpha * 0.5, edge.mu, edge.sigma)

            # elif edge.distribution() == "uniform" or None:
            #     p_ij = invcdf_uniform(1.0 - alpha * 0.5, edge.dist_lb,
            #                         edge.dist_ub)
            #     p_ji = -invcdf_uniform(alpha * 0.5, edge.dist_lb, edge.dist_ub)
            p_ij = invcdf_uniform(1.0 - alpha * 0.5, -edge.Cij, edge.Cij)
            p_ji = -invcdf_uniform(alpha * 0.5, -edge.Cij, edge.Cij)
            stncopy.updateEdge(i, j, p_ij)
            print("pij",p_ij)
            stncopy.updateEdge(j, i, p_ji)
            print("pji",p_ji)
        return stncopy



stn = loadSTNfromJSONfile("small_examples/dynamic1.json")  
newstn = maxgain(stn, debug = True)
print(stn)
print(newstn)
