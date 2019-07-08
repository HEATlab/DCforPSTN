from stn import STN, loadSTNfromJSONfile, invcdf_norm, invcdf_uniform
from util import STNtoDCSTN, PriorityQueue
from dc_stn import DC_STN
import empirical
import random
import glob
import json
import os
from algorithm import *
from math import floor, ceil
from dispatch import *


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
    result = None
    while len(tedges) > 0:
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
            else:
                lower = (upper + lower) // 2
            # finished our search, load the smallest alpha decoupling
            if upper - lower <= 1:
                if result is not None:
                    loweststn = alphaUpdate(stncopy, tedges, lower/1000)
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
                            tedges.pop((edge.i, edge.j))
                else:   
                    if debug:
                        print('could not produce dynamically controllable STNU.')
                    return None
        lower = ceil(lb * 1000) - 1
    return stncopy

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
            p_ij = invcdf_uniform(1.0 - alpha * 0.5, -edge.Cji, edge.Cij)

            p_ji = -invcdf_uniform(alpha * 0.5, -edge.Cji, edge.Cij)
            stncopy.modifyEdge(i, j, p_ij)
            stncopy.modifyEdge(j, i, p_ji)
        return stncopy

def simulate_maxgain(network, shrinked_network, size=200, verbose=False, gauss=False):
    # Collect useful data from the original network
    contingent_pairs = network.contingentEdges.keys()
    contingents = {src: sink for (src, sink) in contingent_pairs}
    uncontrollables = set(contingents.values())

    total_victories = 0
    dc_network = STNtoDCSTN(shrinked_network)
    dc_network.addVertex(ZERO_ID)

    controllability = dc_network.is_DC()
    if verbose:
        print("Finished checking DC...")

    # Detect if the network has an inconsistency in a fixed edge
    verts = dc_network.verts.keys()
    for vert in verts:
        if (vert, vert) in dc_network.edges:
            if verbose:
                print("Checking", vert)
            edge = dc_network.edges[vert, vert][0]
            if edge.weight < 0:
                dc_network.edges[(vert, vert)].remove(edge)
                dc_network.verts[vert].outgoing_normal.remove(edge)
                dc_network.verts[vert].incoming_normal.remove(edge)
                del dc_network.normal_edges[(vert, vert)]

    # Run the simulation
    for j in range(size):
        realization = generate_realization(network, gauss)
        copy = dc_network.copy()
        result = dispatch(network, copy, realization, contingents,
                          uncontrollables, verbose)
        if verbose:
            print("Completed a simulation.")
        if result:
            total_victories += 1

    goodie = float(total_victories / size)
    if verbose:
        print(f"Worked {100*goodie}% of the time.")

    return goodie


if __name__ == "__main__":
    directory = "dataset/uncontrollable"
    data_list = glob.glob(os.path.join(directory, '*.json'))
    comparison = []
    for data in data_list:
        print(data)
        stn = loadSTNfromJSONfile(data)
        newstn = maxgain(stn, debug = False)
        result = simulate_maxgain(stn, newstn, size = 10)
        oldresult = simulation(stn,10)
        print(result, oldresult)
        comparison += [(result, oldresult)]
    




