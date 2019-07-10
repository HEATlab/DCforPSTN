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

    # list of untightened edges
    tedges = stncopy.contingentEdges
    result = None

    rounds = 0

    while len(tedges) > 0:
        rounds += 1

        while upper - lower > 1:
            alpha = alphas[(upper + lower) // 2]

            if debug:
                print('trying alpha = {}'.format(alpha))
            stncopy = alphaUpdate(stncopy, tedges, alpha)
            
            if alpha == .999:
                return stncopy

            dc, conflicts, bounds, weight = DC_Checker(stncopy)
            if dc:
                upper = (upper + lower) // 2
                result = (alpha, stncopy)
            else:
                lower = (upper + lower) // 2
            if debug:
                print("lower", lower/1000, "upper", upper/1000)
            # finished our search, load the smallest alpha decoupling

            if upper - lower <= 1:
                if result is not None:
                    if debug:
                        print('finished binary search, with lower', lower/1000, 'alpha', result[0], 'upper', upper/1000, 'rounds', rounds)
                    stncopy = alphaUpdate(stncopy, tedges, result[0])
                    loweststn = alphaUpdate(stncopy, tedges, result[0]-.001)
                    dc, conflicts, bounds, weight = DC_Checker(loweststn)
                    if dc:
                        return loweststn
                    else:
                        # find the tightest contingent edges
                        tightest = bounds['contingent']
                        print(list(tightest.items()), rounds)
                        for i,j in list(tightest.keys()):
                            edge, bound = tightest[i,j]
                            tedges.pop((edge.i, edge.j))
                else:   
                    if debug:
                        print('could not produce dynamically controllable STNU.')
                    return None
        lower = ceil(lb * 1000) - 1
        upper = result[0]*1000
    return stncopy

def alphaUpdate(inputstn, tedges, alpha):
    stncopy = inputstn.copy()
    # update the edges based on the list of target edges
    if alpha <= 0:
        return stncopy
    else:
        for (i, j), edge in list(tedges.items()):
            if edge.type == "stcu":
                p_ij = invcdf_uniform(1.0 - alpha * 0.5, -edge.Cji, edge.Cij)
                p_ji = -invcdf_uniform(alpha * 0.5, -edge.Cji, edge.Cij)
                stncopy.modifyEdge(i, j, p_ij)
                stncopy.modifyEdge(j, i, p_ji)
            else:
                assert edge.dtype != None
                dist = edge.distribution

                if edge.dtype() == "gaussian":
                    p_ij = invcdf_norm(1.0 - alpha * 0.5, edge.mu, edge.sigma)
                    p_ji = -invcdf_norm(alpha * 0.5, edge.mu, edge.sigma)
                    limit_ij = invcdf_norm(0.997, edge.mu, edge.sigma)
                    limit_ji = -invcdf_norm(0.003, edge.mu, edge.sigma)
                    stncopy.modifyEdge(i, j, p_ij)
                    stncopy.modifyEdge(j, i, p_ji)

                elif edge.dtype() == "uniform":
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
    # directory = "dataset/uncontrollable"
    # data_list = glob.glob(os.path.join(directory, '*.json'))
    # comparison = []
    # for data in data_list:
    #     print(data)
    #     stn = loadSTNfromJSONfile(data)
    #     newstn = maxgain(stn, debug = True)
    #     result = simulate_maxgain(stn, newstn, size = 500)
    #     oldresult = simulation(stn,500)
    #     comparison += [(result, oldresult)]
    # sumNew = 0
    # sumOld = 0
    # improvementCount = 0
    # equalCount = 0
    # for i in range(len(comparison)):
    #     sumNew += comparison[i][0]
    #     sumOld += comparison[i][1]
    #     if comparison[i][0] > comparison[i][1]:
    #         improvementCount += 1
    #     elif comparison[i][0] == comparison[i][1]:
    #         equalCount += 1
    # print(sumNew, sumOld, improvementCount, equalCount, len(comparison))
    # print(comparison)
    stn = loadSTNfromJSONfile('dataset/dreamdata/original_6.json')
    newstn = maxgain(stn, debug = True)
    result = simulate_maxgain(stn, newstn, size = 500)
    print(result)
    




