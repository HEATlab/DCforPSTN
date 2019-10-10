from util import *
from pulp import *
import math
import glob
import json
import os
from stn import STN, loadSTNfromJSONfile
##
# \brief A global variable that stores the max float that will be used to deal
#        with infinite edges.
MAX_FLOAT = sys.float_info.max


## \file LP.py
#  \brief Convert an input STNU to LP form and computes the degree of
#         controllability


##
#  \fn addConstraint(constraint,problem)
#  \brief Adds an LP constraint to the given LP
#
#  @param constraint A constraint that need to be added to the LP problem
#  @param problem    An input LP problem
#
#  @post LP problem with new constraint added
def addConstraint(constraint,problem):
    problem += constraint


##
# \fn setUpCentering(STN)
# \brief Initializes the LP problem and the LP variables
#
# @param STN            An input STNU
# 
# @return   A tuple (bounds, deltas, prob) where bounds and deltas are
#           dictionaries of LP variables, and prob is the LP problem instance
def setUpCentering(STN, proportion=False):
    bounds = {}
    deltas = {}
    offsets = {}

    prob = LpProblem('Expected Value LP', LpMinimize)

    STN = STN.minimal()

    # ##
    # NOTE: Our LP requires each event to occur within a finite interval. If
    #       the input LP does not have finite interval specified for all
    #       events, we want to set the setMakespan to MAX_FLOAT (infinity)
    #       so the LP works
    #
    #       We do not want to run minimal network first because we are going to
    #       modify the contingent edges in LP, while some constraints in
    #       minimal network are obtained through contingent edges
    #
    #       There might be better way to deal with this problem.
    # ##
    for i in STN.verts:
        if STN.getEdgeWeight(0,i) == float('inf'):
            STN.setMakespan(MAX_FLOAT)
            break

    # ##
    # Store Original STN edges and objective variables for easy access.
    # Not part of LP yet
    # ##
    for i in STN.verts:
        bounds[(i,'+')] = LpVariable('t_%i_hi'%i, lowBound=0,
                                            upBound=STN.getEdgeWeight(0,i))

        lowbound = 0 if STN.getEdgeWeight(i,0) == float('inf') else\
                            -STN.getEdgeWeight(i,0)
        bounds[(i,'-')] = LpVariable('t_%i_lo'%i, lowBound=lowbound,
                                            upBound=None)

        addConstraint( bounds[(i,'-')] <= bounds[(i,'+')], prob)

        # TODO: FIND OUT IF NEEDED
        # if i == 0:
        #     addConstraint(bounds[(i,'-')] == 0, prob)
        #     addConstraint(bounds[(i,'+')] == 0, prob)

        if i in STN.uncontrollables:
            addConstraint( bounds[(i,'-')] == bounds[(i,'+')], prob)


    for i,j in STN.edges:
        if (i,j) in STN.contingentEdges:
            offsets[(i,j)] = LpVariable('eps_%i'%j, lowBound=0,
                                                            upBound=None)
            deltas[(i,j)] = LpVariable('delta_%i'%j, lowBound= -offsets[(i,j)],
                                                            upBound=offsets[(i,j)])

            cur_edge = STN.getEdge(i,j)
            if cur_edge.dtype() == "gaussian":
                center = cur_edge.mu
                print('center', center)
            else:
                center = (cur_edge.Cij - cur_edge.Cji)/2


            addConstraint(bounds[(j,'+')] - bounds[(i,'+')] ==
                    center + deltas[(i,j)], prob)
            addConstraint(bounds[(j,'-')] - bounds[(i,'-')] ==
                    center + deltas[(i,j)], prob)

        else:
            # NOTE: We need to handle the infinite weight edges. Otherwise
            #       the LP would be infeasible
            upbound = MAX_FLOAT if STN.getEdgeWeight(i,j) == float('inf') \
                                            else STN.getEdgeWeight(i,j)


            #TODO - shouldnt it be negative?
            lowbound = -MAX_FLOAT if STN.getEdgeWeight(j,i) == -float('inf') \
                                            else -STN.getEdgeWeight(j,i)               
            print('lowbound',lowbound)
            addConstraint(bounds[(j,'+')]-bounds[(i,'-')] <= upbound, prob)
            addConstraint(bounds[(i,'+')]-bounds[(j,'-')] <= lowbound, prob)

    return (bounds, offsets, prob)

def centeringLP(STN, debug=False):
    bounds, offsets, prob = setUpCentering(STN)
    
    Obj = sum([offsets[(i,j)] for i,j in offsets])


    prob += Obj, "Minimize the sum of distances from the expected value for the STN"

    try:
        prob.solve()
    except Exception:
        print("The model is invalid.")
        return 'Invalid', None, None

    # Report status message
    status = LpStatus[prob.status]
    if debug:
        print("Status: ", status)

        for v in prob.variables():
            print(v.name, '=', v.varValue)

    if status != 'Optimal':
        print("The solution for LP is not optimal")
        return status, None, None

    return status, bounds, offsets


def setUpScheduling(STN, debug=False):
    timepoints = {}

    prob = LpProblem('Robust Scheduling LP', LpMaximize)
    radius = LpVariable('rad', lowBound=0, upBound=None)

    # ##
    # NOTE: Our LP requires each event to occur within a finite interval. If
    #       the input LP does not have finite interval specified for all
    #       events, we want to set the setMakespan to MAX_FLOAT (infinity)
    #       so the LP works
    #
    #       We do not want to run minimal network first because we are going to
    #       modify the contingent edges in LP, while some constraints in
    #       minimal network are obtained through contingent edges
    #
    #       There might be better way to deal with this problem.
    # ##
    for i in STN.verts:
        if STN.getEdgeWeight(0,i) == float('inf'):
            STN.setMakespan(MAX_FLOAT)
            break

    # ##
    # Store Original STN edges and objective variables for easy access.
    # Not part of LP yet
    # ##
    for i in STN.verts:
        timepoints[i] = LpVariable('t_%i'%i, lowBound=0, upBound=None)


    for i,j in STN.edges:
        if (i,j) in STN.contingentEdges:
            if i==0 or j==0:
                addConstraint(timepoints[j] - timepoints[i] <= \
                    STN.getEdgeWeight(i,j) + radius, prob)

            else:
                addConstraint(timepoints[j] - timepoints[i] <= \
                    STN.getEdgeWeight(i,j) + math.sqrt(2)*radius, prob)

        else:
            # NOTE: We need to handle the infinite weight edges. Otherwise
            #       the LP would be infeasible
            upbound = MAX_FLOAT if STN.getEdgeWeight(i,j) == float('inf') \
                                            else STN.getEdgeWeight(i,j)
            lowbound = MAX_FLOAT if STN.getEdgeWeight(j,i) == float('inf') \
                                            else STN.getEdgeWeight(j,i)

            addConstraint(timepoints[j] - timepoints[i] <= upbound, prob)
            addConstraint(timepoints[i] - timepoints[j] <= lowbound, prob)

    return (timepoints, radius, prob)

def scheduleLP(STN, debug=False):
    timepoints, radius, prob = setUpScheduling(STN)
    
    Obj = radius

    prob += Obj, "Maximize the robustness of the schedule"

    try:
        prob.solve()
    except Exception:
        print("The model is invalid.")
        return 'Invalid', None, None

    # Report status message
    status = LpStatus[prob.status]
    if debug:
        print("Status: ", status)

        for v in prob.variables():
            print(v.name, '=', v.varValue)

    if status != 'Optimal':
        print("The solution for LP is not optimal")
        return status, None, None

    return status, timepoints, radius


if __name__ == "__main__":
    data = 'dataset/mr_x.json'
    stn = loadSTNfromJSONfile(data)
    print(stn)
    hello = scheduleLP(stn)
    centered = centeringLP(stn, True)

