"""Contains the core of the SREA algorithm.

This file is mostly unchanged from the original RobotBrunch code.

Authors: Jordan R Abrahams, Kyle Lund, Sam Dietrich, Michael Gao
"""

from math import floor, ceil
import pulp

from stn import STN, loadSTNfromJSONfile, invcdf_norm, invcdf_uniform
import random
import glob
import json
import os
from algorithm import *


# \file SREA.py
#
#  \brief Runs the SREA algorithm on an input STN and computes the robustness
#
#  \details To run this file on an STN, use the command:
#  \code{.unparsed}
#  python SREA.py JSONFILE ALPHA
#  \endcode
#
#
#  \note You'll need pulp to be able to run this. On a linux machine run
#  `sudo pip install pulp` or `sudo easy_install -U pulp`.
#  THEN RUN `sudo pulptest`, otherwise it won't work.

# \fn addConstraint(constraint,problem)
#  \brief Adds an LP constraint to the given LP


def addConstraint(constraint, problem):
    problem += constraint
    # print 'adding constraint', constraint

# \fn setUpLP(stn)
#  \brief initializes the LP problem and the LP variables that will not change
#      with alpha
#  \returns A tuple (bounds, deltas, prob) where bounds and deltas are
#      dictionaries of LP variables, and prob is the LP problem instance


def setUpLP(stn, decouple):
    bounds = {}
    deltas = {}

    prob = pulp.LpProblem('PSTN Robust Execution LP', pulp.LpMaximize)

    # ##
    # Store Original STN edges and objective variables for easy access.
    # Not part of LP yet
    # ##

    for i in stn.verts:
        bounds[(i, '+')] = pulp.LpVariable('t_%d_hi' % i,
                                           lowBound=-stn.getEdgeWeight(i, 0),
                                           upBound=stn.getEdgeWeight(0, i))
        bounds[(i, '-')] = pulp.LpVariable('t_%d_lo' % i,
                                           lowBound=-stn.getEdgeWeight(i, 0),
                                           upBound=stn.getEdgeWeight(0, i))
        condition = bounds[(i, '+')] >= bounds[(i, '-')]
        addConstraint(condition, prob) 

    for i, j in stn.edges:
        if (i, j) in stn.contingentEdges:
            deltas[(i, j)] = pulp.LpVariable('delta_%d_%d' %
                                             (i, j), lowBound=0, upBound=None)
            deltas[(j, i)] = pulp.LpVariable('delta_%d_%d' %
                                             (j, i), lowBound=0, upBound=None)

        else:
            # ignore edges from z. these edges are implicitly handled
            # with the bounds on the LP variables
            if i != 0 and not decouple:
                addConstraint(bounds[(j, '+')] - bounds[(i, '-')]
                              <= stn.getEdgeWeight(i, j), prob)
                addConstraint(bounds[(i, '+')] - bounds[(j, '-')]
                              <= stn.getEdgeWeight(j, i), prob)

            elif i != 0 and (i, j) in stn.interagentEdges:
                addConstraint(bounds[(j, '+')] - bounds[(i, '-')]
                              <= stn.getEdgeWeight(i, j), prob)
                addConstraint(bounds[(i, '+')] - bounds[(j, '-')]
                              <= stn.getEdgeWeight(j, i), prob)
    return (bounds, deltas, prob)


##
# \fn srea(inputstn,debug=False,debugLP=False,lb=0.0,ub=0.999)
# \brief Runs the SREA algorithm on an input STN
#
# @param inputstn The STN that we are running SREA on
# @param invCDF_map A dictionary of dictionaries generated from a
#      distgenlib.invcdf_map call. See documentation there for more info.
# @param debug Print optional status messages about alpha levels
# @param debugLP Print optional status messages about each run of the LP
# @param lb The starting lower bound on alpha for the binary search
# @param ub The starting upper bound on alpha for the binary search
#
# @returns a tuple (alpha, outputstn) if there is a solution, or None if there
#     is no solution
def srea(inputstn,
         debug=False,
         debugLP=False,
         returnAlpha=True,
         decouple=False,
         lb=0.0,
         ub=0.999):
    inputstn = inputstn.copy()
    # dictionary of alphas for binary search
    alphas = {i: i / 1000.0 for i in range(1001)}

    # bounds for binary search
    lower = ceil(lb * 1000) - 1
    upper = floor(ub * 1000) + 1

    result = None

    # set up LP
    if not decouple:
        # TODO: Change to faster algorithm?
        inputstn.floyd_warshall()
        print(inputstn)

    bounds, deltas, probBase = setUpLP(inputstn, decouple)

    # First run binary search on alpha
    while upper - lower > 1:
        alpha = alphas[(upper + lower) // 2]
        if debug:
            print('trying alpha = {}'.format(alpha))

        # run the LP
        probContainer = (bounds, deltas, probBase.copy())
        LPbounds = srea_LP(inputstn.copy(),
                           alpha,
                           decouple,
                           debug=debugLP,
                           probContainer=probContainer)

        # LP was feasible, try lower alpha
        if LPbounds is not None:
            upper = (upper + lower) // 2
            result = (alpha, LPbounds)
        # LP was infeasable, try higher alpha
        else:
            lower = (upper + lower) // 2

        # finished our search, load the smallest alpha decoupling
        if upper - lower <= 1:
            if result is not None:
                alpha, LPbounds = result
                if debug:
                    print(
                        'modifying STN with lowest good alpha, {}'.format(alpha))
                for i, sign in LPbounds:
                    if sign == '+':
                        inputstn.updateEdge(
                            0, i, ceil(bounds[(i, '+')].varValue))
                    else:
                        inputstn.updateEdge(
                            i, 0, ceil(-bounds[(i, '-')].varValue))

                if returnAlpha:
                    print(inputstn) #TODO- remove this
                    return alpha, inputstn
                else:
                    return inputstn
    # skip the rest if there was no decoupling at all
    if result is None:
        if debug:
            print('could not produce feasible LP.')
        return 1,None

    # Fail here
    assert(False)


# \fn srea_LP(inputstn,alpha,debug=False,probContainer=None)
#  \brief Runs the robust execution LP on the input STN at the given alpha
#  level
#
#  @param inputSTN The STN used as an input to the LP
#  @param invCDF_map A dictionary of dictionaries generated from a
#       distgenlib.invcdf_map call. See documentation there for more info.
#  @param alpha The risk level (between 0 and 1) that we are using for the LP
#  @param decouple originally was meant to indicate if we wanted decoupling or not
#   but then we discovered that this already decouples the STN
#  @param debug Print optional status messages
#  @param probContainer Optional tuple of LP variables and the LP problem
#       instance, returned from setUpLP
#
#  \returns A dictionary of the LP_variables for the bounds on timepoints.
def srea_LP(inputstn,
            alpha,
            decouple,
            debug=False,
            probContainer=None
            ):

    # Check some types to make sure everything is the correct type
    if not isinstance(inputstn, STN):
        raise TypeError("inputstn is not of type STN")

    alpha = round(float(alpha), 3)

    if probContainer is None:
        if debug:
            print('No saved LP variables, generating all LP variables from current STN')
        bounds, deltas, prob = setUpLP(inputstn, decouple)
    else:
        bounds, deltas, prob = probContainer

    for (i, j), edge in list(inputstn.contingentEdges.items()):
        if edge.dtype() == "gaussian":
            p_ij = invcdf_norm(1.0 - alpha * 0.5, edge.mu, edge.sigma)
            p_ji = -invcdf_norm(alpha * 0.5, edge.mu, edge.sigma)
            limit_ij = invcdf_norm(0.997, edge.mu, edge.sigma)
            limit_ji = -invcdf_norm(0.003, edge.mu, edge.sigma)
        elif edge.dtype() == "uniform":
            p_ij = invcdf_uniform(1.0 - alpha * 0.5, edge.dist_lb,
                                  edge.dist_ub)
            p_ji = -invcdf_uniform(alpha * 0.5, edge.dist_lb, edge.dist_ub)
            limit_ij = invcdf_uniform(0.0, edge.dist_lb, edge.dist_ub)
            limit_ji = -invcdf_uniform(1.0, edge.dist_lb, edge.dist_ub)

        deltas[(i, j)].upBound = limit_ij - p_ij
        deltas[(j, i)].upBound = limit_ji - p_ji

        cons1 = bounds[(j, "+")] - bounds[(i, "+")] == p_ij + deltas[(i, j)]
        cons2 = bounds[(j, "-")] - bounds[(i, "-")] == -p_ji - deltas[(j, i)]
        # Lund et al. LP (3)
        addConstraint(cons1, prob)
        # Lund et al. LP (4)
        addConstraint(cons2, prob)
    # ##
    # Generate the objective function.
    #   Our objective function is SUM delta_ij
    # ##
    deltaSum = sum([deltas[(i, j)] for i, j in deltas])
    prob += deltaSum, 'Maximize time added back to \
        constraints while decoupling'

    if debug:
        prob.writeLP('STN.lp')
        pulp.LpSolverDefault.msg = 10

    # for RPi only
    # mySolver = solvers.GLPK()
    # prob.solve(solver=mySolver)
    # This try except exists because there was a weird bug with Pulp where
    # it was throwing an error instead of just saying that the lp was not
    # resolvable.  I don't know much about the inner workings of pulp and
    # stack overflow suggested I put in this fix so I did.
    # https://stackoverflow.com/questions/27406858/pulp-solver-error
    # try:
    prob.solve()
    # except Exception:
    # return None

    status = pulp.LpStatus[prob.status]
    if debug:
        print('Status:', status)
        # Each of the variables is printed with it's resolved optimum value
        for v in prob.variables():
            print(v.name, '=', v.varValue)
    if status != 'Optimal':
        return None
    return bounds





if __name__ == "__main__":
    comparison = []
    improvement = 0
    tied = 0
    failed = 0
    count = 0
    directory = "dataset/uncontrollable_full"

    data_list = glob.glob(os.path.join(directory, '*.json'))
    # data_list = ['dataset/uncontrollable_full/uncontrollable6.json']
    # data_list = ['dataset/dreamdata/STN_a4_i4_s5_t10000/original_0.json']
    data_list = ['dataset/dreamdata/STN_a2_i4_s5_t5000/original_6.json']



    # testing dream data ##

    directory = 'dataset/dreamdata/'
    folders = os.listdir(directory)
    data_list = []
    result = []
    for folder in folders:
        data = glob.glob(os.path.join(directory, folder, '*.json'))
        data_list += data
    
    for data in data_list:
        print("simulating", data)
        stn = loadSTNfromJSONfile(data)
        inputstn = stn.copy()
        inputstn.floyd_warshall()
        alpha, outputstn = srea(inputstn)
        if outputstn != None:
            a,b,c,d = DC_Checker(outputstn)   
            if not a:
                result += [data, alpha] 
    print(result)
    print(len(result))


        # alpha, newstn = srea(stn, debugLP = True)