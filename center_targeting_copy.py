from util import PriorityQueue, sys
import pulp
import math

from stn import STN, loadSTNfromJSONfile

from dispatch import generate_realization

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

    prob = pulp.LpProblem('Expected Value LP', pulp.LpMinimize)

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
        bounds[(i,'+')] = pulp.LpVariable('t_%i_hi'%i, lowBound=0,
                                            upBound=STN.getEdgeWeight(0,i))

        lowbound = 0 if STN.getEdgeWeight(i,0) == float('inf') else\
                            -STN.getEdgeWeight(i,0)
        bounds[(i,'-')] = pulp.LpVariable('t_%i_lo'%i, lowBound=lowbound,
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
            offsets[(i,j)] = pulp.LpVariable('eps_%i'%j, lowBound=0)
            # deltas[(i,j)] = pulp.LpVariable('delta_%i'%j, lowBound= -offsets[(i,j)],
            #                                                 upBound=offsets[(i,j)])
            deltas[(i,j)] = pulp.LpVariable('delta_%i'%j)
            cur_edge = STN.getEdge(i,j)
            if cur_edge.dtype() == "gaussian":
                center = cur_edge.mu
            else:
                center = (cur_edge.Cij - cur_edge.Cji)/2
            addConstraint(deltas[(i,j)] <= offsets[(i,j)], prob)
            addConstraint(-offsets[(i,j)] <= deltas[(i,j)], prob)
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
            lowbound = MAX_FLOAT if STN.getEdgeWeight(j,i) == float('inf') \
                                            else STN.getEdgeWeight(j,i)

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
    status = pulp.LpStatus[prob.status]
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

    prob = pulp.LpProblem('Robust Scheduling LP', pulp.LpMaximize)
    radius = pulp.LpVariable('rad', lowBound=0, upBound=None)

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
        timepoints[i] = pulp.LpVariable('t_%i'%i, lowBound=0, upBound=STN.getEdgeWeight(0,i))



    for i,j in STN.edges:
        if (i,j) in STN.contingentEdges:
            addConstraint(timepoints[j] - timepoints[i] == STN.getEdgeWeight(i,j), prob)
        else:
            upbound = MAX_FLOAT if STN.getEdgeWeight(i,j) == float('inf') \
                                            else STN.getEdgeWeight(i,j)
            lowbound = MAX_FLOAT if STN.getEdgeWeight(j,i) == float('inf') \
                                            else STN.getEdgeWeight(j,i)
            if i==0 or j==0:
                addConstraint(timepoints[j] - timepoints[i] <= upbound - radius, prob)
                addConstraint(timepoints[i] - timepoints[j] <= lowbound - radius, prob)

            else:
                addConstraint(timepoints[j] - timepoints[i] <= upbound - math.sqrt(2) * radius, prob)
                addConstraint(timepoints[i] - timepoints[j] <= lowbound - math.sqrt(2) * radius, prob)


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
    status = pulp.LpStatus[prob.status]
    if debug:
        print("Status: ", status)

        for v in prob.variables():
            print(v.name, '=', v.varValue)

    if status != 'Optimal':
        print("The solution for LP is not optimal")
        return status, None, None

    return status, timepoints, radius


def simulateCenter(stn, size = 1):
    stncopy = stn.copy()
    failure = 0
    for _ in range(size):
        realization  = generate_realization(stn)
        P = PriorityQueue()
        executed = []
        schedule = {}
        P.push(0,0)
        while not P.isEmpty():
            centerResult = centeringLP(stncopy, False)
            if centerResult[0] == "Optimal":
                stncopy = centerUpdate(stncopy, centerResult[2])
            else:
                failure += 1
                break
            scheduleResult = scheduleLP(stncopy, debug=False)
            if scheduleResult[0] == "Optimal":
                timepoints = scheduleResult[1]
                for i in range(len(timepoints)):
                    vert = int(timepoints[i].name[2:])
                    if vert != 0 and vert not in executed:
                        P.addOrDecKey(vert, timepoints[i].varValue)
            else:
                failure += 1
                break
                
            weight, currVert = P.pop()
            if currVert in stncopy.uncontrollables:
                parentVert = stncopy.parent[currVert]
                duration = realization[(parentVert, currVert)]
                stncopy.modifyEdge(parentVert,currVert,duration)
                stncopy.addEdge(parentVert,currVert,duration,duration)
                executed += [currVert]
                schedule[currVert] = weight

            else:
                stncopy.updateEdge(0, currVert, weight)
                stncopy.updateEdge(currVert, 0, -weight)
                executed += [currVert]
                schedule[currVert] = weight
                while currVert not in stncopy.uncontrollables and not P.isEmpty():
                    weight, currVert = P.pop()
                    stncopy.updateEdge(0, currVert, weight)
                    stncopy.updateEdge(currVert, 0, -weight)
                if currVert in stncopy.uncontrollables:
                    parentVert = stncopy.parent[currVert]
                    duration = realization[(parentVert, currVert)]
                    stncopy.modifyEdge(parentVert,currVert,duration)
                    stncopy.addEdge(parentVert,currVert,duration,duration)
                    executed += [currVert]
                    schedule[currVert] = weight
                else:
                    break
                            

            
            
                


            break
        

def centerUpdate(stn, offsets):
    stncopy = stn.copy()
    for (i,j) in list(offsets.keys()):
        currEdge = stncopy.getEdge(i,j)
        if currEdge.dtype() == "gaussian":
            weight = currEdge.mu + offsets[(i,j)].varValue
            stncopy.updateEdge(i,j,weight)
            stncopy.updateEdge(j,i,-weight)
        else:
            weight = (currEdge.dist_lb + currEdge.dist_ub) / 2 + offsets[(i,j)].varValue
            stncopy.updateEdge(i,j,weight)
            stncopy.updateEdge(j,i,-weight)
    return stncopy


if __name__ == "__main__":
    stn = loadSTNfromJSONfile('dataset/dreamdata/STN_a2_i4_s1_t4000/original_9.json')
    stncopy = stn.copy()
    simulateCenter(stn)