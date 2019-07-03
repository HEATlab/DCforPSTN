from algorithm import *
from pulp import *

## \file relax.py
#  \brief relaxation algorithm for dynamic controllability


##
#  \fn addConstraint(constraint,problem)
#  \brief Adds an LP constraint to the given LP
#
#  @param constraint A constraint that need to be added to the LP problem
#  @param problem    An input LP problem
#
#  @post LP problem with new constraint added
def addConstraint(constraint, problem):
    problem += constraint


##
# \fn getShrinked(STN, bounds, epsilons)
# \brief compute the shrinked contingent intervals and find contingent intervals
#        that are actually shrinked
#
# @param STN            An input STN
# @param bounds         A dictionary of bounds we can relax to resolve conflict
# @param epsilons       A dictionary of variables returned by LP to resolve
#                       the conflict
#
# @return dictionaries of original and shrinked contingent intervals and a list
#         of changed contingent constarints
def getShrinked(STN, bounds, epsilons):
    contingent = bounds['contingent']
    original = {}
    shrinked = {}
    changed = []

    for (i, j) in list(STN.contingentEdges.keys()):
        edge = STN.contingentEdges[(i, j)]
        orig = (-edge.Cji, edge.Cij)
        original[(i, j)] = orig

        if (i, j) not in contingent or epsilons[j].varValue == 0:
            shrinked[(i, j)] = orig
        else:
            eps = epsilons[j].varValue
            _, bound = contingent[(i, j)]
            low, high = orig

            if bound == 'UPPER':
                high -= eps
            else:
                low += eps

            shrinked[(i, j)] = (low, high)
            changed.append((i, j))

    return original, shrinked, changed


##
# \fn relaxNLP(bounds, weight, debug=False)
# \brief run LP to compute amount of uncertainty need to be removed to
#        resolve the negative cycle
#
# @param bounds       A dictionary of bounds we can relax to resolve conflict
# @param weight       Weight of the semi-reducible negative cycle
# @param debug        Flag indicating wehther we want to report information
#
# @return A dictionary Lp variable epsilons (amount of uncertain need to be
#         removed from each contingent interval)
def relaxNLP(bounds, weight, debug=False):
    contingent = bounds['contingent']
    epsilons = {}

    prob = LpProblem('Relaxation LP', LpMinimize)

    eps = []
    for i, j in list(contingent.keys()):
        edge, bound = contingent[(i, j)]
        length = edge.Cij + edge.Cji

        epsilons[j] = LpVariable('eps_%i' % j, lowBound=0, upBound=length)
        eps.append(1.0 / length * epsilons[j])

    s = sum([epsilons[j] for j in epsilons])
    addConstraint(s >= -weight, prob)

    Obj = sum(eps)
    prob += Obj, "Minimize the Uncertainty Removed"

    # write LP into file for debugging (optional)
    if debug:
        prob.writeLP('relax.lp')
        LpSolverDefault.msg = 10

    try:
        prob.solve()
    except Exception:
        print("The model is invalid.")
        return 'Invalid', None

    # Report status message
    status = LpStatus[prob.status]
    if debug:
        print("Status: ", status)

        for v in prob.variables():
            print(v.name, '=', v.varValue)

    if status != 'Optimal':
        print('hi', status)
        print("The solution for LP is not optimal")
        return status, None

    return status, epsilons


##
# \fn relaxDeltaLP(bounds, weight, debug=False)
# \brief run delta LP to compute amount of uncertainty need to be removed to
#        resolve the negative cycle
#
# @param bounds       A dictionary of bounds we can relax to resolve conflict
# @param weight       Weight of the semi-reducible negative cycle
# @param debug        Flag indicating wehther we want to report information
#
# @return A dictionary Lp variable epsilons (amount of uncertain need to be
#         removed from each contingent interval)
def relaxDeltaLP(bounds, weight, debug=False):
    contingent = bounds['contingent']
    epsilons = {}

    prob = LpProblem('Relaxation Delta LP', LpMinimize)
    delta = LpVariable('delta', lowBound=0, upBound=1)

    for i, j in list(contingent.keys()):
        edge, bound = contingent[(i, j)]
        length = edge.Cij + edge.Cji

        epsilons[j] = LpVariable('eps_%i' % j, lowBound=0, upBound=length)
        addConstraint(epsilons[j] == delta * length, prob)

    s = sum([epsilons[j] for j in epsilons])
    addConstraint(s >= -weight, prob)

    Obj = delta
    prob += Obj, "Minimize the Proportion of Uncertainty Removed"

    # write LP into file for debugging (optional)
    if debug:
        prob.writeLP('relax_delta.lp')
        LpSolverDefault.msg = 10

    try:
        prob.solve()
    except Exception:
        print("The model is invalid.")
        return 'Invalid', None

    # Report status message
    status = LpStatus[prob.status]
    if debug:
        print("Status: ", status)

        for v in prob.variables():
            print(v.name, '=', v.varValue)

    if status != 'Optimal':
        print("The solution for LP is not optimal")
        return status, None

    return status, epsilons


##
# \fn optimalRelax(bounds, weight)
# \brief optimal solution for compute relax strategy
#
# @param bounds       A dictionary of bounds we can relax to resolve conflict
# @param weight       Weight of the semi-reducible negative cycle
#
# @return A dictionary of epsilons (amount of uncertain need to be
#         removed from each contingent interval)
def optimalRelax(bounds, weight):
    contingent = [bounds['contingent'][x][0] for x in \
                                list(bounds['contingent'].keys())]
    contingent.sort(key=lambda x: x.Cij + x.Cji, reverse=False)

    length = [e.Cij + e.Cji for e in contingent]
    S = sum(length) + weight
    n = len(contingent)
    if S < 0:
        return None

    m = None
    for i in range(n):
        previous = length[:i]
        test_sum = sum(previous) + (n - i) * length[i]
        if test_sum >= S:
            m = i
            break

    A = (S - sum(length[:m])) / (n - m)
    epsilons = {}
    for e in contingent[m:]:
        epsilons[e.j] = e.Cij + e.Cji - A

    return epsilons


def relaxTotalLength(bounds, weight):
    return null

##
# \fn relaxSearch(STN)
# \brief run relaxation algorithm on an STNU so that it becomes dynamically
#        controllable
#
# @param STN       An STNU we want to relax/process
#
# @return The dynamically controllable relaxed STNU and the number of conflict
#         need to be resolved
def relaxSearch(STN):
    result, conflicts, bounds, weight = DC_Checker(STN.copy(), report=False)

    count = 0
    cycles = []
    weights = []
    while not result:
        edges = [x[0] for x in list(bounds['contingent'].values())]
        cycles.append(edges)
        weights.append(weight)

        epsilons = optimalRelax(bounds, weight)
        if not epsilons:
            print("The STNU cannot resolve the conflict...")
            return None, 0, None

        for (i, j) in list(STN.contingentEdges.keys()):
            if j not in list(epsilons.keys()):
                continue

            edge = STN.contingentEdges[(i, j)]
            if bounds['contingent'][(i, j)][1] == 'UPPER':
                STN.modifyEdge(i, j, edge.Cij - epsilons[j])
            else:
                STN.modifyEdge(j, i, edge.Cji - epsilons[j])

        count += 1
        result, conflicts, bounds, weight = DC_Checker(
            STN.copy(), report=False)

    return STN, count, cycles, weights
