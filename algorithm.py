from util import *

##
# \file algorithm.py
# \brief An implementation of the Conflict Generation Algorithm by William
# \note Read Williams 2017 paper for more details on the algorithm implemented:
#       https://www.ijcai.org/proceedings/2017/0598.pdf

# -------------------------------------------------------------------------
# Helper Functions
# -------------------------------------------------------------------------


##
# \fn extractEdgePath(s, v, labelDist, unlabelDist)
# \brief Extract the edges along the path from vertex v to vertex s
#
# @param s              end node of the path
# @param v              start node of the path
# @param labelDist      a dictionary holding labeled edges
# @param unlabelDist    a dictionary holding unlabeled edges
#
# @return a list of DC STN edge object that are edges along the path from
#         vertex v to s
def extractEdgePath(s, v, labelDist, unlabelDist):
    result = []

    while True:
        # check which edge dictionary we should use
        if v in labelDist and v in unlabelDist:
            if not labelDist[v][1]:
                distArray = unlabelDist
            if not unlabelDist[v][1]:
                distArray = labelDist
        elif v in labelDist:
            distArray = labelDist
        else:
            distArray = unlabelDist

        # add the current edge to the path and continue iterating until
        # we reach the end node
        weight, edge = distArray[v]
        result.append(edge)
        v = edge.j

        if edge.j == s:
            break

    return result


##
# \fn extractConflict(edges, novel, preds)
# \brief Extract edges along the detected semi-reducible negative cycle
#
# @param edges          a list of edges along semi-reducible negative cycle
# @param novel          a list of novel edges
# @param preds          a dictionary of predecessors
#
# @return A list of original edges along the semi-reducible negative cycle
def extractConflict(edges, novel, preds):
    result = []

    for edge in edges:
        entry = (edge.i, edge.j, edge.weight)
        if entry not in novel:
            result.append(edge)
        else:
            result += resolveNovel(edge, novel, preds)

    return result


##
# \fn resolveNovel(e, novel, preds)
# \brief Extract original edges that derive the input novel path
#
# @param e              a novel edge
# @param novel          a list of novel edges
# @param preds          a dictionary of predecessors
#
# @return a list of original edges along the path represented by input novel
#         edge
def resolveNovel(e, novel, preds):
    result = []
    entry = (e.i, e.j, e.weight)
    if entry not in novel:
        result.append(e)
        return result

    labelDist, unlabelDist = preds[e.j]
    distArray = labelDist if e.i in labelDist else unlabelDist

    weight, new_edge = distArray[e.i]

    end = new_edge.j
    while end != e.j:
        add = distArray[end][1]
        if (add.i, add.j, add.weight) in novel:
            result = result + resolveNovel(add, novel, preds)
        else:
            result.append(add)
        end = add.j

    result = result + resolveNovel(new_edge, novel, preds)

    return result


##
# \fn getFinalResult(conflicts, STN, D, report=True)
# \brief extract information about which constraint in the original STNU can
#        be relaxed to resolve the conflict
# \details edited on 7/3/19
#
# @param conflicts      A list of labeled edges along the negative cycle
# @param STN            The original STNU
# @param D              A dictionary containing info about added vertices
# @param C              A list containing all contingent edges
# @param report         Flag indicating whether to print message or not
#
# @return A dictionary containing information about the original constraints
#         in the original STNU that we can relax and whether LOWER or UPPER
#         bound can be relaxed
def getFinalResult(conflicts, STN, D, C, report=True):
    ## Initialize the result dictionary
    result = {}
    result['requirement'] = {}
    result['contingent'] = {}

    ## Loop through all labeled edges in the conflict and find the
    ## corresponding original edge in the STNU
    for edge in conflicts:
        start = edge.i
        end = edge.j

        if start in D:
            continue
        elif end in D:
            e = D[end]
            if start == e.j and edge.type == edgeType.UPPER:
                result['contingent'][(e.i, e.j)] = (e, 'UPPER')
            elif start == e.i:
                result['contingent'][(e.i, e.j)] = (e, 'LOWER')
        elif (start,end) in C:
            e = C[(start, end)]
            if edge.type == edgeType.UPPER:
                result['contingent'][(e.i, e.j)] = (e, 'UPPER')
            else:
                result['contingent'][(e.i, e.j)] = (e, 'LOWER')
        else:
            e = STN.getEdge(start, end)
            if start == e.i:
                result['requirement'][(e.i, e.j)] = (e, 'UPPER')
            else:
                result['requirement'][(e.i, e.j)] = (e, 'LOWER')

    if report:
        print("Reporting Conflicts:")

        print("\nThe requirement edges we can relax are: ")
        for i, j in list(result['requirement'].keys()):
            edge, bound = result['requirement'][(i, j)]
            print(edge, bound)

        print("\nThe contingent edges we can relax are: ")
        for i, j in list(result['contingent'].keys()):
            edge, bound = result['contingent'][(i, j)]
            print(edge, bound)

    return result


# -------------------------------------------------------------------------
# Major Functions: DCDijkstra and DC_Checker
# -------------------------------------------------------------------------


##
# \fn DCDijkstra(G, start, preds, novel, callStack, negNodes)
# \brief Determine if there is any semi-reducible negative cycles in an input
#        labeled graph start with the input vertex and identify the edges
#        along the cycle if there is one
#
# @param G              an input labeled graph to check
# @param start          start node
# @param preds          a dictonary of predecessors edges for each vertex
# @param novel          a list of new edges added
# @param callStack      a list keeping track of the recurrence order
# @param negNodes       a list of negative nodes
#
# @return Return True if there is no semi-reducible negative cycle in the input
#         labeled graph. Otherwise, return False, the edges along the
#         negative cycle and the end node
def DCDijkstra(G, start, preds, novel, callStack, negNodes):
    Q = PriorityQueue()
    labelDist = {}
    unlabelDist = {}
    labelDist[start] = (0, None)
    unlabelDist[start] = (0, None)
    epsilon = 1e-5

    for edge in G.incomingEdges(start):
        if edge.weight < 0:
            Q.push((edge.i, edge.parent), edge.weight)
            if edge.parent == None:
                unlabelDist[edge.i] = (edge.weight, edge)
            else:
                labelDist[edge.i] = (edge.weight, edge)

    if start in callStack[1:]:
        return False, [], start

    preds[start] = (labelDist, unlabelDist)

    while not Q.isEmpty():
        weight, (v, label) = Q.pop()

        if weight >= -epsilon:
            G.addEdge(v, start, weight)
            novel.append((v, start, weight))
            continue

        if v in negNodes:
            newStack = [v] + callStack
            result, edges, end = DCDijkstra(G, v, preds, novel, \
                                                        newStack, negNodes)

            if not result:
                if end != None:
                    edges += extractEdgePath(start, v, labelDist, unlabelDist)
                if end == start:
                    end = None
                return False, edges, end

        for edge in G.incomingEdges(v):
            if edge.weight >= 0 and (edge.type != edgeType.LOWER or \
                                                        edge.parent != label):
                w = edge.weight + weight
                distArray = labelDist if label != None else unlabelDist

                if edge.i not in distArray or w < distArray[edge.i][0]:
                    distArray[edge.i] = (w, edge)
                    Q.addOrDecKey((edge.i, label), w)

    negNodes.remove(start)
    return True, [], None


##
# \fn DC_Checker(STN, report=True)
# \brief Check whether an input STNU is dynamically controllable
#
# \details An STNU is dynamically controllable if there is not semi-reducible
#          negative cycles in its labeled graph.
#
# @param STN    an STN which we want to test
#
# @return Return True if the input STNU is dynamically controllable. Otherwise,
#         return False, conflicts (in labeled graph), conflicts in original STNU,
#         and weights of the negative cycle (conflict).
def DC_Checker(STN, report=True):
    G, C, D = normal(STN.copy())
    negNodes = G.getNegNodes()
    novel = []
    preds = {}

    for v in negNodes:
        result, edges, end = DCDijkstra(G, v, preds, novel, \
                                                    [v], negNodes.copy())

        if not result:
            conflicts = extractConflict(edges, novel, preds)
            bounds = getFinalResult(conflicts, STN, D, C, report=report)

            weight = 0
            for e in edges:
                weight += e.weight

            return False, conflicts, bounds, weight

    return True, [], {}, 0
