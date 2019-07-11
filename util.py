from stn import STN
from stn import loadSTNfromJSONfile
from dc_stn import DC_STN, edgeType
import random
import heapq
import math
import json
import glob
import os
import sys

##
# \file util.py
# \brief helpful functions


##
# \fn STNtoDCSTN(S)
# \brief Convert an STN object to a DC_STN object
#
# @param S  An input STN object to convert
#
# @return A DC_STN object that has the same vertices and edges as the input
def STNtoDCSTN(S):
    new_STN = DC_STN()
    for edge in list(S.edges.values()):
        new_STN.addEdge(edge.i, edge.j, edge.Cij)
        new_STN.addEdge(edge.j, edge.i, edge.Cji)
        if edge.isContingent():
            new_STN.addEdge(
                edge.i,
                edge.j,
                -edge.Cji,
                edge_type=edgeType.LOWER,
                parent=edge.j)
            new_STN.addEdge(
                edge.j,
                edge.i,
                -edge.Cij,
                edge_type=edgeType.UPPER,
                parent=edge.j)
    return new_STN


##
# \fn dc_checking(STN)
# \brief Check if an input STN is dynamically controllable
#
# @param S  An input STN object to check
#
# @return Return True if the input STN is dynamically controllable. Return
#         False otherwise.
def dc_checking(STN):
    return STNtoDCSTN(STN).is_DC()


##
# \fn dc_checking_file(filename)
# \brief Check if an input STN is dynamically controllable
#
# @param filename  An input STN object to check
#
# @return Return True if the input STN is dynamically controllable. Return
#         False otherwise.
def dc_checking_file(filename):
    STN = loadSTNfromJSONfile(filename)
    return dc_checking(STN)


##
# \fn normal(STN)
# \brief convert a given STNU to its normal form labeled graph using DC_STN
#
# \details  a normal form labeled graph is needed for Williams algorithm.
#           The robotbrunch DC_STN class keeps track of the upper, lower case
#           edges, so we need to use it for labeled graph.
#           Changed on 7/3/19
#
# @param STN  an STN to convert
#
# @return Return the DC_STN that is the normal form labeled graph for input STN,
#         a list storing all contingent edges of form [0, x],
#         and a dictionary of just vertices added to create normal form)
def normal(STN):
    new = DC_STN()
    for i in list(STN.verts.keys()):
        new.addVertex(i)

    changed = {}
    contingents = {}
    for e in list(STN.edges.values()):
        if not e.isContingent():
            new.addEdge(e.i, e.j, e.Cij)
            new.addEdge(e.j, e.i, e.Cji)
        else:
            contingents[(e.i,e.j)] = e
            contingents[(e.j,e.i)] = e
            if e.Cji != 0:
                new_vert = len(new.verts)
                changed[new_vert] = e
                new.addVertex(new_vert)
                new.addEdge(e.i, new_vert, -e.Cji)
                new.addEdge(new_vert, e.i, e.Cji)

                upper = e.Cij + e.Cji
                new.addEdge(new_vert, e.j, upper)
                new.addEdge(e.j, new_vert, 0)
                new.addEdge(
                    new_vert, e.j, 0, edge_type=edgeType.LOWER, parent=e.j)
                new.addEdge(
                    e.j,
                    new_vert,
                    -upper,
                    edge_type=edgeType.UPPER,
                    parent=e.j)
            else:
                new.addEdge(e.i, e.j, e.Cij)
                new.addEdge(e.j, e.i, e.Cji)
                new.addEdge(
                    e.i, e.j, -e.Cji, edge_type=edgeType.LOWER, parent=e.j)
                new.addEdge(
                    e.j, e.i, -e.Cij, edge_type=edgeType.UPPER, parent=e.j)
    return new, contingents, changed


##
# \class PriorityQueue
# \brief A simple Priority Queue implementation that pop the item with the
#        lowest priority
# \note  This code is adapted from UC Berkeley AI course project util.py file
class PriorityQueue:

    ##
    # \brief PriorityQueue Constructor
    def __init__(self):
        self.queue = []

    ##
    # \brief Push an element with given priority into the Priority queue
    #
    # @param data         An element to be added
    # @param priority     The priority associated with input data
    #
    # @post A Priority Queue object with input data added
    def push(self, data, priority):
        heapq.heappush(self.queue, (priority, data))

    ##
    # \brief Pop the element with the lowest priority in Priority Queue
    #
    # @return A tuple in the form of (priority, data)
    def pop(self):
        return heapq.heappop(self.queue)

    ##
    # \brief Check if the Priority Queue has element inside
    #
    # @return Return True if the queue is empty and return False otherwise
    def isEmpty(self):
        return len(self.queue) == 0

    ##
    # \brief Add or decrease the priority of an input data
    #
    # \detail If input data is already in the queue, modify its priority if the
    #         input priority is strictly lower than its original priority,
    #         do nothing otherwise. If input data is not in the queue, push it
    #         to the queue with given priority.
    #
    # @param data         An element to be added/modified
    # @param priority     The priority associated with input data
    #
    # @post A Priority Queue object with input data added/modified
    def addOrDecKey(self, data, priority):

        for i, (p, d) in enumerate(self.queue):
            if d == data:
                if priority >= p:
                    break

                del self.queue[i]
                self.push(data, priority)
                break
        else:
            self.push(data, priority)
