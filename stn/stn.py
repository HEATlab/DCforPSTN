import math
import json

# For testing
# import time

## \file stntools.py
#  \brief tools for working with STNs
#  \note this is legacy from the RobotBrunch project


## \class Vertex
#  \brief Represents an STN timepoint
class Vertex(object):

    ## \brief Vertex Constructor
    #  \param nodeID       The unique ID number of the vertex in the STN.
    def __init__(self, nodeID):

        ## The unique ID number of the vertex in the STN.
        self.nodeID = nodeID

    ## \brief Vertex String Representation
    def __repr__(self):
        return "Vertex {}".format(self.nodeID)

    ## \brief Return a copy of this vertex (with identical IDs, so beware).
    def copy(self):
        newVert = Vertex(self.nodeID)
        return newVert

    ## \brief Return a ready-for-json dictionary of this timepoint
    def forJSON(self):
        return {"node_id": self.nodeID}


## \class Edge
#  \brief represents an STN constraint
#  \note distribution is the name of the distribution only
class Edge(object):

    ## \brief Edge Constructor
    #  \param i            The starting node of the edge.
    #  \param j            The ending node of the edge.
    #  \param Tmin         The min time needed to go from node i to node j.
    #  \param Tmax         The max time allowed to go from node i to node j.
    #  \param type         The type of the edge: stc, stcu or pstc
    #  \param distribution The name of the distriution used in the edge.
    def __init__(self, i, j, Tmin, Tmax, type='stc', distribution=None):
        ## The starting node of the edge.
        self.i = i

        ## The ending node of the edge.
        self.j = j

        ## The maximum amount of time allotted.
        self.Cij = Tmax

        ## The negated minimum amount of times allotted.
        #  (In this form for notation)
        self.Cji = -Tmin

        ## The type of the edge
        self.type = type

        ## The string representation of the distribution
        self.distribution = distribution

    ## \brief Return a ready-for-json dictionary of this constraint
    def forJSON(self):
        json = {"first_node": self.i, "second_node": self.j, "type": self.type}

        if self.type == 'pstc' and self.distribution is not None:
            json["distribution"] = {
                "name": self.distribution,
                "type": "Empirical"
            }

        if self.Cji == float('inf'):
            json['min_duration'] = '-inf'
        else:
            json['min_duration'] = -self.Cji

        if self.Cij == float('inf'):
            json['max_duration'] = 'inf'
        else:
            json["max_duration"] = self.Cij

        return json

    ##
    #  \brief gets the weight of the edge between Vertex i and j
    #
    #  \param i the starting node of the edge
    #  \param j the ending node of the edge
    def getWeight(self, i, j):
        if i == self.i:
            return self.Cij
        else:
            return self.Cji

    ## \brief The minimum weight of the edge
    def getWeightMin(self):
        return -self.Cji

    ## \brief The maximum weight of the edge
    def getWeightMax(self):
        return self.Cij

    ## \brief Checks if edge is contingent or not
    def isContingent(self):
        return self.type != 'stc'

    ## \brief The string representation of the edge
    def __repr__(self):
        return "Edge {} => {} [{}, {}], ".format(
            self.i, self.j, -self.Cji, self.Cij) + "type: " + self.type


##
# \class STN
# \brief A representation of an entire STN.
class STN(object):

    ## The Zero Timepoint.
    Z_TIMEPOINT = Vertex(0)

    ## \brief STN constructor
    def __init__(self):
        ## A dictionary of vertices in the STN in the form {NodeID: Node_Object}
        self.verts = {}

        ## A dictionary of edges in the STN in the form
        # {(Node1, Node2): Edge_Object}
        self.edges = {}

        ## a list of uncontrollable events (NodeIDs of vertices with incoming
        # contingencies)
        self.uncontrollables = []

        ## A reverse lookup dictionary of Nodes in contingent edges in the form
        # {NodeID_End: NodeID_Start}
        self.parent = {}

        ## A dictionary of contingent edges in the form
        # {(Node1, Node2): Edge_Object}
        self.contingentEdges = {}

        ## A dictionary of requirement edges in the form
        # {(Node1, Node2): Edge_Object}
        self.requirementEdges = {}

        ## The total amount of time allowed for an STN (implemented in
        #  milliseconds)
        self.makespan = None

    # -------------------------------------------------------------------------
    # Basic functions #
    # -------------------------------------------------------------------------

    ## \brief String representation of the STN
    def __repr__(self):
        toPrint = ""
        for (i, j), edge in sorted(self.edges.items()):

            if edge.i == 0:
                toPrint += "Vertex {}: [{}, {}]".format(edge.j,\
                                                        -edge.Cji,edge.Cij)
            else:
                toPrint += "Edge {} => {}: [{}, {}]".format(
                    edge.i, edge.j, -edge.Cji, edge.Cij)

            toPrint += "\n"
        return toPrint

    ##
    # \brief Returns a copy of the STN
    #
    # @return a copy of the original STN
    def copy(self):
        newSTN = STN()
        for v in self.getAllVerts():
            newSTN.addVertex(v.nodeID)

        for e in self.getAllEdges():
            newSTN.addEdge(e.i, e.j, -e.Cji, e.Cij, e.type, e.distribution)

        newSTN.makespan = self.makespan
        return newSTN

    ##
    # \brief Generates a subSTN from a given vertex list
    #
    # @param vertexList         A List of Vertex objects that are part of the
    #                           STN.
    #
    # @return   a sub STN of the original STN with given vertex list
    def getSubSTN(self, vertexList):
        subSTN = STN()
        vertIDList = [v.nodeID for v in vertexList]

        # Add all edges vertices in stn
        for v in vertexList:
            subSTN.addCreatedVertex(v.copy())

        # Add all edges between vertices in the subSTN
        for e in self.getAllEdges():
            # Check if both ends of the edge are in the new subSTN
            if (e.i in vertIDList) and (e.j in vertIDList):
                subSTN.addCreatedEdge(e.copy())

        return subSTN

    ##
    # \brief Takes in parameters of a vertex and adds the vertex to the STN
    #
    # @param nodeID       The unique ID number of the vertex in the STN.
    #
    # @post a new vertex is added to STN
    def addVertex(self, nodeID):
        assert nodeID not in self.verts
        self.verts[nodeID] = Vertex(nodeID)

    ##
    # \brief Takes in a vertex object and adds it to the STN
    #
    # @param vertex        The vertex to be added to the STN
    #
    # @post a new vertex is added to STN
    def addCreatedVertex(self, vertex):
        nodeID = vertex.nodeID
        assert nodeID not in self.verts
        self.verts[nodeID] = vertex

    ##
    # \brief Takes in the parameters of an edge and adds the edge to the STN
    #
    # @param i            The starting node of the edge.
    # @param j            The ending node of the edge.
    # @param Tmin         The min time needed to go from node i to node j.
    # @param Tmax         The max time allowed to go from node i to node j.
    # @param type         The type of the edge: stc, stcu or pstc
    # @param distribution The name of the distribution used in the edge.
    #
    # @post               A new edge, generated from the arguments is added to
    #                     the STN.
    def addEdge(self, i, j, Tmin, Tmax, type='stc', distribution=None):
        assert i in self.verts and j in self.verts
        assert (i, j) not in self.edges and (j, i) not in self.edges
        newEdge = Edge(i, j, Tmin, Tmax, type, distribution)

        if type == 'stc':
            self.requirementEdges[(i, j)] = newEdge

        else:
            assert j not in self.uncontrollables
            self.contingentEdges[(i, j)] = newEdge
            self.uncontrollables += [j]
            self.parent[j] = i

        self.edges[(i, j)] = newEdge

    ##
    # \brief Takes in a Edge object and adds it to the STN.
    #
    # @param edge    A new Edge object to be added to the STN.
    #
    # @post          A new edge is added to the STN.
    def addCreatedEdge(self, edge):
        i = edge.i
        j = edge.j

        assert (i, j) not in self.edges and (j, i) not in self.edges
        assert i in self.verts and j in self.verts

        if edge.type == 'stc':
            self.requirementEdges[(i, j)] = edge

        else:
            assert j not in self.uncontrollables
            self.contingentEdges[(i, j)] = edge
            self.uncontrollables += [j]
            self.parent[j] = i

        self.edges[(i, j)] = edge

    # -------------------------------------------------------------------------
    # Vertex functions #
    # -------------------------------------------------------------------------

    ##
    # \brief Gets all the Nodes in the STN
    #
    # @return Returns an unordered List of all Node objects in the STN.
    def getAllVerts(self):
        return list(self.verts.values())

    ##
    #  \brief Get a list of edges incident to this node
    #
    #  \param nodeID The ID of the vertex
    #
    #  \return the list of edges incident to the input vertex
    def getEdges(self, nodeID):
        return [
            e for e in self.getAllEdges() if e.i == nodeID or e.j == nodeID
        ]

    ##
    #  \brief Returns the degree (number of edges) of a vertex
    #
    #  \param nodeID The ID of the vertex.
    #
    #  \return the total degree of the given vertex (incoming+outgoing)
    def getDegree(self, nodeID):
        return len(self.getEdges(nodeID))

    ##
    #  \brief Returns a list of nodes adjacent to a given node
    #
    #  \param nodeID The ID of the node
    #
    #  \return a list of vertices adjacent to the given vertex (all neighbors)
    def getAdjacent(self, nodeID):
        adj = []
        edges = self.getEdges(nodeID)
        for e in edges:
            if e.i != nodeID:
                adj.append(e.i)
            if e.j != nodeID:
                adj.append(e.j)
        return adj

    ##
    #  \brief Removes a node from the STP
    #
    #  \param nodeID the ID of the node to be removed
    #
    #  \post a STN with given vertex and all its edges removed
    def removeVertex(self, nodeID):
        if nodeID in self.verts:
            del self.verts[nodeID]

            if nodeID in self.uncontrollables:
                self.uncontrollables.remove(nodeID)

            toRemove = []
            for i, j in self.edges:
                if i == nodeID or j == nodeID:
                    toRemove += [(i, j)]

            for i, j in toRemove:
                del self.edges[(i, j)]

                if (i, j) in self.contingentEdges:
                    self.uncontrollables.remove(j)
                    del self.contingentEdges[(i, j)]

                if (i, j) in self.requirementEdges:
                    del self.requirementEdges[(i, j)]

    ##
    # \brief Gets a node from the STP
    #
    # @param nodeID Integer representing the gloabl ID of the node to get.
    #
    # @return Returns a Vertex with specified global ID. Returns None if it does
    #   not exist.
    def getVertex(self, nodeID):
        if nodeID in self.verts:
            return self.verts[nodeID]
        else:
            return None

    # -------------------------------------------------------------------------
    # Edge functions #
    # -------------------------------------------------------------------------

    ##
    # \brief Gets an Edge from the STP.
    #
    # \details Direction is not accounted for.
    #
    # @param i  The start Node of the edge.
    # @param j  The end Node of the edge.
    #
    # @return Returns an Edge object if one exists between i & j. If not, return
    #   None.
    def getEdge(self, i, j):
        if (i, j) in self.edges:
            return self.edges[(i, j)]
        elif (j, i) in self.edges:
            return self.edges[(j, i)]
        else:
            return None

    ##
    # \brief Gets all edges of the STP
    #
    # @return Returns list of all edges in the STP
    def getAllEdges(self):
        return list(self.edges.values())

    ##
    # \brief Gets a directed edge weight of an edge from the STP
    #
    # \details Direction does matter. If no edge exists between i & j, return
    #   'inf' If i = j, this function returns 0.
    #
    # @param i  The starting Node of the edge.
    # @param j  The ending Node of the edge.
    #
    # @return Returns a float representing the weight from Node i to Node j.
    def getEdgeWeight(self, i, j):
        e = self.getEdge(i, j)

        if e == None:
            if i == j and i in self.verts:
                return 0
            else:
                return float('inf')

        if e.i == i and e.j == j:
            return e.Cij
        else:
            return e.Cji

    ##
    # \brief Checks if an edge exists in the STP, regardless of direction.
    #
    # @param i The first Node of the edge.
    # @param j The second Node of the edge.
    #
    # @return Returns a boolean on whether or not an edge exists between the
    #   inputted nodes. Direction is not accounted for.
    def edgeExists(self, i, j):
        return ((i, j) in self.edges) or ((j, i) in self.edges)

    ##
    # \brief Updates the edge with Node objects i & j.
    #
    # \details Update if the new weight is less than the original weight.
    #    Otherwise, do nothing.
    #
    # @param i          The starting Node of the edge.
    # @param j          The ending Node of the edge.
    # @param w          The new weight to try to update with.
    # @param equality   If set to true, we'll "update" the edge even if the
    #                   original weight is the same as w (the new weight).
    #                   Nothing will actually change, but the function
    #                   will return True.
    #
    # @return Returns a boolean indicating whether or not the update
    #   actually occured.
    def updateEdge(self, i, j, w, equality=False):
        e = self.getEdge(i, j)

        if e == None:
            return False

        if e.i == i and e.j == j:
            if w < e.Cij:
                e.Cij = w
                return True
            else:
                if equality:
                    return w <= e.Cij
                return False

        else:
            if w < e.Cji:
                e.Cji = w
                return True
            else:
                if equality:
                    return w <= e.Cji
                return False

    ##
    # \brief Update the weight of edge with node objects i & j.
    #
    # \details Mofidy the constraint regardless of what new edge weight is.
    #
    # @param i  The starting Node of the edge.
    # @param j  The ending Node of the edge.
    # @param w  The new weight to try to update with.
    #
    # @return Returns boolean whether or not the update actually occured.
    def modifyEdge(self, i, j, w):
        e = self.getEdge(i, j)

        if e == None:
            return False

        else:
            if e.i == i and e.j == j:
                e.Cij = w
            else:
                e.Cji = w

            return True

    ##
    # \brief Get all incoming edges for a given vertex
    #
    # @param nodeID     Integer representing the gloabl ID of the node.
    #
    # @return Returns a list of Edge objects that are incoming edges for
    #         the input vertex
    def getIncoming(self, nodeID):
        return [e for e in self.getAllEdges() if e.j == nodeID]

    ##
    # \brief Get all outgoing edges for a given vertex
    #
    # @param nodeID     Integer representing the gloabl ID of the node.
    #
    # @return Returns a list of Edge objects that are outgoing edges for
    #         the input vertex
    def getOutgoing(self, nodeID):
        return [e for e in self.getAllEdges() if e.i == nodeID]

    ##
    # \brief Get all incoming edges for a given vertex
    #
    # \details the input vertex needs to be an uncontrollable nodes. Otherwise
    #          the function would raise AssertionError. Note, an uncontrollable
    #          event cannot have two incoming contingent edges
    #
    # @param nodeID     The node ID for an uncontrollable event.
    #
    # @return Returns the only incoming contingent edge for the input
    #         uncontrollable vertex
    def getIncomingContingent(self, nodeID):
        assert nodeID in self.uncontrollables
        ctg = [self.contingentEdges[(i,j)] for (i,j) in self.contingentEdges \
                                                            if j == nodeID]

        if len(ctg) != 1:
            print('E: {} incoming contingent edges!\n{}'.format(len(ctg), ctg))
            raise AssertionError

        return ctg[0]

    ##
    # \brief Remove an input egde from the STN
    #
    # \details Direction does not matter.
    #
    # @param i  The node ID for the first node in the edge
    # @param j  The node ID for the second node in the edge
    #
    # @post A STN with edge between node i and j removed
    def removeEdge(self, i, j):
        if not self.edgeExists(i, j):
            raise ValueError("The input edge does not exist")

        to_remove = (i, j) if (i, j) in self.edges else (j, i)
        del self.edges[to_remove]

        if to_remove in self.contingentEdges:
            del self.contingentEdges[to_remove]
            self.uncontrollables.remove(to_remove[1])

        else:
            del self.requirementEdges[to_remove]

    # -------------------------------------------------------------------------
    # Other functions #
    # -------------------------------------------------------------------------

    ##
    # \brief Convert an input STN to distance matrix (list of lists) format
    #
    # @return a distance matrix that represents an input STN
    def toMatrix(self):
        num = len(self.verts) if 0 in self.verts else len(self.verts) + 1
        matrix = []

        for i in range(num):
            row = []
            for j in range(num):
                dis = self.getEdgeWeight(i, j)

                if j == 0 and dis == float('inf'):
                    dis = 0

                row.append(dis)
            matrix.append(row)

        return matrix

    ##
    # \brief set the makespan of the STN
    #
    # @param makespan a number representing the total time allowed for the STN
    #                 in milliseconds
    #
    # @post STN with makespan set to the input value
    def setMakespan(self, makespan):
        self.makespan = makespan

        if 0 not in self.verts:
            self.addVertex(self.Z_TIMEPOINT)

        count = 0
        for vert in self.verts:
            if vert != 0:
                val = self.getEdgeWeight(0, vert)
                if val > makespan:
                    if (0, vert) in self.edges:
                        self.modifyEdge(0, vert, makespan)
                    else:
                        self.addEdge(0, vert, 0, makespan, type='stc')
                else:
                    count += 1

        if count == len(self.verts) - 1:
            for vert in self.verts:
                if vert != 0:
                    self.modifyEdge(0, vert, makespan)

    ##
    # \brief Return a ready-for-json dictionary of this STN
    #
    # @return a json object that is ready to store for future use
    def forJSON(self):
        jsonSTN = {}

        # Add the vertices
        jsonSTN['nodes'] = []
        for v in self.getAllVerts():
            if v.nodeID == 0:
                continue

            json = v.forJSON()
            jsonSTN['nodes'].append(json)

        # Add the edges
        jsonSTN['constraints'] = []
        for c in self.getAllEdges():
            jsonSTN['constraints'].append(c.forJSON())

        return jsonSTN

    ##
    # \brief Write the STN object to a json file for later use
    #
    # @param name  The name of the output json file
    # @param path  The directory we should store the output file to
    #
    # @post a JSON file containing STN info stored to desired directory
    def toJSON(self, name, path):
        jsonstr = self.forJSON()

        with open(path + name, 'w') as outfile:
            json.dump(jsonstr, outfile)

        outfile.close()

    ##
    #  \brief Runs the Floyd-Warshal algorithm on an STN
    #
    #  @return Return a STN object that is the minimal network of the original
    #          STN. If the input STN is not consistent, return None
    def minimal(self):
        minSTN = self.copy()

        verts = list(range(len(minSTN.verts)))
        B = [[minSTN.getEdgeWeight(i, j) for j in verts] for i in verts]

        for k in verts:
            for i in verts:
                for j in verts:
                    B[i][j] = min(B[i][j], B[i][k] + B[k][j])
                    minSTN.updateEdge(i, j, B[i][j])
                    # minSTN.updateEdge()

        for e in minSTN.getAllEdges():
            if e.getWeightMin() > e.getWeightMax():
                return None

        return minSTN

    ##
    # \brief Run minimal STN and check if the STN is consistent or not
    #
    # @return Returns true if the given STN is consistent. Otherwise, returns
    #         False
    def isConsistent(self):
        return self.minimal() != None

    ##
    # \fn altConsistent(self)
    #
    # Uses Floyd Warshall to check for consistency. Is hopefull quicker than the
    # the old isConsistent function.
    def altConsistent(self):
        # Number of vertices
        events = list(self.verts.keys())
        n = len(events)

        # Initialize table
        dp = {}
        for i in events:
            for j in events:
                dp[i, j] = self.getEdgeWeight(i, j)

        # Future entries
        for k in range(1, n + 1):
            for i in events:
                for j in events:
                    mid = events[k - 1]
                    new_weight = dp[i, mid] + dp[mid, j]
                    if new_weight < dp[i, j]:
                        dp[i, j] = new_weight

            for i in events:
                if dp[i, i] < 0:
                    return False

        return True

    ##
    # \brief check whether the input STNU is strongly controllable or not
    #
    # \details   An STN is strongly controllable if any assignment of
    #            values for executable timepoints is guaranteed to be
    #            consistent with all constraints (irrespective osf contingent
    #            edges) (copied from Lund et al. 2017).
    #
    # @param debug      Flag indicating whether want to print message for debug
    # @param returnSTN  Flg indicating wehther want to return the reduced STN
    #                   with controllable events or not
    #
    # @return Returns True if STNU is strongly controllable, and False otherwise
    #         If returnSTN is True, then also return the reduced STN
    #
    # NOTE: This function has not been tested, and might not be reliable!!
    def isStronglyControllable(self, debug=False, returnSTN=False):

        if not self.isConsistent():
            if returnSTN:
                return False, None
            else:
                return False

        newSTN = STN()
        for v in self.getAllVerts():
            if v.nodeID not in self.uncontrollables:
                newSTN.addVertex(v.nodeID)

        for e in self.getAllEdges():
            if e.type == 'stc':
                if debug:
                    print("dealing with edge {}->{}".format(e.i, e.j))
                if e.i in self.uncontrollables:
                    cont_edge = self.getEdge(self.parent[e.i], e.i)
                    i = self.parent[e.i]
                    l_i = cont_edge.getWeightMin()
                    u_i = cont_edge.getWeightMax()
                else:
                    i = e.i
                    l_i = 0
                    u_i = 0

                if e.j in self.uncontrollables:
                    cont_edge = self.getEdge(self.parent[e.j], e.j)
                    j = self.parent[e.j]
                    l_j = cont_edge.getWeightMin()
                    u_j = cont_edge.getWeightMax()
                else:
                    j = e.j
                    l_j = 0
                    u_j = 0

                # This is taken from the 1998 paper by
                # Vidal et al. on controllability
                if debug:
                    print("Cij: {}  Cji: {}  l_i: {}  u_i: {}  l_j: {}  \
                                u_j: {}".format(e.Cij, e.Cji, l_i, u_i, l_j,
                                                u_j))

                lower_bound = -e.Cji + u_i - l_j
                upper_bound = e.Cij + l_i - u_j

                if (i, j) in newSTN.edges or (j, i) in newSTN.edges:
                    newSTN.updateEdge(i, j, upper_bound)
                    newSTN.updateEdge(j, i, -lower_bound)
                    if debug:
                        print("updated edge {}->{}: [{},{}]".format(
                            i, j, lower_bound, upper_bound))
                else:
                    newSTN.addEdge(i, j, lower_bound, upper_bound)
                    if debug:
                        print("added edge {}->{}: [{},{}]".format(
                            i, j, lower_bound, upper_bound))

        if debug:
            print(newSTN)

        if returnSTN:
            if newSTN.isConsistent():
                return True, newSTN
            else:
                return False, None

        return newSTN.isConsistent()
