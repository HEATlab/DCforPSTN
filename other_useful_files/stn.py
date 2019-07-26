"""
File:
    This file holds the core classes needed to construct STNs and PSTNs.
"""

import math

from distempirical import norm_sample, uniform_sample

MAX_FLOAT = 1.7976931348623157e+308


class Vertex(object):
    """Represents a timepoint in an STN"""
    # \brief Vertex Constructor
    #  \param nodeID The unique ID number of the vertex in the STN.
    #  \param localID The ID number of the vertex for only the agent
    #  that owns the vertex.
    #  \param ownerID The ID number of the agent that owns the vertex.
    #  \param location The grid point number indicating the physical
    #         location of the node.
    #  \param executed Flag for the simulator that indicates if the vertex
    #  has
    #                      been executed.

    def __init__(self, nodeID, ownerID, location=None):
        # The unique ID number of the vertex in the STN.
        self.nodeID = nodeID
        # The ID number of the agent that owns the vertex.
        self.ownerID = ownerID
        # The grid point number indicating the physical location of the node.
        self.location = location
        # Flag for the simulator that indicates if the vertex has been
        # executed.
        self.executed = False

    # \brief Vertex String Representation
    def __repr__(self):
        return "Vertex {} owned by agent {}".format(self.nodeID, self.ownerID)

    def __eq__(self, other):
        if other is None:
            return False
        return (self.nodeID == other.nodeID and
                self.ownerID == other.ownerID and
                self.location == other.location and
                self.executed == other.executed)

    # \brief Return a copy of this vertex (with identical IDs and
    #   execution, so beware).
    def copy(self):
        newVert = Vertex(self.nodeID, self.ownerID, self.location)
        if self.executed:
            newVert.execute()
        return newVert

    # \brief Return a ready-for-json dictionary of this timepoint
    def for_json(self):
        try:
            local_ID = self.localID
        except AttributeError:
            local_ID = None
        return {"node_id": self.nodeID, "local_id": local_ID,
                "owner_id": self.ownerID, "location": self.location,
                "executed": self.executed}

    # \brief sets the vertex as executed by the agents
    def execute(self):
        self.executed = True

    def is_executed(self):
        return self.executed


# \class Edge
#  \brief represents an STN constraint
#  \note distribution is the name of the distribution only
class Edge(object):

    # How many resamples should we make until we give up?
    RESAMPLES_UNTIL_QUIT = 100

    # \brief Edge Constructor
    #  \param i            The starting node of the edge.
    #  \param j            The ending node of the edge.
    #  \param Tmin         The min time needed to go from node i to node j.
    #  \param Tmax         The max time allowed to go from node i to node j.
    #  \param distribution The name of the distriution used in the edge.
    #  \param fake         Flag indicating if the edge is fake
    def __init__(self, i, j, Tmin, Tmax, distribution=None, fake=False):
        self.i = i
        """The starting node id for the edge."""
        self.j = j
        """The ending node id for the edge."""
        self.Cij = Tmax
        """The maximum amount of time alloted."""
        self.Cji = -Tmin
        """The negated minimum amount of time allotted."""
        self.distribution = distribution
        """The string representation for this edge's probability
            distribution
        """
        self._sampled_time = 0

    def __eq__(self, other):
        if other is None:
            return False
        return (self.i == other.i and self.j == other.j and
                self.Cij == other.Cij and self.Cji == other.Cji and
                self.distribution == other.distribution and
                self._sampled_time == other._sampled_time)

    # \brief Return a ready-for-json dictionary of this constraint
    def for_json(self):
        json = {"first_node": self.i,
                "second_node": self.j}

        if self.distribution is not None:
            json["distribution"] = {"name": self.distribution,
                                    "type": "Empirical"}

        if self.Cji == float('inf'):
            json['min_duration'] = '-inf'
        else:
            json['min_duration'] = -self.Cji

        if self.Cij == float('inf'):
            json['max_duration'] = 'inf'
        else:
            json["max_duration"] = self.Cij

        return json

    # \brief gets the weight of the edge between Vertex i and j
    #  \param i the starting node of the edge
    #  \param j the ending node of the edge
    def get_weight(self, i, j):
        if i == self.i:
            return self.Cij
        else:
            return self.Cji

    # \brief The minimum weight of the edge
    def get_weight_min(self):
        return -self.Cji

    # \brief The maximum weight of the edge
    def get_weight_max(self):
        return self.Cij

    # \brief Checks if edge has a probability distribution
    def is_contingent(self):
        return self.distribution is not None

    # \brief The string representation of the edge
    def __repr__(self):
        return "Edge {} => {} [{}, {}]".format(self.i, self.j, -self.Cji,
                                               self.Cij)

    def resample(self, random_state):
        """ Retrieve a new sample from this contingent edge.
        Raises an exception if this is a requirement edge.

        Returns:
            A float selected from this edge's contingent distribution.
        """
        sample = None
        if not self.is_contingent():
            raise TypeError("Cannot sample requirement edge")
        if self.distribution[0] == "N":
            sample = norm_sample(self.mu, self.sigma, random_state)
        elif self.distribution[0] == "U":
            sample = uniform_sample(self.dist_lb, self.dist_ub, random_state)
        # We have to use integers because of rounding errors.
        self._sampled_time = round(sample)
        return self._sampled_time

    def sampled_time(self) -> float:
        """Gets the sampled time for this contingent edge.
        If this edge has not been resampled, it will return 0.0 or the last
        sample collected.
        """
        return self._sampled_time

    def copy(self):
        new_edge = Edge(self.i, self.j, -self.Cji, self.Cij, self.distribution)
        new_edge._sampled_time = self._sampled_time
        return new_edge

    def dtype(self):
        """Returns the distribution edge type as a String. If no there is
            distribution for this edge, return None.
        """
        if self.distribution is None:
            return None
        if self.distribution[0] == "U":
            return "uniform"
        elif self.distribution[0] == "N":
            return "gaussian"
        else:
            return "unknown"

    @property
    def mu(self):
        name_split = self.distribution.split("_")
        if len(name_split) != 3 or name_split[0] != "N":
            raise ValueError("No mu for non-normal dist")
        return float(name_split[1]) * 1000

    @property
    def sigma(self):
        name_split = self.distribution.split("_")
        if len(name_split) != 3 or name_split[0] != "N":
            raise ValueError("No sigma for non-normal dist")
        return float(name_split[2]) * 1000

    @property
    def dist_ub(self):
        name_split = self.distribution.split("_")
        if len(name_split) != 3 or name_split[0] != "U":
            raise ValueError("No upper bound for non-uniform dist")
        return float(name_split[2]) * 1000

    @property
    def dist_lb(self):
        name_split = self.distribution.split("_")
        if len(name_split) != 3 or name_split[0] != "U":
            raise ValueError("No lower bound for non-uniform dist")
        return float(name_split[1]) * 1000

    def cap(self):
        """Caps this edge's Cij and Cji properties to a "max" floating point
            value.
        """
        self.Cij = min(MAX_FLOAT, max(-MAX_FLOAT, self.Cij))
        self.Cji = min(MAX_FLOAT, max(-MAX_FLOAT, self.Cji))

##
# \class STN
# \brief A representation of an entire STN.


class STN(object):

    # The Zero Timepoint.
    Z_TIMEPOINT = Vertex(0, None, None)

    # \brief STN constructor
    def __init__(self):
        # A dictionary of vertices in the STN in the form {NodeID: Node_Object}
        self.verts = {}

        # A dictionary of edges in the STN in the form
        # {(Node1, Node2): Edge_Object}
        self.edges = {}

        # A reverse lookup dictionary of Nodes in contingent edges in the form
        # {NodeID_End, NodeID_Start}
        self.parent = {}

        self.received_timepoints = []
        """Holds a list of vertex IDs of received/contingent timepoints"""

        # A dictionary of contingent edges in the form
        # {(Node1, Node2): Edge_Object}
        self.contingent_edges = {}

        # A dictionary of interagent edges in the form
        # {(Node1, Node2): Edge_Object}
        self.interagent_edges = {}

        # A dictionary of requirement edges in the form
        # {(Node1, Node2): Edge Object}
        self.requirement_edges = {}

        # The total amount of time allowed for an STN (implemented in
        # milliseconds)
        self.makespan = None

        # List of agents involved in the STN.
        self.agents = []

        self.name = "Unnamed STN"
        """Identifying name of the STN"""

    # \brief String representation of the STN
    def __str__(self):
        to_print = ""
        for (i, j), edge in sorted(self.edges.items()):
            if edge.i == 0:
                to_print += "Vertex {}: [{}, {}]".format(edge.j,
                                                         -edge.Cji, edge.Cij)
                if edge.distribution is not None:
                    to_print += " ({})".format(edge.distribution)
                if self.get_vertex(j).is_executed():
                    to_print += " Ex"
            else:
                to_print += "Edge {} => {}: [{}, {}]".format(
                    edge.i, edge.j, -edge.Cji, edge.Cij)
                if edge.distribution is not None:
                    to_print += " ({})".format(edge.distribution)
            to_print += "\n"
        return to_print

    ##
    # \fn copy
    # \brief Returns a copy of the STN
    #
    def copy(self):
        new_stn = STN()
        for v in self.get_all_verts():
            new_stn.add_vertex(v.nodeID, v.ownerID, v.location)
            new_stn.verts[v.nodeID].executed = v.executed

        for e in self.get_all_edges():
            ecopy = e.copy()
            new_stn.add_created_edge(ecopy)

        # Copy the agents list over
        new_stn.agents = list(self.agents)
        new_stn.makespan = self.makespan
        return new_stn

    ##
    # \fn getAgentSubSTN
    #
    # \brief Returns a new STN that represents a single agent's subSTN
    #
    # @param agent              int representing the owner ID of the vertices
    #                           we want to
    #                           grab.
    # @param includeZTimepoint  Include the zero timepoint when generating a
    #                           new subSTN.
    #                           Default set to True.
    # @return                   Returns an STN which is a subSTN of this
    #                           current object.
    #                           Every vertex belonging to 'agent' will be part
    #                           of this
    #                           new subSTN, and so will all interagent edges
    #                           and tris.
    def get_agent_substn(self, agent, includeZTimepoint):
        # Just get a subSTN where all the verts belong to one agent.
        return self.get_substn(self.getAgentVerts(agent), includeZTimepoint)

    ##
    # \fn getSubSTN
    #
    # \brief Generates a subSTN from a given vertex list
    #
    # @param vertexList         A List of Vertex objects that are part of the
    #                           STN.
    # @param includeZTimepoint  Boolean, do we want to add the Z-timepoint?
    def get_substn(self, vertexList, includeZTimepoint):
        subSTN = STN()
        vertIDList = [v.nodeID for v in vertexList]
        if includeZTimepoint:
            vertIDList.append(0)

        agentsFound = []
        for v in vertexList:
            subSTN.add_created_vertex(v.copy())
            if not (v.ownerID in agentsFound and v.ownerID is not None):
                agentsFound.append(v.ownerID)

        # Add and execute the zero timepoint if desired.
        if includeZTimepoint:
            subSTN.add_vertex(0, None, None)
        # Add all edges between vertices in the subSTN
        for e in self.get_all_edges():
            # Check if both ends of the edge are in the new subSTN
            if (e.i in vertIDList) and (e.j in vertIDList):
                # TODO: should use copy()
                subSTN.add_created_edge(e.copy())

        # Assign the agents used in the STN
        subSTN.agents = agentsFound

        return subSTN

    ##
    # \brief Takes in parameters of a vertex and adds the vertex to the STN
    #
    # @param nodeID       The unique ID number of the vertex in the STN.
    # @param localID      The ID number of the vertex for only the agent that
    #                     owns the vertex.
    # @param ownerID      The ID number of the agent that owns the vertex.
    # @param location     The grid point number indicating the physical
    #                     location of the node.
    def add_vertex(self, nodeID, ownerID, location=None):
        self.verts[nodeID] = Vertex(nodeID, ownerID, location)

    ##
    # \fn add_created_vertex
    # \brief Takes in a vertex and adds it to the STN
    #
    # @param vertex        The vertex to be added to the STN

    def add_created_vertex(self, vertex):
        nodeID = vertex.nodeID
        self.verts[nodeID] = vertex

    def add_edge(self, i, j, Tmin, Tmax, distribution=None):
        """Takes in the parameters of an edge and adds the edge to the STN

        Args:
            i: The starting node of the edge.
            j: The ending node of the edge.
            Tmin: The min time needed to go from node i to node j.
            Tmax: The max time allowed to go from node i to node j.
            distribution: The name of the distribution used in the edge.

        Post:
            A new edge, generated from the arguments is added to the STN.
        """
        if i not in self.verts or j not in self.verts:
            raise ValueError("Vertex pair does not exist")
        new_edge = Edge(i, j, Tmin, Tmax, distribution)
        self.edges[(i, j)] = new_edge
        if distribution is not None:
            self.contingent_edges[(i, j)] = new_edge
            self.received_timepoints += [j]
            self.parent[j] = i
        elif self.get_vertex(i).ownerID != self.get_vertex(j).ownerID \
                and self.get_vertex(i).ownerID is not None \
                and self.get_vertex(j).ownerID is not None:
            self.interagent_edges[(i, j)] = new_edge
        else:
            self.requirement_edges[(i, j)] = new_edge

    def add_created_edge(self, edge):
        """Takes an existing Edge object and adds it to the STN.

        Direction matters. Will add to the needed contingent, interagent,
        or requirement dicts.

        Args:
            edge (Edge): Edge to add.
        """
        i = edge.i
        j = edge.j

        self.edges[(i, j)] = edge
        if edge.distribution is not None:
            self.contingent_edges[(i, j)] = edge
            self.received_timepoints.append(j)
            self.parent[j] = i
        elif self.get_vertex(i).ownerID != self.get_vertex(j).ownerID and \
                self.get_vertex(i).ownerID is not None:
            self.interagent_edges[(i, j)] = edge
        else:
            self.requirement_edges[(i, j)] = edge

    # -------------------------------------------------------------------------
    # Agent functions #
    # -------------------------------------------------------------------------

    ##
    # \fn getAgentVerts
    # \brief Gets the vertices an agent owns
    #
    # @param agentID Integer representing the ID of the agent.
    # @return Returns a list of Nodes that belong to the specified agent.

    def getAgentVerts(self, agentID):
        return [v for v in self.get_all_verts() if v.ownerID == agentID]

    ##
    # \fn get_all_verts
    # \brief Gets all the Nodes in the STN
    #
    # @return Returns an unordered List of all Node objects in the STN.
    def get_all_verts(self):
        return list(self.verts.values())

    # \brief Return a list of edges incident to this node
    #  \param nodeID The ID of the vertex
    def get_edges_incident(self, nodeID):
        return [e for e in self.get_all_edges() if e.i == nodeID
                or e.j == nodeID]

    # \brief Returns the degree (number of edges) of a vertex
    #  \param nodeID The ID of the vertex.

    def get_degree(self, node_id):
        return len(self.get_edges_incident(node_id))

    # \brief Returns a list of nodes adjacent to a given node
    #  \param nodeID The ID of the node

    def get_adjacent(self, nodeID):
        adj = []
        edges = self.get_edges_incident(nodeID)
        for e in edges:
            if e.i != nodeID:
                adj.append(e.i)
            if e.j != nodeID:
                adj.append(e.j)
        return adj

    def is_interagent_edge(self, from_vert_id, to_vert_id,
                           bidirectional=False):
        """ Returns a boolean indicating whether the edge between two vertices
        is an interagent edge.

        Args:
            from_vert_id (int): From vertex id.
            to_vert_id (int): To vertex id.
            bidirectional (:obj:`bool`, optional): Should we search
                bidirectionally?

        Returns:
            Returns True if the edge is an interagent edge. Returns False if
            the edge is not, or does not exist.
        """
        if not bidirectional:
            return ((from_vert_id, to_vert_id) in self.interagent_edges)
        else:
            return ((from_vert_id, to_vert_id) in self.interagent_edges or
                    (to_vert_id, from_vert_id) in self.interagent_edges)

    # \brief Removes a node from the STP
    #  \param nodeID the ID of the node to be removed

    def remove_vertex(self, nodeID):
        if nodeID in self.verts:
            del self.verts[nodeID]

            if nodeID in self.received_timepoints:
                self.received_timepoints.remove(nodeID)
            # Clear edges
            toRemove = []
            for i, j in self.edges:
                if i == nodeID or j == nodeID:
                    toRemove += [(i, j)]
            for i, j in toRemove:
                del self.edges[(i, j)]
                if (i, j) in self.contingent_edges:
                    del self.contingent_edges[(i, j)]
                if (i, j) in self.interagent_edges:
                    del self.interagent_edges[(i, j)]
                if (i, j) in self.requirement_edges:
                    del self.requirement_edges[(i, j)]

            # self.tris = [t for t in self.tris
            #             if t.i != nodeID and t.j != nodeID and t.k != nodeID]

    ##
    # \fn get_vertex
    # \brief Gets a node from the STP
    #
    # @param nodeID Integer representing the gloabl ID of the node to get.
    # @return Returns a Vertex with specified global ID. Returns None if it
    # does not exist.
    def get_vertex(self, nodeID):
        if nodeID in self.verts:
            return self.verts[nodeID]
        else:
            return None

    # -------------------------------------------------------------------------
    # Edge functions #
    # -------------------------------------------------------------------------

    ##
    # \fn get_edge
    # \brief Gets an Edge from the STP.
    #
    # \details Direction is not accounted for.
    #
    # @param i The first Node of the edge.
    # @param j The second Node of the edge.
    # @return Returns an Edge object if one exists between i & j. If not,
    #   return None.
    def get_edge(self, i, j):
        if (i, j) in self.edges:
            return self.edges[(i, j)]
        elif (j, i) in self.edges:
            return self.edges[(j, i)]
        else:
            return None

    # \brief Ges all edges of the STP

    def get_all_edges(self):
        return list(self.edges.values())

    def get_edge_weight(self, i, j):
        """Gets the directed edge weight of an edge from the STN

        Direction does matter. If no edge exists between i & j, returns
        infinity.
        If i = j, this function returns 0.

        Args:
            i (int): From node Id
            j (int): To node Id

        Return:
            Returns distance graph weight from i to j as float.
        """
        e = self.get_edge(i, j)
        if e is None:
            if i == j and i in self.verts:
                return 0
            else:
                return float('inf')
        if e.i == i and e.j == j:
            return e.Cij
        else:
            return e.Cji

    ##
    # \fn edgeExists
    # \brief Checks if an edge exists in the STP, regardless of direction.
    #
    # @param i The first Node of the edge.
    # @param j The second Node of the edge.
    # @return Returns a boolean on whether or not an edge exists between the
    #   inputted nodes. Direction is not accounted for.
    def edge_exists(self, i, j):
        return ((i, j) in self.edges) or ((j, i) in self.edges)

    def update_edge(self, i, j, w, equality=False, force=False, create=False):
        """Updates the edge with node ids i & j.

        If the new weight is less than the original weight. Otherwise, do
            nothing. Return whether the updates were successful.

        Args:
            i (int): The starting Node of the edge.
            j (int): The ending Node of the edge.
            w (float): The new weight to try to update with.
            equality (bool, optional): If set to true, we'll "update" the edge
                even if the original weight is the same as w (the new weight).
                Nothing will actually change, but this function will return
                True. Default is false.
            force (bool, optional): If set to true, will ignore previous
                constraints and overwrite them.
            create (bool, optional): If set to true, will create an edge if
                none existed previously.

        Returns:
            Returns boolean whether or not the update actually occured.
        """
        e = self.get_edge(i, j)
        if e is None:
            if not create:
                return False
            self.add_edge(i, j, -float("inf"), w)
            return True
        if e.i == i and e.j == j:
            if w < e.Cij or force:
                e.Cij = w
                return True
            else:
                if equality:
                    return w <= e.Cij
                return False
        else:
            if w < e.Cji or force:
                e.Cji = w
                return True
            else:
                if equality:
                    return w <= e.Cji
                return False

    def get_assigned_time(self, node_id):
        """Retrieves the assigned time of a node.

        Args:
            node_id (int): Node id to check.

        Returns:
            If the node is assigned (has an exact, specific time), return that
                time. Otherwise, return None.
        """
        if node_id == 0:
            # This is the Z-node.
            return 0.0
        if (self.get_edge_weight(0, node_id) !=
                -self.get_edge_weight(node_id, 0)):
            return None
        return self.get_edge_weight(0, node_id)

    # returns True if the time point has been executed or if the previous
    # time point has been executed, meaning that if there is an edge leading
    # to nodeID then the first timepoint of that edge has been executed
    def incoming_executed(self, node_id):
        if self.verts[node_id].is_executed():
            return True

        if node_id in self.received_timepoints:
            ctg_e = self.get_incoming_contingent(node_id)
            return self.verts[ctg_e.i].executed

        ex = [self.verts[e.i].is_executed() for e in self.get_all_edges()
              if e.j == node_id]

        return all(ex)

    def outgoing_executed(self, nodeID):
        if not self.verts[nodeID].executed:
            return False
        ex = [self.verts[e.j].executed for e in self.get_all_edges()
              if e.i == nodeID]
        return all(ex)

    def get_incoming(self, node_id):
        return [e for e in self.get_all_edges() if e.j == node_id]

    def get_outgoing(self, node_id):
        return [e for e in self.get_all_edges() if e.i == node_id]

    def get_incoming_contingent(self, nodeID):
        ctg = [self.contingent_edges[(i, j)] for (i, j)
               in self.contingent_edges if j == nodeID]
        if len(ctg) > 1:
            print('[Error]: {} incoming contingent edges!\n{}'.format(
                len(ctg), ctg))
            raise AssertionError
        if len(ctg) <= 0:
            return None
        return ctg[0]

    def print_executed(self, nodeID):
        nodes = [e.i for e in self.get_all_edges() if e.j ==
                 nodeID and not e.fake]
        for n in nodes:
            for node in self.get_all_verts():
                print("I have: %d" % (node.nodeID))
                if node.nodeID == n:
                    print("%d: %d" % (node.nodeID, node.executed))

    def execute(self, nodeID):
        if nodeID in self.verts:
            self.verts[nodeID].execute()

    def set_makespan(self, makespan):
        self.makespan = makespan
        currentMakespan = 0
        for vert in self.verts:
            if vert != 0:
                val = self.edges[(0, vert)].Cij
                if val > currentMakespan:
                    currentMakespan = val
        for vert in self.verts:
            if vert != 0:
                if self.edges[(0, vert)].Cij == currentMakespan:
                    self.edges[(0, vert)].Cij = makespan

    def for_json(self):
        jsonSTN = {}

        # Add the vertices
        jsonSTN['nodes'] = []
        for v in self.get_all_verts():
            if v.nodeID == 0:
                continue
            json = v.for_json()
            json['min_domain'] = -self.get_edge_weight(v.nodeID, 0)
            json['max_domain'] = self.get_edge_weight(0, v.nodeID)
            jsonSTN['nodes'].append(json)

        # Add the edges
        jsonSTN['constraints'] = []
        for c in self.get_all_edges():
            if c.i == 0:
                continue
            jsonSTN['constraints'].append(c.for_json())

        jsonSTN['num_agents'] = len(self.agents)

        return jsonSTN

    def add_all_edges(self):
        verts = list(range(len(self.verts)))
        for i in verts:
            for j in verts:
                if not self.get_edge(i, j) and i != j:
                    self.add_edge(i, j, -float('inf'), float('inf'))

    # initial attempt (Summer 2017) to fix the cycle problem in running the
    # timeline simulation flips edges and adjusts weights accordingly
    def flip_edges(self):
        new_edges = {}
        for key in self.edges:
            value = self.edges[key]
            max_val = value.Cij
            min_val = -value.Cji
            val1 = math.copysign(1.0, max_val)
            val2 = math.copysign(1.0, min_val)
            if val1 == -1.0 and val2 == -1.0:
                max_val = -1.0 * max_val
                min_val = -1.0 * min_val
                new_edges[(key[1], key[0])] = Edge(
                    value.j,
                    value.i,
                    max_val,
                    min_val,
                    distribution=value.distribution
                )
            else:
                new_edges[key] = self.edges[key]
        self.edges = new_edges

    # \fn FloydWarshall()
    # \brief Runs the Floyd-Warshal algorithm on an STN

    def floyd_warshall(self, create=False):
        verts = self.verts
        B = {}
        for u in self.verts.keys():
            for v in self.verts.keys():
                B[(u, v)] = self.get_edge_weight(u, v)
        for k in verts.keys():
            for i in verts.keys():
                for j in verts.keys():
                    B[(i, j)] = min(B[(i, j)], B[(i, k)] + B[(k, j)])
                    self.update_edge(i, j, B[(i, j)], create=create)

        for e in self.get_all_edges():
            if e.get_weight_min() > e.get_weight_max():
                return False
        return True

    def cap_edges(self):
        """Removes any excessively large edges, and replaces them with a
            very, very large floating point number.
            This floating point number can be found in MAX_FLOAT.
        """
        for k in self.edges.keys():
            self.edges[k].cap()

    # \brief minimizes the stn if possible
    #  and returns true if consistent
    # def minimize(self):
    #    # Triangulate the STN
    #    stnCopy = self.copy()
    #    agentWait = []
    #    # Set up the wait timer for agent load balancing
    #    for a in range(len(self.agents)):
    #        agentWait.append(0)

    #    # Perform the reduction - we want a graph that has existing or
    #    # inferred edges
    #    # between every pair of vertices (i.e. triangulated)
    #    while len(stnCopy.verts) > 1:
    #        for a in range(len(self.agents)):
    #            if agentWait[a] == 0:
    #                # number of triangles created by triagulating using the
    #                # minimum degree vertex
    #                created = stnreduce(stnCopy, self.agents[a], self)
    #                agentWait[a] = created
    #            else:
    #                agentWait[a] -= 1

    #    # Add all triangles initially
    #    updateQueue = Queue()
    #    triangleQueue = Queue()
    #    for a in self.agents:
    #        for tri in self.getAgentTriangles(a):
    #            triangleQueue.put(tri)

    #    # Get PPC
    #    return DItriSTP(self, updateQueue, triangleQueue)

    # \brief minimizes the stn if possible
    #  and returns true if consistent
    # def minimize_alpha(self, alpha):
    #    # Triangulate the STN
    #    stnCopy = self.copy()
    #    agentWait = []
    #    # Set up the wait timer for agent load balancing
    #    for a in range(len(self.agents)):
    #        agentWait.append(0)

    #    # Perform the reduction
    #    while len(stnCopy.verts) > 1:
    #        for a in range(len(self.agents)):
    #            if agentWait[a] == 0:
    #                created = stnreduce(stnCopy, self.agents[a], self)
    #                agentWait[a] = created
    #            else:
    #                agentWait[a] -= 1

    #    # Add all triangles initially
    #    updateQueue = Queue()
    #    triangleQueue = Queue()
    #    for a in self.agents:
    #        for tri in self.getAgentTriangles(a):
    #            triangleQueue.put(tri)

    #    # Get PPC
    #    return DItriSTP_alpha(self, updateQueue, triangleQueue,alpha)

    # This functionality was used in Lund et al. 2017 to expand the amount of
    #   probability mass captured by the algorithm given a risk budget.
    # def minimize_binary_search(self, lb = 0.0, ub = 0.999):
    #    bestSTN = None
    #    alphas = { i : str(round(i/1000.0,3)) for i in range(1001)}
    #    # bounds for binary search
    #    lower = ceil(lb * 1000) - 1
    #    upper = floor(ub * 1000) + 1
    #    # Run binary search on alpha
    #    while upper - lower > 1:
    #        alpha = alphas[(upper + lower) // 2]
    #        tempSTN = self.copy()
    #        # try lower alpha
    #        if tempSTN.minimize_alpha(alpha):
    #            bestSTN = tempSTN.copy()
    #            upper = (upper + lower) // 2
    #        # try higher alpha
    #        else:
    #            lower = (upper + lower) // 2

    #    if bestSTN == None:
    #        return None, None
    #    else:
    #        return bestSTN, alpha

    # An STN is strongly controllable if any assignment of values for
    #   executable timepoints is guaranteed to be consistent with all
    #   constraints (irrespective of contingent edges)
    #   (copied from Lund et al. 2017).

    def isStronglyControllable(self, debug=False, returnSTN=False):
        if not self.possiblyStronglyControllable:
            if returnSTN:
                return False, None
            else:
                return False

        stn_copy = self.copy()
        if not stn_copy.floydWarshall():
            if returnSTN:
                return False, None
            else:
                return False

        new_stn = STN()
        for v in self.get_all_verts():
            if v.nodeID not in self.received_timepoints:
                new_stn.add_vertex(v.nodeID, v.localID, v.ownerID, v.location)

        for e in self.get_all_edges():
            if not e.fake and not e.distribution:
                if debug:
                    print("dealing with edge {}->{}".format(e.i, e.j))
                if e.i in self.received_timepoints:
                    cont_edge = self.get_edge(self.parent[e.i], e.i)
                    i = self.parent[e.i]
                    l_i = cont_edge.get_weight_min()
                    u_i = cont_edge.get_weight_max()
                else:
                    i = e.i
                    l_i = 0
                    u_i = 0

                if e.j in self.received_timepoints:
                    cont_edge = self.get_edge(self.parent[e.j], e.j)
                    j = self.parent[e.j]
                    l_j = cont_edge.get_weight_min()
                    u_j = cont_edge.get_weight_max()
                else:
                    j = e.j
                    l_j = 0
                    u_j = 0

                # This is taken from the 1998 paper by Vidal et al. on
                # controllability
                if debug:
                    st = "Cij: {}  Cji: {}  l_i: {}  u_i: {}  l_j: {}  u_j: {}"
                    print(st.format(e.Cij, e.Cji, l_i, u_i, l_j, u_j))
                lower_bound = -e.Cji + u_i - l_j
                upper_bound = e.Cij + l_i - u_j

                if (i, j) in new_stn.edges:
                    new_stn.update_edge(i, j, upper_bound)
                    new_stn.update_edge(j, i, -lower_bound)
                    if debug:
                        print(
                            "updated edge {}->{}: [{},{}]"
                            .format(i, j, lower_bound, upper_bound)
                        )
                else:
                    new_stn.add_edge(i, j, lower_bound, upper_bound)
                    if debug:
                        print(
                            "added edge {}->{}: [{},{}]"
                            .format(i, j, lower_bound, upper_bound)
                        )

        new_stn.agents = list(self.agents)

        if returnSTN:
            if new_stn.floydWarshall():
                for edge in list(new_stn.edges.values()):
                    stn_copy.update_edge(edge.i, edge.j, edge.Cij)
                    stn_copy.update_edge(edge.j, edge.i, edge.Cji)
                return True, stn_copy
            else:
                return False, None

        return new_stn.floydWarshall()
