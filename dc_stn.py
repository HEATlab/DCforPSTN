## \file dc_stn.py
#  \brief An implementation of an STNU with dynamic controllability checking
#  \note Read the Morris 2003 paper on dynamic controllability for more info on
#    the algorithm used
# (https://pdfs.semanticscholar.org/6fb1/b64231b23924bd9e2a1c49f3b282a99e2f15.pdf).

# This file consists almost entirely of legacy code, inherited from the earlier
# project RobotBrunch (in particular, documentation is sparse). The only
# additions are in the "New functions" section at the bottom of the file. The
# main purpose of this file was to originally check dynamic controllability -
# a quicker check is provided in algorithm.py, which we prefer.

## \fn enum(*sequential, **named)
#  \brief Implements enumerations
#  \details Found on the internet at:
# http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
def enum(*sequential, **named):
    enums = dict(list(zip(sequential, list(range(len(sequential))))), **named)
    return type('Enum', (), enums)

edgeTypes = ["NORMAL","LOWER","UPPER"]
edgeType = enum(*edgeTypes)


class DC_Vertex:
    def __init__(self,nodeID):
        # nodeIDs should be unique to vertices.
        self.nodeID = nodeID
        self.outgoing_normal = []
        self.incoming_normal = []
        self.outgoing_upper = []
        self.incoming_upper = []
        self.outgoing_lower = []
        self.incoming_lower = []

    def __repr__(self):
        return f"Node {self.nodeID}"

    def isSpecial(self):
        return len(self.incoming_upper) != 0 or len(self.incoming_lower) != 0

    ## Equality only tests for nodeID.
    def __eq__(self, other):
        same_name = self.nodeID == other.nodeID
        return same_name

    def copy(self):
        vertex = DC_Vertex(self.nodeID)

        for edge in self.outgoing_normal:
            vertex.outgoing_normal.append(edge.copy())
        for edge in self.incoming_normal:
            vertex.incoming_normal.append(edge.copy())
        for edge in self.outgoing_upper:
            vertex.outgoing_upper.append(edge.copy())
        for edge in self.incoming_upper:
            vertex.incoming_upper.append(edge.copy())
        for edge in self.outgoing_lower:
            vertex.outgoing_lower.append(edge.copy())
        for edge in self.incoming_lower:
            vertex.incoming_lower.append(edge.copy())

        return vertex

## Directed, weighted edges.
class DC_Edge:
    def __init__(self,i,j,weight,type=edgeType.NORMAL,parent=None,fake=False):
        self.i = i
        self.j = j
        self.weight = weight
        self.type = type
        self.parent = parent
        self.fake = fake

    def __repr__(self):
        s = ""
        return "Edge {} => {}, {}, Type: {}, Label: {}".format(self.i, self.j,
                            self.weight, self.type, self.parent)

    def copy(self):
        edge = DC_Edge(self.i, self.j, self.weight,
                self.type, self.parent, self.fake)
        return edge

    def __eq__(self, other):
        src   = self.i == other.i
        snk   = self.j == other.j
        wght  = self.weight == other.weight
        types = self.type == other.type
        prnt  = self.parent == self.parent
        fake  = self.fake == other.fake

        return src and snk and wght and types and prnt and fake

## \class DC_STN
#  \brief an implementation of an STN for determining dynamic controllability.
class DC_STN(object):
    def __init__(self):
        self.verts = {}

        self.edges = {}
        self.normal_edges = {}
        self.upper_case_edges = {}
        self.lower_case_edges = {}

    def copy(self):
        new_dc_stn = DC_STN()
        for k, v in self.verts.items():
            new_dc_stn.verts[k] = v.copy()
        for k, v in self.edges.items():
            new_dc_stn.edges[k] = [edge.copy() for edge in v]
        for k, v in self.normal_edges.items():
            new_dc_stn.normal_edges[k] = v.copy()
        for k, v in self.upper_case_edges.items():
            new_dc_stn.upper_case_edges[k] = v.copy()
        for k, v in self.lower_case_edges.items():
            new_dc_stn.lower_case_edges[k] = v.copy()
        return new_dc_stn

    def __repr__(self):
        s = ""
        vert_ids = sorted(self.verts.keys())
        edge_keys = [(i,j) for i in vert_ids for j in vert_ids if i<j]
        for (i,j) in edge_keys:
            edge_list1 = self.edges.get( (i,j), [])
            if i != j:
                edge_list2 = self.edges.get( (j,i), [])
            else:
                edge_list2 = []
            if len(edge_list1) == 0 and len(edge_list2) == 0:
                continue

            for_edge,rev_edge,for_up_edge,rev_up_edge,for_low_edge,rev_low_edge = (None,)*6

            for edge in edge_list1:
                if edge.type == edgeType.NORMAL:
                    for_edge = edge
                if edge.type == edgeType.UPPER:
                    for_up_edge = edge
                if edge.type == edgeType.LOWER:
                    for_low_edge = edge
            for edge in edge_list2:
                if edge.type == edgeType.NORMAL:
                    rev_edge = edge
                if edge.type == edgeType.UPPER:
                    rev_up_edge = edge
                if edge.type == edgeType.LOWER:
                    rev_low_edge = edge

            if for_edge is None:
                for_edge = DC_Edge(None,None,float('inf'))
            if rev_edge is None:
                rev_edge = DC_Edge(None,None,float('inf'))

            spec_edges = [for_up_edge,rev_up_edge,for_low_edge,rev_low_edge]
            if (for_edge.fake or rev_edge.fake) and all([e is None for e in spec_edges]):
                continue

            #forward contingent edge
            if (not for_low_edge is None) and (not rev_up_edge is None):
                s += "{} => {} : [{:.0f},{:.0f}] Contingent\n"\
                .format(i,j,-rev_edge.weight,for_edge.weight)
            #reverse contingent edge
            elif (not for_up_edge is None) and (not rev_low_edge is None):
                s += "{} => {} : [{:.0f},{:.0f}] Contingent\n"\
                .format(j,i,-for_edge.weight,rev_edge.weight)
            #forward wait
            elif (not rev_up_edge is None):
                s += "{} => {} : [{:.0f},{:.0f}] <Wait ({}): {:.0f}>\n"\
                .format(i,j,-rev_edge.weight,for_edge.weight,rev_up_edge.parent,
                              -rev_up_edge.weight)
            #reverse wait
            elif (not for_up_edge is None):
                s += "{} => {} : [{:.0f},{:.0f}] <Wait ({}): {:.0f}>\n"\
                .format(j,i,-for_edge.weight,rev_edge.weight,for_up_edge.parent,
                              -for_up_edge.weight)
            #normal edge
            else:
                s += "{} => {} : [{:.0f},{:.0f}]\n"\
                .format(i,j,-rev_edge.weight,for_edge.weight)

        return s

    ## \fn addVertex(self,nodeID)
    #  \brief adds a vertex
    def addVertex(self,nodeID):
        if nodeID not in self.verts:
            self.verts[nodeID] = DC_Vertex(nodeID)


    ## #\fn addEdge(self,i,j,weight,edge_type = edgeType.NORMAL,parent=None,debug=False)
    #  Attempts to add an edge from i to j. If there is already an edge there, we
    #  will attempt to update that edge. Returns a boolean indicating whether or
    #  not the STN is changed
    def addEdge(self,i,j,weight,edge_type = edgeType.NORMAL,parent=None,debug=False,fake=False):
        if edge_type == edgeType.NORMAL and (i,j) in self.normal_edges:
            old_edge = self.normal_edges[(i,j)]
            if weight < old_edge.weight:
                old_edge.weight = weight
                if debug:
                    print("\nUpdated edge {}--->{}".format(i,j))
                if (i,j) in self.upper_case_edges:
                    if weight <= self.upper_case_edges[(i,j)].weight:
                        self.remove_upper_edge(i,j)
                return True
            return False
        if edge_type == edgeType.UPPER and (i,j) in self.upper_case_edges:
            old_edge = self.upper_case_edges[(i,j)]
            if parent != old_edge.parent:
                raise ValueError('Cannot have two upper-case edges with different'+\
                                 'parents\nOn edge %d->%d, parents %d and %d'\
                                 %(i,j,old_edge.parent,parent) )
            else:
                if weight < old_edge.weight:
                    old_edge.weight = weight
                    if debug:
                        print("\nUpdated edge {}-U->{}".format(i,j))
                    return True
                return False
        if edge_type == edgeType.LOWER and (i,j) in self.lower_case_edges:
            raise ValueError('Cannot update lower-case edges.\nOn edge %d->%d'%(i,j) )
        else: # if there was no edge present in the STN, add one.

            #don't add upper case edges corresponding to waits that will never be used
            if edge_type == edgeType.UPPER and (i,j) in self.normal_edges:
                if weight >= self.normal_edges[(i,j)].weight:
                    return False

            if i not in self.verts:
                self.addVertex(i)
            if j not in self.verts:
                self.addVertex(j)

            newEdge = DC_Edge(i,j,weight,edge_type,parent,fake)
            try:
                self.edges[(i,j)] += [newEdge]
            except KeyError:
                self.edges[(i,j)] = [newEdge]

            start_vert = self.verts[i]
            end_vert = self.verts[j]

            if newEdge.type == edgeType.NORMAL:
                self.normal_edges[(i,j)] = newEdge
                start_vert.outgoing_normal.append(newEdge)
                end_vert.incoming_normal.append(newEdge)
                if debug:
                    print("\nAdded edge {}--->{}".format(i,j))
                if (i,j) in self.upper_case_edges:
                    if weight <= self.upper_case_edges[(i,j)].weight:
                        self.remove_upper_edge(i,j)
            elif newEdge.type == edgeType.UPPER:
                self.upper_case_edges[(i,j)] = newEdge
                start_vert.outgoing_upper.append(newEdge)
                end_vert.incoming_upper.append(newEdge)
                if debug:
                    print("\nAdded edge {}-U->{}".format(i,j))
            elif newEdge.type == edgeType.LOWER:
                self.lower_case_edges[(i,j)] = newEdge
                start_vert.outgoing_lower.append(newEdge)
                end_vert.incoming_lower.append(newEdge)
                if debug:
                    print("\nAdded edge {}-l->{}".format(i,j))

            return True

    def remove_upper_edge(self,i,j):
        if (i,j) not in self.upper_case_edges:
            return False
        else:
            edge = self.upper_case_edges[(i,j)]
            self.verts[i].outgoing_upper.remove(edge)
            self.verts[j].incoming_upper.remove(edge)
            self.edges[(i,j)].remove(edge)
            del self.upper_case_edges[(i,j)]


    def is_DC(self,debug_flag=False):
        i = 0
        while i < len(self.edges):
            if not self.all_max_consistent():
                return False

            num_reductions = 0
            num_reductions += self.no_case_reductions(debug_flag)
            num_reductions += self.upper_case_reductions(debug_flag)
            num_reductions += self.cross_case_reductions(debug_flag)
            num_reductions += self.lower_case_reductions(debug_flag)
            num_reductions += self.label_removal_reductions(debug_flag)
            if num_reductions == 0:
                return True
            i += 1
        return False

    ## \fn all_max_consistent(self)
    #  \brief determines if the all-max projection of this STN is consistent
    #  \note does not change the input STN
    def all_max_consistent(self):
        verts = list(self.verts.keys())
        num_verts = len(verts)
        B = [ [self.max_projection(i,j) for j in verts] for i in verts ]
        #check for starting negative cycles
        for i in range(num_verts):
            if B[i][i] < 0:
                return False
        #run Floyd Warshal
        for k in range(num_verts):
            for i in range(num_verts):
                for j in range(num_verts):
                    #update edge i->j using the path i->k->j
                    B[i][j] = min(B[i][j], B[i][k] + B[k][j])
                    #check for a negative self-loop
                    if i == j and B[i][j] < 0:
                        return False
        return True

    ## \fn max_projection(self,i,j)
    #  \brief computes the edge weight  between i and j in the all-max projection
    #    of the STNU
    def max_projection(self,i,j):
        if i == j:
            max_proj = 0
        else:
            max_proj = float('inf')
        for edge in self.edges.get( (i,j), []):
            if edge.type != edgeType.LOWER and edge.weight < max_proj:
                max_proj = edge.weight
        return max_proj

    ## \fn no_case_reductions(self)
    #  \brief performs needed no-case reductions on the input STN, and returns the
    #     number of reductions performed
    def no_case_reductions(self,debug):
        num_reductions = 0
        # we only perform no-case reductions when the resultant edge starts or ends
        # at a special node
        for first_edge in list(self.normal_edges.values()):
            start_vert = self.verts[first_edge.i]
            mid_vert = self.verts[first_edge.j]
            for second_edge in mid_vert.outgoing_normal:
                if second_edge.j == first_edge.i:
                    continue
                end_vert = self.verts[second_edge.j]
                #original algorithm makes sure that this path starts or ends with a
                #special node we don't, because ommitting that constraint results in a
                #minimal network which is easier to compare to other results
                # if you want this constraint back, just use the following line:
                #if start_vert.isSpecial() or end_vert.isSpecial():
                new_weight = first_edge.weight + second_edge.weight
                if self.addEdge(first_edge.i,second_edge.j,new_weight,debug=debug,fake=True):
                    num_reductions += 1
                    if debug:
                        print(("no-case reduction on    {}--->{}--->{}  "+\
                              "adds {}--->{}, weight {}").format(first_edge.i,
                                first_edge.j,second_edge.j,first_edge.i,second_edge.j,
                                new_weight))
        return num_reductions

    ## \fn no_upper_reductions(self)
    #  \brief performs needed upper-case reductions on the input STN, and returns
    #    the number of reductions performed
    def upper_case_reductions(self,debug):
        num_reductions = 0
        # we only perform upper-case reductions with upper-case edges between two
        # special nodes
        for uc_edge in list(self.upper_case_edges.values()):
            start_vert = self.verts[uc_edge.i]
            end_vert = self.verts[uc_edge.j]
            if start_vert.isSpecial() and end_vert.isSpecial():
                for edge in self.verts[uc_edge.i].incoming_normal:
                    new_weight = uc_edge.weight + edge.weight
                    if self.addEdge(edge.i,uc_edge.j,new_weight,
                                                   edge_type = edgeType.UPPER,
                                                   parent = uc_edge.parent,
                                                   debug = debug):
                        num_reductions += 1
                        if debug:
                            print(("upper-case reduction on {}--->{}-U->{}  "+\
                                  "adds {}-U->{}, weight {}").format(edge.i,uc_edge.i,
                                                                   uc_edge.j,edge.i,uc_edge.j,
                                                                   new_weight))

        return num_reductions

    ## \fn cross_case_reductions(self)
    #  \brief performs all possible cross-case reductions on the input STN, and
    #    returns the number of reductions performed
    def cross_case_reductions(self,debug):
        num_reductions = 0
        for lc_edge in list(self.lower_case_edges.values()):
            for uc_edge in self.verts[lc_edge.j].outgoing_upper:
                if lc_edge.parent != uc_edge.parent and uc_edge.weight < 0:
                    new_weight = lc_edge.weight + uc_edge.weight
                    if self.addEdge(lc_edge.i,uc_edge.j,new_weight,
                                                   edge_type = edgeType.UPPER,
                                                   parent = uc_edge.parent,
                                                   debug=debug):
                        num_reductions += 1
                        if debug:
                            print(("cross-case reduction on {}-l->{}-U->{}  "+\
                                  "adds {}-U->{}, weight {}").format(lc_edge.i,uc_edge.i,
                                                                uc_edge.j,lc_edge.i,uc_edge.j,
                                                                new_weight))
        return num_reductions

    ## \fn lower_case_reductions(self)
    #  \brief performs all possible lower-case reductions on the input STN, and
    #    returns the number of reductions performed
    def lower_case_reductions(self,debug):
        num_reductions = 0
        for lc_edge in list(self.lower_case_edges.values()):
            for edge in self.verts[lc_edge.j].outgoing_normal:
                if edge.weight < 0:
                    new_weight = lc_edge.weight + edge.weight
                    if self.addEdge(lc_edge.i,edge.j,new_weight,
                                                edge_type = edgeType.NORMAL,debug=debug):
                        num_reductions += 1
                        if debug:
                            print(("lower-case reduction on {}-l->{}--->{}  "+\
                                  "adds {}--->{}, weight {}").format(lc_edge.i,edge.i,
                                                                   edge.j,lc_edge.i,edge.j,
                                                                   new_weight))
        return num_reductions

    ## \fn label_removal_reductions(self)
    #  \brief performs all possible label removal reductions on the input STN, and
    #    returns the number of reductions performed
    def label_removal_reductions(self,debug):
        num_reductions = 0
        for uc_edge in list(self.upper_case_edges.values()):
            for lc_edge in self.verts[uc_edge.j].outgoing_lower:
                if lc_edge.parent == uc_edge.parent:
                    if uc_edge.weight  >= -lc_edge.weight:
                        if self.addEdge(uc_edge.i,uc_edge.j,uc_edge.weight,
                                                  edge_type = edgeType.NORMAL,debug=debug):
                            num_reductions += 1
                            if debug:
                                print(("label removal reduction on {}-l->{}-U->{}  "+\
                                      "adds {}--->{}, weight {}").format(lc_edge.i,uc_edge.i,
                                                                  uc_edge.j,uc_edge.i,uc_edge.j,
                                                                  uc_edge.weight))
        return num_reductions


    # -------------------------------------------------------------------------
    # New functions
    # -------------------------------------------------------------------------

    ##
    # \fn incomingEdges(self, start)
    # \brief Get all incoming edges for a given vertex
    #
    # @param start    nodeID for the vertex we want to get the incoming edges
    #
    # @return A list of DC Edges object that are incoming edges for start
    def incomingEdges(self, start):
        vertex = self.verts[start]

        result = []
        result += vertex.incoming_lower
        result += vertex.incoming_upper
        result += vertex.incoming_normal

        return result

    ##
    # \fn getNegNodes(self)
    # \brief Get the set of all vertices with incoming negative edges
    #
    # @return A list of nodeID with incoming negative edges
    def getNegNodes(self):
        negNodes = []
        for v in list(self.verts.keys()):
            incoming = self.incomingEdges(v)

            for edge in incoming:
                if edge.weight < 0:
                    negNodes.append(v)
                    break

        return negNodes
