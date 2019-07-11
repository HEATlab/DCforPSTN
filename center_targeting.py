from stn import STN, loadSTNfromJSONfile
from util import STNtoDCSTN, PriorityQueue
from dc_stn import DC_STN
from relax import relaxSearch, relaxNextExecutable
import random
import json
import simulation as sim
import collections as ct



def makeExpectedSTN(network: STN):
    centerSTN = network.copy()
    for edge in centerSTN.edges:
        newVal = (network.getEdgeWeight(edge[0], edge[1]) - network.getEdgeWeight(edge[1], edge[0]))/2
        centerSTN.modifyEdge(edge[0], edge[1], newVal)
        centerSTN.modifyEdge(edge[1], edge[0], - newVal)
    return centerSTN
        

def dispatchCenterTargeting(network: STN,
             dc_network: DC_STN,
             realization: dict,
             contingent_map: dict,
             uncontrollable_events,
             verbose=False) -> bool:

    # Dispatch the modified network and assume we have a zero reference point
    enabled = {0}
    not_executed = set(network.verts.keys())
    executed = set()
    current_time = 0.0

    schedule = {}

    print("not executed are", not_executed)
    time_windows = {event: [0, float('inf')] for event in not_executed}
    current_event = None

    centerSTN = makeExpectedSTN(network)

    #pre-processing step
    vertex_edge_map = ct.defaultdict(list)
    for (vert_1, vert_2) in network.edges.keys():
        vertex_edge_map[vert_2].append(network.edges[(vert_1,vert_2)])
        vertex_edge_map[vert_1].append(network.edges[(vert_1,vert_2)])

    while len(not_executed) > 0:
        # Find next event to execute
        min_time = float('inf')

        # add newly enabled events to enabled
        for event in not_executed:
            ready = True
            if event in enabled:
                continue
            #check if enabled
            for edge in vertex_edge_map[event]:
                if edge.i == event and edge.Cij <= 0:
                    if edge.j not in executed:
                        ready = False
                        break
                elif edge.j == event and edge.Cji <= 0:
                    if edge.i not in executed:
                        ready = False
                        break
            if ready:
                enabled.add(event) 
        # choose event to execute
        print(enabled) 
        print(schedule)
        for ready_event in enabled:
            lower_bound = time_windows[ready_event][0]
            if lower_bound < min_time:
                min_time = lower_bound
                current_event = ready_event
            
        # schedule the current event
        is_uncontrollable = current_event in uncontrollable_events

        prev_event = None
        if is_uncontrollable:
            delay = realization[current_event][0]
            prev_event = realization[current_event][1]
        else:
            try:
                range_freedom = float('inf')
                for edge in vertex_edge_map[current_event]:
                    if edge.Cij - edge.Cji < range_freedom:
                        if edge.i == current_event and edge.Cij <= 0:
                            prev_event = edge.j
                            range_freedom = edge.Cij - edge.Cji
                        elif edge.j == current_event and edge.Cji <= 0:
                            prev_event = edge.i
                            range_freedom = edge.Cij - edge.Cji
                delay = relaxNextExecutable(network, centerSTN, current_event, prev_event)
            except KeyError:
                delay = 0

        current_time = current_time + delay
        schedule[current_event] = current_time

        # update the centered STN
        centerSTN.modifyEdge(prev_event,current_event, delay)
        centerSTN.modifyEdge(current_event, prev_event, -delay)


        not_executed.remove(current_event)
        enabled.remove(current_event)
        executed.add(current_event)

    good = empirical.scheduleIsValid(network, schedule)
    return good


























    
    
    
    # enabled = {}
    # not_executed = set(network.verts.keys())
    # executed = set()
    # current_time = 0.0

    # centerSTN = makeExpectedSTN(network)
    # schedule = {}

    # time_windows = {event: [0, float('inf')] for event in not_executed}
    # current_event = None

    # while len(not_executed) > 0:
    #     min_time = float('inf')
    #     # Pick an event to schedule
    #     for event in enabled:
    #         lower_bound = time_windows[event][0]
    #         if event in uncontrollable_events:
    #             if lower_bound < min_time:
    #                 min_time = lower_bound
    #                 current_event = event
    #         else:
    #             if lower_bound < min_time:
    #                 min_time = lower_bound
    #                 current_event = event

    #     is_uncontrollable = current_event in uncontrollable_events
    #     current_time = min_time
    #     schedule[current_event] = current_time

    #     # Propagate the constraints
    #     for nodes, edge in dc_network.normal_edges.items():
    #         if edge.i == current_event:
    #             new_upper_bound = edge.weight + current_time
    #             if new_upper_bound < time_windows[edge.j][1]:
    #                 time_windows[edge.j][1] = new_upper_bound
    #         if edge.j == current_event:
    #             new_lower_bound = current_time - edge.weight
    #             if new_lower_bound > time_windows[edge.i][0]:
    #                 time_windows[edge.i][0] = new_lower_bound

    #     # Add newly enabled events
    #     for event in not_executed:
    #         if verbose:
    #             print("***")
    #             print("Checking event", event)
    #         if (event not in enabled) and (event not in uncontrollable_events):
    #             ready = True
    #             outgoing_reqs = dc_network.verts[event].outgoing_normal
    #             # Check required constraints
    #             for edge in outgoing_reqs:
    #                 # For required
    #                 if edge.weight < 0:
    #                     if edge.j not in executed:
    #                         if verbose:
    #                             print(event, "was not enabled because of",
    #                                   edge)
    #                         ready = False
    #                         break

