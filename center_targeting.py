from stn import STN, loadSTNfromJSONfile
from util import STNtoDCSTN, PriorityQueue
from dc_stn import DC_STN
from relax import relaxSearch
import empirical
import random
import json


def makeExpectedSTN(network: STN):
    centerSTN = network.copy()
    for edge in centerSTN.edges:
        newVal = (network.getEdgeWeight(edge.i, edge.j) + network.getEdgeWeight(edge.j, edge.i))/2
        centerSTN.modifyEdge(edge.i, edge.j, newVal)
    return centerSTN
        

def dispatchCenterTargeting(network: STN,
             dc_network: DC_STN,
             realization: dict,
             contingent_map: dict,
             uncontrollable_events,
             verbose=False) -> bool:

    # Dispatch the modified network and assume we have a zero reference point
    enabled = {}
    not_executed = set(network.verts.keys())
    executed = set()
    current_time = 0.0

    centerSTN = makeExpectedSTN(network)
    schedule = {}

    time_windows = {event: [0, float('inf')] for event in not_executed}
    current_event = null

    while len(not_executed) > 0:
        min_time = float('inf')
        # Pick an event to schedule
        for event in enabled:
            lower_bound = time_windows[event][0]
            if event in uncontrollable_events:
                if lower_bound < min_time:
                    min_time = lower_bound
                    current_event = event
            else:
                if lower_bound < min_time:
                    min_time = lower_bound
                    current_event = event

        is_uncontrollable = current_event in uncontrollable_events
        current_time = min_time
        schedule[current_event] = current_time

        # Propagate the constraints
        for nodes, edge in dc_network.normal_edges.items():
            if edge.i == current_event:
                new_upper_bound = edge.weight + current_time
                if new_upper_bound < time_windows[edge.j][1]:
                    time_windows[edge.j][1] = new_upper_bound
            if edge.j == current_event:
                new_lower_bound = current_time - edge.weight
                if new_lower_bound > time_windows[edge.i][0]:
                    time_windows[edge.i][0] = new_lower_bound

        # Add newly enabled events
        for event in not_executed:
            if verbose:
                print("***")
                print("Checking event", event)
            if (event not in enabled) and (event not in uncontrollable_events):
                ready = True
                outgoing_reqs = dc_network.verts[event].outgoing_normal
                # Check required constraints
                for edge in outgoing_reqs:
                    # For required
                    if edge.weight < 0:
                        if edge.j not in executed:
                            if verbose:
                                print(event, "was not enabled because of",
                                      edge)
                            ready = False
                            break

