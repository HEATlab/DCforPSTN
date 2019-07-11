from stn import STN, loadSTNfromJSONfile
from util import STNtoDCSTN, PriorityQueue
from dc_stn import DC_STN
from relax import relaxSearch, relaxNextExecutable
import empirical
import random
import json
import collections as ct

# For faster checking in safely_scheduled
import simulation as sim

##
# \file dispatch.py
# \brief Hosts method to dispatch STNUs that were modified by the
#        old dynamic checking algorithm
# \note More detailed explanation about the dispatch algorithm can be found in:
#       https://pdfs.semanticscholar.org/0313/af826f45d090a63fd5d787c92321666115c8.pd

ZERO_ID = 0


##
# \fn simulate_and_save(file_names, size, out_name)
# \brief Keep track of dispatch results on networks
def simulate_and_save(file_names: list, size: int, out_name: str):
    rates = {}
    # Loop through files and record the dispatch success rates and
    # approximated probabilities
    for name in file_names:
        success_rate = simulate_file(name, size)
        rates[name] = success_rate

    # Save the results
    with open(out_name, 'w') as out_json:
        out_json.dump(rates)
    print("Results saved to", out_name)


##
# \fn simulate_file(file_name, size)
# \brief Record dispatch result for single file
def simulate_file(file_name, size, verbose=False, gauss=False, relaxed=False, ct=False) -> float:
    network = loadSTNfromJSONfile(file_name)
    result = simulation(network, size, verbose, gauss, relaxed, ct)
    if verbose:
        print(f"{file_name} worked {100*result}% of the time.")
    return result


##
# \fn simulation(network, size)
def simulation(network: STN, size: int, verbose=False, gauss=False, relaxed=False, ct=False) -> float:
    # Collect useful data from the original network
    contingent_pairs = network.contingentEdges.keys()
    contingents = {src: sink for (src, sink) in contingent_pairs}
    uncontrollables = set(contingents.values())

    if relaxed:
        dispatching_network, count, cycles, weights = relaxSearch(network.copy())
    else:
        dispatching_network = network

    total_victories = 0
    dc_network = STNtoDCSTN(dispatching_network)
    dc_network.addVertex(ZERO_ID)

    controllability = dc_network.is_DC()
    if verbose:
        print("Finished checking DC...")

    # Detect if the network has an inconsistency in a fixed edge
    verts = dc_network.verts.keys()
    for vert in verts:
        if (vert, vert) in dc_network.edges:
            if verbose:
                print("Checking", vert)
            edge = dc_network.edges[vert, vert][0]
            if edge.weight < 0:
                dc_network.edges[(vert, vert)].remove(edge)
                dc_network.verts[vert].outgoing_normal.remove(edge)
                dc_network.verts[vert].incoming_normal.remove(edge)
                del dc_network.normal_edges[(vert, vert)]

    # Run the simulation
    for j in range(size):
        realization = generate_realization(network, gauss)
        copy = dc_network.copy()
        if ct:
            result = dispatchCenterTargeting(dispatching_network, copy, realization, contingents,
                          uncontrollables, verbose)
        else:
            result = dispatch(dispatching_network, copy, realization, contingents,
                          uncontrollables, verbose)
        if verbose:
            print("Completed a simulation.")
        if result:
            total_victories += 1

    goodie = float(total_victories / size)
    if verbose:
        print(f"Worked {100*goodie}% of the time.")

    return goodie


##
# \fn dispatch(network, dc_network, realization, contingent_map,
#           uncontrollable_events, verbose)
# \brief Run an early-first scheduling algorithm on a network
#
# @param network                The original STNU we are scheduling on
# @param dc_network             The modified STNU with inferred constraints
# @param realization            An assignment of values for contingent edges
# @param contingent_map         A dictionary for contingent edges
# @param uncontrollable_events  A collection of uncontrollables
# @param verbose                Prints extra statements when set to True
#
# @post A flag which is True precisely when dispatch is succeeds
def dispatch(network: STN,
             dc_network: DC_STN,
             realization: dict,
             contingent_map: dict,
             uncontrollable_events,
             verbose=False) -> bool:

    # Dispatch the modified network and assume we have a zero reference point
    enabled = {ZERO_ID}
    not_executed = set(dc_network.verts.keys())
    executed = set()
    current_time = 0.0

    schedule = {}

    time_windows = {event: [0, float('inf')] for event in not_executed}
    current_event = ZERO_ID
    if verbose:
        print("Beginning dispatch...")
    while len(not_executed) > 0:
        # Find next event to execute
        min_time = float('inf')

        if verbose:
            print("\n\nNetwork looks like: ")
            print(dc_network)

            print("Current time windows: ", time_windows)
            print("Currently enabled: ", enabled)
            print("Already executed: ", executed)
            print("Still needs to be executed: ", not_executed)

        # Pick an event to schedule
        for event in enabled:
            lower_bound = time_windows[event][0]
            if event in uncontrollable_events:
                if lower_bound < min_time:
                    min_time = lower_bound
                    current_event = event
            else:
                # Check that the wait constraints on the event are satisfied
                waits = dc_network.verts[event].outgoing_upper
                lower_bound = time_windows[event][0]

                for edge in waits:
                    if edge.parent != event:
                        if (edge.parent not in executed):
                            if edge.j not in executed:
                                continue
                            lower_bound = max(lower_bound,
                                              schedule[edge.j] - edge.weight)

                if lower_bound < min_time:
                    min_time = lower_bound
                    current_event = event

        print("current event is", current_event)

        is_uncontrollable = current_event in uncontrollable_events

        if verbose:
            if is_uncontrollable:
                print("This event is uncontrollable!")
        current_time = min_time
        schedule[current_event] = current_time

        # Quicker check for scheduling errors
        if not sim.safely_scheduled(network, schedule, current_event):
            if verbose:
                print("------------------------------------------------------")
                print("Failed -- event", current_event,
                      "violated a constraint.")
                print(
                    f"At this time, we still had {len(not_executed)} "
                    f"out of {len(dc_network.verts)} events left to schedule")
                verbose = False
                print("------------------------------------------------------")
            return False

        # If the executed event was a contingent source
        if current_event in contingent_map:
            uncontrollable = contingent_map[current_event]
            delay = realization[uncontrollable][0]
            set_time = current_time + delay
            enabled.add(uncontrollable)
            time_windows[uncontrollable] = [set_time, set_time]

        if is_uncontrollable:
            # Remove waits
            original_edges = list(dc_network.upper_case_edges.items())
            for nodes, edge in original_edges:
                if edge.parent == current_event:
                    if (current_event != edge.i) and (current_event != edge.j):
                        # Modifying the network
                        dc_network.remove_upper_edge(edge.i, edge.j)

        not_executed.remove(current_event)
        enabled.remove(current_event)
        executed.add(current_event)
        print("not executed are", not_executed)

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

                # Check wait constraints
                outgoing_upper = dc_network.verts[event].outgoing_upper
                for edge in outgoing_upper:
                    if edge.weight < 0:
                        label_wait = (edge.parent not in executed)
                        main_wait = (edge.j not in executed)
                        if label_wait and main_wait:
                            ready = False
                            if verbose:
                                print(event, "was not enabled because of",
                                      edge)
                            break

                if ready:
                    if verbose:
                        print("Looks like we enabled", event)
                    enabled.add(event)

    # The realization should be preserved for src, sink in contingent_map
    if verbose:
        print("\n\nFinal schedule is: ")
        print(schedule)
        print("Network is: ")
        print(network)

    good = empirical.scheduleIsValid(network, schedule)
    if verbose:
        msg = "We're safe!" if good else "We failed!"
        print(msg)
    return good


##
# \fn generate_realization(network)
# \brief Uniformly at random pick values for contingent edges in STNU
def generate_realization(network: STN, gauss=False) -> dict:
    realization = {}
    for nodes, edge in network.contingentEdges.items():
        if gauss:
            mu = (edge.Cji + edge.Cij)/2 
            sd = (edge.Cij - edge.Cji)/4
            realization[nodes[1]] = (random.normalvariate(mu, sd), nodes[0])
        else:
            realization[nodes[1]] = (random.uniform(-edge.Cji, edge.Cij), nodes[0])
    return realization


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
            current_time = current_time + delay
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
                current_time = schedule[prev_event] + delay
            except KeyError:
                delay = 0

        schedule[current_event] = current_time

        # update the centered STN
        centerSTN.modifyEdge(prev_event,current_event, delay)
        centerSTN.modifyEdge(current_event, prev_event, -delay)


        not_executed.remove(current_event)
        enabled.remove(current_event)
        executed.add(current_event)

    good = empirical.scheduleIsValid(network, schedule)
    return good
