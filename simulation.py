from util import PriorityQueue
from stn import STN, loadSTNfromJSONfile

import empirical as emp
import random

##
# \file simulation.py
# \brief Hosts some cheap implementations of dynamic execution strategies on STNUs
# \note This file assumes that any event with id zero is necessarily the zero timepoint

MAX_FLOAT = 100000000000000000.0
ZERO_ID = 0


##
# \fn early_execution(network, realization)
# \brief Runs an STNU simulation where the agent follows an early execution strategy
#
# @param network       STNU we will run simulation on
# @param realization   Dictionary from uncontrollables to contingent edge values
#
# @return A bool, which is True if and only if the execution is successful
def early_execution(network: STN, realization: dict) -> bool:
    ## Bookkeeping for events
    all_uncontrollables = set(network.uncontrollables)
    unused_events = set(network.verts.keys())
    not_scheduled = PriorityQueue()
    final_schedule = {}

    # Mapping from contingent sources to uncontrollables
    contingent_pairs = network.contingentEdges.keys()
    disabled_uncontrollables = {src: sink for (src, sink) in contingent_pairs}

    # Initialize bounds for simulation - starts off with just controllables
    # and zero time point
    controllable_bounds = find_bounds(network)
    true_weight = {}
    for event in controllable_bounds:
        not_scheduled.push(event, controllable_bounds[0])
        true_weight[event] = controllable_bounds[0]
    not_scheduled.addOrDecKey(ZERO_ID, 0)

    # Run simulation
    old_time = 0
    while len(unused_events) > 0:
        current_time, activated_event = not_scheduled.pop()

        # This check ensures that we popped out a valid time_point
        # A better way to deal with this would be to just figure out a way to
        # increase priorities of elements in a heap
        if activated_event in true_weight:
            if true_weight[activated_event] > current_time:
                continue

        unused_events.remove(activated_event)
        final_schedule[activated_event] = current_time

        assert old_time < current_time, "Chronology violated!"

        if activated_event in disabled_uncontrollables:
            # If this is a contingent source, we add the associated uncontrollable sink
            # to the queue
            uncontrollable = disabled_uncontrollables[activated_event]
            delay = realization[uncontrollable]
            not_scheduled.push(uncontrollable, current_time + delay)

        # Update the bounds for all other timepoints
        # We only care about events being moved later in time
        relevant_edges = network.getEdges(activated_event)
        for edge in relevant_edges:
            if (edge.j == activated_event) and (
                    edge.i not in all_uncontrollables):
                if needs_early_update(edge, activated_event, current_time,
                                      true_weight):
                    lower_bound = current_time - edge.Cij
                    true_weight[edge.i] = lower_bound

        # Keep track of this for next iteration of loop
        old_time = current_time
    # Check if we dispatched succesfully
    return emp.scheduleIsValid(network, final_schedule)


##
# \fn dispatch(network, sample_size)
# \brief Simulates dispatch on the network a sample_size number of times,
#        and then reports the success rate of dispatch
#
# @param network      STNU we schedule with
# @param sample_size  Number of times we schedule on the network
def dispatch(network: STN, sample_size: int) -> float:
    successes = 0
    for sample in range(sample_size):
        if simulate_once(network):
            successes += 1
    success_rate = float(successes / sample_size)
    print(f"Dispatch was succesful {100*success_rate}% of the time.")
    return float(successes / sample_size)


##
# \fn simulate_once(network)
# \brief Generate a realization randomly, and simulate execution of the network for
#        that realization
#
# @param network      STNU we run the simulation on
#
# @return A bool which is True if and only if execution is successful
def simulate_once(network: STN) -> bool:
    # Generate the realization
    realization = {}
    for nodes, edge in network.contingentEdges.items():
        realization[nodes[1]] = random.uniform(-edge.Cji, edge.Cij)
    # Run the simulation
    return early_execution(network, realization)


# -------------------------------------------------------------------------
#  Simulation Helpers
# -------------------------------------------------------------------------
##
# \fn needs_early_update()
# \brief Checks if the endpoint of a particular edge needs to have its
#        planned time updated in early execution
#
# @param edge
# @param fixed_event
# @param fixed_value
# @param planned_times
#
# @return True if and only if we should modify the edge's source vertex
def needs_early_update(edge, fixed_event, fixed_value, planned_times):
    new_time = fixed_value - edge.Cij
    if new_time > planned_times[edge.i]:
        return True


##
# \fn safely_scheduled(network, partial, event)
#
# @param network      An input STNU
# @param partial      A partial schedule, inputed as a dictionary from
#                     event IDs to values
#
# @return True if and only if the assigned value of event is consistent in
#         the partial schedule
def safely_scheduled(network: STN, partial: dict, event) -> bool:
    assert event in partial, "Event not in partial schedule!"
    epsilon = 0.001
    edges = network.getEdges(event)
    for edge in edges:
        if (edge.i in partial) and (edge.j in partial):
            start, end = (edge.i, edge.j)
            lBound, uBound = (-edge.Cji, edge.Cij)

            boundedAbove = (partial[end] - partial[start]) <= uBound + epsilon
            boundedBelow = (partial[end] - partial[start]) >= lBound - epsilon

            if ((not boundedAbove) or (not boundedBelow)):
                #print("Violated constraint", edge)
                return False
    return True


# -------------------------------------------------------------------------
#  Modify Networks
# -------------------------------------------------------------------------
##
# \fn find_bounds(network)
# \brief
#
# @param network      The STNU to compute bounds for early execution.
#
# @return A dictionary from controllable events to the implied lower and
#         upper bounds (lb, ub) relative to the zero time point
def find_bounds(network: STN) -> dict:
    # Add zero timepoint
    if ZERO_ID not in network.verts:
        network.addVertex(ZERO_ID)
    # Add bounds relative to zero timepoint
    adjacent_to_zero = set(network.getAdjacent(ZERO_ID))
    events = network.verts.keys()
    bounds = {}
    for event in events:
        if (event != ZERO_ID) and (event not in adjacent_to_zero):
            return 0
        else:
            return 0

    # To make sure zero timepoint starts first
    bounds[ZERO_ID] = (-1.0, 0.0)
    return bounds


##
# \fn set_dynamic_zeropoint(network)
# \brief
#
# @param
# @param
def set_dynamic_zeropoint(network: STN):
    network = network.copy()
    largish = 1000000.0

    if ZERO_ID not in network.verts:
        network.addVertex(ZERO_ID)

    adjacent_events = set(network.getAdjacent(ZERO_ID))
    for event in network.verts:
        if (event not in adjacent_events) and (event != ZERO_ID):
            network.addEdge(ZERO_ID, event, 0.0, largish)

    return network


##
# \fn def make_graph(network)
# \brief
#
#
#
def make_graph(network: STN):
    network = set_dynamic_zeropoint(network)

    # Create the graph in adjacency list format
    events = network.verts.keys()
    graph = {event: {} for event in events}
    for nodes, edge in network.edges.items():
        i, j = nodes
        graph[i][j] = [edge.Cij]
        graph[j][i] = [edge.Cji]
        if edge.type != 'stc':
            graph[i][j].append(-edge.Cji)
            graph[j][i].append(-edge.Cij)

    return graph


##
# \fn get_weight(graph, event_1, event_2)
def get_weight(graph, event_1, event_2) -> float:
    weights = graph[event_1]
    if event_2 not in weights:
        return MAX_FLOAT
    else:
        if type(weights[event_2]) == list:
            return min(weights[event_2])
        else:
            return weights[event_2]


##
# \fn minimize_stnu(graph)
# \brief Uses Floyd Warshall to find minimum graph
#
#
def minimize_stnu(graph):
    # We'll use the list here to order the events
    events = list(graph.keys())
    num_events = len(events)

    # DP Table -- index as dist_table[k][event_1][event_2]
    dist_table = [{
        event_1: {event_2: 0
                  for event_2 in events if event_2 != event_1}
        for event_1 in events
    } for k in range(num_events + 1)]
    # Initialize
    for event_1 in events:
        for event_2 in events:
            if event_2 != event_1:
                dist_table[0][event_1][event_2] = get_weight(
                    graph, event_1, event_2)

    for k in range(1, num_events + 1):
        for event_1 in events:
            for event_2 in events:
                if event_2 != event_1:
                    old_weights = dist_table[k - 1]
                    new_weights = dist_table[k]
                    new_event = events[k - 1]

                    old_path = old_weights[event_1][event_2]

                    fst_half = get_weight(old_weights, event_1, new_event)
                    snd_half = get_weight(old_weights, new_event, event_2)

                    new_path = fst_half + snd_half

                    new_weights[event_1][event_2] = min(old_path, new_path)

    minimized = dist_table[-1]

    return minimized
