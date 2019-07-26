##
# \file stnjsontools.py
#
# \brief Contains library functions for reading json files that represent STNs.
#

import json
import os.path
from .stn import STN

##
# \fn loadSTNfromjson
#
# \brief Wrapper for loadSTNfromjsonobj that loads an STN from a json string.
#
# \deprecated Notice, this function was quite horribly implemented the first
#   time around. A correction was made that breaks the entirely inconsistent
#   interface that was originally there.
#
# @param json_str   String representation of a full json object.
# @param reduction  Not sure what this does.
# @return           Returns the dictionary of the form...
#                   {'stn' : STN, 'agent_count' : numAgents}


def load_stn_from_json(json_str, using_pstn=True):
    jsonstn = json.loads(json_str)
    return load_stn_from_json_obj(jsonstn,
                                  using_pstn=using_pstn)


##
# \fn loadSTNfromjsonfile
# \brief Wrapper function for loadSTNfromjson to allow reading from files.
#
# @param filepath Path of file to read in
# @param reduction Triangulate the STN, I guess?
# @return Returns a dictionary that has the format
#   {'stn' : STN, 'agent_count' : numAgents}
def load_stn_from_json_file(filepath, using_pstn=True, from_millis=False):
    with open(filepath, 'r') as f:
        output_dict = load_stn_from_json(f.read(),
                                         using_pstn=using_pstn)
    output_dict["stn"].name = os.path.basename(filepath)
    return output_dict


##
# \fn loadstnfromjsonobj
# \brief Returns a dictionary that gives us an STN from a json file, with
#   a few extra details.
#
# @param jsonstn json object to use that represents an STN.
# @param reduction Triangulate the STN, I guess?
# @return Returns a dictionary that has the format
#   {'stn' : STN, 'agent_count' : numAgents}
def load_stn_from_json_obj(jsonstn, using_pstn=True,
                           from_millis=False):
    stn = STN()

    # Add the root vertex and put it in the T_x set
    stn.add_vertex(0, None, None)
    agents = []

    # specify conversion factor.
    if from_millis:
        conv_fact = 0.001
    else:
        conv_fact = 1

    # Add the vertices
    for v in jsonstn['nodes']:
        # Accumulate a list of all the owners to retrieve the agents.
        if not v['owner_id'] in agents:
            agents.append(v['owner_id'])

        # We don't necessarily need a location, just set it to None if we
        # don't.
        if not ('location' in v):
            v['location'] = None

        stn.add_vertex(v['node_id'], v['owner_id'],
                       v['location'])

        stn.add_edge(0, v['node_id'], float(v['min_domain']) * conv_fact,
                     float(v['max_domain']) * conv_fact)
        if 'executed' in v:
            if v['executed']:
                stn.execute(v['node_id'])

    # Add the edges
    for e in jsonstn['constraints']:
        if 'distribution' in e and using_pstn:
            stn.add_edge(e['first_node'],
                         e['second_node'],
                         float(e['min_duration']) * conv_fact,
                         float(e['max_duration']) * conv_fact,
                         e['distribution']['name'])
        else:
            stn.add_edge(e['first_node'],
                         e['second_node'],
                         float(e['min_duration']) * conv_fact,
                         float(e['max_duration']) * conv_fact)

    stn.agents = agents

    # if reduction:
    #    # Triangulate the STN
    #    stnCopy = stn.copy()
    #    agentWait = []
    #    # Set up the wait timer for agent load balancing
    #    for a in stn.agents:
    #        agentWait.append(0)

    #    # Perform the reduction
    #    while len(stnCopy.verts) > 1:
    #        for i, a in enumerate(stn.agents):
    #            if agentWait[i] == 0:
    #                created = stnreduce(stnCopy, a, stn)
    #                agentWait[i] = created
    #            else:
    #                agentWait[i] -= 1

    # Return back an dictionary for easy labelling of return types.
    # This is done due to legacy support, but it's not too bad here.
    output_dict = {'stn': stn}

    # Ideally, this would return a class, however, this is already butchering
    # quite a bit of legacy support, and a class would end up being even more
    # damaging.
    return output_dict
