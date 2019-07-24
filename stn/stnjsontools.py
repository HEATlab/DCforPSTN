##
# \file stnjsontools.py
#
# \brief Contains library functions for reading JSON files that represent STNs.
# \note This is legacy from the RobotBrunch project.

import json
from stn import STN


##
# \fn loadSTNfromJSONfile
# \brief Wrapper function for loadSTNfromJSON to allow reading from files.
#
# @param filepath       Path of file to read in
# @param using_PSTN     Flag indicating whether the input STN is a PSTN
#
# @return Returns a STN object loaded from the json object
def loadSTNfromJSONfile(filepath, using_PSTN=True):
    with open(filepath, 'r') as f:
        stn = loadSTNfromJSON(f.read(), using_PSTN=using_PSTN)
    return stn


##
# \fn loadSTNfromJSON
#
# \brief Wrapper for loadSTNfromJSONobj that loads an STN from a JSON string.
#
# @param json_str       String representation of a full JSON object.
# @param using_PSTN     Flag indicating whether the input STN is a PSTN
#
# @return Returns a STN object loaded from the json string
def loadSTNfromJSON(json_str, using_PSTN=True):
    jsonSTN = json.loads(json_str)
    return loadSTNfromJSONobj(jsonSTN, using_PSTN=using_PSTN)


##
# \fn loadSTNfromJSONobj
# \brief Returns a dictionary that gives us an STN from a JSON file, with
#   a few extra details.
#
# @param jsonSTN        json object to use that represents an STN.
# @param using_PSTN     Flag indicating whether the input STN is a PSTN
#
# @return Returns a STN object loaded from the json object
def loadSTNfromJSONobj(jsonSTN, using_PSTN=True):
    stn = STN()

    # Add the root vertex and put it in the T_x set
    stn.addVertex(0)

    # Add the vertices
    for v in jsonSTN['nodes']:
        stn.addVertex(v['node_id'])
        if 'min_domain' in v:       
            stn.addEdge(0, v['node_id'], float(v['min_domain']),
                float(v['max_domain']))
        else:
            if not stn.edgeExists(0, v['node_id']):
                stn.addEdge(0,v['node_id'], float(0), float('inf'))


    # Add the edges
    for e in jsonSTN['constraints']:
        if stn.edgeExists(e['first_node'], e['second_node']):
            stn.updateEdge(e['first_node'], e['second_node'],float(e['max_duration']))
            stn.updateEdge(e['second_node'], e['first_node'],float(e['min_duration']))
        else:
            if using_PSTN and 'distribution' in e:
                stn.addEdge(e['first_node'], e['second_node'],
                            float(max(0,e['min_duration'])), float(e['max_duration']),
                            e['distribution']['type'], e['distribution']['name'])
            elif 'type' in e:
                if e['type'] == 'stcu':
                    dist = "U_"+str(e['min_duration']) + "_" + str(e['max_duration'])
                    stn.addEdge(e['first_node'], e['second_node'],
                        float(max(0,e['min_duration'])), float(e['max_duration']),
                        e['type'], dist)
                else:
                    stn.addEdge(e['first_node'], e['second_node'],
                                float(e['min_duration']), float(e['max_duration']),
                                e['type'])
            else:
                stn.addEdge(e['first_node'], e['second_node'],
                            float(e['min_duration']), float(e['max_duration']))

    return stn
