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
def loadSTNfromJSONfile(filepath, using_PSTN=False):
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
def loadSTNfromJSON(json_str, using_PSTN=False):
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
def loadSTNfromJSONobj(jsonSTN, using_PSTN=False):
    stn = STN()

    # Add the root vertex and put it in the T_x set
    stn.addVertex(0)

    # Add the vertices
    for v in jsonSTN['nodes']:
        stn.addVertex(v['node_id'])

    # Add the edges
    for e in jsonSTN['constraints']:
        if using_PSTN and 'distribution' in e:
            stn.addEdge(e['first_node'], e['second_node'],
                        float(e['min_duration']), float(e['max_duration']),
                        e['type'], e['distribution']['name'])
        else:
            stn.addEdge(e['first_node'], e['second_node'],
                        float(e['min_duration']), float(e['max_duration']),
                        e['type'])

    return stn
