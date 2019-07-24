import json
from stn import STN, loadSTNfromJSONfile

NO_AGENT = 0


def mit2stn(fp: str) -> list:
    """Convert MIT JSON file to a PSTN via the STN class.
    Args:
        fp: File path for the JSON file.
        add_z (boolean, optional): Add an extra z timepoint to the STN.
        connect_origin (boolean, optional): Connect every event to the z
            timepoint.
        cap (boolean, optional):
    Returns:
        Returns an STN object made from the passed in file.
    """
    with open(fp, "r") as f:
        jo = json.loads(f.read())
    stnlist = []
    for inst in jo["instances"]:
        # inst is an STP dict
        for arr_name, arr in inst.items():
            stn = _make_stn(arr)
            stn.name = arr_name
            stnlist.append(stn)
    return stnlist
    



def _make_stn(arr):
    """Construct STN from MIT edges"""
    # To store name to id mappings.
    name_to_id = {}
    # Number of events encountered so far.
    event_count = 0
    # Create the STN.
    stn = STN()
    stn.agents = [NO_AGENT]

    for mit_edge in arr:
        start_name = mit_edge["start_event_name"]
        end_name = mit_edge["end_event_name"]
        if start_name not in name_to_id:
            name_to_id[start_name] = event_count
            stn.addVertex(event_count)
            event_count += 1
        if end_name not in name_to_id:
            name_to_id[end_name] = event_count
            stn.addVertex(event_count)
            event_count += 1

        dist = _get_dist(mit_edge)
        if dist is None:
            stn.addEdge(name_to_id[start_name],
                         name_to_id[end_name],
                         float(mit_edge["properties"]["lb"] * 1000),
                         float(mit_edge["properties"]["ub"] * 1000))
        else:
            stn.addEdge(name_to_id[start_name],
                         name_to_id[end_name],
                         -float("inf"),
                         float("inf"),
                         distribution=dist)

    return stn

def _get_dist(mit_edge):
    """Extract distribution information from an edge."""
    edge_type = mit_edge["type"]
    dist = None
    if edge_type == "uncontrollable_probabilistic":
        dist_type = mit_edge["properties"]["distribution"]["type"]
        if dist_type == "gaussian":
            mean = mit_edge["properties"]["distribution"]["mean"]
            sd = mit_edge["properties"]["distribution"]["variance"]**0.5
            dist = "N_{}_{}".format(mean, sd)
        elif dist_type == "uniform":
            lb = mit_edge["properties"]["distribution"]["lb"]
            ub = mit_edge["properties"]["distribution"]["ub"]
            dist = "U_{}_{}".format(lb, ub)
        else:
            print("Unaccounted for distribution: {}".format(dist_type))
    elif edge_type == "controllable":
        return None
    elif edge_type == "uncontrollable_bounded":
        lb = mit_edge["properties"]["lb"]
        ub = mit_edge["properties"]["ub"]
        dist = "U_{}_{}".format(lb, ub)
    else:
        raise ValueError("Edge type not found: {}".format(edge_type))
    return dist
