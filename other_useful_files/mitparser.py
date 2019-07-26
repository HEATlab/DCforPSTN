import json

from stn import STN


NO_AGENT = 0


def mit2stn(fp: str, add_z=False, connect_origin=False, cap=False) -> list:
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
            verbose = False
            if arr_name == '2_rovers_1':
                verbose = True
            stn = _make_stn(arr, add_z, connect_origin, verbose)
            stn.name = arr_name
            if cap:
                stn.cap_edges()
            stnlist.append(stn)
    return stnlist


def _make_stn(arr, add_z, connect_origin, verbose):
    """Construct STN from MIT edges"""
    # To store name to id mappings.
    name_to_id = {}
    # Number of events encountered so far.
    event_count = 0
    # Create the STN.
    stn = STN()
    stn.agents = [NO_AGENT]
    if add_z:
        stn.add_vertex(0, None)
        event_count += 1

    for mit_edge in arr:
        start_name = mit_edge["start_event_name"]
        end_name = mit_edge["end_event_name"]
        if start_name not in name_to_id:
            name_to_id[start_name] = event_count
            stn.add_vertex(event_count, NO_AGENT)
            event_count += 1
        if end_name not in name_to_id:
            name_to_id[end_name] = event_count
            stn.add_vertex(event_count, NO_AGENT)
            event_count += 1

        dist = _get_dist(mit_edge)
        if dist is None:
            stn.add_edge(name_to_id[start_name],
                         name_to_id[end_name],
                         float(mit_edge["properties"]["lb"]) * 1000,
                         float(mit_edge["properties"]["ub"]) * 1000)
            # print(start_name, end_name)
            # print(mit_edge["properties"]["lb"])
            # print(mit_edge["properties"]["ub"])
            # print(float(mit_edge["properties"]["ub"]) * 1000)
            # print(float(mit_edge["properties"]["lb"]) * 1000)
            # print("_____________")

        elif dist[0] == "U":
            try:
                stn.add_edge(name_to_id[start_name],
                         name_to_id[end_name],
                         float(mit_edge["properties"]["distribution"]["lb"]) * 1000,
                         float(mit_edge["properties"]["distribution"]["ub"]) * 1000,
                         distribution=dist)
            except KeyError:
                stn.add_edge(name_to_id[start_name],
                         name_to_id[end_name],
                         float(mit_edge["properties"]["lb"]) * 1000,
                         float(mit_edge["properties"]["ub"]) * 1000,
                         distribution=dist)

        else:
            stn.add_edge(name_to_id[start_name],
                         name_to_id[end_name],
                         -float("inf"),
                         float("inf"),
                         distribution=dist)

        if connect_origin:
            for i in stn.verts.keys():
                if not stn.edge_exists(0, i) and i != 0:
                    stn.add_edge(0, i, 0, float("inf"))
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



if __name__ == '__main__':
    stnlist = mit2stn("libheat/stntools/car_sharing.json")
    for stn in stnlist:
        outname = 'libheat/stntools/car_data/' + stn.name + ".json"
        with open(outname, 'w') as f:
            json.dump(stn.for_json(), f)