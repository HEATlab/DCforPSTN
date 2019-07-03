import json
import glob
import os
from dispatch import simulate_file
from probability import prob_of_DC_file
from relax import relaxSearch
from stn import loadSTNfromJSONfile


def generate_DDC_result(data_path, sim_num, out_name):
    data_list = glob.glob(os.path.join(data_path, '*.json'))

    result = {}

    for data in data_list:
        dispatch = simulate_file(data, sim_num)
        ddc = prob_of_DC_file(data)

        path, name = os.path.split(data)
        result[name] = [ddc, dispatch]

    # Save the results
    with open(out_name, 'w') as f:
        json.dump(result, f)

    print("Results saved to", out_name)


def generate_result_relax(data_path, out_name):
    data_list = glob.glob(os.path.join(data_path, '*.json'))

    result = {}

    for data in data_list:
        STN = loadSTNfromJSONfile(data)
        _, count, _, _ = relaxSearch(STN)

        path, name = os.path.split(data)
        result[name] = count

    # Save the results
    with open(out_name, 'w') as f:
        json.dump(result, f)

    print("Results saved to", out_name)


if __name__ == "__main__":
    generate_result_relax('dataset/uncontrollable', 'result_relax.json')
    generate_DDC_result('dataset/uncontrollable', 5000, 'result_dynamic.json')
