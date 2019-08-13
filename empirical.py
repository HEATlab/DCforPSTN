from stn import STN, loadSTNfromJSONfile
from relax import *
from util import *
from dispatch import *
from probability import *
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import glob
import json
import os
import random
import math
import csv

##
# \file empirical.py
# \brief perform empirical analysis on degree of controllability metrics.


# -------------------------------------------------------------------------
# Generating Results
# -------------------------------------------------------------------------

def vary_risk(data_path, risk_levels, sim_num, out_name):
    nested_folders = False
    data_list = glob.glob(os.path.join(data_path, '*.json'))
    if len(data_list) == 0:
        folders = os.listdir(data_path)
        data_list = []
        for folder in folders:
            data = glob.glob(os.path.join(data_path, folder, '*.json'))
            data_list += data

    result = {}

    for data in data_list:
        result[data] = [data]
        for risk in risk_levels:
            print("executing:", data)
            dispatch_ml, times_ml = simulate_file(data, sim_num, False, True, True, risk)
            dispatch_reg, times_reg = simulate_file(data, sim_num, False, True, False, risk)

            durations_ml = [times_ml[2*i + 1] - times_ml[2*i] for i in range(int(len(times_ml)/2))]
            durations_reg = [times_reg[2*i + 1] - times_reg[2*i] for i in range(int(len(times_reg)/2))]
            relax_time = durations_ml[0]
            dc_time = durations_ml[1]
            time_ml = sum(durations_ml[2:])/len(durations_ml[2:])
            time_reg = sum(durations_reg[2:])/len(durations_reg[2:])
            result[data] += [dispatch_ml, dispatch_reg, relax_time, dc_time, time_ml, time_reg]
            with open(out_name, 'w') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(data, dispatch_ml, dispatch_reg, relax_time, dc_time, time_ml, time_reg)

    #return result

    #Save the results
    # with open(out_name, 'w') as f:
    #     json.dump(result, f)
    # with open(out_name, 'w') as csv_file:
    #     writer = csv.writer(csv_file)
    #     for key, value in result.items():
    #         writer.writerow(value)


def check_list(data_list, sim_num, out_name):

    if len(data_list) == 0:
        folders = os.listdir(data_path)
        data_list = []
        for folder in folders:
            data = glob.glob(os.path.join(data_path, folder, '*.json'))
            data_list += data

    result = {}

    for data in data_list:
        print("executing:", data)
        dispatch_ml, times_ml = simulate_file(data, sim_num, False, False, True)

        result[data] = [dispatch_ml]
    
    return result




def compare_ml_reg(data_path, sim_num, out_name, gauss, risk):

    nested_folders = False
    data_list = glob.glob(os.path.join(data_path, '*.json'))
    if len(data_list) == 0:
        folders = os.listdir(data_path)
        data_list = []
        for folder in folders:
            data = glob.glob(os.path.join(data_path, folder, '*.json'))
            data_list += data

    result = {}
    for data in data_list:
        print("executing:", data)
        dispatch_ml, times_ml = simulate_file(data, sim_num, False, gauss, True, risk)
        dispatch_reg, times_reg = simulate_file(data, sim_num, False, gauss, False, risk)

        durations_ml = [times_ml[2*i + 1] - times_ml[2*i] for i in range(int(len(times_ml)/2))]
        durations_reg = [times_reg[2*i + 1] - times_reg[2*i] for i in range(int(len(times_reg)/2))]
        relax_time = durations_ml[0]
        dc_time = durations_ml[1]

        reg_dc = durations_reg[0]

        time_ml = sum(durations_ml[2:])/len(durations_ml[2:])
        time_reg = sum(durations_reg[1:])/len(durations_reg[1:])
        result[data] = [dispatch_ml, dispatch_reg, relax_time, dc_time, time_ml, 0, reg_dc, time_reg]
    
    #return result

    #Save the results
    # with open(out_name, 'w') as f:
    #     json.dump(result, f)
    with open(out_name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in result.items():
            writer.writerow([key, value[0], value[1], value[2], value[3], value[4], value[5], value[6], value[7]])


def generate_DDC_result(data_path, sim_num, out_name, gauss, relaxed):
    data_list = glob.glob(os.path.join(data_path, '*.json'))

    result = {}

    for data in data_list:
        dispatch = simulate_file(data, sim_num, False, gauss, relaxed)
        ddc = prob_of_DC_file(data, gauss)

        path, name = os.path.split(data)
        result[name] = [ddc, dispatch]
    
    return result

    # Save the results
    # with open(out_name, 'w') as f:
    #     json.dump(result, f)



def generate_result_relax(data_path, out_name):
    data_list = glob.glob(os.path.join(data_path, '*.json'))

    result = {}

    for data in data_list:
        STN = loadSTNfromJSONfile(data)
        _, count, _, _ = relaxSearch(STN)

        path, name = os.path.split(data)
        result[name] = count

    # return result
    # Save the results
    with open(out_name, 'w') as f:
        json.dump(result, f)


def compute_corr(values:dict):
    pairs = list(values.items())
    x_vals = [x[1][0] for x in pairs]
    y_vals = [x[1][1] for x in pairs]
    return np.corrcoef(x_vals,y_vals)

def identify_outliers(values:dict, threshold):
    pairs = list(values.items())
    x_vals = [x[1][0] for x in pairs]
    y_vals = [x[1][1] for x in pairs]
    outliers = []
    for i in range(len(x_vals)):
        if abs(x_vals[i] - y_vals[i]) > threshold:
            outliers.append(pairs[i])
    return outliers

def plot_from_dict(values:dict):
    pairs = list(values.items())
    x_vals = [x[1][0] for x in pairs]
    y_vals = [x[1][1] for x in pairs]
    plt.scatter(x_vals,y_vals)
    plt.xlabel("Degree of Controllability")
    plt.ylabel("Success Rate")
    plt.xticks(np.arange(0, 1, step=0.2))
    plt.show()


def sample_from_folder(folder_name, gauss=False, success='default', LP='original'):
    results = {}
    for filename in os.listdir(folder_name):
        print(filename)
        STN = loadSTNfromJSONfile(folder_name + filename)
        degree, success_rate = sample(STN, success, LP, gauss)
        results[filename] = (degree, success_rate)
    return results



##
# \fn scheduleIsValid(network: STN, schedule: dict) -> STN
# \brief Given an STNU and schedule, checks if the schedule is valid or not.
#
# @param network       An input STNU
#
# @param schedule      A dictionary whose keys are vertex IDs and whose values
#                      are the times selected for those events.
#
# @return              True/False if the schedule is valid/invalid.
def scheduleIsValid(network: STN, schedule: dict) -> STN:
    ## Check that the schedule is actually defined on all relevant vertices
    # This number is arbitrary - any sufficiently small, positive constant works
    epsilon = 0.001
    vertices = network.getAllVerts()
    for vertex in vertices:
        vertexID = vertex.nodeID
        try:
            assert vertexID in schedule
        except AssertionError:
            return False

    # Check that the schedule is valid
    edges = network.getAllEdges()
    for edge in edges:
        if edge.type == 'stc' or edge.type == None:
        # Loop through the constraints
            start = edge.i
            fin   = edge.j
            uBound = edge.Cij
            lBound = -edge.Cji

            boundedAbove = (schedule[fin] - schedule[start]) <= uBound + epsilon
            boundedBelow = (schedule[fin] - schedule[start]) >= lBound - epsilon

            # Check if constraint is not satisfied
            if ((not boundedAbove) or (not boundedBelow)):
                return False

    return True




# -------------------------------------------------------------------------
#  Main function
# -------------------------------------------------------------------------

if __name__ == '__main__':
    #plot()
    compare_ml_reg("dataset/rover_data", 200, "result/rover_ml_timing_knuth.csv", False, 0.25)
