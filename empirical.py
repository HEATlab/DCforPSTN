from stn import STN, loadSTNfromJSONfile
from LP import *
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

##
# \file empirical.py
# \brief perform empirical analysis on degree of controllability metrics.


# -------------------------------------------------------------------------
# Generating Results
# -------------------------------------------------------------------------

def generate_DDC_result(data_path, sim_num, out_name, gauss):
    data_list = glob.glob(os.path.join(data_path, '*.json'))

    result = {}

    for data in data_list:
        dispatch = simulate_file(data, sim_num, False, gauss)
        ddc = prob_of_DC_file(data, gauss)

        path, name = os.path.split(data)
        result[name] = [ddc, dispatch]
    
    #return result

    # Save the results
    with open(out_name, 'w') as f:
        json.dump(result, f)



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


# -------------------------------------------------------------------------
# Strong controllability
# -------------------------------------------------------------------------


##
# \fn newInterval(STN, epsilons)
# \brief compute shrinked contingent intervals
#
# @param STN            an input STN
# @param epsilons       a dictionary of epsilons returned by our LP
#
# @return Return a list of original contingent intervals and a list of shrinked
#         contingent intervals
def newInterval(STN, epsilons):
    original = []
    shrinked = []
    for edge in list(STN.contingentEdges.values()):
        orig = (-edge.Cji, edge.Cij)
        original.append(orig)

        low, high = epsilons[(edge.j, '-')].varValue, epsilons[(edge.j, '+')].varValue
        new = (-edge.Cji+low, edge.Cij-high)
        shrinked.append(new)

    return original, shrinked


##
# \fn calculateMetric(original, shrinked)
# \brief Compute our degree of strong controllability
#
# @param original       A list of original contingent intervals
# @param shrinked       A list of shrinked contingent intervals
#
# @return the value of degree of strong controllability
def calculateMetric(original, shrinked, gauss=False):
    if not gauss:
        orig = 1
        new = 1
        for i in range(len(original)):
            x, y = original[i]
            orig *= y-x

            a, b = shrinked[i]
            new *= b-a
        return new, orig, float(new/orig)
    else:
        total = 1
        for i in range(len(original)):
            x, y = original[i]
            mean = (x + y)/2
            sd = (y - mean)/2
            a, b = shrinked[i]
            prob_contained = norm.cdf(b, mean, sd) - norm.cdf(a, mean, sd)
            total *= prob_contained
            
        return total


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
        assert vertexID in schedule

    # Check that the schedule is valid
    edges = network.getAllEdges()
    for edge in edges:
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
#  Sample to get success rate
# -------------------------------------------------------------------------

##
# \fn sampleOnce(original, shrinked)
# \brief Check whether a randomly generated realization is inside strong
#        controllable region
#
# @param original       A list of original contingent intervals
# @param shrinked       A list of shrinked contingent intervals
#
# @return Return True if the random realization falls into the strongly
#         controllable region. Return False otherwise
def sampleOnce(original, shrinked, gauss=False):
    for i in range(len(original)):
        x,y = original[i]
        a,b = shrinked[i]
        if not gauss:
            real = random.uniform(x, y)
        else:
            real = random.normalvariate((x+y)/2, (y-x)/4)
        if real < a or real > b:
            return False

    return True


##
# \fn getSchedule(STN, schedule)
# \brief Construct a possible schedule for an STN given a fixed decision
#
# @param STN            An STNU we want to test
# @param schedule       A dictionary with the fixed decision
#
# @return a schedule for the given STNU
def getSchedule(STN, schedule, gauss):
    for edge in list(STN.contingentEdges.values()):
        start_time = schedule[edge.i]
        x = -edge.Cji
        y = edge.Cij
        if gauss:
            real = random.normalvariate((x+y)/2, (y-x)/4)
        else:
            real = random.uniform(x, y)
        time = start_time + real
        schedule[edge.j] = time
    return schedule


##
# \fn altSampleOnce(STN, schedule)
# \brief Another strategy of sampling to test strong controllability
#
# @param STN            An STNU we want to test
# @param schedule       A dictionary with the fixed decision
#
# @return Return True if the schedule generate is valid. Return False otherwise
def altSampleOnce(STN, schedule, gauss):
    s = getSchedule(STN, schedule, gauss)
    if scheduleIsValid(STN, s):
        return True
    return False


##
# \fn sample(STN, success='default', LP='original')
# \brief Compute the success rate of an STNU by randomly sample 50000 times
#
# \note There are three kinds of LPs we can use to compute the amount of
#       uncertainty removed from each contingent interval
#
# @param STN      An STN to test
# @param LP       The type of LP we want to use
#
# @return The degree of controllability and the success rate for input STN
def sample(STN, success='default', LP='original', gauss=False):
    if LP == 'original':
        _, bounds, epsilons = originalLP(STN.copy(), naiveObj=False)
    elif LP == 'proportion':
        _, _, bounds, epsilons = proportionLP(STN.copy())
    else:
        _, _, bounds, epsilons = maxminLP(STN.copy())

    original, shrinked = newInterval(STN, epsilons)
    if not gauss:
        degree = calculateMetric(original, shrinked)[2]
    else: 
        degree = calculateMetric(original, shrinked, gauss)

    schedule = {}
    for i in list(STN.verts.keys()):
        if i not in STN.uncontrollables:
            time = (bounds[(i, '-')].varValue + bounds[(i, '+')].varValue)/2
            schedule[i] = time

    # Collect the sample data.
    count = 0
    for i in range(10000):
        if success=='default':
            result = sampleOnce(original, shrinked, gauss)
        else:
            result = altSampleOnce(STN, schedule.copy(), gauss)
        if result:
            count += 1

    success = float(count/10000)

    return degree, success


##
# \fn sampleAll(listOfFile, success='default', LP='original')
# \brief Compute the success rate for a list of STNUs
#
# @param STN      An STN to test
# @param LP       The type of LP we want to use
#
# @return a list of (degree, success) tuple for STNUs in the list
def sampleAll(listOfFile, success='default', LP='original'):
    result = {}
    for fname in listOfFile:
        p, f = os.path.split(fname)
        print("Processing file: ", f)
        STN = loadSTNfromJSONfile(fname)
        degree, success = sample(STN, success=success, LP=LP)
        result[f] = (degree, success)


    return result



# ---------------------------------
#  Analyze result from the solver
# ---------------------------------


##
# \fn actual_vol(result_name)
# \brief Calculate the actual volume from solver's result
#
# @param result_name        The path to json file that contains solver's result
#
# @return A dictionary in which keys are the name of the network and values are
#         the maximized volume that guarantees strong controllability
def actual_vol(result_name):
    with open(result_name, 'r') as f:
        result = json.loads(f.read())

    actual_Dict = {}
    for x in list(result.keys()):
        v = result[x]
        actual = math.exp(v)
        actual_Dict[x] = actual

    return actual_Dict


##
# \fn compare(actual_Dict)
# \brief Compute the actual and approximated degree of strong controllability
#
# @param actual_Dict        A dictionary containing actual volume
#
# @return A dictionary in which keys are the name of the network and values are
#         (approximation, actual degree)
def compare(actual_Dict):
    dynamic_folder = input("Please input directory with DC STNUs:\n")
    uncertain_folder = input("Please input directory with uncertain STNUs:\n")

    compare_Dict = {}
    for x in list(actual_Dict.keys()):
        actual_volume = actual_Dict[x]

        if x[:7] == 'dynamic':
            fname = os.path.join(dynamic_folder, x + '.json')
        else:
            fname = os.path.join(uncertain_folder, x + '.json')

        STN = loadSTNfromJSONfile(fname)

        _, _, epsilons = originalLP(STN.copy())
        original, shrinked = newInterval(STN, epsilons)

        old, new, degree = calculateMetric(original, shrinked)
        actual = float(actual_volume / old)
        compare_Dict[x] = (degree, actual)

    return compare_Dict



def plot():
    # Plot actual vs approximated
    result_name = input("Please input path to result json file: \n")
    actual_Dict = actual_vol(result_name)
    compare_Dict = compare(actual_Dict)

    L = list(compare_Dict.values())
    x = [d[0] for d in L]
    y = [d[1] for d in L]

    plt.plot(x,y,'o')
    plt.xlim(-0.04, 1.04)
    plt.ylim(-0.04, 1.04)
    plt.xlabel('Approximated degree of strong controllability')
    plt.ylabel('Actual degree of strong controllability')
    plt.title('Accuracy of Approximation Using DSC LP')

    out_folder = input("Please input the output directory:\n")
    fname = os.path.join(out_folder, 'accuracy.png')
    plt.savefig(fname, format='png')
    plt.close()

    # Plot success rate
    dynamic_folder = input("Please input directory with DC STNUs:\n")
    uncertain_folder = input("Please input directory with uncertain STNUs:\n")

    listOfFile = []
    listOfFile += glob.glob(os.path.join(dynamic_folder, '*.json'))
    listOfFile += glob.glob(os.path.join(uncertain_folder, '*.json'))

    resultD_1= sampleAll(listOfFile, success='new')
    result_1 = list(resultD_1.values())
    x_1 = [d[0] for d in result_1]
    y_1 = [d[1] for d in result_1]
    plt.plot(x_1, y_1, 'o')
    plt.xlim(-0.04, 1.04)
    plt.ylim(-0.04, 1.04)
    plt.xlabel("Approximated degree of strong controllability")
    plt.ylabel("Probabiliy of success")
    plt.title("Success rate of ...")

    out_folder = input("Please input the output directory:\n")
    fname = os.path.join(out_folder, 'success_rate.png')
    plt.savefig(fname, format='png')
    plt.close()

# -------------------------------------------------------------------------
# Dynamic controllability
# -------------------------------------------------------------------------


##
# \fn dynamicMetric(STN, new_STN)
# \brief compute the degree of controllability
#
# @param STN        An input STNU
# @param new_STN    Original STNU with contingent intervals shrinked to be DC
#
# @return degree of dynamic controllability
def dynamicMetric(STN, new_STN):
    original = [(-e.Cji, e.Cij) for e in list(STN.contingentEdges.values())]
    shrinked = [(-e.Cji, e.Cij) for e in list(new_STN.contingentEdges.values())]
    return calculateMetric(original, shrinked)


##
# \fn computeDynamic(nlp=True)
# \brief compute degree of controllability for all uncontrollable STNUs we have
#
# @param nlp        Flag indicating whether we want to use NLP
#
# @return A dictionary in which keys are names of the STNU json file and value
#         is the degree of controllability
def computeDynamic(nlp=False):
    uncertain_folder = input("Please input uncertain STNUs folder:\n")
    chain_folder = input("Please input chain STNUs folde:\n")

    listOfFile = []
    listOfFile += glob.glob(os.path.join(uncertain_folder, '*.json'))
    listOfFile += glob.glob(os.path.join(chain_folder, '*.json'))

    degree = {}
    for fname in listOfFile:
        p, f = os.path.split(fname)
        print("Processing: ", f)

        STN = loadSTNfromJSONfile(fname)
        new_STN, count = relaxSearch(STN.copy(), nlp=nlp)

        if not new_STN:
            degree[f] = 0
        else:
            degree[f] = dynamicMetric(STN.copy(), new_STN.copy())[2]

    return degree



##
# \fn generateData(num)
# \brief generate uncontrollable STNUs with decent degree of dynamic
#        controllability
#
# @param num    number of STNUs we want to generate
def generateData(num):
    data_folder = input("Please input destination directory:\n")
    while num != 0:
        new = generateChain(50, 2500)
        result, conflicts, bounds, weight = DC_Checker(new.copy(), report=False)

        if result:
            print("Failed. Dynamically controllable...")
            continue

        new_STN, count = relaxSearch(new.copy(), nlp=False)
        if not new_STN:
            print("Failed. Not able to resolve conflict...")
            continue

        degree = dynamicMetric(new.copy(), new_STN.copy())[2]
        if degree >= 0.2:
            fname = 'new' + str(num) + '.json'
            print("\nGENERATED ONE SUCCESSFUL CHAIN!!!!!!\n")
            new.toJSON(fname, data_folder)
            num -= 1
        else:
            print("Failed. Degree is too small....")



##
# \fn readNeos(filename, json_folder)
# \brief compute degree of dynamic controllability from the output file of
#        Neos Server
#
# @param filename       The name of the input txt file
# @json_folder          The directory containing STNUs we generated
#
# @return an STNU's volume of Omega and Omega' and the computed degree of DC
def readNeos(filename, json_folder):
    f = open(filename, 'r')
    for i in range(3):
        line = f.readline()

    obj_value = float(line[17:])
    actual = math.exp(obj_value)

    p, f = os.path.split(filename)
    fname = f[:-4] + '.json'
    json_file = os.path.join(json_folder, fname)

    STN = loadSTNfromJSONfile(json_file)
    result, conflicts, bounds, weight = DC_Checker(STN.copy(), report=False)
    contingent = bounds['contingent']

    total = 1
    for (i,j) in list(STN.contingentEdges.keys()):
        edge = STN.contingentEdges[(i,j)]
        length = edge.Cij + edge.Cji
        total *= length

        if (i,j) not in contingent:
            actual *= length

    return actual, total, float(actual/total)


##
# \fn processNeos()
# \brief compute degree of dynamic controllability from the output file of
#        Neos Server for all new chains we generated
#
# @return A dictionary of dictionary with information about an STNU's
#         volume of Omega and Omega', and the computed degree of DC
def processNeos():
    txt_folder = input("Please input folder with txt Neos file:\n")
    json_folder = input("Please input folder with json file:\n")

    result = {}
    text_L = glob.glob(os.path.join(txt_folder, '*.txt'))
    for filename in text_L:
        p, f = os.path.split(filename)
        fname = f[:-4] + '.json'
        print("Processing: ", fname)

        new, orig, degree = readNeos(filename, json_folder)
        result[fname] = {}
        result[fname]['shrinked'] = new
        result[fname]['original'] = orig
        result[fname]['degree'] = degree

    output_folder = input("Please input output folder:\n")
    filename = os.path.join(output_folder, 'result_neos.json')

    with open(filename, 'w') as f:
        json.dump(result, f)

    return result



##
# \fn processOptimal()
# \brief compute degree of dynamic controllability using the optimal solution
#        for all new chains we generated
#
# @return A dictionary of dictionary with information about an STNU's
#         volume of Omega and Omega', and the computed degree of DC
def processOptimal():
    json_folder = input("Please input folder with json file:\n")
    json_list = glob.glob(os.path.join(json_folder, '*.json'))

    result = {}
    for fname in json_list:
        p, f = os.path.split(fname)
        print("Processing: ", f)

        STN = loadSTNfromJSONfile(fname)
        new_STN, count = relaxSearch(STN.copy())
        new, orig, degree = dynamicMetric(STN.copy(), new_STN.copy())

        result[f] = {}
        result[f]['shrinked'] = new
        result[f]['original'] = orig
        result[f]['degree'] = degree

    output_folder = input("Please input output folder:\n")
    filename = os.path.join(output_folder, 'result_optimal.json')

    with open(filename, 'w') as f:
        json.dump(result, f)

    return result



# -------------------------------------------------------------------------
# Generate STNU
# -------------------------------------------------------------------------

##
# \fn generateChain(task, free)
# \brief generate a consistent STNUs in a chainlike structure
#
# \details The chainlike STNU is very common in real life application, such as
#          AUV need to drive to different sites and complete task at each site.
#          Driving to different cite is contingent, but the agent can decide
#          how long it takes to complete the task.
#
# @param task  The number of tasks need to be completed
# @param free  The total length of the free constraint intervals we want
#              in the generated STNU
#
# @return Return the generated STNU
def generateChain(task, free):
    totalEvent = 2 * (task+1)

    while True:
        new = STN()
        for i in range(totalEvent):
            new.addVertex(i)

        L = [random.randint(0, 100) for i in range(task)]
        s = sum(L)
        L = [int(x/s*free) for x in L]
        diff = free - sum(L)
        L[-1] += diff

        bounds = []
        for i in range(totalEvent-1):
            type = 'stcu' if i % 2==0 else 'stc'
            if type == 'stcu':
                lowBound = random.randint(0,50)
                length = random.randint(1,50)
                bounds.append((lowBound, lowBound+length))
                new.addEdge(i, i+1, lowBound, lowBound+length, type='stcu')
            else:
                lowBound = random.randint(0,100)
                length = L[int((i-1)/2)]
                bounds.append((lowBound, lowBound+length))
                new.addEdge(i, i+1, lowBound, lowBound+length)

        low = sum([x[0] for x in bounds])
        high = sum([x[1] for x in bounds])
        S = sum([e.Cij+e.Cji for e in list(new.contingentEdges.values())])
        # makespan = random.randint(int(0.5*low), low)
        makespan = low + int(0.6*S)
        print(low, makespan, high)
        new.addEdge(0,task*2+1, 0, makespan)

        if new.isConsistent():
            return new


def generateParallelChain(agent, task):
    total_event = ((2 * task) + 1) * agent + 1

    while True:
        new = STN()
        new.addVertex(0)

        for i in range(total_event):
            new.addVertex(i+1)

        contingent = True
        for i in range(agent):
            start = ((2 * task) + 1) * i + 1
            end = ((2 * task) + 1) * (i + 1)
            new.addEdge(0, start, 0, 15)

            for j in range(start, end):
                type = 'stcu' if contingent else 'stc'
                contingent = not contingent

                if type == 'stcu':
                    # low = round(random.uniform(10, 20), 2)
                    # high = round(random.uniform(30, 40), 2)
                    low = random.randint(10, 20)
                    high = random.randint(30, 40)
                    new.addEdge(j, j+1, low, high, type='stcu')
                else:
                    # low = round(random.uniform(5, 10), 2)
                    # high = round(random.uniform(30, 35), 2)
                    low = random.randint(5, 10)
                    high = random.randint(35, 40)
                    new.addEdge(j, j+1, low, high)

            new.addEdge(end, total_event, -10, 10)

        num_activity = (2 * task) + 1
        max_length = max([e.Cij + e.Cji for e in list(new.edges.values())])
        up_bound = max_length * num_activity

        # low = round(random.uniform(0.35*up_bound, 0.45*up_bound), 2)
        # high = round(random.uniform(0.5*up_bound, 0.6*up_bound), 2)
        low = random.randint(int(0.45*up_bound), int(0.53*up_bound))
        high = random.randint(int(0.55*up_bound), int(0.65*up_bound))
        new.addEdge(0, total_event, low, high)

        print("\n\nChecking consistensy...")
        if not new.isConsistent():
            continue

        print("Checking Dynamic Controllability...")
        try:
            result, conflicts, bounds, weight = DC_Checker(new.copy(), report=False)
        except Exception:
            continue

        if result:
            return new

# -------------------------------------------------------------------------
#  Main function
# -------------------------------------------------------------------------

if __name__ == '__main__':
    #plot()
    print("ran")
