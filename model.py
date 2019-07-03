from stn import STN
from stn import loadSTNfromJSONfile
from algorithm import *
import os
import glob

## \file model.py
#  \brief convert STNU to optimization problem in an AMPL format file

MAX_FLOAT = 10000000000

# -------------------------------------------------------------------------
# Strong controllability
# -------------------------------------------------------------------------

##
# \fn prepare(STN)
# \brief extract all constraints, variables, and Objective function for the
#        optimization problem of an STNU
#
# @param STN        An STN to be processed
#
# @return A dictionary of bound variables and a dictionary of epsilon variable
#         in the form of (var, lowbound, upbound). A dictionary of constraints
#         and the objective function in string form.
def prepare(STN):
    bounds = {}
    epsilons = {}
    constraints = {}

    for i in STN.verts:
        if STN.getEdgeWeight(0,i) == float('inf'):
            STN.setMakespan(MAX_FLOAT)
            break

    for i in list(STN.verts.keys()):
        bounds[(i, '+')] = ('T%iHI'%i, 0, STN.getEdgeWeight(0,i))

        low = 0 if STN.getEdgeWeight(i,0) == float('inf') else\
                            -STN.getEdgeWeight(i,0)

        bounds[(i, '-')] = ('T%iLO'%i, low, None)

        if i in STN.uncontrollables:
            epsilons[(i, '-')] = ('EPS%iLO'%i, 0, None)
            epsilons[(i, '+')] = ('EPS%iHI'%i, 0, None)
            constraints[i] = bounds[(i, '-')][0] + " <= " + bounds[(i, '+')][0]
        elif i == 0:
            constraints[(i, '-')] = bounds[(i, '-')][0] + " == 0"
            constraints[(i, '+')] = bounds[(i, '+')][0] + " == 0"
        else:
            constraints[i] = bounds[(i, '-')][0] + " == " + bounds[(i, '+')][0]


    for i,j in STN.edges:
        if (i,j) in STN.contingentEdges:

            constraints[(i, j, '+')] = bounds[(j, '+')][0] + ' - ' + \
                bounds[(i, '+')][0] + ' == ' + str(STN.getEdgeWeight(i,j)) \
                + ' - ' + epsilons[(j, '+')][0]
            constraints[(i, j, '-')] = bounds[(j, '-')][0] + ' - ' + \
                bounds[(i, '-')][0] + ' == ' + str(-STN.getEdgeWeight(j,i))\
                + ' + ' + epsilons[(j, '-')][0]


        else:
            # NOTE: We need to handle the infinite weight edges. Otherwise
            #       the LP would be infeasible
            high = MAX_FLOAT if STN.getEdgeWeight(i,j) == float('inf') \
                                            else STN.getEdgeWeight(i,j)
            low = MAX_FLOAT if STN.getEdgeWeight(j,i) == float('inf') \
                                            else STN.getEdgeWeight(j,i)

            constraints[(i, j, '+')] = bounds[(j, '+')][0] + ' - ' + \
                                        bounds[(i, '-')][0] + ' <= ' + str(high)
            constraints[(i, j, '-')] = bounds[(i, '+')][0] + ' - ' + \
                                        bounds[(j, '-')][0] + ' <= ' + str(low)


    Obj = ''
    for i in STN.uncontrollables:
        Obj += 'log(' + bounds[(i, '+')][0] + ' - ' + bounds[(i, '-')][0] \
                                                                        + ') + '
    Obj = Obj[:-3]

    return bounds, epsilons, constraints, Obj


##
# \fn modelObj(STN, fname)
# \brief convert an STN obj to AMPL format model file to use in the solver
#
# @param STN        An STN object to be processed
# @param fname      The filename for the output model file
#
# @post An AMPL format model file that contains optimization problem for the
#       input STN
def modelObj(STN, fname):
    bounds, epsilons, constraints, Obj = prepare(STN)

    f = open(fname, 'w')
    for v, l, h in list(bounds.values()):
        line = 'var ' + v + ' >= ' + str(l)

        if h != None:
            line += ', <= ' + str(h)

        line += ';'
        f.write(line +'\n')

    for v, l, h in list(epsilons.values()):
        line = 'var ' + v + ' >= ' + str(l)

        if h != None:
            line += ', <= ' + str(h)

        line += ';'
        f.write(line +'\n')

    Obj_line = '\nmaximize VOLUME: ' + Obj + ';\n\n'
    f.write(Obj_line)

    for i in list(STN.verts.keys()):
        line = 'subject to '
        if i == 0:
            upper = 'ZERO_U: ' + constraints[(i, '-')] + ';\n'
            lower = 'ZERO_L: ' + constraints[(i, '+')] + ';\n'
            f.write(line + upper + line + lower)
            del constraints[(i, '-')]
            del constraints[(i, '+')]
        else:
            line += 'LIMIT' + str(i) + ': ' + constraints[i] + ';\n'
            f.write(line)
            del constraints[i]


    for i,j,x in list(constraints.keys()):
        con = constraints[(i, j, x)]
        L = 'L' if x == '-' else 'U'
        line = 'subject to C' + str(i) + '_' + str(j) + L + ': ' + con + ';\n'
        f.write(line)

    f.close()


##
# \fn modelFile(path)
# \brief extract model file from an STN json file
#
# @param path       The path to an STN json file
#
# @post An AMPL model file
def modelFile(path):
    STN = loadSTNfromJSONfile(path)
    p, f = os.path.split(path)
    fname = f[:-5] + '.mod'
    fname = os.path.join('../../../model', fname)
    modelObj(STN, fname)


# -------------------------------------------------------------------------
# Dynamic controllability
# -------------------------------------------------------------------------

##
# \fn prepareDynamic(STN)
# \brief extract all constraints, variables, and Objective function for the
#        optimization problem of an STNU (degree of DC)
#
# @param STN        An STN to be processed
#
# @return A dictionary of epsilons variables and constraint, objective
#         function string
def prepareDynamic(STN):
    epsilons = {}
    result, conflicts, bounds, weight = DC_Checker(STN.copy())

    contingent = bounds['contingent']

    constraint = ''
    Obj = ''
    for i,j in list(contingent.keys()):
        edge, bound = contingent[(i,j)]
        length = edge.Cij + edge.Cji
        epsilons[j] = ('EPS_%i'%j, 0, length)

        constraint += epsilons[j][0] + ' + '
        Obj += 'log(' + str(length) + ' - ' + epsilons[j][0] + ') + '

    constraint = constraint[:-3] + ' >= ' + str(-weight)
    Obj = Obj[:-3]

    return epsilons, constraint, Obj



##
# \fn modelObjDynamic(STN, fname)
# \brief convert an STN obj to AMPL format model file to use in the solver
#        for compute degree of DC
#
# @param STN        An STN object to be processed
# @param fname      The filename for the output model file
#
# @post An AMPL format model file that contains optimization problem for the
#       input STN for computing degree of DC
def modelObjDynamic(STN, fname):
    epsilons, constraint, Obj = prepareDynamic(STN.copy())

    f = open(fname, 'w')
    for v, l, h in list(epsilons.values()):
        line = 'var ' + v + ' >= ' + str(l)

        if h != None:
            line += ', <= ' + str(h)

        line += ';'
        f.write(line +'\n')

    Obj_line = '\nmaximize VOLUME: ' + Obj + ';\n\n'
    f.write(Obj_line)

    constraint_line = 'subject to WEIGHT: ' + constraint + ';\n'
    f.write(constraint_line)

    f.close()



##
# \fn modelDynamicFile(path)
# \brief extract model file from an STN json file
# \note This function is used to compute degree of dynamic controllability
#
# @param path       The path to an STN json file
#
# @post An AMPL model file
def modelDynamicFile(path):
    STN = loadSTNfromJSONfile(path)
    p, f = os.path.split(path)
    fname = f[:-5] + '.mod'
    fname = os.path.join('../../../model/dynamic/model', fname)
    modelObjDynamic(STN, fname)

# -------------------------------------------------------------------------
#  Main function
# -------------------------------------------------------------------------


def main():
    directory = input("Please input the folder containing STN json file:\n")
    listOfFile = glob.glob(os.path.join(directory, '*.json'))

    for path in listOfFile:
        p, f = os.path.split(path)
        print("Processing: ", f)
        modelDynamicFile(path)

if __name__ == '__main__':
    main()
