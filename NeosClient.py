#!/usr/bin/env python

# Copyright (c) 2017 NEOS-Server
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Python source code - Python XML-RPC client for NEOS Server
"""
# Added
import re

import argparse
import os
import sys
import time
import glob
import json
try:
    import xmlrpc.client as xmlrpclib
except ImportError:
    import xmlrpclib

##
# \file NeosClient.py
# \brief NEOS client, adapted from https://neos-server.org/neos/downloads.html,
#        that submits optimization problems to NEOS and reports the values
#        achieved by the optimal points.


##
# \fn getObjValue(xml_name, outfolder, username, user_password, output)
# \brief Submit job to neos server and let it report the final result
#
# \details This function is modified from the original NeoClient.py file so
#          that it can report the final objective value of the job submitted
#          and allow us to process large amount of file at once. If you submit
#          a job with username and password, you can view the job status in
#          neos server account
#
# @param xml_name           The filename/path to the xml file we want to process
# @param username           The username of your account to the neos server
# @param user_password      The password of your account to the neos server
# @param output             Flag indicating whether we want to report
#                           final result
#
# @return if output, return the objective value from the final result
def getObjValue(xml_name, outfolder, username=None, \
                                    user_password=None, output=False):
    neos = xmlrpclib.ServerProxy("https://neos-server.org:3333")

    alive = neos.ping()
    if alive != "NeosServer is alive\n":
        sys.stderr.write("Could not make connection to NEOS Server\n")
        sys.exit(1)

    if xml_name == "queue":
        msg = neos.printQueue()
        sys.stdout.write(msg)
    else:
        xml = ""
        try:
            xmlfile = open(xml_name, "r")
            buffer = 1
            while buffer:
                buffer = xmlfile.read()
                xml += buffer
            xmlfile.close()
        except IOError as e:
            sys.stderr.write("I/O error(%d): %s\n" % (e.errno, e.strerror))
            sys.exit(1)
        if username and user_password:
            (jobNumber, password) = neos.authenticatedSubmitJob(
                xml, username, user_password)
        else:
            (jobNumber, password) = neos.submitJob(xml)

        start = time.time()

        p, f = os.path.split(xml_name)
        fname = os.path.join(outfolder, f[:-4] + '.txt')
        file = open(fname, 'w')
        file.write(
            "Job number = %d\nJob password = %s\n" % (jobNumber, password))

        sys.stdout.write(
            "Job number = %d\nJob password = %s\n" % (jobNumber, password))
        if jobNumber == 0:
            sys.stderr.write("NEOS Server error: %s\n" % password)
            sys.exit(1)
        else:
            offset = 0
            status = ""
            while status != "Done":
                end = time.time()
                if end - start >= 70:
                    file.write('\nJob still processing...\n')
                    file.close()
                    print("Job hasn't finished processing yet...")
                    return 'waiting'

                time.sleep(1)
                status = neos.getJobStatus(jobNumber, password)

            msg = neos.getFinalResults(jobNumber, password)

            decoded_msg = msg.data.decode()
            objective_pattern = re.compile("\nObjective ([\+0-9.e-]+)\n")

            if output:
                out = objective_pattern.findall(decoded_msg)
                print(out)
                if len(out) == 0:
                    file.write("\nPossibly unbounded...\n")
                    file.close()
                    return None

                file.write("Objective Value: {}".format(float(out[0])))
                file.close()
                return float(out[0])
            else:
                sys.stdout.write(decoded_msg)


##
# \fn main()
# \brief I/O for submitting specific sequence of jobs.
def main():
    xml_folder = input("Please input directory for xml files:\n")
    xml_L = glob.glob(os.path.join(xml_folder, '*.xml'))
    username = input("Please input Neos username:\n")
    password = input("Please input Neos password:\n")
    outfolder = input("Please enter output folder for job file:\n")

    result = {}
    result['normal'] = {}
    result['unbounded'] = []
    result['waiting'] = []
    for fname in xml_L:
        p, f = os.path.split(fname)
        print("Processing: ", f)
        obj = getObjValue(fname, outfolder, username, password, True)

        print('Timer started. ')
        time.sleep(60)
        print('One minute passed...\n')

        if obj == 'waiting':
            result['waiting'].append(f)
        elif not obj:
            result['unbounded'].append(f)
        else:
            # Ignore the file's extension, and save its results.
            f = f[:-4]
            result['normal'][f] = obj

    with open('result.json', 'w') as f:
        json.dump(result, f)


if __name__ == '__main__':
    main()
