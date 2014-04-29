#!/usr/bin/env python
"""
------------------------------------------------
File Information:

@file = example.py
@author: Yael Elmatad
@email: yael.elmatad@tapad.com
@date: Aug 21, 2013
-------------------------------------------------
Copyright Information:

-------------------------------------------------
Some general description of what this script/module does:

Example useage of modularityClustering.  To run pleae try:
    python example.py karate.csv
    [Automatic assumption, verbose=False (change for useful comments)]
    [Automatic assumption, symmetric=True ie edge 0->1 is connectd to 1->0 with same strength, change if not true]
"""
#================================
#Imports
import sys
sys.path.append('../src')
from modularity_clustering import modularity_cluster
#================================


#================================
#Main

def main():
    numargs = 2#EDIT ME
    if len(sys.argv) != (numargs+1): #EDIT ME!!
        print "please give me",numargs,"argument(s) [CSV file with 3 columns -- Edge data] and an output file root name"
        sys.exit(1)

    inputfile = sys.argv[1]
    mc = modularity_cluster(verbose=True)
    mc.loadEdges(inputfile, ignoreHeader=False)

    mc.findCommunities(stopAtFirstNegativeDeltaQ=False)
   
    rootFN = sys.argv[2]

    jsonFN = rootFN + ".json"
    gnuFN = rootFN + ".tsv"
    MRFN = rootFN + ".MR"

    mc.printClustersJSON(outputFile=jsonFN)
    mc.printgnu(filename=gnuFN)
    mc.printMR(outputFile=MRFN)

if __name__ == '__main__':
    main()

#================================


