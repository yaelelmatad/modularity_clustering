#!/usr/bin/env python
"""
------------------------------------------------
File Information:

@file = modularityClustering.py
@author: Yael Elmatad
@email: yael.elmatad@tapad.com
@date: <date>
-------------------------------------------------
Copyright Information:

-------------------------------------------------
Some general description of what this script/module does:

This is a code that takes a closed network and spits out the best clusters.

It is based on an algorithm devloped by M.E.J. Newman and described in Phys Rev E 69, 066133 (2004)

It clusters based on maximizing modularity.

"""
#================================
#Imports
import sys
import math
from copy import deepcopy as dcopy

from operator import itemgetter

import json
import csv
import pprint

from collections import defaultdict

from sets import Set

#================================
#Classes
class modCluster(object):
    '''modality clustering for simple weighted matrices'''

    def __init__(self, verbose=True):
        self.verbose = verbose                  #if you want it to give helpful comments
        self.nNodes = 0                         #number of nodes in graph
        self.nodes = Set()                      #node names
        self.edgeData = {}                      #edge data (weights divided by total weight)
        self.comparableCommunities = set()         #communities which still have outgoing edges
        self.deltaQs = {}

        self.currCommunityData = {}             #current grouping data
        self.currQValue = 0                     #current value of Q, modularity
        self.maxQValueCommData = {}             #current found break up that maximizes Q
        self.maxQValue = 0                      #current max value of Q
        self.header = ""                        #if hte CSV has a header, it will be stored here
        self.isDone = False                     #checks if the algorithm is done doing it's thing
        self.originalData = {}                  #the original, unreweighted data
        self.symmetric=False                    #is the data 'symmetric' ie the weight of edge i->j same as j->i
        self.totalEdges = 0                     #sum of all edges and their weights.
        self.pp = pprint.PrettyPrinter(indent=4) #no one likes to look at gross dictionaries, make them look pretty
        if self.verbose:
            print "All done initiating, please call modCluster.loadEdges(filename(CSV!), ignoreHeader=BOOL, symmetric=BOOL) now to load a CSV"
        self.currMembership = {}                #holds current round membership data
        self.membership={} #                    #holds max("true") membership data
        self.QPath = []                         #holds the history of Q values
        self.tiny = 0.000000000000000000000000000001
        self.stopAtNegativeDeltaQ = False

    def updateMaxMembership(self):
        for comm in self.maxQValueCommData:
            for node in self.maxQValueCommData[comm]["members"]:
                self.membership[node]=comm

    def updateCurrMembershipAll(self):
        '''currMembership is an inversion of community data for fast checking of which community a member belongs to'''
        for comm in self.currCommunityData:
            for node in self.currCommunityData[comm]["members"]:
                self.currMembership[node]=comm

    def updateCurrMembership(self,node,comm):
        '''handles single node updates on currMembership'''
        self.currMembership[node]=comm

    def getMembership(self,node):
        '''returns the index of the community to which node belongs to -- for current max!'''
        return self.membership[node]

    def areInSameComm(self,node1,node2):
        '''tests if two nodes are in the same community -- for current max!'''
        if self.membership[node1] == self.membership[node2]:
            return True
        else:
            return False

    def loadEdges(self,fn,ignoreHeader = False, symmetric=False):
        '''takes a filename and reads it in'''
        '''file format: CSV: node1, node2, weight'''
        '''send ignoreHeader=True if want to ignore first line)'''
        '''send symmetric=True if the edge strenghts are not unidirectional'''
        self.symmetric=symmetric
        loadFile=open(fn,'rb')
        inputCSV=csv.reader(loadFile, skipinitialspace=True, dialect = "excel")
        if (ignoreHeader):
            self.header = inputCSV.next()

        for row in inputCSV:
            if len(row)==3:
                n1 = row[0]
                n2 = row[1]
                weight = float(row[2]) 
                self.edgeData[(n1,n2)]=weight
                self.originalData[(n1,n2)]=weight
                self.edgeData[(n2,n1)]=weight
                self.originalData[(n2,n1)]=weight #let's keep the original data around, mmkay?
                self.nodes.add(n1)
                self.nodes.add(n2)

            else:
                #data should be 3 columns, otherwise, please reformat and try again!
                print 'Something has gone wrong, row of wrong length'
                print 'Length =',len(row)
                print 'Row text = ', row
                sys.exit(1)

        loadFile.close()
        #make my life easy and make sure all keys used exist.
        factor = 0.5
        for n1 in self.nodes:
            for n2 in self.nodes:
                if n1 != n2:
                    if self.edgeData.has_key((n1,n2)):
                        self.totalEdges += factor*self.edgeData[(n1,n2)]

        for n1,n2 in self.edgeData:
            print n1,n2,self.edgeData[(n1,n2)]


        for key in self.edgeData: #reweight by total value
            self.edgeData[key] = self.edgeData[key]/self.totalEdges

        self.nNodes = len(self.nodes)
        #now we set up the current communities with nNodes number of communities
        self.setUpCommunity()
        if self.verbose:
            print "All done setting up, now call findCommunities(stopAtFirstNegativeDeltaQ=BOOL)"

    def printClustersJSON(self,outputFile="communities.json"):
        '''prints community data in pretty json format'''
        jsonFile = open(outputFile,'wb')
        json.dump(self.maxQValueCommData,jsonFile,indent=4)
        jsonFile.close()

    def removeLargerClusters(self, minClusterSize = 4):
        '''gets rid of small clusters'''
        if self.verbose:
            print "Removing clusters smaller than", minClusterSize
        clustersRemoved = 0
        for comm in self.maxQValueCommData.keys():
            if len(self.maxQValueCommData[comm]["members"]) < minClusterSize:
                del self.maxQValueCommData[comm]
                clustersRemoved +=1
        if self.verbose:
            print "Removed", clustersRemoved, "clusters"
            print len(self.maxQValueCommData), "clusters remaining"
        return

    def printMR(self, outputFile="MR_output.dat"):
        oF = open(outputFile,"wb")
        for comm in self.maxQValueCommData:
            for member in self.maxQValueCommData[comm]["members"]:
                toPrint = "\""+member+"\" -> " + str(comm) + ",\n"
                oF.write(toPrint)

        oF.close()

        

    def findJoinAndUpdateQ(self):
        '''this method goes through one pass of algorithm and combines two communities together'''
        #this is the bread and butter, it passes through the community, first checking if it's done or not, and finds the join that maxes DeltaQ and then does the join

        if len(self.currCommunityData)==1 or len(self.comparableCommunities)==0:
            #"All done!"
            self.isDone = True
            self.updateMaxMembership()
            if self.verbose:
                print "DONE: Max Q value found:", self.maxQValue            
                print "If you want to print the output please call: printClustersJSON(outputFile=\"communities.json\")"
                print "For convenient printing with gnuplot please call: printgnu(filename=\"gnuOut.dat\") (plot w/ 'u 3:4:5')"
                print "For using with network analysis (MR) call: printMR(outputFile=\"MR_output.dat\")"
                print "This is the Q path starting from the first Q computed (all separated)"
                for elem in self.QPath:
                    if self.maxQValue > elem - self.tiny and self.maxQValue < elem + self.tiny:
                        print elem, "MAX"
                    else:
                        print elem
            return

        deltaQ, (i,j) = self.findNextPair()
        self.joinNextPair(i,j)
        self.currQValue+=deltaQ
        #self.computeQ()
        if self.currQValue > self.maxQValue:
            self.maxQValue = self.currQValue
            self.maxQValueCommData = dcopy(self.currCommunityData)
        if deltaQ < 0 and self.stopAtNegativeDeltaQ:
            self.isDone = True
            self.updateMaxMembership()
            if self.verbose:
                print "DONE: Max Q value found:", self.maxQValue            
                print "If you want to print the output please call: printClustersJSON(outputFile=\"communities.json\")"
                print "For convenient printing with gnuplot please call: printgnu(filename=\"gnuOut.dat\") (plot w/ 'u 3:4:5')"
                print "For using with network analysis (MR) call: printMR(outputFile=\"MR_output.dat\")"
                print "This is the Q path starting from the first Q computed (all separated)"
                for elem in self.QPath:
                    if self.maxQValue > elem - self.tiny and self.maxQValue < elem + self.tiny:
                        print elem, "MAX"
                    else:
                        print elem
            return



    def findCommunities(self,stopAtFirstNegativeDeltaQ = True):
        self.stopAtNegativeDeltaQ = stopAtFirstNegativeDeltaQ
        '''call this routine to actually do the loop'''
        i = 0
        while not self.isDone:
            i+=1
            if i%50 == 0 and self.verbose:
                print "Clustering Pass", i
            self.findJoinAndUpdateQ()
            self.QPath.append(self.currQValue)


    def joinNextPair(self,i,j):
        '''this routine joins two community clusters i and j'''
        #add the members of j to i

        if i > j:
            print "I SHOULD NEVER BE BIGGER THAN J!"

        for member in self.currCommunityData[j]["members"]:
            self.updateCurrMembership(member,i) #HERE: Check this
        
        self.currCommunityData[i]["members"].extend(self.currCommunityData[j]["members"])

        #combine values of eii
        self.currCommunityData[i]["e"][i]+=self.currCommunityData[j]["e"][j] + self.currCommunityData[j]["e"][i] + self.currCommunityData[i]["e"][j]

        #fix values of eik for curr comm, combine values of eij for other comm, fix values of a
        self.currCommunityData[i]["a"]+=self.currCommunityData[j]["a"]

        keysToUpdate = set()
        for comm in self.comparableCommunities: #HERE: Is a good place for imporovemnt.  When we combine memberships maybe we only need to do it for nonzero values?

            if comm != i and comm !=j:
                if self.currCommunityData[i]["e"].has_key(comm) or self.currCommunityData[j]["e"].has_key(comm) or self.currCommunityData[comm]["e"].has_key(j) or self.currCommunityData[comm]["e"].has_key(i):
                    if comm < i:
                        keysToUpdate.add((comm,i))
                    else:
                        keysToUpdate.add((i,comm))
                    self.currCommunityData[comm]["e"][i]+=self.currCommunityData[comm]["e"][j]
                    self.currCommunityData[i]["e"][comm]+=self.currCommunityData[j]["e"][comm]
                    #remove all refs to j
                    del self.currCommunityData[comm]["e"][j]
                    if comm in self.deltaQs and j in self.deltaQs[comm]:
                        del self.deltaQs[comm][j]

        if (i,j) in keysToUpdate:
            keysToUpdate.remove((i,j))

        #remove cluster j
        del self.currCommunityData[i]["e"][j]
        del self.currCommunityData[j]
        self.comparableCommunities.remove(j)
        print "cluster",j,"has been merged with cluster", i

        removeIFlag = False
        if len(self.currCommunityData[i]["e"]) == 1:
            if i in self.comparableCommunities:
                removeIFlag =True
                self.comparableCommunities.remove(i)
                print "removing from comparison cluster", i

        if not removeIFlag:
            for (c1,c2) in keysToUpdate:
                if c1 < c2:
                    self.deltaQs[c1][c2] = self.currCommunityData[c2]["e"][c1] + self.currCommunityData[c1]["e"][c2]-2*self.currCommunityData[c1]["a"]*self.currCommunityData[c2]["a"]
                elif c2 < c1:
                    print "c2 should never be smaller than c1 once you get here!!", c1, c2
                    self.deltaQs[c2][c1] = self.currCommunityData[c2]["e"][c1] + self.currCommunityData[c1]["e"][c2]-2*self.currCommunityData[c1]["a"]*self.currCommunityData[c2]["a"]
            del self.deltaQs[j]
            del self.deltaQs[i][j]
        else:
            del self.deltaQs[j]
            del self.deltaQs[i]

            for key1 in self.deltaQs:
                if i in self.deltaQs[key1]:
                    del self.deltaQs[key1][i]
                if j in self.deltaQs[key1]:
                    del self.deltaQs[key1][j]









    

    def findNextPair(self):
        '''goes through the data and finds the next pair to join that maximizes deltaQ'''
        if len(self.currCommunityData)==1 or len(self.comparableCommunities)==0:
            #only one community, nothing to do here.
            return

        maxs = {}
        for key1 in self.deltaQs:
            if len(self.deltaQs[key1]) > 0:
                keyMax, valMax = max(self.deltaQs[key1].items(), key=lambda x: x[1])
                maxs[(key1,keyMax)] = valMax


        if len(maxs.items()) == 0:
            print self.comparableCommunities
            print "maxItems:", maxs.items()
            self.printClustersJSON(outputFile="something_wrong_comm.json")
            return

        (key1,key2),deltaQ = max(maxs.items(), key=lambda x: x[1])
        currMax = deltaQ
        return deltaQ, (key1,key2)


        #
        #for c1 in self.comparableCommunities:
        #    for c2 in self.comparableCommunities:
        #        if c1 != c2:
        #            if (c1,c2) not in seenPairs:
        #                seenPairs.add((c1,c2))
        #                seenPairs.add((c2,c1))
        #                #deltaQ =eij+eji-2*aj*ai
        #                self.deltaQs[c1][c2] = deltaQ
        #                self.deltaQs[c2][c1] = deltaQ
        #                if currMax < deltaQ or firstPass:
        #                    currMax = deltaQ
        #                    currJoin = (c1,c2)
        #                    firstPass = False

        #deltaQ = currMax #for transparancy
        #return deltaQ, currJoin
        
    def setUpCommunity(self):
        '''initialize every node to be it's own community island and empty out all values'''
        '''then call computeEs on these island'''

        #Do stuff with es here!
        for i,node in enumerate(self.nodes):
            self.currCommunityData[i]={}
            self.currCommunityData[i]["members"]=[]
            self.currCommunityData[i]["members"].append(node)
            self.currMembership[node]=i
            self.currCommunityData[i]["e"] = defaultdict(float)
            self.currCommunityData[i]["e"][i] =0
            self.currCommunityData[i]["a"] = 0
            self.comparableCommunities.add(i)
        
        self.setUpEs()
        self.setUpQ()


    def computeQ(self):
        #computes value of Q
        self.currQValue = 0
        for c1 in self.currCommunityData:
            self.currQValue+=(self.currCommunityData[c1]["e"][c1]-self.currCommunityData[c1]["a"]**2)


    def setUpQ(self):
        '''for first pass, compute Q'''
        for c1 in self.currCommunityData:
            self.currQValue+=(self.currCommunityData[c1]["e"][c1]-self.currCommunityData[c1]["a"]**2)
        self.maxQValue= self.currQValue
        self.maxQValueCommData = dcopy(self.currCommunityData)
        #self.pp.pprint(self.maxQValueCommData)



    def setUpEs(self):
        '''loops over all edges to set up nNodes number of single component communities to initialize algorithm'''
        for c1 in self.currCommunityData:
            for c2 in self.currCommunityData:
                if c1 != c2:
                    #self.currCommunityData[c1]["e"][c2] = 0
                    for m1 in self.currCommunityData[c1]["members"]:
                        for m2 in self.currCommunityData[c2]["members"]:
                            if m1 == m2: #CAREFUL SHOULD ONLY HAPPEN IF COMM=COMM2!!
                                #do nothing, should never add weight of self with self
                                print "something has gone horribly wrong, something is a member of two communities at once!"
                                print "c1: ", c1, "c2:", c2, "member: ", m1
                            elif (m1,m2) in self.edgeData:
                                self.currCommunityData[c1]["e"][c2]+=self.edgeData[(m1,m2)]*0.5
                                self.currCommunityData[c1]["a"]+=self.currCommunityData[c1]["e"][c2]

        for c1 in self.currCommunityData:
            self.deltaQs.setdefault(c1,{})
            for c2 in self.currCommunityData[c1]["e"].keys():
                if c1 < c2:
                    currDeltaQ = self.currCommunityData[c2]["e"][c1] + self.currCommunityData[c1]["e"][c2]-2*self.currCommunityData[c1]["a"]*self.currCommunityData[c2]["a"]
                    if currDeltaQ > 0:
                        self.deltaQs[c1][c2]=currDeltaQ

        
        #if you want to see all went well
        #self.pp.pprint(self.currCommunityData)


    def printAllMembers(self,filename="allRemainingMembers.dat"):
        '''useful after self.removeLargerClusters to print remaining members in a list'''
        oF = open(filename,"w")
        if self.verbose:
            print "Printing members in clusters"
        for comm in self.maxQValueCommData.keys():
            for member in self.maxQValueCommData[comm]["members"]:
                oF.write(member+"\n")

        oF.close()
        return

    def printgnu(self,filename="gnuOut.dat"):
        '''nice routine to print for use for gnuplot''' 
        oF = open(filename,"w")
        members = []
        for comm in self.maxQValueCommData:
            members.extend(self.maxQValueCommData[comm]["members"])
        oF.write("#Plotting Commands: set pm3d map; plot \"gnuOut.dat\" u 3:4:5 w image\n")
        oF.write("#Node1 Node2 OrderIndex1 OrderIndex2 OriginalWeight\n")

        for i,m1 in enumerate(members):
            for j,m2 in enumerate(members):
                if self.originalData.has_key((m1,m2)):
                    toPrint = m1 + " " + m2 + " " + str(i) + " " + str(j) + " " + str(self.originalData[(m1,m2)]) + "\n"
                else:
                    toPrint = m1 + " " + m2 + " " + str(i) + " " + str(j) + " " + str(0.0) + "\n"
               
                oF.write(toPrint)
            oF.write("\n")
        oF.close()

#END Classes
#================================

#================================
#Main
#================================

def main():
    #dummy main, don't use
    numargs = 0 
    if len(sys.argv) != (numargs+1): 
        print "please give me",numargs,"argument(s)"
        sys.exit(1)


if __name__ == '__main__':
    main()



