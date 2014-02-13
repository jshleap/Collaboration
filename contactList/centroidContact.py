#!/usr/bin/python
'''
Utility scripts for contacts
Copyright (C) 2012 Christian Blouin, Jose Sergio Hleap
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: cblouin@cs.dal.ca

'''

import sys
import getContactList

class Contacts:
    def __init__(self, filename, thres):

        load = getContactList.pdbFile(filename)
        self.thres = thres # Represents distance by which to compare to to find contacts.
        self.getContacts(load.data)
     
    def getContacts(self, protein):
        
        self.results = {} # Contains residue indices.
        self.listLens = {} # Lengths of all chain contact lists.
        
        # Go through all chains in the protein.
        
        for ch in protein.chains:
            
            print "Getting contact list for chain " + ch
            chain = protein.chains[ch] # A chain in the protein.
            centroids = protein.GetAllCentroid(ch) # Gets centroids for the chain.
            self.results[ch] = self.checkAll(centroids)
            
        # Go through all results and make one list.
        
        adjList = set(self.results[(self.results.keys())[0]])
        print "Forming a final adjacency list..."
        
        for ch in self.results:
            self.listLens[ch] = len(self.results[ch])
            adjList = adjList.union(set(self.results[ch]))
        
        self.adjList = adjList # The final adjacency list.
           
    def checkAll(self, cen):
        
        # Compares within a SINGLE chain (it does not compare to
        # any residues outside that chain).
        
        # Uses centroid list of PDB atom structures 
        # output from PDBnet.
        
        contactList = []
        
        for i in cen:
            
            # Get a centroid for a residue in chain.
            centroid1 = i
            resid1 = centroid1.parent
            
            for s in cen:
                
                if s is i: # Do not want to test against self.
                    continue
                
                # Compare to another residue.
                centroid2 = s
                resid2 = centroid2.parent
                
                avg = centroid1.DistanceTo(centroid2) # Calculate distance.
                
                if (avg < self.thres):
                    if not (int(resid2.index),int(resid1.index)) in contactList:
                        # Appends only indices to save on processing time.
                        contactList.append((int(resid1.index),int(resid2.index))) 
    
        return contactList  
    
def processFile(filename, thres):
    # Allows for non-command line access to this file
    # in order to get data.
    data = Contacts(filename,thres)
    return data.adjList
                    
if __name__ == "__main__":
    data = Contacts(sys.argv[2],float(sys.argv[1]))
    print "------------------------------"
    print "Length of final adjacency list: " + str(len(data.adjList))
    for ch in data.listLens:
        print "Length of list for chain " + ch + ": " + str(data.listLens[ch])
    out = data.adjList
    fout = open('adjList-centroids.txt','w')
    for i in out:
        fout.write(str(i) + "\n")
    print "Adjacency list saved to ./adjList-centroids.txt"
    fout.close()
        