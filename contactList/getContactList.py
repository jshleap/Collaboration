#!/usr/bin/python
'''
Utility scripts for contacts
Copyright (C) 2012 Christian Blouin, Jose Sergio Hleap, Alex Safatli

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
from utils import PDBnet

class pdbFile:
    def __init__(self, fileName):
        
        self.fileName = fileName
        
        # Read PDB file using PDBnet.
        
        self.data = PDBnet.PDBstructure(fileName)
        

class Contacts:
    def __init__(self, filename, thres=1.5):

        self.thres = thres

        if filename:
            # Process PDB file.
            proteinFile = pdbFile(filename)
            # Get processed PDB file data.
            self.getContacts(proteinFile.data)
    
    @staticmethod        
    def getIndices(listin):
        
        listout = []
        
        # Returns a new list that displays the residue
        # names for output purposes.
        
        for i in listin:
            listout.append((int(i[0].index),int(i[1].index)))
            
        return listout
     
    def getContacts(self, protein):
        
        self.maps = {} # Contains residue classes.
        self.results = {} # Contains residue indices.
        self.listLens = {} # Lengths of all chain contact lists.
        self.avgLen = 0 # Average length of chains.        
        
        # Go through all chains in the protein.
        
        for ch in protein.chains:
            
            print "Getting contact list for chain " + ch
            chain = protein.chains[ch] # A chain in the protein.
            self.maps[ch] = self.checkWithin(chain) 
            self.results[ch] = self.getIndices(self.maps[ch])
            
        # Go through all results and make one list.
        
        adjList = set(self.results[(self.results.keys())[0]])
        print "Forming a final adjacency list..."
        
        for ch in self.results:
            self.listLens[ch] = len(self.results[ch])
            self.avgLen += self.listLens[ch]
            adjList = adjList.union(set(self.results[ch]))
        
        self.avgLen /= len(self.listLens) # Final average chain length.
        self.adjList = adjList # The final adjacency list.
           
    def checkWithin(self, chain):
        
        # Compares within a SINGLE chain (it does not compare to
        # any residues outside that chain).
        
        contactList = []
        
        for r in chain:
            
            # Get a residue in the chain.
            resid1 = chain[r]
            
            # Create a dict with all entries but this one.
            others = self.excludeItem(chain, r)
        
            for a in resid1.atoms:
                
                # Get an atom in the residue.
                atom1 = resid1.atoms[a]
                # If backbone atom, continue.
                if atom1.name in ["C","O","N"]:
                    continue
                
                for s in others:
                    
                    # Compare to another residue.
                    resid2 = others[s]
                    
                    for b in resid2.atoms:
                        
                        atom2 = resid2.atoms[b]
                       
                        # If backbone atom, continue.
                        if atom2.name in ["C","O","N"]:
                            continue
                        
                        # Compare distance.
                        # If within 1.5 Angstroms, go ahead
                        # and check to see if the contact has
                        # already been listed.
                        if (atom1.DistanceTo(atom2) < self.thres):
                            if not (resid2,resid1) in contactList:
                                contactList.append((resid1,resid2))
                            break
                     
        return contactList  

    def excludeItem(self, listin, itemexc):
        
        # Outputs a structure excluding a
        # given item.
        
        listNew = listin.copy()
        del listNew[itemexc]
        return listNew
    
def processFile(filename,thres=1.5):
    # Allows for non-command line access to this file
    # in order to get data.
    return Contacts(filename,thres)
                    
if __name__ == "__main__":
    data = Contacts(sys.argv[1])
    print "------------------------------"
    print "Length of final adjacency list: " + str(len(data.adjList))
    for ch in data.listLens:
        print "Length of list for chain " + ch + ": " + str(data.listLens[ch])
    out = data.adjList
    fout = open('adjList.txt','w')
    for i in out:
        fout.write(str(i) + "\n")
    print "Adjacency list saved to ./adjList.txt"
    fout.close()
        