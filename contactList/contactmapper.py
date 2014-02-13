#!/usr/bin/python
'''
This Script will check if the residues of a given protein (PDB) are in contact

Contactmapper Copyright (C) 2012 Christian Blouin
Later versions: Copyright (C) 2013 Christian Blouin, Jose Sergio Hleap

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

Requirements:
1) Python:
   a) PDBnet.py in folder utils, availale at https://github.com/AlexSafatli/LabBlouinTools
   #######################################
   #PDBnet.py is a python script and have#
   #to be provided with this code        #
   #######################################

Also, SET THE PATH in .bashrc to the folder containing this script and and the PYTHONPATH to PDBnet!!

To set the PYTHON PATH in UBUNTU:
1. Go to your home directory in the terminal
2. type: nano .bashrc
3. scroll down the document 
4. at the end of the document write:
   export PYTHONPATH=$PYTHONPATH:<path to PDBnet.py>
5. Re-start your terminal
'''

import sys,os
from utils import PDBnet

class pdbFile:
    def __init__(self, fileName):
        
        self.fileName = fileName
        
        # Read PDB file using PDBnet.
        
        self.data = PDBnet.PDBstructure(fileName)
        

class Contacts:
    def __init__(self, prefix, multiple, corrected, thres=4.5):
        self.prefix = prefix
        self.thres = thres
        self.pdbs = []       
        # Load GM index equivalency
        self.datalength = 0
        self.indexmap = self.GetIndexLandmarks(prefix+'.landmarks')
        # Load Protein
        if multiple:
            pdb = self.indexmap.keys()
            chas=[]
            for p in pdb:
                if not corrected:
                    self.pdbs.append(PDBnet.PDBstructure(p+'.pdb'))
                    p = ''.join(e for e in p if e.isalnum())
                    chas.append(p)
                else:
                    self.pdbs.append(PDBnet.PDBstructure(p+'.pdb'))
                    ch=''.join(e for e in p if e.isalnum())[4]
                    chas.append(ch)
            #set an initial file
            self.proteinFile = PDBnet.PDBstructure(pdb[0]+'.pdb')
            self.proteinFile.chains={}
            self.proteinFile.chainsOrder={}
            self.proteinFile.orderofchains=[]
            for pr in range(len(self.pdbs)):
                if not corrected:
                    chain = pdb[pr][-1]
                else:
                    chain=chas[pr]
                self.proteinFile.chains[pdb[pr]] = self.pdbs[pr].chains[chain]
                self.proteinFile.chainsOrder[pdb[pr]] = self.pdbs[pr].chainsOrder[chain]
                self.proteinFile.orderofchains.append(pdb[pr])
        else:
            self.proteinFile = PDBnet.PDBstructure(prefix+'.pdb')
        
        # Do the deed.
        self.contactmap = []
        for i in range(self.datalength):
            self.contactmap.append([])
            
        self.getContacts(self.proteinFile)
    
    def GetIndexLandmarks(self, filename):
        '''
            Load the index equivalencies between chains and the GM landmarks from 
            a landmark file.
        '''
        fin = open(filename)
        out = {}
        
        for line in fin:
            if line[0] == '>':
                # Dictionary of residues indices indexed by GM indices 
                chain = line[1:-1]
                out[chain] = []
                
            elif len(line) > 2:
                out[chain].append(line.split()[1])
            
        self.datalength = len(out[chain])
        return out
        
     
    def getContacts(self, protein):       
        # Go through all chains in the protein.
        
        for ch in protein.chains:
            print "Getting contact list for chain " + ch
            self.checkWithin(ch) 
            
           
    def checkWithin(self, chainname):
        
        # Compares within a SINGLE chain (it does not compare to
        # any residues outside that chain).
        
        contactList = []
        chain = self.proteinFile.chains[chainname]
        
        for a in range(len(self.indexmap[chainname])):
            resA = chain[self.indexmap[chainname][a]]
            for b in range(a+1,len(self.indexmap[chainname])):
                resB = chain[self.indexmap[chainname][b]]
                
                foundmatch = False
                for atomA in resA.atoms:
                    if atomA in ['C', 'N', 'O']:
                        continue
                    atomA = resA.atoms[atomA]
                    
                    for atomB in resB.atoms:
                        if atomB in ['C', 'N', 'O']:
                            continue   
                        atomB = resB.atoms[atomB]

                        if atomA.DistanceTo(atomB) <= self.thres:
                            foundmatch = True
                            
                        if foundmatch:
                            break
                    if foundmatch: 
                        break
                    
                if foundmatch:
                    if not b in self.contactmap[a]:
                        self.contactmap[a].append(b)
                    if not a in self.contactmap[b]:
                        self.contactmap[b].append(a)
                        
    def WriteToFile(self, prefix):
        fout = open(prefix+ '.contacts', 'w')
        
        for a in range(len(self.contactmap)):
            for b in self.contactmap[a]:
                fout.write('(%d,%d)\n'%(a,b))
        
        fout.close()
                    
if __name__ == "__main__":
    prefix = sys.argv[1]
    multiple = False
    corrected=False
    for arg in sys.argv[1:]:
        if arg == '-m' or arg == '-multiple':
            multiple = True
            print 'Make sure that all the PDBs are in the current working directory'
        if arg == '-c':
            corrected=True
            print 'Allowing corrected (by Modeller) PDBs'
    if not multiple:
        if os.system('find %s.pdb'%(prefix)) == 256:
            print 'PDB file not in directory'
            sys.exit(-1)
    data = Contacts(prefix,multiple,corrected)
    data.WriteToFile(prefix)
        
