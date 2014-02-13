#!/usr/bin/python
'''
Utility scripts for contacts
Copyright (C) 2012 Alex Safatli, Christian Blouin, Jose Sergio Hleap

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

E-mail: iltafas@gmail.com

'''

import centroidContact
import getContactList
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Generates an atom-based adjacency list using threshold amounts ranging from
# 1 to 3 (Angstroms), with 0.1 increments, and plots the amount of
# resulting contacts vs. this threshold.

# If run from command line: $ python contacts-atomDistance.py pdbFile.pdb

fIn = sys.argv[1]
numContacts = [] # List to hold the number of resulting contacts.
cutoffs = [x*0.1 for x in xrange(10,31)] # Cutoffs ranging from 1 to 3, 0.1 increments.

# Get atom-based adjacency list.
print "\nLoading file: " + fIn
print "Generating atom-based contact list for each threshold. This will take a while.\n"

for x in cutoffs:
    print "\nThreshold = " + str(x) + "\n"
    atomBased = getContactList.processFile(fIn,x)
    REF = atomBased.adjList # Adjacency list.
    numContacts.append(len(set(REF))) # Number of contacts.

# Plot the data.  
plt.plot(cutoffs,numContacts)
plt.title("Threshold vs. Number of Contacts for Atom-Based Method")
plt.ylabel("Number of Contacts")
plt.xlabel("Threshold (Angstroms)")
pp = PdfPages('contact-thres-plot.pdf')
plt.savefig(pp, format='pdf')
pp.close()

# Save raw data to CSV file.
fout = open('contact-threshold.csv','w')
fout.write("Threshold (Angstroms)" + "\t" + "Number of Contacts" + "\n")
for x in xrange(0,len(cutoffs)):
    fout.write(str(cutoffs[x]) + "\t" + str(numContacts[x]) + 
               "\n")
fout.close()