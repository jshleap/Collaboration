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

# Generates an adjacency list using the all-atom method
# (getContactList.py) and using the centroid method (centroidContact.py).
# Symmetric difference is then found using the following:
#
# REF <-- the set of contacts using atom-based method, as a list of tuples.
# Mc <-- the adjacency list for a given cutoff.
#
# symmetric difference = 0.5 * (len(set(REF).difference(set(Mc))) 
# + len(set(Mc).difference(set(REF))))

# If run from command line: $ python contacts-symmDiff.py pdbFile.pdb

fIn = sys.argv[1]
aLists = [] # List to hold adjacency lists for cutoffs ranging from 3 to 20, 0.5 increments.
symmDiffs = [] # List to hold the symmetric differences for these cutoffs.
cutoffs = [x*0.5 for x in xrange(6,41)] # Cutoffs ranging from 3 to 20, 0.5 increments.

# Get atom-based adjacency list.
print "\nLoading file: " + fIn
print "Will first generate atom-based contact list. This will take up to a few minutes.\n"
atomBased = getContactList.processFile(fIn)
REF = atomBased.adjList # Adjacency list.
avgLen = atomBased.avgLen # Average number of amino acids.
print "\nAverage Number of Amino Acids: " + str(avgLen)

# Get centroid-based adjacency lists. Calculate symm difference.
print "\nNow, will generate centroid-based adjacency list. This will take a little while.\n"
for x in cutoffs:
    print "\nCutoff = " + str(x) + "\n"
    Mc = centroidContact.processFile(fIn,x)
    aLists.append(Mc)
    symmDiff = 0.5 * (len(set(REF).difference(set(Mc)))+
                      len(set(Mc).difference(set(REF))))
    symmDiff /= avgLen # Scale symmetric difference by dividing by average number of amino acids.
    symmDiffs.append(symmDiff)

# Plot the data.  
plt.plot(cutoffs,symmDiffs)
plt.title("Symmetric Differences for Centroid vs. Atom-based Contact Methods")
plt.ylabel("Symmetric Difference / Avg Number of Amino Acids")
plt.xlabel("Cutoff (Angstroms)")
pp = PdfPages('symmDiff-plot.pdf')
plt.savefig(pp, format='pdf')
pp.close()

# Save raw data to text file.
fout = open('symmDiffs.txt','w')
fout.write("Average Number of Amino Acids: " + str(avgLen) + "\n")
fout.write("Data for Plot (Symmetric Difference, Cutoff):\n")
t = 0
for x in symmDiffs:
    fout.write(str(x) + "\t" + str(cutoffs[t]) + "\n")
    t += 1
fout.close()