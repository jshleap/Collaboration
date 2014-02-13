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
# Data for three plots are then found as follows:
#
# True Positive (TP): Number of contacts at a given threshold (also found with atom method).
# False Positive (FP): Number of contacts at a given threshold (not found in atom method).
# False Negative (FN): Number of contacts from atom method not predicted at a given threshold.
#
# specificity Sp = TP / (TP+FP)
# sensitivity Sv = TP / (TP+FN)
# F score = 2 * (Sp*Sv)/(Sp+Sv)

# If run from command line: $ python contacts-classification.py pdbFile.pdb

fIn = sys.argv[1]
TPs = [] # List to hold True Positives.
FPs = [] # List to hold False Positives.
FNs = [] # List to hold False Negatives.
specificities = [] # List to hold the specificities for these cutoffs.
sensitivities = [] # List to hold the sensitivities for these cutoffs.
fScores = [] # List to hold the F Scores for these cutoffs.
cutoffs = [x*0.5 for x in xrange(6,41)] # Cutoffs ranging from 3 to 20, 0.5 increments.

# Get atom-based adjacency list.
print "\nLoading file: " + fIn
print "Will first generate atom-based contact list. This will take up to a few minutes.\n"
atomBased = getContactList.processFile(fIn)
REF = atomBased.adjList # Adjacency list.

# Get centroid-based adjacency lists. Calculate appropriately.
print "\nNow, will generate centroid-based adjacency lists. This will take a little while.\n"
for x in cutoffs:
    print "\nCutoff = " + str(x) + "\n"
    c = centroidContact.processFile(fIn,x)
    
    TP = len(set(REF).intersection(set(c)))
    FP = len(set(c).difference(set(REF)))
    FN = len(set(REF).difference(set(c)))
    
    TPs.append(TP)
    FPs.append(FP)
    FNs.append(FN)
    
    Sp = float(TP)/(TP+FP)
    Sv = float(TP)/(TP+FN)
    
    specificities.append(Sp)
    sensitivities.append(Sv)
    
    # Avoid division by zero.
    fScore = 0 if ((Sp+Sv) == 0) else (2.0*((Sp*Sv)/(Sp+Sv)))
    
    fScores.append(fScore)

# Plot the data.  
plt.plot(cutoffs,specificities)
plt.title("Specificities for Contact Determination Methods")
plt.ylabel("Specificity")
plt.xlabel("Cutoff (Angstroms)")
pp = PdfPages('contact-Sp-plot.pdf')
plt.savefig(pp, format='pdf')
pp.close()
plt.clf()
plt.plot(cutoffs,sensitivities)
plt.title("Sensitivities for Contact Determination Methods")
plt.ylabel("Sensitivity")
plt.xlabel("Cutoff (Angstroms)")
pp = PdfPages('contact-Sv-plot.pdf')
plt.savefig(pp, format='pdf')
plt.clf()
pp.close()
plt.plot(cutoffs,fScores)
plt.title("F Scores for Contact Determination Methods")
plt.ylabel("F Score")
plt.xlabel("Cutoff (Angstroms)")
pp = PdfPages('contact-Fscore-plot.pdf')
plt.savefig(pp, format='pdf')
pp.close()

# Save raw data to CSV file.
fout = open('classifications.csv','w')
fout.write("Cutoff (Angstroms)" + "\t" + "Specificity" + "\t" 
           + "Sensitivity" + "\t" + "F Score" +
           "\t" + "TP" + "\t" + "FP" + "\t" + "FN" + "\n")
for x in xrange(0,len(cutoffs)):
    fout.write(str(cutoffs[x]) + "\t" + str(specificities[x]) + 
               "\t" + str(sensitivities[x]) + "\t" +  str(fScores[x]) 
               + "\t" +  str(TPs[x]) + "\t" +  str(FPs[x]) 
               + "\t" +  str(FNs[x]) + "\n")
fout.close()