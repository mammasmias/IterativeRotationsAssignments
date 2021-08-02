#!/usr/bin/env python2.7

import sys
import numpy as np
import hungarian
from collections import Counter 
import operator
import argparse


def kabsch(A, B):
   """
   Kabsch Algorithm as implemented by Jimmy Charnley Kromann

   Calculate RMSD between two XYZ files

   by: Jimmy Charnley Kromann <jimmy@charnley.dk> and 
   Lars Andersen Bratholm <larsbratholm@gmail.com>
   project: https://github.com/charnley/rmsd
   license: https://github.com/charnley/rmsd/blob/master/LICENSE

   A - set of coordinates
   B - set of coordinates

   Performs the kabsch algorithm to calculate the RMSD between A and B

   Returns an RMSD
   """
   A_new = np.array(A)
   A_new = A_new - sum(A_new) / len(A_new)
   A = A_new
   B_new = np.array(B)
   B_new = B_new - sum(B_new) / len(B_new)
   B = B_new

   # Compute covariance matrix
   C = np.dot(np.transpose(A), B)

   # Compute singular value decomposition (SVD)
   V, S, W = np.linalg.svd(C)
   d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

   if d:
      S[-1] = -S[-1]
      V[:, -1] = -V[:, -1]

   # Compute rotation matrix
   U = np.dot(V, W)

   # Rotate A
   A = np.dot(A, U)

   return rmsd(A, B)

def rmsd(V, W):
    """
    V - set of coordinates
    W - set of coordinates

    Returns root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def read_xyz(filename, noHydrogens):
   """
   Reads an xyz file into lists
   filename - name of xyz file
   noHydrogens - if true, hydrogens are ignored; if not, hydrogens are included

   Returns a tuple (a, b)
    where a is a list of coordinate labels and b is a set of coordinates
   (i.e) a = ["O", "H", "H"], b = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]
   """
   xyz = open(filename, "r")
   num_atoms = int(xyz.readline().strip())
   xyz.readline()

   unsorted_labels = []
   unsorted_coords = []
   for i in range(num_atoms):
      line = xyz.readline().strip().split()
      #if noHydrogens and line[0].upper() == "H" : 
      if noHydrogens and line[0].upper().startswith("H") : 
         continue
      else:
         unsorted_labels.append(line[0].upper())
         unsorted_coords.append([float(line[1]), float(line[2]), float(line[3])])
   xyz.close()
   NA = len(unsorted_labels)
   return unsorted_labels, unsorted_coords, NA

def sorted_xyz(filename, noHydrogens):
   """
   Reads an xyz file into lists
   filename - name of xyz file
   noHydrogens - if true, hydrogens are ignored; if not, hydrogens are included

   Returns a tuple (a, b) sorted by atom labels and coordinates
    where a is a list of coordinate labels and b is a set of coordinates
   (i.e) a = ["O", "H", "H"], b = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]
   Sorts the file by atom labels first and coordinates second 
   such that atoms of the same label/type are grouped together
   """
   xyz = open(filename, "r")
   num_atoms = int(xyz.readline().strip())
   xyz.readline()

   sortedlabels = []
   sortedcoords = []
   sortedorder = []
   sortedlines = []
   atomcount = 0
   for i in range(num_atoms):
      line = xyz.readline().strip().split()
      if noHydrogens and line[0].upper().startswith("H") : 
         continue
      else:
         sortedlines.append([line[0].upper(),float(line[1]), float(line[2]), float(line[3]), atomcount])
         atomcount += 1
   xyz.close()

   # sort by element followed by first coordinate (x) then second coordinate (y)
   sortedlines.sort(key=lambda x: (x[0],x[1],x[2]))
   for i in range(atomcount):
      sortedlabels.append(sortedlines[i][0])
      sortedcoords.append([float(sortedlines[i][1]), float(sortedlines[i][2]), float(sortedlines[i][3])])
      sortedorder.append(sortedlines[i][4])
   xyz.close()
   #print sortedlabels
   NA = len(sortedlabels)
   return sortedlabels, sortedcoords, NA, sortedorder

def parse_for_atom(labels, coords, atom):
   """
   labels - a list of coordinate labels
   coords - a set of coordinates
   atom - the atom that is to be parsed

   Returns a set of coordinates corresponding to parsed atom
   """
   atom_coords = []
   for i in range(len(labels)):
      if labels[i] == atom:
           atom_coords.append(coords[i])
   return atom_coords

def transform_coords(coords, swap, reflect):
   """
   coords - a set of coordinates
   swap - the swap transformation (i.e. (0, 2, 1) --> (x, z, y))
   reflect - the reflection transformation (i.e. (1, -1, 1) --> (x, -y, z))

   Returns the transformed coordinates
   """
   new_coords = []
   for i in range(len(coords)):
      new_coords.append([coords[i][swap[0]]*reflect[0], \
                          coords[i][swap[1]]*reflect[1], \
                          coords[i][swap[2]]*reflect[2]])
   return new_coords

def transform_atoms(coords, swap, reflect, atom_indices):
   """
   coords - a set of coordinates
   swap - the swap transformation (i.e. (0, 2, 1) --> (x, z, y))
   reflect - the reflection transformation (i.e. (1, -1, 1) --> (x, -y, z))
   atom_indices - indices of all desired atoms in [coords]

   Returns coordinates after transforming specific atoms
   """
   new_coords = [x[:] for x in coords]
   for i in atom_indices:
      new_coords[i][0] = coords[i][swap[0]]*reflect[0]
      new_coords[i][1] = coords[i][swap[1]]*reflect[1]
      new_coords[i][2] = coords[i][swap[2]]*reflect[2]
   return new_coords

def permute_coords(coords, permutation):
   """
   UNUSED at the moment 

   coords - a set of coordinates
   permutation - permutation of atoms (i.e. [0, 2, 3, 1])

   Returns the permuted coordinates
   """
   new_coords = []
   for i in permutation:
      new_coords.append(coords[i])
   return new_coords

def permute_atoms(coords, permutation, atom_indices):
   """
   coords - a set of coordinates
   permuation - a permutation of atoms
   atom_indices - indices of all desired atoms in [coords]

   Returns the coordinates after permuting just the specified atom
   """
   new_coords = coords[:]
   for i in range(len(permutation)):
      j = atom_indices[permutation[i]]
      k = atom_indices[i]
      new_coords[k] = coords[j] 
   return new_coords

def permute_all_atoms(labels, coords, permutation):
   """
   labels - atom labels 
   coords - a set of coordinates
   permuation - a permutation of atoms

   Returns the permuted labels and coordinates
   """
   new_coords = coords[:]
   new_labels = labels[:]
   for i in range(len(permutation)):
      new_coords[permutation[i]] = coords[i]
      new_labels[permutation[i]] = labels[i]
   return new_labels, new_coords

def get_atom_indices(labels, atom):
   """
   labels - a list of coordinate labels ("Elements")
   atom - the atom whose indices in labels are sought
   Returns a list of all location of [atom] in [labels]
   """
   indices = []
   for i in range(len(labels)):
      if labels[i] == atom:
         indices.append(i)
   return indices

def coords_to_xyz(labels, coords):
   """
   Displays coordinates
   """
   s = ""
   for i in range(len(coords)):
      s += labels[i] + " " + str(coords[i][0]) + " " +  \
             str(coords[i][1]) + " " + str(coords[i][2]) + "\n"
   return s

def write_to_xyz(num_atoms, name, labels, coords):
   """
   num_atoms - number of atoms 
   name - name of file to write coordinates to
   labels - a list of coordinate labels ("Elements")
   coords - a list of XYZ coordinates

   Writes the Cartesian coordinates to a file called 'name'
   """
   xyz = open(name, "w")
   xyz.write(str(num_atoms) + "\n")
   xyz.write(name + "\n")
   xyz.write(coords_to_xyz(labels, coords))
   xyz.close()

def main():

   description = """
This code uses the Kuhn-Munkres or Hungarian algorithm to optimally align two
arbitrarily ordered isomers. Given two isomers A and B whose Cartesian
coordinates are given in XYZ format, it will optimally align B on A to minimize
the Kabsch root-mean-square deviation (RMSD) between structure A and B after

1) a Kuhn-Munkres assignment/reordering (quick) 
2) a Kuhn-Munkres assignment/reordering factoring in axes swaps and reflections (~48x slower)

We recommend the second method although the first one would still be better
than RMSD calculations without atom reorderings.  

A web server with this implementation is available at http://www.arbalign.org

While this script is kept as minimal as possible in order to ensure ease of use
and portability, it does require these two Python packages beyond what's
included in standard python installations.

1) Python Numpy module 
2) Python Hungarian module by Harold Cooper 
   (Hungarian: Munkres' Algorithm for the Linear Assignment Problem in Python.
   https://github.com/Hrldcpr/Hungarian) 
   This is a wrapper to a fast C++ implementation of the Kuhn-Munkres algorithm.
   The installation instructions are described at https://github.com/Hrldcpr/Hungarian

Other optional tools are:

1) PrinCoords.py - using principal coordinates generally yields better
   alignment (lower RMSDs).  A Python script to convert molecules from arbitrary
   to principal coordinate system is included.

2) In cases where one wants to use atom types including connectivity and
   hybridization information, it is necessary to use OpenBabel to convert the
   Cartesian coordinates to SYBYL Mol2 (sy2) and MNA (mna) formats.  

The best way to take advantage of these two optional tools is probably to use
the attached driver script (ArbAlign-driver.py) The syntax looks like

   Usage: ArbAlign-driver.py -<flag> <filename_1.xyz> <filename_2.xyz>"
        : where the <flag> is "
        : -l   match by atom or element label "
        : -t   match by SYBYL atom type"
        : -c   match by NMA atom connectivity type"
        "
     Eg.: ArbAlign-driver.py -b -N cluster1.xyz cluster2.xyz"
        : ArbAlign-driver.py -T cluster1.xyz cluster2.xyz"
        : ArbAlign-driver.py -C cluster1.xyz cluster2.xyz"
        "
   This matches the Cartesian coordinates of the file1 and file2 using the \
         Kuhn-Munkres algorithm based on atom labels (-l), type (-t) or \
         connectivity (-t). "
   It produces s-file1.xyz and s-file2-matched.xyz which are the sorted and \
         matched file1 and file2.xyz, respectively."

"""

   epilog = """
The code will provide the following:
1) The initial Kabsch RMSD
2) The final Kabsch RMSD after the application of the Kuhn-Munkres algorithm 
3) The coordinates corresponding to the best alignment of B on A to a file called B-aligned_to-A.xyz

If you find this script useful for any publishable work, please cite the corresponding paper:
  Berhane Temelso, Joel M. Mabey, Toshiro Kubota, Nana Appiah-padi, George C. Shields
  J. Chem. Info. Model. X(Y), 2016
"""

   # Parse arguments and provide usage information when necessary
   parser = argparse.ArgumentParser( description=description, formatter_class=argparse.RawDescriptionHelpFormatter,epilog=epilog)

   parser.add_argument('xyz1', metavar='A.xyz', type=str)
   parser.add_argument('xyz2', metavar='B.xyz', type=str)
   parser.add_argument('-s', '--simple', action='store_true', help='Perform Kuhn-Munkres \
           assignment reordering without axes swaps and reflections; \
           the default is to perform axes swaps and reflections')
   parser.add_argument('-n', '--noHydrogens', action='store_true', help='Ignore hydrogens; \
           the default is to include all atoms ')
   parser.add_argument('-z', '--zero', action='store_false', help='does nothing')

   args = parser.parse_args()

   # Read in original coordinates and labels of xyz1 and xyz2
   a_labels, a_coords, NA_a = read_xyz(args.xyz1, args.noHydrogens)
   b_labels, b_coords, NA_b = read_xyz(args.xyz2, args.noHydrogens)
   b_init_labels = b_labels 
   b_init_coords = b_coords
   
   #Calculate the initial unsorted all-atom RMSD as a baseline
   A_all = np.array(a_coords)
   B_all = np.array(b_coords)

   #If the two molecules are of the same size, get 
   if NA_a == NA_b:
      InitRMSD_unsorted = kabsch(A_all,B_all)
      print "InitRMSD_unosrted",InitRMSD_unsorted
   else:
      print "Error: unequal number of atoms. " + str(NA_a) + " is not equal to " + str(NA_b) 
      sys.exit()

   """
   If the initial RMSD is zero (<0.001), then the structured are deemed identical already and 
   we don't need to do any reordering, swapping, or reflections
   """
   if InitRMSD_unsorted < 0.001:
      print "The structures are identical. No reordering, swapping or reflection needed."
      print "All-atom RMSD: %2.3f" % float(InitRMSD_unsorted)
      name = str(args.xyz2.split(".xyz")[0]) + "-aligned_to-" + str(args.xyz1)
      write_to_xyz(num_atoms, name, b_labels, b_coords)
      print "Best alignment of " + str(args.xyz2) + " on " + str(args.xyz1) + " is written to " + str(name)
      sys.exit()

   #If ignoring hydrogens, the coordinates are written to a file with "noHydrogens.xyz" ending
   if args.noHydrogens:
      name = str(args.xyz1.split(".xyz")[0]) + "-noHydrogens" + ".xyz"
      write_to_xyz(NA_a, name, a_labels, a_coords)
      print "Coordinates of " + str(args.xyz1) + " without hydrogens is written to " + str(name)

   """
   Read in the original coordinates and labels of xyz1 and xyz2, 
   and sort them by atom labels so that atoms of the same label/name are grouped together

   Then, count how many types of atoms, and determine their numerical frequency
   """
   a_labels, a_coords, NA_a, order = sorted_xyz(args.xyz1, args.noHydrogens)
   Uniq_a = list(set(a_labels))
   list.sort(Uniq_a)
   N_uniq_a = len(Uniq_a)
   Atom_freq_a = dict(Counter(a_labels))

   b_labels, b_coords, NA_b, junk = sorted_xyz(args.xyz2, args.noHydrogens)
   Uniq_b = list(set(b_labels))
   list.sort(Uniq_b)
   N_uniq_b = len(Uniq_b)
   Atom_freq_b = dict(Counter(b_labels))

   """
   If the number and type of atoms in the two structures are not equal, exit with 
   an error message
   """
   if (NA_a == NA_b) & (Uniq_a == Uniq_b) & (Atom_freq_a == Atom_freq_b) :
      num_atoms = NA_a
      num_uniq = N_uniq_a
      Uniq = Uniq_a  #list(set(a_labels))
      Atom_freq = Atom_freq_a
      #del Atom_freq['H']
      #print Atom_freq
      Sorted_Atom_freq = sorted(Atom_freq.items(), key=operator.itemgetter(1), reverse=True)
      print Sorted_Atom_freq
      """
      Atom = sorted(Uniq, key=operator.itemgetter(0), reverse=True)
      print Atom
      print num_uniq
      """
   else:
      print "Unequal number or type of atoms. Exiting ... "
      print "Atoms in 1st molecule" +str(Atom_freq_a)
      print "Atoms in 2nd molecule" +str(Atom_freq_b)
      sys.exit()

   A_all = np.array(a_coords)
   A_all = A_all - sum(A_all) / len(A_all)
   B_all = np.array(b_coords)
   B_all = B_all - sum(B_all) / len(B_all)
   InitRMSD_sorted = kabsch(A_all,B_all)
   
   """
   Dynamically generate hashes of coordinates and atom indices for every atom type
   """
   a_Coords = {}
   a_Indices = {}
   b_Coords = {}
   b_Indices = {}
   Perm = {}
   for i in range(len(Uniq)):
      a_Coords[Uniq[i]]  = 'a_' + Uniq[i] + 'coords' 
      b_Coords[Uniq[i]]  = 'b_' + Uniq[i] + 'coords' 
      a_Indices[Uniq[i]] = 'a_' + Uniq[i] + 'indices' 
      b_Indices[Uniq[i]] = 'b_' + Uniq[i] + 'indices' 
      Perm[Uniq[i]]      = 'perm_' + Uniq[i]
      #print "Atom_freq is " + str(Atom_freq[Uniq[i]]) 
      vars()[Perm[Uniq[i]]] = []
      #vars()[Perm[Uniq[i]]] = build_perm(Atom_freq[Uniq[i]])
      #for n in range(Atom_freq[Uniq[i]]): 
      #   vars()[Perm[Uniq[i]]] += [(n,n)]
      vars()[a_Coords[Uniq[i]]]  = parse_for_atom(a_labels, a_coords, str(Uniq[i]))
      vars()[a_Indices[Uniq[i]]] = get_atom_indices(a_labels, str(Uniq[i]))
      #print Uniq[i]
      #print vars()[a_Coords[Uniq[i]]]
      vars()[b_Coords[Uniq[i]]]  = parse_for_atom(b_labels, b_coords, str(Uniq[i]))
      vars()[b_Indices[Uniq[i]]] = get_atom_indices(b_labels, str(Uniq[i]))
      #print vars()[b_Indices[Uniq[i]]]

   l = 0 
   A = np.array(vars()[a_Coords[Uniq[l]]])
   A = A - sum(A) / len(A)
   B = np.array(vars()[b_Coords[Uniq[l]]])
   B = B - sum(B) / len(B)

   '''
   For each atom type, we can do a Kuhn-Munkres assignment in the initial 
   coordinates or the many swaps and reflections thereof

   If a single Kuhn-Munkres assignment is requested with a -s or --simple flag,
   no swaps and reflections are considered. Otherwise, the default is to perform 
   a combination of 6 axes swaps and 8 reflections and do Kuhn-Munkres assignment 
   on all 48 combinations. 
   '''
   if args.simple:      # will do nothing
      swaps = [(0, 1, 2)]
      reflects = [(1, 1, 1)]
   else:                # will perform swaps and reflections
      swaps = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
      reflects = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1), \
               (-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1)]
   B_t = []
   for i in swaps:
      for j in reflects:
         B_t.append([transform_coords(B, i, j), i, j])

   rmsds = []
   # Performs the munkres algorithm on each set of transformed coordinates
   for i in range(len(B_t)):
      l = 0
      cost_matrix = np.array([[np.linalg.norm(a - b) \
                                  for b in B_t[i][0]] for a in A])
      LAP = hungarian.lap(cost_matrix)
      vars()[Perm[Uniq[l]]] = []
      for j in range(len(LAP[0])):
         vars()[Perm[Uniq[l]]] += [(j,LAP[0][j])]
      vars()[Perm[Uniq[l]]] = sorted( vars()[Perm[Uniq[l]]], key = lambda x: x[0])
      vars()[Perm[Uniq[l]]] = [x[1] for x in vars()[Perm[Uniq[l]]]]

      # If there's more than one atom type, loop through each unique atom type 
      if num_uniq == 1:
         #print str(vars()[b_Indices[Uniq[l]]])
         b_perm = permute_atoms(b_coords, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
         b_final = transform_coords(b_perm, B_t[i][1], B_t[i][2])
         print str(Uniq[l]) + " Swap: " + str(B_t[i][1]) + " Refl: " + str(B_t[i][2]) + " RMSD: " + str(kabsch(a_coords, b_final)) + " " + str(vars()[Perm[Uniq[l]]])
         rmsds.append([kabsch(a_coords, b_final), B_t[i][1], B_t[i][2], b_final, vars()[Perm[Uniq[l]]]])
         rmsds = sorted(rmsds, key = lambda x: x[0])
      else: 
         #print str(vars()[b_Indices[Uniq[l]]])
         b_perm = permute_atoms(b_coords, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
         b_trans = transform_coords(b_perm, B_t[i][1], B_t[i][2])
         #print str(b_trans)
         #vars()[b_Coords[Uniq[l+1]]] = parse_for_atom(b_labels, b_trans, Uniq[l+1])
         while l < num_uniq:
            if l > 0:
               vars()[b_Coords[Uniq[l]]] = parse_for_atom(b_labels, b_final, Uniq[l])
            else:
               vars()[b_Coords[Uniq[l]]] = parse_for_atom(b_labels, b_trans, Uniq[l])
            cost_matrix = np.array([[np.linalg.norm(a- b) \
                                         for b in np.array(vars()[b_Coords[Uniq[l]]])] \
                                         for a in np.array(vars()[a_Coords[Uniq[l]]])])
            LAP = hungarian.lap(cost_matrix)
            vars()[Perm[Uniq[l]]] = []
            for k in range(len(LAP[0])):
               vars()[Perm[Uniq[l]]] += [(k,LAP[0][k])]
            vars()[Perm[Uniq[l]]] = sorted( vars()[Perm[Uniq[l]]], key = lambda x: x[0])
            vars()[Perm[Uniq[l]]] = [x[1] for x in vars()[Perm[Uniq[l]]]]
            #print str(vars()[b_Indices[Uniq[l]]])
            b_final = permute_atoms(b_trans, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
            b_trans = b_final
            l += 1
            q = l - 1 
            print str(Uniq[q]) + " Swap: " + str(B_t[i][1]) + " Refl: " + str(B_t[i][2]) + " RMSD: " + str(kabsch(a_coords, b_final)) + " " + str(vars()[Perm[Uniq[q]]])
            rmsds.append([kabsch(a_coords, b_final), B_t[i][1], B_t[i][2], b_final])
            rmsds = sorted(rmsds, key = lambda x: x[0])
            #print "Permutation: " + str(vars()[Perm[Uniq[q]]])
   
   if not args.simple:
      print "Swap Transform: " + str(rmsds[0][1])
      print "Reflection Transform: " + str(rmsds[0][2])

   #print "Permutation: " + str(rmsds[0][4])
   FinalRMSD = float(rmsds[0][0])
   if FinalRMSD < float(InitRMSD_unsorted): 
      print "Initial unsorted RMSD: %2.3f" % float(InitRMSD_unsorted)
      print "Initial   sorted RMSD: %2.3f" % float(InitRMSD_sorted)
      print "Best             RMSD: %2.3f" % float(rmsds[0][0])
   else:
      print "The initial alignment is already optimal."
      print "Initial and final RMSD: %2.3f" % float(InitRMSD_unsorted)

   if args.noHydrogens:
      name = str(args.xyz2.split(".xyz")[0]) + "-aligned_to-" + \
             str(args.xyz1.split(".xyz")[0]) + "-noHydrogens.xyz" 
   else:
      name = str(args.xyz2.split(".xyz")[0]) + "-aligned_to-" + str(args.xyz1)

   if FinalRMSD < float(InitRMSD_unsorted): 
      b_final_labels, b_final_coords = permute_all_atoms(b_labels, rmsds[0][3], order)
   else:
      b_final_labels = b_init_labels
      b_final_coords = b_init_coords

   write_to_xyz(num_atoms, name, b_final_labels, b_final_coords)
   print "Best alignment of " + str(args.xyz2) + " with " + str(args.xyz1) + " is written to " + str(name)

   #miha
   # permute B
#   bb = b_final_coords
#   for i in range(num_atoms):
#       ind = order[i]
#       b_final_coords[i] = bb[ind]

#   print order
#   print len(a_coords)
#   print ''
#   for i in a_coords:
#       print 1, i[0],i[1],i[2]
#
#   print len(a_coords)
#   print ''
#   for i in range(num_atoms):
#       ind = order[i]
##       ind = i
#       r = np.array(b_final_coords[ind])
#       print 1, r[0],r[1],r[2]
#
#   rmsd1=0.0
#   for i in range(num_atoms):
#       ind = order[i]
##       ind = i
#       rij = np.array(a_coords[i])-np.array(b_final_coords[ind])
#       rmsd1+=np.dot(rij,rij)
#
#   print 'rmsd:',(1.0/num_atoms*rmsd1)**0.5
#   print rmsd( a_coords, b_final_coords)
#   print kabsch( a_coords, b_final_coords)

if __name__ == "__main__":
   main()
