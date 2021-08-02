#!/usr/bin/env python2.7

import sys
import subprocess
import argparse
#import numpy as np
#import hungarian
#from collections import Counter 
#import operator

def read_and_write_xyz(InFile, OutFile, AlignBy):
   """
   Reads an xyz file into lists
   InFile - name of XYZ file to be read to
   OutFile - name of XYZ file to be written to
   AlignBy - Alignment based on label (l), type (t) or connectivity (c)

   Writes OutFile
   """
   in_xyz = open(InFile, "r")
   num_atoms = int(in_xyz.readline().strip())
   in_xyz.readline()
   #print AlignBy
   unsorted_labels = []
   unsorted_coords = []
   for i in range(num_atoms):
      line = in_xyz.readline().strip().split()
      unsorted_labels.append(line[0].upper())
      if AlignBy == 't':
         unsorted_labels[i] = str(unsorted_labels[i].split(".")[0])
      elif AlignBy == 'c':
         unsorted_labels[i] = str(unsorted_labels[i].split("-")[0])
      else:
         AlignBy = 'l'
         #str(unsorted_labels[i])
      unsorted_coords.append([float(line[1]), float(line[2]), float(line[3])])
   in_xyz.close()

   out_xyz = open(OutFile, "w")
   out_xyz.write(str(num_atoms) + "\n")
   out_xyz.write(OutFile + "\n")
   out_xyz.write(coords_to_xyz(unsorted_labels, unsorted_coords))
   out_xyz.close()

def coords_to_xyz(labels, coords):
   """
   Displays coordinates
   """
   s = ""
   for i in range(len(coords)):
      s += labels[i] + " " + str(coords[i][0]) + " " +  \
             str(coords[i][1]) + " " + str(coords[i][2]) + "\n"
   return s

def main():

   scrdir = sys.path[0]
   description = """
This script is a driver to the Kuhn-Munkres alignment script to optimally align two
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

Other optional tools if you plan to align molecules by atom type or connectivity instead of atom
label alone are:

1) PrinCoords.py - using principal coordinates generally yields better
   alignment (lower RMSDs).  A Python script to convert molecules from arbitrary
   to principal coordinate system is included.

2) genType.csh and genConn.csh -- in cases where one wants to use atom types including connectivity and
   hybridization information, it is necessary to use OpenBabel to convert the
   Cartesian coordinates to SYBYL Mol2 (sy2) and MNA (mna) formats. These two shell scripts will
   take the SYBYL Mol2 (sy2) and MNA (mna) information and cast the atom name in the XYZ file to
   reflect that information. 

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
   parser.add_argument('-s', '--simple', action='store_true', \
           help='Perform Kuhn-Munkres \
           assignment reordering without axes swaps and reflections; \
           the default is to perform axes swaps and reflections')
   parser.add_argument('-n', '--noHydrogens', action='store_true', \
           help='Ignore hydrogens; the default is to include all atoms ')
   parser.add_argument('-b', '--by', type=str, choices=['l', 't', 'c'], default='l', \
           help='match atoms by l-label, SYBYL t-type, or NMA connectivity (-c). \
                The default is by atom label (-l)')
   args = parser.parse_args()
   # See if noHydrogens or simple options are selected
   if args.noHydrogens :
       if args.simple:
         sargs = str('-sn') 
       else:
         sargs = str('-n')
       in_xyz1 = str(args.xyz1.split(".xyz")[0]) + "-noHydrogens" + ".xyz"
       in_xyz2 = str(args.xyz2.split(".xyz")[0]) + "-noHydrogens" + ".xyz"
       ref_file = in_xyz1
       t_ref_file = "t-"+ str(in_xyz1)
       c_ref_file = "c-"+ str(in_xyz1)
       OutFile = str(args.xyz2.split(".xyz")[0]) + "-aligned_to-" + str(args.xyz1.split(".xyz")[0])  + "-noHydrogens.xyz"
       t_OutFile = "t-"+ str(args.xyz2.split(".xyz")[0])  + "-aligned_to-t-" + str(args.xyz1.split(".xyz")[0])  + "-noHydrogens.xyz"
       c_OutFile = "c-"+ str(args.xyz2.split(".xyz")[0])  + "-aligned_to-c-" + str(args.xyz1.split(".xyz")[0])  + "-noHydrogens.xyz"
   elif args.simple:
       sargs = str('-s')
       in_xyz1 = str(args.xyz1)
       in_xyz2 = str(args.xyz2)
       ref_file = in_xyz1
       t_ref_file = "t-"+ str(in_xyz1)
       c_ref_file = "c-"+ str(in_xyz1)
       OutFile = str(in_xyz2.split(".xyz")[0]) + "-aligned_to-" + str(in_xyz1)
       t_OutFile = "t-" + str(in_xyz2.split(".xyz")[0]) + "-aligned_to-t-" + str(in_xyz1)
       c_OutFile = "c-" + str(in_xyz2.split(".xyz")[0]) + "-aligned_to-c-" + str(in_xyz1)
   else:
       sargs = str('-z') 
       in_xyz1 = str(args.xyz1)
       in_xyz2 = str(args.xyz2)
       ref_file = in_xyz1
       t_ref_file = "t-"+ str(in_xyz1)
       c_ref_file = "c-"+ str(in_xyz1)
       OutFile = str(args.xyz2.split(".xyz")[0]) + "-aligned_to-" + str(args.xyz1)
       t_OutFile = "t-" + str(args.xyz2.split(".xyz")[0]) + "-aligned_to-t-" + str(args.xyz1)
       c_OutFile = "c-" + str(args.xyz2.split(".xyz")[0]) + "-aligned_to-c-" + str(args.xyz1)

   if args.by == 'l':
      subprocess.call([str(scrdir)+"/ArbAlign.py", sargs, args.xyz1, args.xyz2])
      read_and_write_xyz(ref_file, 'reference.xyz', 'l')
      read_and_write_xyz(OutFile, 'molecule-aligned_to-reference.xyz', 'l')
   elif args.by == 't':
      tprin_file1 = str("t-") + str(args.xyz1)
      tprin_file2 = str("t-") + str(args.xyz2)
      subprocess.call([str(scrdir)+"/genTypes.csh", args.xyz1])
      subprocess.call([str(scrdir)+"/genTypes.csh", args.xyz2])
      subprocess.call([str(scrdir)+"/ArbAlign.py", sargs, tprin_file1, tprin_file2])
      OutFile = t_OutFile
      read_and_write_xyz(t_ref_file, 'reference.xyz', 't')
      read_and_write_xyz(OutFile, 'molecule-aligned_to-reference.xyz', 't')
   elif args.by == 'c':
      cprin_file1 = str("c-") + str(args.xyz1)
      cprin_file2 = str("c-") + str(args.xyz2)
      subprocess.call([str(scrdir)+"/genConn.csh", args.xyz1])
      subprocess.call([str(scrdir)+"/genConn.csh", args.xyz2])
      subprocess.call([str(scrdir)+"/ArbAlign.py", sargs, cprin_file1, cprin_file2])
      OutFile = c_OutFile
      read_and_write_xyz(c_ref_file, 'reference.xyz', 'c')
      read_and_write_xyz(OutFile, 'molecule-aligned_to-reference.xyz', 'c')
   else:
       print "Error: unsure how to perform the alignment"


if __name__ == "__main__":
   main()
