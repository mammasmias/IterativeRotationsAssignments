#!/usr/bin/env python2.7

"""
moi.py

(c) 2008, James Stroud
This module is released under the GNU Public License 3.0.
See <http://www.gnu.org/licenses/gpl.html>.

A module to find the moment of inertia.

See <http://en.wikipedia.org/wiki/Moment_of_inertia> for
a discussion.

Anywhere "ref" is optional defaults to the center of mass
except where indicated.
The "sel" argument is always a list of Atoms.
""" 

import sys
import math
from textwrap import dedent
from ConfigParser import ConfigParser
from operator import mul
from itertools import izip
from numpy import *
from numpy.linalg import *

PLANCK = 6.626068e-34
PLANCK_SQR = PLANCK**2
C = 2.99792458e8 
BOLTZ = 1.3806503e-23
NA = 6.022e23
SQRT_PI = math.sqrt(math.pi)
PI_SQR = math.pi**2

elements = {
             'H': 1.00782503207,
             'D': 2.0141017778,
             'T': 3.0160492777,
             'HE': 4.002602, 
             'C12': 12.000000,
             'C13': 13.0033548378,
             'O16': 16.00,
             'O18': 18.00,
             'C' : 12.0107, 
             'N' : 14.00674, 
             'O' : 15.9994, 
             'F' : 18.9984032, 
             'NE' : 20.1797, 
             'NA' : 22.98977, 
             'MG' : 24.305, 
             'AL' : 26.981538, 
             'SI' : 28.0855, 
             'P' : 30.973761, 
             'S' : 32.066, 
             'CL' : 35.4527, 
             'AR' : 39.948, 
             'K' : 39.0983, 
             'CA' : 40.078, 
             'SC' : 44.95591, 
             'TI' : 47.867, 
             'V' : 50.9415, 
             'CR' : 51.9961, 
             'Mn' : 54.938049, 
             'FE' : 55.845, 
             'CO' : 58.9332, 
             'NI' : 58.6934, 
             'CU' : 63.546, 
             'ZN' : 65.39, 
             'GA' : 69.723, 
             'GE' : 72.61, 
             'AS' : 74.9216, 
             'SE' : 78.96, 
             'BR' : 79.904, 
             'KR' : 83.8, 
             'RB' : 85.4678, 
             'SR' : 87.62, 
             'Y' : 88.90585, 
             'ZR' : 91.224, 
             'NB' : 92.90638, 
             'MO' : 95.94, 
             'TC' : -98, 
             'RU' : 101.07, 
             'RH' : 102.9055, 
             'PD' : 106.42, 
             'AG' : 107.8682, 
             'CD' : 112.411, 
             'IN' : 114.818, 
             'SN' : 118.71, 
             'SB' : 121.76, 
             'TE' : 127.6, 
             'I' : 126.90447, 
             'XE' : 131.29, 
             'CS' : 132.90545, 
             'BA' : 137.327, 
             'HF' : 178.49, 
             'TA' : 180.9479, 
             'W' : 183.84, 
             'RE' : 186.207, 
             'OS' : 190.23, 
             'IS' : 192.217, 
             'PT' : 195.078, 
             'AU' : 196.96655, 
             'HG' : 200.59, 
             'TL' : 204.3833, 
             'PB' : 207.2, 
             'BI' : 208.98038
           }

class Atom(object):
  def __init__(self, x, y, z, mass):
    self.x = float(x)
    self.y = float(y)
    self.z = float(z)
    self.mass = float(mass)
  def __str__(self):
    return "%s @ (%s, %s, %s)" % (self.x, self.y, self.z, self.mass)

class Point(object):
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
  def __str__(self):
    return "(%s, %s, %s)" % (self.x, self.y, self.z)

def read_xyz(filename):
  s = open(filename).readlines()
  atoms = []
  for aline in s[2:]:
    e, x, y, z = aline.upper().split()
    mass = elements[e.upper()]
    x = float(x)
    y = float(y)
    z = float(z)
    atom = Atom(x, y, z, mass)
    atoms.append(atom) 
  return atoms 

def readXYZLabels(filename):
  s = open(filename).readlines()
  return [k.split()[0].upper() for k in s[2:]]

def read_config(cfgfile):
  cfg = ConfigParser()
  cfg.read(cfgfile)
  files = []
  for section in sorted(cfg.sections()):
    files.append((cfg.get(section, 'file'), float(cfg.get(section, 'symm'))))
  return files

def prod(alist):
  return reduce(mul, alist)

def molecular_mass(sel):
  """
  Calculates the molecular mass 
  """
  mol_mass= sum(atom.mass for atom in sel)
  return mol_mass 

def com(sel):
  """
  Calculates the center of mass of a list of atoms.
  """
  m = sum(atom.mass for atom in sel)
  c = []
  for i in 'xyz':
    mr = []
    for atom in sel:
      mr.append(atom.mass * getattr(atom, i))
    c.append(sum(mr)/m)
  return Point(*c)

def dist(atom, ref=None):
  """
  Calculates distance between an Atom (or Point) and a ref.
  If ref is not given, then the origin is assumed.
  """
  if ref is None:
    d = math.sqrt(atom.x**2 + atom.y**2 + atom.z**2)
  else:
    x = atom.x - ref.x
    y = atom.y - ref.y
    z = atom.z - ref.z
    d = math.sqrt(x**2 + y**2 + z**2)
  return d

def kd(i, j):
  """
  The Kronaker delta.
"""
  return int(i==j)

def moi_element(i, j, sel, ref):
  """
  Generalized moment of inertia element--both diagonal and off
  diagonal. Coordinate axes i and j should be 'x', 'y', or 'z'.
  The sel argument is a list of Atoms.
  The ref arguemnt is a Point.
  """
  delta = kd(i, j)
  refi = getattr(ref, i)
  refj = getattr(ref, j)
  el = []
  for atom in sel:
    r2d = (dist(atom, ref)**2) * delta
    ri = getattr(atom, i) - refi
    rj = getattr(atom, j) - refj
    el.append(atom.mass*(r2d - ri*rj))
  return sum(el)

def tensor_moi(sel, ref=None):
  """
  Calculates the cartesian moment of inertia tensor.
  """
  if ref is None:
    ref = com(sel)
  indices = [(i,j) for i in 'xyz' for j in 'xyz']
  tensor = [moi_element(i, j, sel, ref) for (i,j) in indices]
  return [tensor[0:3], tensor[3:6], tensor[6:9]]

def principal_moi(sel, ref=None):
  moi = tensor_moi(sel, ref)
  pmoi = eigvalsh(moi)*1.6605e-47
  return [PLANCK/(8*PI_SQR*pmoi[i]*1e9) for i in xrange(3)]

def rot_const(sel, ref=None):
  moi = tensor_moi(sel, ref)
  pmoi = eigvalsh(moi)*1.6605e-47
  eigenVectors = eig(moi)[1]

def eigenVectors(sel, ref=None):
  moi = tensor_moi(sel, ref)
  #eigen_vectors = eig(moi)[1]
  evals,eigen_vectors = eig(moi)
  index = evals.argsort()[::-1]
  eigen_vectors = eigen_vectors[:,index]
  return eigen_vectors

def scalar_moi(sel, ref=None):
  """
  Calculates the scalar moment of inertia.
  """
  if ref is None:
    ref = com(sel)
  i = []
  for atom in sel:
    i.append(atom.mass * dist(atom, ref)**2)
  return sum(i)

def rotational_entropy(I, symm, N, T):
  """
  from http://www.codessa-pro.com/descriptors/thermodynamic/rem.htm
  """
  k1 = N * BOLTZ
  k2 = SQRT_PI / symm
  k3 = 8 * PI_SQR * BOLTZ * T / PLANCK_SQR
  return k1 * math.log(k2 * prod(math.sqrt(k3 * Ij) for Ij in I))

def rotationMatrix(A):
  '''
  calculates the rotation matrix to get to A
  A - eigenvectors of MOI matrix
  '''
  I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
  cosines = [dot(I[k], A[k]) for k in range(3)]
  angles = [math.acos(cosines[k]) for k in range(3)]
  sines = [math.sin(angles[k]) for k in range(3)]
  Rx = [[1, 0, 0], [0, cosines[0], -sines[0]], [0, sines[0], cosines[0]]]
  Ry = [[cosines[1], 0, sines[1]], [0, 1, 0], [-sines[1], 0, cosines[1]]]
  Rz = [[cosines[2], -sines[2], 0], [sines[2], cosines[2], 0], [0, 0, 1]]
  R = dot(dot(Rx, Rz), Ry)
  return R
 
def printMatrix(A):
  for i in range(len(A)):
    for j in range(len(A[0])):
      print '%8f' % A[i][j],
    print                    

def generateXYZFile(N, L, C):
  num_atoms = len(C)
  filename = N + '-prin.xyz'
  f = open(filename, 'w')
  f.write(str(num_atoms)+'\n')
  f.write(N+'\n')
  for i in range(len(L)):
    #f.write(str(L[i]) + ' ' + str(C[i][0]) + ' ' + str(C[i][1]) + ' ' + str(C[i][2])+'\n')
    f.write('%-2s  %12.6f  %12.6f  %12.6f \n' % (L[i], C[i][0], C[i][1], C[i][2]))
  f.close() 
      
def usage():
  u = """
      usage: moi.py configfile.cfg

        NOTE: configfile.cfg is an ini type file.

        EXAMPLE

          [sa-dimer]
          file = SA-dimer.xyz
          symm = 2

          [some-tetramer]
          file = Some-tetramer.xyz
          symm = 4
      """
  print dedent(u)
  sys.exit()

def main():
  '''
  try:
    cfgfile = sys.argv[1]
  except:
    usage()

  files = read_config(cfgfile)
  '''
  atomXYZ = sys.argv[1]
  #symmetry = sys.argv[2]
  atoms = read_xyz(atomXYZ)
  
  orig_coords = [[atoms[k].x, atoms[k].y, atoms[k].z] for k in range(len(atoms))]
  #print '\nOriginal Coordinates:'
  #printMatrix(orig_coords)
  com_atoms = com(atoms)
  #print '\nCOM'
  #print 'x: ', com_atoms.x, 'y: ', com_atoms.y, 'z: ', com_atoms.z
  com_coords = [[orig_coords[k][0] - com_atoms.x, orig_coords[k][1] - com_atoms.y, orig_coords[k][2] - com_atoms.z] \
                 for k in range(len(orig_coords))]
  #print '\nCOM Coordinates'
  #printMatrix(com_coords)
  eigen_vectors = eigenVectors(atoms)
  #print '\nEigen Vectors'
  #printMatrix(eigen_vectors)
  #print '\nRotation Matrix'
  #printMatrix(rotationMatrix(eigen_vectors))
  new_coords = dot(com_coords, eigen_vectors)
  #print '\nNew Coordinates'
  #printMatrix(new_coords)
  
  labels = readXYZLabels(atomXYZ)
  name = atomXYZ[:-4]
  generateXYZFile(name, labels, new_coords)

  B = rot_const(atoms)
  M = molecular_mass(atoms)
  #print 'Mass ', M
  #print 'GHZ ', B
  """
  rot_ent = rotational_entropy(I, symm, 1.0, 298.15)
  print
  print "=="
  print "file:", filename
  print "I:", I
  print "rot_ent:", rot_ent
  print "=="
  print
  print "=="
  print "file:", filename
    print "=="
  """
  '''
  for (filename, symm) in files:
    atoms = read_xyz(filename)
    B = rot_const(atoms)
    M = molecular_mass(atoms)
    print "MASS  ", M 
    print "GHZ   ", B
  '''

if __name__ == "__main__":
  main()
