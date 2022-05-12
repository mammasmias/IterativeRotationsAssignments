import numpy as np

#
# import the IRA module
#
import ira_mod

#
# auxiliary function for reading .xyz files
#
def read_xyz(fname):
  typ  = np.loadtxt( fname, skiprows = 2, usecols = [0], dtype = int )
  coords = np.loadtxt( fname, skiprows = 2, usecols = [1,2,3] )
  nat = len( typ )
  return nat, typ, coords



#
# read .xyz files "s1.xyz" and "s2.xyz"
#

nat1, typ1, coords1 = read_xyz( "example_inputs/s1.xyz" )
nat2, typ2, coords2 = read_xyz( "example_inputs/s2.xyz" )




#
# initialize rotation matrix and translation vector
#

r=np.zeros(( 3, 3 ))
t=np.zeros( 3 )




#
# call ira + svd in single call.
# Attention!: when structures have different number of atoms, the structure with
# lower number of atoms should belabelled as structure 1!
#
# input variables:
#   nat1    -> integer, number of atoms in structure 1
#   typ1    -> array of integers size(nat1), atomic types of atoms in structure 1
#   coords1 -> array of floats size(3,nat1), coordinates of atoms in structure 1
#   ---||---   similarly for structure 2
#   kmax_factor -> float, multiplicative factor for the radial distance of the basis check, should be >= 1. Default is 1.2
# output variables:
#   r   -> float size(3,3) rotation matrix
#   t   -> float size(3), translation vector
#   p   -> integer array size(nat2), permutation order
#   hd  -> float, Hausdorff distance computed after the matching
#   rmsd -> float, root-mean-square-distance computed after the matching
#

kmax_factor = 1.2

r, t, p, hd, rmsd = ira_mod.ira_svd( nat1, typ1, np.transpose(coords1), nat2, typ2, np.transpose(coords2), kmax_factor )



#
# apply found transformation:
#

#
# permutation order p is always the permutation of structure 2.
# permute with [p-1] because fortran starts arrays at 1 and python at 0
#
typ2=typ2[p-1]
coords2=coords2[p-1]
#
# rotate and translate
#
for i in range(nat2):
    coords2[i] = np.matmul( r, coords2[i] ) + t
#
# Alternatively, r and t can be applied to structure 1, in the following way:
#
#for i in range(nat1):
#    coords1[i] = np.matmul( np.transpose(r), coords1[i] ) - np.matmul( np.transpose(r), t )



#print(r)
#print(t)
#print(p)

print( 'matched structures' )
print ( nat1 )
print()
print(coords1)

print(nat2)
print()
print( coords2 )

print( 'computed dH, rmsd:' )
print(hd, rmsd)
