import numpy as np

#
# auxiliary function for reading .xyz files
#
def read_xyz(fname):
    typ = np.genfromtxt( fname, skip_header=2, usecols=[0], dtype=None, encoding=None )
    coords = np.genfromtxt( fname, skip_header = 2, usecols = [1,2,3], dtype=np.float64 )
    nat = len( typ )
    cc=np.ndarray( (nat,3), dtype=np.float64, order="C" )
    cc=coords
    return nat, typ, cc


#
# read .xyz files "s1.xyz" and "s2.xyz" containing some near-congruent structures
#

nat1, typ1, coords1 = read_xyz( "example_inputs/s1.xyz" )
nat2, typ2, coords2 = read_xyz( "example_inputs/s2.xyz" )




#
# import the IRA module
#
import ira_mod
#
#=========================================================================
# For more information about the functions in the ira_mod, you can use the
# `help()` function of python, e.g. `help(ira_mod)`
#=========================================================================
#



#
# initialize IRA
#
ira = ira_mod.IRA()

#=================================================================
# Demonstration how to match structures with equal number of atoms:
#=================================================================
print("====================================================")
print("=  Demonstration Python program to call IRA for    =")
print("=  matching structures with equal number of atoms  =")
print("====================================================")

print( nat1 )
print( " original structure 1" )
for i in range( nat1 ):
    print( "%i %.4f %.4f %.4f" %(typ1[i], coords1[i][0], coords1[i][1], coords1[i][2]) )
print( nat2 )
print( " original structure 2" )
for i in range( nat2 ):
    print( "%i %.4f %.4f %.4f" %(typ2[i], coords2[i][0], coords2[i][1], coords2[i][2]) )
print("")
#
# compute Hausdorff distance between the two structures as-read from the file,
# using the CShDA algorithm
#
perm, dists = ira.cshda( nat1, typ1, coords1, nat2, typ2, coords2 )
#
# The `dists` array contains distances between the assigned atom pairs.
# The Hausdorff distance is the maximal distance, `np.max( dists )`
#
print( "initial Hausdorff distance:", np.max(dists) )


#
# perform the IRA shape-matching algorithm:
#
kmax_factor = 1.8
rmat, tr, perm, dh = ira.match( nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor )

print( "Hausdorff distance after matching:", dh )

print( "Rotation matrix:" )
print( rmat )
print( "translation vector" )
print( tr )
print( "permutation of atoms:" )
print( perm )
print("")

#
# apply found transformation:
#
for i in range( nat2 ):
    coords2[i] = np.matmul( rmat, coords2[i] ) + tr
#
# apply permutation
#
coords2[:] = coords2[perm]
typ2[:] = typ2[perm]




#
# structures 1 and 2 should now be matched:
#
print( nat2 )
print( " matched structure 2" )
for i in range( nat2 ):
    print( "%i %.4f %.4f %.4f" %(typ2[i], coords2[i][0], coords2[i][1], coords2[i][2]) )





# ==============================
# Demonstration for structures with different number of atoms:
# ==============================
print("")
print("=====================================================")
print(" Matching structures with different number of atoms  ")
print("=====================================================")

#
# re-read the original structures
#
nat1, typ1, coords1 = read_xyz( "example_inputs/s1.xyz" )
nat2, typ2, coords2 = read_xyz( "example_inputs/s2.xyz" )
#
# for purpose of demonstration:
#    remove last three atoms from structure1, such that it is smaller
#    than structure2
#
nat1 = nat1 - 3
typ1=typ1[:-3]
coords1=coords1[:-3]


print( nat1 )
print( " original structure 1 (smaller by 3 atoms)" )
for i in range( nat1 ):
    print( "%i %.4f %.4f %.4f" %(typ1[i], coords1[i][0], coords1[i][1], coords1[i][2]) )

print( nat2 )
print( " original structure 2" )
for i in range( nat2 ):
    print( "%i %.4f %.4f %.4f" %(typ2[i], coords2[i][0], coords2[i][1], coords2[i][2]) )
print("")

# compute cshda of the two structures
perm_p, dists_p = ira.cshda( nat1, typ1, coords1, nat2, typ2, coords2 )

# when different number of atoms, the `dists_p` array has size `nat2`, such that the first
# `nat1` number of values contains distances ebtween assigned pairs of atoms, and the values
# beyond this index are "very large", i.e. "not assigned".
# The Hausdorff distance is thus the maximal value among the
# first `nat1` values: np.max( dists_p[:nat1] )
print("initial Hausdorff of partial:", np.max( dists_p[:nat1] ) )


#
# find matching
#
rmat_p, tr_p, perm_p, dh_p = ira.match( nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor )

print("Hausdorff distance of matched:", dh_p )
print( "Rotation matrix:" )
print( rmat_p )
print( "translation vector" )
print( tr_p )
print( "permutation of atoms:" )
print( perm_p )
print("")

#
# apply found transformation:
#
for i in range( nat2 ):
    coords2[i] = np.matmul( rmat_p, coords2[i] ) + tr_p
#
# apply permutation
#
coords2[:] = coords2[perm_p]
typ2[:] = typ2[perm_p]

print( nat2 )
print( " matched structure 2 (the first 10 atoms should match the smaller structure1)" )
for i in range( nat2 ):
    print( "%i %.4f %.4f %.4f" %(typ2[i], coords2[i][0], coords2[i][1], coords2[i][2]) )
