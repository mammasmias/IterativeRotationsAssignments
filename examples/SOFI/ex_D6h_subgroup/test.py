import ira_mod

sofi=ira_mod.SOFI()

import numpy as np


# define the axes
ax = np.ndarray( (7,3), dtype = float )
ax[0] = np.array([  0.9896, -0.0941, -0.1090 ])
ax[1] = np.array([ -0.0518, -0.9391,  0.3397 ])
ax[2] = np.array([ -0.1343, -0.3305, -0.9342 ])
ax[3] = np.array([ -0.0222,  0.6480, -0.7613 ])
ax[4] = np.array([ -0.1120, -0.9785, -0.1729 ])
ax[5] = np.array([ -0.1422, -0.7558, -0.6392 ])
ax[6] = np.array([  0.0904, -0.1833,  0.9789 ])


## construct the D6h object
d6h_matrix = np.zeros( (24, 3, 3), dtype = float )

# fill the theta operations

# identity
d6h_matrix[0] = sofi.construct_operation( "E", ax[0], 0.0 )
# inversion
d6h_matrix[1] = sofi.construct_operation( "I", ax[0], 0.0 )

m=2
for i in range(7):
   # construct S0 on all axes
   d6h_matrix[m] = sofi.construct_operation( "S", ax[i], 0.0 )
   m += 1
   # construct C2 on all axes
   d6h_matrix[m] = sofi.construct_operation( "C", ax[i], 0.5 )
   m += 1

# fill the rest on ax[0]
# C3
d6h_matrix[m] = sofi.construct_operation( "C", ax[0], 1/3.0 )
m+=1
d6h_matrix[m] = sofi.construct_operation( "C", ax[0], -1/3.0 )
m+=1

# C6
d6h_matrix[m] = sofi.construct_operation( "C", ax[0], 1/6.0 )
m+=1
d6h_matrix[m] = sofi.construct_operation( "C", ax[0], -1/6.0 )
m+=1

# S3
d6h_matrix[m] = sofi.construct_operation( "S", ax[0], 1/3.0 )
m+=1
d6h_matrix[m] = sofi.construct_operation( "S", ax[0], -1/3.0 )
m+=1

# S6
d6h_matrix[m] = sofi.construct_operation( "S", ax[0], 1/6.0 )
m+=1
d6h_matrix[m] = sofi.construct_operation( "S", ax[0], -1/6.0 )

# confirm we have D6h
d6h_pg = sofi.get_pg( 24, d6h_matrix )
print( "initial PG is:", d6h_pg )
op0, ax0, angle0 = sofi.get_unique_ax_angle( 24, d6h_matrix )
print( "The operations are:")
print( f"Op   angle    axis")
for i in range( 24 ):
   print( f"{op0[i]:3s} {angle0[i]:6.3f} {ax0[i][0]:9.4f} {ax0[i][1]:9.4f} {ax0[i][2]:9.4f}")

print()



# take all operations from axes 5 and 7 (which is 4 and 6 counting from 0)
four_ops1 = np.zeros( (4,3,3), dtype = float)
four_ops1[0] = sofi.construct_operation( "S", ax[4], 0.0)
four_ops1[1] = sofi.construct_operation( "C", ax[4], 0.5)
four_ops1[2] = sofi.construct_operation( "S", ax[6], 0.0)
four_ops1[3] = sofi.construct_operation( "C", ax[6], 0.5)

# construct matrix combinations
n_combo1, matrix_combo1 = sofi.mat_combos(4, four_ops1 )
print( "combinations of matrices from axes 5 and 7 are:", n_combo1 )

# get the op, ax, angle data from matrices
ops1, ax1, angle1 = sofi.get_unique_ax_angle( n_combo1, matrix_combo1 )
# print that info
print( "The operations are:")
print( f"Op   angle    axis")
for i in range( n_combo1 ):
   print( f"{ops1[i]:3s} {angle1[i]:6.3f} {ax1[i][0]:9.4f} {ax1[i][1]:9.4f} {ax1[i][2]:9.4f}")

# find the PG
pg_combo1 = sofi.get_pg( n_combo1, matrix_combo1 )
print("Corresponding PG:", pg_combo1)
print()




# take all operations from axes 6 and 7 (which is 5 and 6 counting from 0)
four_ops2 = np.zeros( (4,3,3), dtype = float)
four_ops2[0] = sofi.construct_operation( "S", ax[5], 0.0)
four_ops2[1] = sofi.construct_operation( "C", ax[5], 0.5)
four_ops2[2] = sofi.construct_operation( "S", ax[6], 0.0)
four_ops2[3] = sofi.construct_operation( "C", ax[6], 0.5)

# construct matrix combinations
n_combo2, matrix_combo2 = sofi.mat_combos(4, four_ops2)
print( "combinations of matrices from axes 6 and 7 are:", n_combo2 )

# get the op, ax, angle data from matrices
ops2, ax2, angle2 = sofi.get_unique_ax_angle( n_combo2, matrix_combo2 )
# print that
print( "The operations are:", n_combo2)
print( f"Op   angle    axis")
for i in range( n_combo2 ):
   print( f"{ops2[i]:3s} {angle2[i]:6.3f} {ax2[i][0]:9.4f} {ax2[i][1]:9.4f} {ax2[i][2]:9.4f}")

# get the PG
pg_combo2 = sofi.get_pg( n_combo2, matrix_combo2 )
print("Corresponding PG:", pg_combo2)
print()




# take all operations from axes 5, 6, and 7 (4, 5, 6 counting from 0)
# this should contruct back the D6h from beginning
six_ops = np.zeros( (6,3,3), dtype = float )
six_ops[0] = sofi.construct_operation( "C", ax[4], 0.5)
six_ops[1] = sofi.construct_operation( "S", ax[4], 0.0)
six_ops[2] = sofi.construct_operation( "C", ax[5], 0.5)
six_ops[3] = sofi.construct_operation( "S", ax[5], 0.0)
six_ops[4] = sofi.construct_operation( "C", ax[6], 0.5)
six_ops[5] = sofi.construct_operation( "S", ax[6], 0.0)
n_combo6, matrix_combo6 = sofi.mat_combos(6, six_ops)
pg6 = sofi.get_pg( n_combo6, matrix_combo6 )
print("combos of all operations from axes 5,6,7 give PG:",pg6)

