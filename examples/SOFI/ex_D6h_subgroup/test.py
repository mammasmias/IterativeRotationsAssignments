import ira_mod

sofi=ira_mod.SOFI('../../../src/libira.so')

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
d6h=ira_mod.sym_data()
d6h.n_sym = 24
d6h.matrix = np.zeros( (24, 3, 3), dtype = float )

# fill the theta operations

# identity
d6h.matrix[0] = sofi.construct_operation( "Id", ax[0], 0.0 )
# inversion
d6h.matrix[1] = sofi.construct_operation( "I", ax[0], 0.0 )

m=2
for i in range(7):
   # construct S0 on all axes
   theta = sofi.construct_operation( "S", ax[i], 0.0 )
   d6h.matrix[m] = theta
   m += 1
   # construct C2 on all axes
   theta = sofi.construct_operation( "C", ax[i], 0.5 )
   d6h.matrix[m] = theta
   m += 1

# fill the rest on ax[0]
d6h.matrix[m] = sofi.construct_operation( "C", ax[0], 1/3.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "C", ax[0], -1/3.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "C", ax[0], 1/6.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "C", ax[0], -1/6.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "S", ax[0], 1/3.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "S", ax[0], -1/3.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "S", ax[0], 1/6.0 )
m+=1
d6h.matrix[m] = sofi.construct_operation( "S", ax[0], -1/6.0 )

# confirm we have D6h
d6h.pg = sofi.get_pg( d6h.n_sym, d6h.matrix )
print( "initial PG is:", d6h.pg )


# take all operations from axes 5 and 7 (which is 4 and 6 counting from 0)
four_ops1 = np.zeros( (4,3,3), dtype = float)
four_ops1[0] = sofi.construct_operation( "S", ax[4], 0.0)
four_ops1[1] = sofi.construct_operation( "C", ax[4], 0.5)
four_ops1[2] = sofi.construct_operation( "S", ax[6], 0.0)
four_ops1[3] = sofi.construct_operation( "C", ax[6], 0.5)

# construct matrix combinations
n_combo1, matrix_combo1 = sofi.mat_combos(4, four_ops1 )
# see which pg we get
pg_combo1 = sofi.get_pg( n_combo1, matrix_combo1 )
print( "combinations of matrices from axes 5 and 7 give PG:", pg_combo1 )


# take all operations from axes 6 and 7 (which is 5 and 6 counting from 0)
four_ops2 = np.zeros( (4,3,3), dtype = float)
four_ops2[0] = sofi.construct_operation( "S", ax[5], 0.0)
four_ops2[1] = sofi.construct_operation( "C", ax[5], 0.5)
four_ops2[2] = sofi.construct_operation( "S", ax[6], 0.0)
four_ops2[3] = sofi.construct_operation( "C", ax[6], 0.5)

# construct matrix combinations
n_combo2, matrix_combo2 = sofi.mat_combos(4, four_ops2)
# see which pg we get
pg_combo2 = sofi.get_pg( n_combo2, matrix_combo2 )
print( "combinations of matrices from axes 6 and 7 give PG:", pg_combo2 )


six_ops = np.zeros( (6,3,3), dtype = float )

six_ops[0] = sofi.construct_operation( "C", ax[4], 0.5)
six_ops[1] = sofi.construct_operation( "S", ax[4], 0.0)
six_ops[2] = sofi.construct_operation( "C", ax[5], 0.5)
six_ops[3] = sofi.construct_operation( "S", ax[5], 0.0)
six_ops[4] = sofi.construct_operation( "C", ax[6], 0.5)
six_ops[5] = sofi.construct_operation( "S", ax[6], 0.0)
n_combo6, matrix_combo6 = sofi.mat_combos(6, six_ops)
pg6 = sofi.get_pg( n_combo6, matrix_combo6 )
print("combos of axes 5,6,7:",pg6)

