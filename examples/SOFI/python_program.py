import numpy as np
import sys

m=len(sys.argv)
if m==1:
    print("================================")
    print("usage: ./test_py.py filename.xyz")
    print("================================")
    quit()

fname=sys.argv[1]


def read_xyz(fname):
  import numpy as np
  typ = np.genfromtxt( fname, skip_header=2, usecols=[0], dtype=None, encoding=None )
  coords = np.genfromtxt( fname, skip_header = 2, usecols = [1,2,3], dtype=np.float64 )
  nat = len( typ )
  cc=np.ndarray( (nat,3), dtype=np.float64, order="C" )
  cc=coords
  return nat, typ, cc


import ira_mod

nat, typ, coords = read_xyz(fname)

## recenter to chosen origin;
# in this case take the geometric center
gc = np.mean( coords, axis=0 )
# gc=coords[:][3]
coords = coords - gc

## initialize SOFI
sofi = ira_mod.SOFI()

## threshold for symmetries
sym_thr = 0.2

## compute all sofi data and store into "sym" object
sym = sofi.compute( nat, typ, coords, sym_thr )


sym.print()

