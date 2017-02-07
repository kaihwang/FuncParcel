import sys
import FuncParcel
import numpy as np

# script to iterate modularity detection by calling function in FuncParcel

subject = sys.stdin.read().strip('\n')
print(subject)
print(type(subject))
ci, q = FuncParcel.iterate_modularity_partition(subject, 100)
fn = '/home/despoB/kaihwang/bin/FuncParcel/Data/Subject_Partition/%s_ci' %subject
ci = ci +1
np.savetxt(fn, ci, fmt='%3.d') 