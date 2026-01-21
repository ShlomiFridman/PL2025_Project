import time
import numpy as np
from scipy.linalg import expm

mat = np.loadtxt("exp_data.csv", delimiter=",")

start_time = time.time()
res_mat = expm(mat)
total_time = time.time() - start_time

print(f"Input sub-matrix (5x5):\n{mat[:5,:5]}\n")
print(f"Total time: {total_time:.3f} seconds\n")
print(f"Result sub-matrix (5x5) exponential:\n{res_mat[:5,:5]}")