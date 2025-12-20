import numpy as np
import math
import time

def get_pade_coefficients(q):
    """
    Computes coefficients for diagonal Pade approximant of order q.
    """
    coeffs = []
    factorial_2q = math.factorial(2 * q)
    factorial_q = math.factorial(q)
    
    for j in range(q + 1):
        num = math.factorial(2 * q - j) * factorial_q
        denom = factorial_2q * math.factorial(j) * math.factorial(q - j)
        coeffs.append(num / denom)
    return coeffs

def gauss_jordan_solve(A, B):
    """
    Solves AX = B using Gauss-Jordan elimination.
    Replaces np.linalg.solve.
    """
    # Work on copies to avoid modifying originals
    A_sys = A.astype(float).copy()
    B_sys = B.astype(float).copy()
    n = A_sys.shape[0]
    
    # 1. Construct Augmented Matrix [A | B]
    # We stack them horizontally
    aug = np.hstack([A_sys, B_sys])
    
    # 2. Gaussian Elimination
    for i in range(n):
        # Pivot: Find row with largest value in current column i
        pivot_row = np.argmax(np.abs(aug[i:, i])) + i
        
        # Swap current row with pivot row
        if i != pivot_row:
            aug[[i, pivot_row]] = aug[[pivot_row, i]]
            
        pivot_val = aug[i, i]
        
        # Check for singularity
        if abs(pivot_val) < 1e-15:
            raise ValueError("Matrix is singular, cannot invert.")
            
        # Normalize the pivot row
        aug[i] = aug[i] / pivot_val
        
        # Eliminate all other rows
        for j in range(n):
            if i != j:
                factor = aug[j, i]
                aug[j] -= factor * aug[i]
                
    # 3. Extract Solution
    # The right half of the augmented matrix is now A^(-1)B, which is X
    X = aug[:, n:]
    return X

def pade_approximant(A, q=6):
    n = A.shape[0]
    N = np.zeros((n, n), dtype=float) 
    D = np.zeros((n, n), dtype=float)
    I = np.eye(n, dtype=float)
    
    coeffs = get_pade_coefficients(q)
    
    A_pow = I.copy()
    for j in range(q + 1):
        c_j = coeffs[j]
        term = c_j * A_pow
        N += term
        # Denominator terms alternate signs
        D += term * ((-1)**j)
        
        if j < q:
            A_pow = np.dot(A_pow, A)
            
    # Solve D * R = N
    return gauss_jordan_solve(D, N)

def matrix_exponential(file_path, epsilon, q=6):
    # Load the matrix
    try:
        A = np.loadtxt(file_path, delimiter=',')
    except Exception as e:
        print(f"Error: {e}")
        return None

    # Calculate Infinity Norm (max row sum)
    A_norm = np.max(np.sum(np.abs(A), axis=1))
    
    # --- Step 1: Scaling ---
    # We scale until norm <= epsilon
    if A_norm <= epsilon:
        s = 0
    else:
        # 2^s >= norm / epsilon  ->  s >= log2(norm/epsilon)
        s = int(np.ceil(np.log2(A_norm / epsilon)))
    
    A_scaled = A / (2**s)
    
    # --- Step 2: Pade Approximation ---
    R = pade_approximant(A_scaled, q)
    
    # --- Step 3: Squaring ---
    F = R
    for _ in range(s):
        F = np.dot(F, F)
        
    return F

# --- Example Usage: Quantum Hamiltonian Simulation ---

# 3. Compute using our custom Pade function
start_time = time.perf_counter()
result = matrix_exponential('exp_data.csv', epsilon=1e5)
end_time = time.perf_counter()
duration = end_time - start_time

print(f"\nFinal Result ({duration:.4f} seconds):\n", result)