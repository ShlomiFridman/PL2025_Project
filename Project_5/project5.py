import numpy as np
import time

epsilon = 1e-6
scale_pow = 16

mat = np.loadtxt("exp_data.csv", delimiter=",")

n = mat.shape[0]

def calc_vec_norm(vec):
    return np.sqrt(np.dot(vec, vec))

def find_bk(A):
    b_curr = np.zeros(n)
    b_next = np.ones(n)
    while calc_vec_norm(b_next - b_curr) >= epsilon:
        b_curr = b_next
        b_next = A @ b_curr
        b_next /= calc_vec_norm(b_next)
    return b_curr

def calc_mat_norm(A, bk):
    bk_t = bk.T
    mu_org = (bk_t @ A @ bk) / (bk_t @ bk)
    return mu_org

def find_m(A, bk):
    mu_org = calc_mat_norm(A, bk)
    mu_acc = mu_org
    m = 1
    m_fact = 1

    while (mu_acc / m_fact) >= epsilon:
         mu_acc *= mu_org
         m += 1
         m_fact *= m

    return m

def calc_exp_mat(A, m):
    A_acc = np.eye(n, dtype=A.dtype)
    res_mat = np.eye(n, dtype=A.dtype)
    k_fact = 1

    for i in range(1, m+1):
        A_acc = A_acc @ A
        k_fact *= i
        res_mat += A_acc/k_fact
    return res_mat

def run_prog(A):
    startTime = time.time()
    
    A_scaled = A/(2**scale_pow)
    
    bk = find_bk(A_scaled)
    m = find_m(A_scaled, bk)
    res_mat = calc_exp_mat(A_scaled, m)
    
    for i in range(scale_pow):
        res_mat = res_mat @ res_mat
    
    endTime = time.time()
    return res_mat, (endTime-startTime)


res_mat, total_time = run_prog(mat)
print(f"Before:\n{mat[:5,:5]}\n")
print(f"Total time: {total_time:.3f} seconds\n")
print(res_mat[:5,:5])
