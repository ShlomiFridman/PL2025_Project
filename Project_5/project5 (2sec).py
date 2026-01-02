import numpy as np
import time

epsilon = 1e-8

mat = np.loadtxt("exp_data.csv", delimiter=",")

n = mat.shape[0]
    
def calc_vec_norm(vec: np.array) -> float:
    return np.sqrt(np.sum(np.square(vec)))

def find_k(A: np.array) -> tuple[int, np.array]:
    k=1
    b_curr = np.ones(n)
    b_next = np.zeros(n)
    while True:
        b_next = A @ b_curr
        b_next /= calc_vec_norm(b_next)
        
        if calc_vec_norm(b_next - b_curr) < epsilon:
            return k, b_curr
        else:
            k += 1
            b_curr = b_next
        
def find_m(A: np.array, bk: np.array) -> int:
    bk_t = bk.T
    mu_org = (bk_t @ A @ bk) / (np.dot(bk_t, bk))
    mu_acc = mu_org
    m = 1
    m_fact = 1
    while (mu_acc / m_fact) >= epsilon:
         mu_acc *= mu_org
         m += 1
         m_fact *= m
         
    return m

def calc_exp_mat(A: np.array, m: int) -> np.array:
    A_acc = np.eye(n)
    k_fact = 1
    res_mat = np.eye(n)
    for i in range(1, m+1):
        A_acc @= A
        k_fact *= i
        res_mat += A_acc/k_fact
    return res_mat

def run_prog(A: np.array) -> tuple[np.array, float]:
    startTime = time.time()
    
    scale_pow = 16
    A_scaled = A/(2**scale_pow)
    
    k,bk = find_k(A_scaled)
    m = find_m(A_scaled, bk)
    res_mat = calc_exp_mat(A_scaled, m)
    
    for i in range(scale_pow):
        res_mat = res_mat @ res_mat
    
    endTime = time.time()
    return res_mat, (endTime-startTime)


res_mat, total_time = run_prog(mat)
print(f"Total time: {total_time:.3f} seconds")
print("Result sub-matrix (5x5):")
print(res_mat[:5,:5])
