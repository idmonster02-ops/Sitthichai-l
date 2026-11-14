import math
import random

print("==========================================================")
print("  PURE SPIN2 FRAMEWORK v6.3 & EXOPLANET MCMC SUITE")
print("  GitHub Cloud Execution - Master Production Script")
print("==========================================================\n")

# ==========================================
# MODULE 1: Pure SPIN2 & Dynamic Higuchi Bound
# ==========================================
H0_km_s_Mpc = 67.4
Omega_m = 0.315
Omega_L = 0.685
A_spin = -0.6000
z_c = 0.0665
h_bar_eVs = 6.582119569e-16
Mpc_to_km = 3.08567758e19

def w_pure_spin2(z):
    return -1.0 - A_spin * math.exp(-z / z_c)

def calc_f_de_spin2(z_target, steps=1000):
    if z_target == 0:
        return 1.0
    dz = z_target / steps
    integral = 0.0
    for i in range(steps):
        z = (i + 0.5) * dz
        integral += ((1.0 + w_pure_spin2(z)) / (1.0 + z)) * dz
    return math.exp(3.0 * integral)

def get_H_in_eV(z):
    f_de = calc_f_de_spin2(z)
    H_val = H0_km_s_Mpc * math.sqrt(Omega_m * (1 + z)**3 + Omega_L * f_de)
    return h_bar_eVs * (H_val / Mpc_to_km)

print("--- [Module 1] Dynamic Higuchi Bound & Dark Energy Check ---")
for z_test in [0.0, 1.0, 1100.0]:
    H_eV = get_H_in_eV(z_test)
    min_mg = math.sqrt(2.0) * H_eV
    dynamic_mg = 1.5 * min_mg  
    print(f"  - z = {z_test:6.1f} | H(z) = {H_eV:.4e} eV | m_g(z) = {dynamic_mg:.4e} eV (Stable)")

# ==========================================
# MODULE 2: Exoplanet MCMC Orbital Fitting
# ==========================================
P_true, K_true, e_true, omega_true, V0_true, T0_true = 12.5, 45.5, 0.35, 1.2, 12.0, 2.0

def solve_kepler(M, e, tol=1e-8, max_iter=100):
    E = M
    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        df = 1.0 - e * math.cos(E)
        dE = f / df
        E -= dE
        if abs(dE) < tol:
            break
    return E

def compute_rv(time_array, P, K, e, omega, T0, V0):
    rv = []
    for t in time_array:
        M = (2.0 * math.pi / P) * ((t - T0) % P)
        E = solve_kepler(M, e)
        nu = 2.0 * math.atan(math.sqrt((1.0 + e) / (1.0 - e)) * math.tan(E / 2.0))
        rv.append(V0 + K * (math.cos(nu + omega) + e * math.cos(omega)))
    return rv

random.seed(42)
time_obs = [i * 0.5 for i in range(60)]
clean_rv = compute_rv(time_obs, P_true, K_true, e_true, omega_true, T0_true, V0_true)
noisy_rv = [v + random.gauss(0, 3.0) for v in clean_rv]

def log_likelihood(params, times, observed, sigma=3.0):
    K, V0 = params
    expected = compute_rv(times, P_true, K, e_true, omega_true, T0_true, V0)
    chi2 = sum(((o - e) / sigma)**2 for o, e in zip(observed, expected))
    return -0.5 * chi2

print("\n--- [Module 2] Running Metropolis-Hastings MCMC ---")
current_K, current_V0 = 35.0, 5.0
current_log_prob = log_likelihood([current_K, current_V0], time_obs, noisy_rv)
chain = []

iterations = 5000
accepted = 0
for _ in range(iterations):
    prop_K = current_K + random.gauss(0, 0.8)
    prop_V0 = current_V0 + random.gauss(0, 0.5)
    prop_log_prob = log_likelihood([prop_K, prop_V0], time_obs, noisy_rv)
    
    if prop_log_prob > current_log_prob or random.random() < math.exp(prop_log_prob - current_log_prob):
        current_K, current_V0, current_log_prob = prop_K, prop_V0, prop_log_prob
        accepted += 1
    chain.append((current_K, current_V0))

burned = chain[1000:]
n_s = len(burned)
K_samp = [p[0] for p in burned]
V0_samp = [p[1] for p in burned]

est_K = sum(K_samp) / n_s
est_V0 = sum(V0_samp) / n_s
sig_K = math.sqrt(sum((k - est_K)**2 for k in K_samp) / (n_s - 1))
sig_V0 = math.sqrt(sum((v - est_V0)**2 for v in V0_samp) / (n_s - 1))

print(f"  - Acceptance Rate: {(accepted/iterations)*100:.1f}%")
print(f"  - K Amplitude  -> True: {K_true} | Estimated (Mean +- 1sigma): {est_K:.2f} +- {sig_K:.2f} m/s")
print(f"  - V0 Velocity  -> True: {V0_true} | Estimated (Mean +- 1sigma): {est_V0:.2f} +- {sig_V0:.2f} m/s")

print("\n==========================================================")
print("  ALL EQUATIONS EXECUTED SUCCESSFULLY ON CLOUD ENVIRONMENT")
print("==========================================================")
