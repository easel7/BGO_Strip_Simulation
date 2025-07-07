import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ==== Biexponential ====
def biexp_model(x, b, p1, delta):
    p2 = p1 + delta
    c = 10000 - b
    return 25.5 * (b * (-p1) * np.exp(p1 * x) + c * (-p2) * np.exp(p2 * x))

def objective(params, x, y, yerr):
    b, p1, delta = params
    y_fit = biexp_model(x, b, p1, delta)
    residual = (y - y_fit) / yerr
    return np.sum(residual**2)

# ==== Data ====
xdata = np.array([ 12.75 , 38.25 , 63.75 , 89.25 , 114.75 , 140.25 , 165.75 , 191.25 , 216.75 , 242.25 , 267.75 , 293.25 , 318.75 , 344.25 ]) 
ydata = np.array([ 1278.5 , 1118 , 976.5 , 807.5 , 751.5 , 635.5 , 556 , 474 , 423 , 368.5 , 324 , 272 , 263 , 229 ]) 
yerr = np.array([ 25.2834 , 23.6432 , 22.0964 , 20.0935 , 19.3843 , 17.8255 , 16.6733 , 15.3948 , 14.543 , 13.5739 , 12.7279 , 11.6619 , 11.4673 , 10.7005 ])

# ==== Fit ====
# b, p1, delta
initial_guess = [5000, -1/200, 1/200-1/170]
bounds = [(0, 10000), (-1e-1, -1e-5), (-1e-1, -1e-10)]
result = minimize(objective, initial_guess, args=(xdata, ydata, yerr), method='L-BFGS-B', bounds=bounds)

# ==== Fit Result ====
if not result.success:
    print("❌ Fail：", result.message)
    exit()

b_fit, p1_fit, delta_fit = result.x
p2_fit = p1_fit + delta_fit
lambda1 = -1 / p1_fit
lambda2 = -1 / p2_fit
chi2_min = objective([b_fit, p1_fit, delta_fit], xdata, ydata, yerr)
frac_fit = b_fit / 10000

print(f"✅ Success：b = {b_fit:.2f} (Ratio {frac_fit:.3f}), λp = {lambda1:.2f}, λd = {lambda2:.2f}, χ^2 = {chi2_min:.2f}")

# ==== Prepare ====
x_fit = np.linspace(min(xdata), max(xdata), 300)
y_fit = biexp_model(x_fit, b_fit, p1_fit, delta_fit)

# ==== 1. λp vs λd  ====
λ1_vals = np.linspace(lambda1 * 0.7, lambda1 * 1.3, 100)
λ2_vals = np.linspace(lambda2 * 0.7, lambda2 * 1.3, 100)
L1, L2 = np.meshgrid(λ1_vals, λ2_vals)

chi2_l1l2 = np.full_like(L1, np.nan)
for i in range(L1.shape[0]):
    for j in range(L1.shape[1]):
        l1 = L1[i, j]
        l2 = L2[i, j]
        # if l2 <= l1:
        #     continue
        p1 = -1 / l1
        p2 = -1 / l2
        delta = p2 - p1
        chi2_l1l2[i, j] = objective([b_fit, p1, delta], xdata, ydata, yerr)

# ==== 2. λp vs b/10000  ====
frac_vals = np.linspace(0.3, 0.7, 100)   # b/10000 的Ratio
# frac_vals = np.linspace(1e-3, 1, 100)
L1_, F = np.meshgrid(λ1_vals, frac_vals)

chi2_l1f = np.full_like(L1_, np.nan)
for i in range(L1_.shape[0]):
    for j in range(L1_.shape[1]):
        l1 = L1_[i, j]
        frac = F[i, j]
        b = frac * 10000
        if b <= 0 or b >= 10000:
            continue
        p1 = -1 / l1
        delta = p2_fit - p1  # Fix p2
        # if delta <= 0:
        #     continue
        chi2_l1f[i, j] = objective([b, p1, delta], xdata, ydata, yerr)

# ==== Plot ====
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

# 1️⃣ Data & Fit
axs[0].errorbar(xdata, ydata, yerr=yerr, fmt='o', label='MC', capsize=3)
axs[0].plot(x_fit, y_fit, 'r-', label='Fit')
axs[0].set_title("Data + Fit")
axs[0].set_xlabel("Depth (mm)")
axs[0].set_ylabel("$N_{int}$")
axs[0].legend()
axs[0].grid(True)

# 2️⃣ λp vs λd Contour
cs1 = axs[1].contourf(L1, L2, chi2_l1l2, levels=30, cmap='viridis')
axs[1].contour(L1, L2, chi2_l1l2, levels=[chi2_min + 2.3, chi2_min + 6.18], colors='red', linestyles='--')
axs[1].plot(lambda1, lambda2, 'r*', label='Best Fit')
axs[1].set_xlabel("λp (mm)")
axs[1].set_ylabel("λd (mm)")
axs[1].set_title("χ^2 vs λp, λd")
axs[1].legend()
fig.colorbar(cs1, ax=axs[1], label='χ^2')

# 3️⃣ λp vs b/10000 Contour
cs2 = axs[2].contourf(L1_, F, chi2_l1f, levels=30, cmap='plasma')
axs[2].contour(L1_, F, chi2_l1f, levels=[chi2_min + 2.3, chi2_min + 6.18], colors='blue', linestyles='--')
axs[2].plot(lambda1, frac_fit, 'r*', label='Best Fit')
axs[2].set_xlabel("λp (mm)")
# axs[2].set_yscale('log')
axs[2].set_ylabel("Np/10000")
axs[2].set_title("χ^2 vs λp, Np / 10000 = r_p")
axs[2].legend()
fig.colorbar(cs2, ax=axs[2], label='χ^2')

plt.tight_layout()
plt.show()