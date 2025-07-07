import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ==== fix parameters ====
total_strength = 10000
p1_fixed = -1 / 207.321  # λ₁ ≈ 200 mm

# ==== Biexponential ====
def biexp_fixed(x, b, delta):
    p2 = p1_fixed + delta
    c = total_strength - b
    return b * np.exp(p1_fixed * x) + c * np.exp(p2 * x)

# ==== chi 2 function  ====
def objective_fixed(params, x, y, yerr):
    b, delta = params
    y_fit = biexp_fixed(x, b, delta)
    residual = (y - y_fit) / yerr
    return np.sum(residual**2)

# ==== Data ====
# 1000 GeV 0.5
# xdata = np.array([25.5, 51, 76.5, 102, 127.5, 153, 178.5, 204,
#                   229.5, 255, 280.5, 306, 331.5, 357])
# ydata = np.array([8721.5, 7603.5, 6627, 5819.5, 5068, 4432.5,
#                   3876.5, 3402.5, 2979.5, 2611, 2287, 2015,
#                   1752, 1523])
# yerr = np.array([93.389, 87.1981, 81.4064, 76.2856, 71.1899,
#                  66.577, 62.2615, 58.331, 54.5848, 51.0979,
#                  47.8226, 44.8888, 41.8569, 39.0256])

# 1000 GeV 0.04
xdata = np.array([
    25.5, 51, 76.5, 102, 127.5, 153, 178.5, 204,
    229.5, 255, 280.5, 306, 331.5, 357
])

ydata = np.array([
    8840.64, 7764.04, 6859.76, 6069.28, 5346.76, 4719.08,
    4183.32, 3716.68, 3285.4, 2905.4, 2563.92,
    2302.96, 2022.48, 1768.64
])

yerr = np.array([
    94.0247, 88.1138, 82.8237, 77.9056, 73.1215, 68.6956,
    64.6786, 60.9646, 57.3184, 53.9018, 50.6352,
    47.9892, 44.9720, 42.0552
])

# ====Inital &  Boundary ====
initial_guess = [4000,  1/200-1/170]
bounds = [(0, total_strength), (-1e-1, -1e-6)]

# ==== Fit ====
result = minimize(objective_fixed, initial_guess,
                  args=(xdata, ydata, yerr),
                  method='L-BFGS-B', bounds=bounds)

# ==== Fit Result ====
if not result.success:
    print("❌ Fail：", result.message)
    exit()

b_fit, delta_fit = result.x
c_fit = total_strength - b_fit
p2_fit = p1_fixed + delta_fit
lambda1 = -1 / p1_fixed
lambda2 = -1 / p2_fit
frac_fit = b_fit / total_strength
chi2_min = objective_fixed([b_fit, delta_fit], xdata, ydata, yerr)

print(f"✅ Success：b = {b_fit:.2f} (Ratio {frac_fit:.3f})")
print(f"λ₁ (Fix) = {lambda1:.2f} mm")
print(f"λ₂ = {lambda2:.2f} mm")
print(f"χ² = {chi2_min:.2f}")

# ==== Fit ====
x_fit = np.linspace(min(xdata), max(xdata), 300)
y_fit = biexp_fixed(x_fit, *result.x)

# ==== λ₂ vs b/10000 contour  ====
lambda2_vals = np.linspace(lambda2 * 0.8, lambda2 * 1.3, 100)
# frac_vals = np.linspace(0.3, 0.7, 100)
frac_vals = np.linspace(1e-3, 0.9, 100)
L2, F = np.meshgrid(lambda2_vals, frac_vals)

chi2_map = np.full_like(L2, np.nan)
for i in range(L2.shape[0]):
    for j in range(L2.shape[1]):
        lam2 = L2[i, j]
        frac = F[i, j]
        b = frac * total_strength
        p2 = -1 / lam2
        delta = p2 - p1_fixed
        # if delta <= 0:
        #     continue
        chi2_map[i, j] = objective_fixed([b, delta], xdata, ydata, yerr)

# ==== Plot ====
fig, axs = plt.subplots(1, 2, figsize=(14, 5))

# 🎯 Subplot 1：Data & Fit
axs[0].errorbar(xdata, ydata, yerr=yerr, fmt='o', label='MC data', capsize=3)
axs[0].plot(x_fit, y_fit, 'r-', label='Fitted')
axs[0].set_xlabel("Depth (mm)")
axs[0].set_ylabel("Counts")
axs[0].set_title("Fit with Fixed λ₁ and Total Strength")
axs[0].legend()
axs[0].grid(True)

# 🎯 Subplot 2：Contour  (λ₂ vs b/10000)
cs = axs[1].contourf(L2, F, chi2_map, levels=30, cmap='viridis', extend='both')
axs[1].contour(L2, F, chi2_map, levels=[chi2_min + 2.3, chi2_min + 6.18], colors='red', linestyles='--')
axs[1].plot(lambda2, frac_fit, 'r*', label='Best Fit')
axs[1].set_xlabel("λ₂ (mm)")
axs[1].set_yscale('log')
axs[1].set_ylabel("Np/10000 (log-scale)")
axs[1].set_title("χ² Contour: λ₂ vs Np/10000= r_p")
axs[1].legend()
fig.colorbar(cs, ax=axs[1], label='χ²')

plt.tight_layout()
plt.show()