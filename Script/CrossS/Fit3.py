import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ==== 固定参数 ====
total_strength = 10000
p1_fixed = -1 / 207.321  # λ₁ ≈ 200 mm

# ==== 模型 ====
def biexp_fixed(x, b, delta):
    p2 = p1_fixed + delta
    c = total_strength - b
    return b * np.exp(p1_fixed * x) + c * np.exp(p2 * x)

# ==== 目标函数 ====
def objective_fixed(params, x, y, yerr):
    b, delta = params
    y_fit = biexp_fixed(x, b, delta)
    residual = (y - y_fit) / yerr
    return np.sum(residual**2)

# ==== 数据 ====
xdata = np.array([25.5, 51, 76.5, 102, 127.5, 153, 178.5, 204,
                  229.5, 255, 280.5, 306, 331.5, 357])
ydata = np.array([8721.5, 7603.5, 6627, 5819.5, 5068, 4432.5,
                  3876.5, 3402.5, 2979.5, 2611, 2287, 2015,
                  1752, 1523])
yerr = np.array([93.389, 87.1981, 81.4064, 76.2856, 71.1899,
                 66.577, 62.2615, 58.331, 54.5848, 51.0979,
                 47.8226, 44.8888, 41.8569, 39.0256])

# ==== 初始值 & 边界 ====
initial_guess = [4000,  1/170 - 1/200]
bounds = [(0, total_strength), (1e-6, 1e-1)]

# ==== 拟合 ====
result = minimize(objective_fixed, initial_guess,
                  args=(xdata, ydata, yerr),
                  method='L-BFGS-B', bounds=bounds)

# ==== 提取拟合结果 ====
if not result.success:
    print("❌ 拟合失败：", result.message)
    exit()

b_fit, delta_fit = result.x
c_fit = total_strength - b_fit
p2_fit = p1_fixed + delta_fit
lambda1 = -1 / p1_fixed
lambda2 = -1 / p2_fit
frac_fit = b_fit / total_strength
chi2_min = objective_fixed([b_fit, delta_fit], xdata, ydata, yerr)

print(f"✅ 拟合成功：b = {b_fit:.2f} (比例 {frac_fit:.3f})")
print(f"λ₁ (固定) = {lambda1:.2f} mm")
print(f"λ₂ = {lambda2:.2f} mm")
print(f"χ² = {chi2_min:.2f}")

# ==== 拟合曲线 ====
x_fit = np.linspace(min(xdata), max(xdata), 300)
y_fit = biexp_fixed(x_fit, *result.x)

# ==== λ₂ vs b/10000 contour 图 ====
lambda2_vals = np.linspace(lambda2 * 0.8, lambda2 * 1.3, 100)
frac_vals = np.linspace(0.3, 0.7, 100)
L2, F = np.meshgrid(lambda2_vals, frac_vals)

chi2_map = np.full_like(L2, np.nan)
for i in range(L2.shape[0]):
    for j in range(L2.shape[1]):
        lam2 = L2[i, j]
        frac = F[i, j]
        b = frac * total_strength
        p2 = -1 / lam2
        delta = p2 - p1_fixed
        if delta <= 0:
            continue
        chi2_map[i, j] = objective_fixed([b, delta], xdata, ydata, yerr)

# ==== 画图 ====
fig, axs = plt.subplots(1, 2, figsize=(14, 5))

# 🎯 子图1：数据点和拟合曲线
axs[0].errorbar(xdata, ydata, yerr=yerr, fmt='o', label='MC data', capsize=3)
axs[0].plot(x_fit, y_fit, 'r-', label='Fitted')
axs[0].set_xlabel("Depth (mm)")
axs[0].set_ylabel("Counts")
axs[0].set_title("Fit with Fixed λ₁ and Total Strength")
axs[0].legend()
axs[0].grid(True)

# 🎯 子图2：Contour 图 (λ₂ vs b/10000)
vmin = chi2_min
vmax = chi2_min + 20
levels = np.linspace(vmin, vmax, 30)

cs = axs[1].contourf(L2, F, chi2_map, levels=levels, cmap='viridis', extend='both')
axs[1].contour(L2, F, chi2_map, levels=[chi2_min + 2.3, chi2_min + 6.18], colors='red', linestyles='--')
axs[1].plot(lambda2, frac_fit, 'r*', label='Best Fit')
axs[1].set_xlabel("λ₂ (mm)")
axs[1].set_ylabel("b / 10000")
axs[1].set_title("χ² Contour: λ₂ vs b/10000")
axs[1].legend()
fig.colorbar(cs, ax=axs[1], label='χ²')

plt.tight_layout()
plt.show()