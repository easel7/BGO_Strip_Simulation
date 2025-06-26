import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ==== 固定参数 ====
total_strength = 10000
p1_fixed = -1 / 205.76  # λ₁ ≈ 200 mm

# ==== 模型 ====
def biexp_fixed(x, b, delta):
    p2 = p1_fixed + delta
    c = total_strength - b
    return  25.5*(- b * p1_fixed * np.exp(p1_fixed * x)  - c * p2 * np.exp(p2 * x))

# ==== 目标函数 ====
def objective_fixed(params, x, y, yerr):
    b, delta = params
    y_fit = biexp_fixed(x, b, delta)
    residual = (y - y_fit) / yerr
    return np.sum(residual**2)

# ==== 数据 ====
# # 1000 GeV 0.5
# xdata = np.array([ 12.75 , 38.25 , 63.75 , 89.25 , 114.75 , 140.25 , 165.75 , 191.25 , 216.75 , 242.25 , 267.75 , 293.25 , 318.75 , 344.25 ]) 
# ydata = np.array([ 1278.5 , 1118 , 976.5 , 807.5 , 751.5 , 635.5 , 556 , 474 , 423 , 368.5 , 324 , 272 , 263 , 229 ]) 
# yerr = np.array([ 25.2834 , 23.6432 , 22.0964 , 20.0935 , 19.3843 , 17.8255 , 16.6733 , 15.3948 , 14.543 , 13.5739 , 12.7279 , 11.6619 , 11.4673 , 10.7005 ])

# 1000 GeV 0.04
xdata = np.array([ 12.75 , 38.25 , 63.75 , 89.25 , 114.75 , 140.25 , 165.75 , 191.25 , 216.75 , 242.25 , 267.75 , 293.25 , 318.75 , 344.25 ]) 
ydata = np.array([ 1159.36 , 1076.6 , 904.28 , 790.48 , 722.52 , 627.68 , 535.76 , 466.64 , 431.28 , 380 , 341.48 , 260.96 , 280.48 , 253.84 ]) 
yerr = np.array([ 32.5756 , 31.476 , 28.7973 , 26.9901 , 25.7838 , 24.0598 , 22.2049 , 20.7422 , 19.9698 , 18.7536 , 17.7932 , 15.4942 , 16.1332 , 15.3705  ])
# ==== 初始值 & 边界 ====
initial_guess = [6000,  1/200-1/170]
bounds = [(0, total_strength), (-1e-1, -1e-6)]

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

print(f"Delta = {delta_fit:.2f}")
print(f"✅ 拟合成功：b = {b_fit:.2f} (比例 {frac_fit:.3f})")
print(f"λ₁ (固定) = {lambda1:.2f} mm")
print(f"λ₂ = {lambda2:.2f} mm")
print(f"χ² = {chi2_min:.2f}")

# ==== 拟合曲线 ====
x_fit = np.linspace(min(xdata), max(xdata), 300)
y_fit = biexp_fixed(x_fit, *result.x)

# ==== λ₂ vs b/10000 contour 图 ====
lambda2_vals = np.linspace(lambda2 * 0.8, lambda2 * 1.3, 100)
# frac_vals = np.linspace(0.3, 0.7, 100)
frac_vals = np.linspace(1e-3, 1, 100)
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

# ==== 画图 ====
fig, axs = plt.subplots(1, 2, figsize=(13, 6))

# 🎯 子图1：数据点和拟合曲线
axs[0].errorbar(xdata, ydata, yerr=yerr, fmt='o', label='MC data', capsize=3)
axs[0].plot(x_fit, y_fit, 'r-', label='Fitted')
axs[0].set_xlabel("Depth (mm)")
axs[0].set_ylabel("N_{int}")
axs[0].set_title("Fit with Fixed λ₁ and Total Strength")
axs[0].legend()
axs[0].grid(True)

# 修改 contourf：限制 χ² 范围
# vmin = chi2_min
# vmax = chi2_min + 10  # 只关注低 χ² 区域
# levels = np.linspace(vmin, vmax, 30)
# 🎯 子图2：Contour 图 (λ₂ vs b/10000)
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