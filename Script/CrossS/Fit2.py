import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# 模型：无 a 偏移量
def biexp_model(x, b, p1, delta):
    p2 = p1 + delta
    c = 10000 - b
    return b * np.exp(p1 * x) + c * np.exp(p2 * x)

# 加权残差平方和
def objective(params, x, y, yerr):
    b, p1, delta = params
    y_fit = biexp_model(x, b, p1, delta)
    residual = (y - y_fit) / yerr
    return np.sum(residual**2)

# 数据
xdata = np.array([25.5, 51, 76.5, 102, 127.5, 153, 178.5, 204, 229.5, 255, 280.5, 306, 331.5, 357])
ydata = np.array([8721.5, 7603.5, 6627, 5819.5, 5068, 4432.5, 3876.5, 3402.5, 2979.5, 2611, 2287, 2015, 1752, 1523])
yerr  = np.array([93.389, 87.1981, 81.4064, 76.2856, 71.1899, 66.577, 62.2615, 58.331, 54.5848, 51.0979, 47.8226, 44.8888, 41.8569, 39.0256])

# 初始值（不再有 a）
initial_guess = [5000, -1/200, 1/114 - 1/200]  # b, p1, delta

# 边界
bounds = [
    (0, 10000),       # b
    (-1e-1, -1e-5),   # p1 < 0
    (1e-6, 1e-1),     # delta > 0 ⇒ p2 > p1 ⇒ λ1 < λ2
]

# 拟合
result = minimize(objective, initial_guess, args=(xdata, ydata, yerr), method='L-BFGS-B', bounds=bounds)

# 输出
if result.success:
    b_fit, p1_fit, delta_fit = result.x
    c_fit = 10000 - b_fit
    p2_fit = p1_fit + delta_fit
    print("✅ 拟合成功（无偏移量 a）：")
    print(f"b = {b_fit:.2f}")
    print(f"c = {c_fit:.2f}")
    print(f"p1 = {p1_fit:.6f}    ⇒ λ1 = {-1/p1_fit:.2f}")
    print(f"p2 = {p2_fit:.6f}    ⇒ λ2 = {-1/p2_fit:.2f}")
else:
    print("❌ 拟合失败：", result.message)

# 画图
x_fit = np.linspace(min(xdata), max(xdata), 300)
y_fit = biexp_model(x_fit, *result.x)

plt.errorbar(xdata, ydata, yerr=yerr, fmt='o', label='MC with error', capsize=3)
plt.plot(x_fit, y_fit, 'r-', label='Fitted (no offset)')
plt.xlabel("Depth (mm)")
plt.ylabel("Counts")
plt.title("Double Exponential Fit (No Offset Term)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()