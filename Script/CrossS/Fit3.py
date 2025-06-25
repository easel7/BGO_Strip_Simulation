import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# 固定参数
total_strength = 10000
p1_fixed = -1 / 207.321  # lambda_1 = 200 mm

# 模型（只有两个自由参数 b, delta）
def biexp_fixed(x, b, delta):
    p2 = p1_fixed + delta
    c = total_strength - b
    return b * np.exp(p1_fixed * x) + c * np.exp(p2 * x)

# 目标函数（最小化 chi2）
def objective_fixed(params, x, y, yerr):
    b, delta = params
    y_fit = biexp_fixed(x, b, delta)
    residual = (y - y_fit) / yerr
    return np.sum(residual ** 2)

# 数据
xdata = np.array([25.5, 51, 76.5, 102, 127.5, 153, 178.5, 204,
                  229.5, 255, 280.5, 306, 331.5, 357])
ydata = np.array([8721.5, 7603.5, 6627, 5819.5, 5068, 4432.5,
                  3876.5, 3402.5, 2979.5, 2611, 2287, 2015,
                  1752, 1523])
yerr = np.array([93.389, 87.1981, 81.4064, 76.2856, 71.1899,
                 66.577, 62.2615, 58.331, 54.5848, 51.0979,
                 47.8226, 44.8888, 41.8569, 39.0256])

# 初始猜测
initial_guess = [5000, 1/114 - 1/200]  # b, delta

# 边界
bounds = [
    (0, total_strength),     # b
    (1e-6, 1e-1),             # delta > 0 ⇒ p2 > p1 ⇒ λ1 < λ2
]

# 拟合
result = minimize(objective_fixed, initial_guess,
                  args=(xdata, ydata, yerr),
                  method='L-BFGS-B', bounds=bounds)

# 输出结果
if result.success:
    b_fit, delta_fit = result.x
    c_fit = total_strength - b_fit
    p2_fit = p1_fixed + delta_fit
    print("✅ 拟合成功（固定总强度=10000，p1=-1/200）：")
    print(f"b = {b_fit:.2f}")
    print(f"c = {c_fit:.2f}")
    print(f"p1 (固定) = {p1_fixed:.6f}   ⇒ λ1 = {-1/p1_fixed:.2f}")
    print(f"p2 = {p2_fit:.6f}            ⇒ λ2 = {-1/p2_fit:.2f}")
else:
    print("❌ 拟合失败：", result.message)

# 绘图
x_fit = np.linspace(min(xdata), max(xdata), 300)
y_fit = biexp_fixed(x_fit, *result.x)

plt.errorbar(xdata, ydata, yerr=yerr, fmt='o', label='MC + error', capsize=3)
plt.plot(x_fit, y_fit, 'r-', label='Fitted (fixed total & p1)')
plt.xlabel("Depth (mm)")
plt.ylabel("Counts")
plt.title("Double Exponential Fit (Fixed Total & λ₁)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()