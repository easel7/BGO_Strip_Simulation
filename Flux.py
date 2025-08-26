import numpy as np
import matplotlib.pyplot as plt

# 定义函数
def f1(E):
    return 7.58e-5 * (E/1000.)**(-2.772) * (1 + (E/480.)**5)**(0.173/5.) * E**2.7

def f2(R):
    return 0.4544 * (R/45.)**(-2.849) * (1 + (R/340.)**(0.1333/0.024))**0.024 * R**2.7

def E2R(E):
    mp = 0.938  # proton mass in GeV
    return np.sqrt((E + mp)**2 - mp**2)

# 能量范围
E = np.logspace(1, 4, 500)  # 从1e1到1e5，500个点

# 画图
plt.figure(figsize=(7,5))
plt.semilogx(E, f1(E), label="DAMPE")
plt.semilogx(E, f2(E2R(E)), label="AMS02")
plt.xlabel("E (GeV) ")
plt.ylabel("Flux (GeV * m * sr * s)^-1 * E^2.6")
plt.ylim(2000, 20000)  # 设置 y 轴范围
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.show()