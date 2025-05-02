import ROOT
import numpy as np
import matplotlib.pyplot as plt

# 载入 ROOT 文件
filename = "/Users/xiongzheng/software/B4/B4e/ROOT/Proton_1000GeV.root"
file = ROOT.open(filename)
tree = file["B4"]

# 读取分支
energy_vecs = tree["BarEnergyVector"].arrays(library="np")["BarEnergyVector"]
rms_vecs    = tree["RMS"].arrays(library="np")["RMS"]
fh_depth    = tree["First_Had_Depth"].array(library="np")
fh_layer    = tree["First_Had_Layer"].array(library="np")
fh_type     = tree["First_Had_Type"].array(library="np")
nhits       = tree["Nhits"].array(library="np")

# 直方图容器
h_max_min0 = []
h_max_min1 = []
h_max_min2 = []
h_max_min3 = []

# 用于 scatter plot 的容器
sum_len0 = []
sum_len1 = []
sum_len2 = []
sum_len3 = []

# 遍历事件进行分析
for i in range(len(energy_vecs)):
    energy = energy_vecs[i]
    if len(energy) < 1:
        continue

    bar_change = np.diff(energy)  # 模拟“增长段”
    positive_indices = np.where(bar_change > 0)[0]

    if len(positive_indices) == 0:
        continue

    growth_start = positive_indices[0]
    growth_sum = np.sum(bar_change[positive_indices])
    growth_len = len(positive_indices)

    max_val = np.max(bar_change[positive_indices])
    min_val = np.min(bar_change[positive_indices])
    delta_max_min = max_val - min_val

    # 分类填入
    t = fh_type[i]
    if t == 0:
        h_max_min0.append(delta_max_min)
        sum_len0.append((growth_len, growth_sum))
    elif t == 1:
        h_max_min1.append(delta_max_min)
        sum_len1.append((growth_len, growth_sum))
    elif t == 2:
        h_max_min2.append(delta_max_min)
        sum_len2.append((growth_len, growth_sum))
    else:
        h_max_min3.append(delta_max_min)
        sum_len3.append((growth_len, growth_sum))

# 示例绘图
plt.hist(h_max_min0, bins=60, range=(-1, 5), alpha=0.5, label="Type 0")
plt.hist(h_max_min1, bins=60, range=(-1, 5), alpha=0.5, label="Inelastic")
plt.hist(h_max_min2, bins=60, range=(-1, 5), alpha=0.5, label="Elastic")
plt.hist(h_max_min3, bins=60, range=(-1, 5), alpha=0.5, label="Pass")
plt.legend()
plt.title("Δ(Max - Min) in Growth Region")
plt.show()