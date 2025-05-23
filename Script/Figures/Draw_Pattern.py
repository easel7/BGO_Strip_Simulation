import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec

import math  
import os


# Find the most longest increasing segment, find where it start and it value, and how long it is 
def find_max_positive_bin_segment(bar_info):
    max_sum = 0
    max_len = 0
    max_start_bin = -1

    curr_sum = 0
    curr_len = 0
    curr_start_bin = -1

    for i in range(len(bar_info)):
        content = bar_info[i]

        if content > 0:
            if curr_len == 0:
                curr_start_bin = i  

            curr_sum += content
            curr_len += 1

            if curr_sum > max_sum or (curr_sum == max_sum and curr_len > max_len):
                max_sum = curr_sum
                max_len = curr_len
                max_start_bin = curr_start_bin
        else:
            curr_sum = 0
            curr_len = 0
            curr_start_bin = -1

    return max_sum, max_len, max_start_bin

# Use ar_Change_info Array to Find the largest value and it index
def find_max_value_in_positive_segment(bar_info, start_bin, length):
    max_value = -1e9 
    max_bin = -1

    for i in range(start_bin, start_bin + length):
        content = bar_info[i]
        if content > max_value:
            max_value = content
            max_bin = i

    return max_value, max_bin

# Modified Sigmoid function
def modified_sigmoid(x, E0, Ymin, Ymax, Xmid, Slope):
    return Ymin + E0 * x + (Ymax - Ymin - E0 * x) / (1 + np.exp(-(x - Xmid) / Slope))

def Mod_Sigmoid_Percentile(x, Xmid, Slope):
    return 1 / (1 + np.exp(-(x - Xmid) / Slope))

# Import ROOT
string2 = "Proton_1000GeV"
file_path = f"/Users/xiongzheng/software/B4/B4e/Root/{string2}.root" ### Change your path!!!
with uproot.open(file_path) as f:
    tree = f["B4"]
    # Import Brancch
    data = tree.arrays([
        "BarEnergyVector", 
        "RMS", 
        "First_Had_Depth", 
        "First_Had_Layer", 
        "First_Ine_Depth", 
        "First_Ine_Layer", 
        "First_Had_Type", 
        "Nhits"
    ])

# Read Entry
entry = 5                                     ### Change your entry, you can make this entry as a loop !!!
for entry in range(10):  # from 0 - 10
    print(f"Processing entry: {entry}")
    Energy_vec  = data["BarEnergyVector"][entry]
    RMS_vec     = data["RMS"][entry]
    FH_Dep      = data["First_Had_Depth"][entry]
    FH_Lay      = data["First_Had_Layer"][entry]
    FH_Type     = data["First_Had_Type"][entry]
    FI_Dep      = data["First_Ine_Depth"][entry]
    FI_Lay      = data["First_Ine_Layer"][entry]
    Nhits       = data["Nhits"][entry]

    color_map = {
        -1:'orange',
        1: 'red',
        2: 'magenta',
    }

    if FH_Type == 1:
        string1 = "Inelastic"
    elif FH_Type == 2:
        string1 = "Elastic"
    else:
        string1 = "Pass"
    line_color = color_map.get(FH_Type, 'black')

    # Create 2D histogram
    hXZ = np.full((14, 22), -5.0)
    hYZ = np.full((14, 22), -5.0)

    # Fill the 2D histogram
    for i, energy in enumerate(Energy_vec):
        layer = i // 22
        bar = i % 22
        logE = math.log10(energy) if energy >= 1e-2 else -5.0

        if layer % 2 == 0:
            hXZ[layer, bar] = logE
        else:
            hYZ[layer, bar] = logE

    layer_start = 4
    RMS_threshold = 15.0
    bar_info = [0, 0]
    bar_info_assigned = False

    # ---------- Reconstruction ----------

    # ---------- Step1. find the thin trajory ----------
    for k in range(layer_start, 13):
        if RMS_vec[k] <= RMS_threshold and RMS_vec[k+1] <= RMS_threshold:
            max_index1 = np.argmax(Energy_vec[k*22:(k+1)*22])
            bar1 = max_index1
            max_index2 = np.argmax(Energy_vec[(k+1)*22:(k+2)*22])
            bar2 = max_index2
            if k % 2 == 0:
                bar_info = [bar1, bar2]
            else:
                bar_info = [bar2, bar1]
            bar_info_assigned = True
            break

    # ---------- Fit the track of shower axis----------
    if not bar_info_assigned:
        print("No bar_info assigned, start fitting...")

        def fit_func(x, bars, energies):
            return np.sum((bars - x[0])**2 * energies)

        # Fit the odd layer
        bars = np.arange(22)
        energies = Energy_vec[layer_start*22:(layer_start+1)*22]
        res = minimize(fit_func, [11], args=(bars, energies))
        bar_info[0] = round(res.x[0])

        # Fit the even layer
        energies = Energy_vec[(layer_start+1)*22:(layer_start+2)*22]
        res = minimize(fit_func, [11], args=(bars, energies))
        bar_info[1] = round(res.x[0])


    print("bar_info[0]=",bar_info[0])
    print("bar_info[1]=",bar_info[1])

    # Create histogram to store the information
    bar_Energy_info = np.zeros(14)
    bar_Change_info = np.zeros(14)
    bar_Accumu_info = np.zeros(14)
    bar_Accumu_error = np.zeros(14)

    for layer in range(14):
        print(f"RMS {layer} , {RMS_vec[layer]}")
        center_bar = bar_info[0] if layer % 2 == 0 else bar_info[1]
        for k in range(center_bar - 1, center_bar + 2):
            if 0 <= k < 22:
                bar_Energy_info[layer] += Energy_vec[layer * 22 + k]
        
        if layer == 0:
            bar_Change_info[0] = np.log10(bar_Energy_info[0] / 0.023)
            bar_Accumu_info[0] = bar_Energy_info[0]
            bar_Accumu_error[0] = 0.3 * bar_Accumu_info[0]
        else:
            prev = bar_Energy_info[layer-1]
            curr = bar_Energy_info[layer]
            bar_Accumu_info[layer] = bar_Accumu_info[layer-1] + bar_Energy_info[layer]
            bar_Accumu_error[layer] = 0.3 * bar_Accumu_info[layer]

            if prev == 0:
                bar_Change_info[layer-1] = -5
            elif curr == 0:
                bar_Change_info[layer-1] = -4
            else:
                bar_Change_info[layer] = np.log10(curr / prev)

    # Use bar_Change_info Array to Find the most longest increasing segment, find where it start and it value, and how long it is 
    max_sum, max_len, max_start_bin = find_max_positive_bin_segment(bar_Change_info)
    print(f"Positive Increase Sum: {max_sum}, Positive Increase Length: {max_len}, Start at Bin: {max_start_bin}")

    # Use ar_Change_info Array to Find the largest value and it index
    max_value, max_bin = find_max_value_in_positive_segment(bar_Change_info, max_start_bin, max_len)
    print(f"Max Changed Value: {max_value} at Bin: {max_bin}")

    # --------------------------- PLOT -------------------------------------


    fig, axs = plt.subplots(2, 3, figsize=(18, 12))

    masked_hXZ = np.ma.masked_less(hXZ, -2)  # mask the value < -2 
    masked_hYZ = np.ma.masked_less(hYZ, -2)
    # Create colormap：
    cmap = plt.cm.viridis
    cmap.set_bad(color='white')  # make the value < -2 is white

    # 1st SubPlot XZ plane
    axs[0, 0].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    im1 = axs[0, 0].imshow(masked_hXZ, cmap='viridis', extent=[0, 22, 14, 0], aspect='auto', vmin=-2, vmax=2)
    axs[0, 0].set_title("XZ plane")
    axs[0, 0].set_xticks(np.arange(0, 23, 1))           
    axs[0, 0].set_xticklabels(np.arange(-11, 12, 1))    
    fig.colorbar(im1, ax=axs[0, 0], label='log10(Energy)')
    box0 = patches.Rectangle((bar_info[0]-1, layer_start), width=3, height=14-layer_start, linewidth=1, edgecolor=line_color, facecolor='none')
    axs[0, 0].add_patch(box0)
    axs[0, 0].plot([0, 22], [FH_Lay]*2, linestyle='--', color=line_color, linewidth=1)
    axs[0, 0].scatter(bar_info[0]+0.5, max_start_bin ,  marker='*', color=line_color, s=100, label='Begin to Increase Point')
    axs[0, 0].scatter(bar_info[0]+0.5, max_bin ,  marker='^', color=line_color, s=100, label='Most Changed Point')


    # 2nd SubPlot YZ plane
    axs[0, 1].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    im2 = axs[0, 1].imshow(masked_hYZ, cmap='viridis', extent=[0, 22, 14, 0], aspect='auto', vmin=-2, vmax=2)
    axs[0, 1].set_title("YZ plane")
    axs[0, 1].set_xticks(np.arange(0, 23, 1))
    axs[0, 1].set_xticklabels(np.arange(-11, 12, 1))
    fig.colorbar(im2, ax=axs[0, 1], label='log10(Energy)')

    box1 = patches.Rectangle((bar_info[1]-1, layer_start), width=3, height=14-layer_start, linewidth=1, edgecolor=line_color, facecolor='none')
    axs[0, 1].add_patch(box1)
    axs[0, 1].plot([0, 22], [FH_Lay]*2, linestyle='--', color=line_color, linewidth=1)
    axs[0, 1].scatter(bar_info[1]+0.5, max_start_bin ,  marker='*', color=line_color, s=100, label='Begin to Increase Point')
    axs[0, 1].scatter(bar_info[1]+0.5, max_bin ,  marker='^', color=line_color, s=100, label='Most Changed Point')


    # 3rd Subplot Energy Array
    axs[1, 0].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    axs[1, 0].hist(
        np.arange(14),  
        bins=np.arange(15),  
        weights=bar_Energy_info, 
        histtype='step',  
        linewidth=1.5,
        color='black'
    )
    axs[1, 0].set_xlim(0, 14)
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_title("BGO Core Axis Energy Deposit")
    axs[1, 0].set_xlabel("BGO Layer")
    axs[1, 0].set_ylabel("3 Bars Energy / GeV")
    hist_E_min = np.min(bar_Energy_info)  
    hist_E_max = np.max(bar_Energy_info)  
    axs[1, 0].plot([FH_Lay, FH_Lay], [hist_E_min, hist_E_max], linestyle='--', color=line_color, linewidth=1)
    axs[1, 0].scatter(max_start_bin+0.5, bar_Energy_info[max_start_bin], marker='*', color=line_color, s=100, label='Begin to Increase Point')
    axs[1, 0].scatter(max_bin+0.5, bar_Energy_info[max_bin], marker='^', color=line_color, s=100, label='Most Changed Point')
    axs[1, 0].scatter(FI_Dep/25.5, bar_Energy_info[int(FI_Dep/25.5)], marker='X', color=line_color, s=100, label='Inelastic Point')

    # 4th Subplot Change Rate
    axs[1, 1].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    axs[1, 1].hist(
        np.arange(14),  
        bins=np.arange(15) , 
        weights=bar_Change_info,  
        histtype='step',  
        linewidth=1.5,
        color='black'
    )
    axs[1, 1].set_xlim(0, 14)
    axs[1, 1].set_ylim([-0.5, 2.5])
    axs[1, 1].set_title("Energy Deposit Change Rate")
    # axs[1, 1].set_xlabel("BGO Layer")
    axs[1, 1].set_ylabel(r"$\log_{10}(E_i / E_{i-1})$")
    hist_C_min = np.min(bar_Change_info)  
    hist_C_max = np.max(bar_Change_info)  
    axs[1, 1].plot([FH_Lay, FH_Lay], [hist_C_min, hist_C_max], linestyle='--', color=line_color, linewidth=1)
    axs[1, 1].scatter(max_start_bin+0.5, bar_Change_info[max_start_bin], marker='*', color=line_color, s=100, label='Begin to Increase Point')
    axs[1, 1].scatter(max_bin+0.5, bar_Change_info[max_bin], marker='^', color=line_color, s=100, label='Most Changed Point')
    axs[1, 1].scatter(FI_Dep/25.5, bar_Change_info[int(FI_Dep/25.5)], marker='X', color=line_color, s=100, label='Inelastic Point')

    # 5th Subplot Change Rate
    axs[1, 2].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    axs[1, 2].hist(
        np.arange(14),  
        bins=np.arange(15),  
        weights=bar_Accumu_info, 
        histtype='step',  
        linewidth=1.5,
        color='black'
    )
    axs[1, 2].set_xlim(0, 14)
    axs[1, 2].set_yscale('log')
    axs[1, 2].set_title("Accumulated Energy Deposit")
    axs[1, 2].set_xlabel("BGO Layer")
    axs[1, 2].set_ylabel("3 Bars Energy / GeV")
    hist_E_min = np.min(bar_Accumu_info)  
    hist_E_max = np.max(bar_Accumu_info)  
    axs[1, 2].plot([FH_Lay, FH_Lay], [hist_E_min, hist_E_max], linestyle='--', color=line_color, linewidth=1)
    axs[1, 2].scatter(max_start_bin+0.5, bar_Accumu_info[max_start_bin], marker='*', color=line_color, s=100, label='Begin to Increase Point')
    axs[1, 2].scatter(max_bin+0.5, bar_Accumu_info[max_bin], marker='^', color=line_color, s=100, label='Most Changed Point')
    axs[1, 2].scatter(FI_Dep/25.5, bar_Accumu_info[int(FI_Dep/25.5)], marker='X', color=line_color, s=100, label='Inelastic Point')


    # 1. Prepare fitting Array
    x_vals = np.arange(14) + 0.5
    y_vals = np.array(bar_Accumu_info)
    y_errs = np.array(bar_Accumu_error)

    # 2. Mod sigmoid function (done, above)
    # 3. Initial Value for fitting
    p0 = [0.02 ,y_vals[0], y_vals[-1], 6.0, 1.0]  # Ymin, Ymax, Xmid, Slope

    # 4. Fit
    popt, pcov = curve_fit(modified_sigmoid, x_vals, y_vals, sigma=y_errs, p0=p0, absolute_sigma=True)
    E0, Ymin, Ymax, Xmid, Slope = popt
    fit_vals = modified_sigmoid(x_vals, *popt)
    percentile2 = Mod_Sigmoid_Percentile(FI_Dep/25.5, Xmid, Slope)
    axs[1, 2].text(
        0.6, 0.85,  # x, y in axes fraction (0~1)
        f"FI Layer Percentile: {percentile2*100:.2f}%",
        transform=axs[1, 2].transAxes,
        fontsize=10,
        color='blue',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="blue")
    )
    # 5. Residual: (Data - Fit)/Err
    residuals = (y_vals - fit_vals) / y_errs

    # 6. Add fitting line to the histogram
    axs[1, 2].plot(x_vals, fit_vals, label='Sigmoid Fit', color='red', linestyle='-')
    axs[1, 2].legend()

    # 7. Add Residual Subplot
    axs[0, 2].bar(x_vals, residuals, width=0.8, color='gray', edgecolor='black')
    axs[0, 2].axhline(0, color='red', linestyle='--', linewidth=1)
    axs[0, 2].set_xlim(0, 14)
    axs[0, 2].set_ylim(-max(abs(residuals)) * 1.2, max(abs(residuals)) * 1.2)
    axs[0, 2].set_title("Residual (Data - Fit) / Err")
    axs[0, 2].set_xlabel("BGO Layer")
    axs[0, 2].set_ylabel("Residual")
    axs[0, 2].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)


    if entry < 100:
        save_path = f"/Users/xiongzheng/software/B4/B4e/Script/Figures/{string1}/{string2}/{entry}_PythonVer.png" ### Change your path!!!
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path)

    plt.tight_layout()
    # plt.show()

    entry += 1