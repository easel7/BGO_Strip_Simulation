std::vector<int> g_fit_bars;
std::vector<double> g_fit_energies;
std::vector<double> g_fit_errors;
double g_fit_total_energy = 0;

//--------------------------------- Reconstruction -------------------------------------------

int FindMaxMiddleIndex(const std::vector<double>* p_EnergyVec, int layer) 
{
    auto p_start = p_EnergyVec->begin() + layer * 22;
    auto p_end   = p_EnergyVec->begin() + (layer + 1) * 22;

    double max_sum = -1e9;
    int best_index = -1;

    for (auto it = p_start + 1; it < p_end - 1; ++it) {
        double left = *(it - 1);
        double mid  = *it;
        double right = *(it + 1);

        double sum = left + mid + right;

        if (mid >= left && mid >= right) { // 中间最大
            if (sum > max_sum) {
                max_sum = sum;
                best_index = std::distance(p_EnergyVec->begin(), it); // 得到整体中的 index
            }
        }
    }
    return best_index; // 返回的是中间值的 index
}

bool AssignBarInfoFromRMS(const std::vector<double>* p_RMSVec,
    const std::vector<double>* p_EnergyVec,
    const std::vector<double>* p_L_EnergyVec,
    int bar_info[2],
    int layer_start = 4,
    double RMS_threshold = 15.0) 
{
    for (int k = layer_start; k <= 10; ++k) 
    {
        if ((*p_RMSVec)[k] <= RMS_threshold && (*p_L_EnergyVec)[k]     > 0 && 
        (*p_RMSVec)[k + 1] <= RMS_threshold && (*p_L_EnergyVec)[k + 1] > 0 && 
        (*p_RMSVec)[k + 2] <= RMS_threshold && (*p_L_EnergyVec)[k + 2] > 0 && 
        (*p_RMSVec)[k + 3] <= RMS_threshold && (*p_L_EnergyVec)[k + 3] > 0) 
        {
            int max_index1 = FindMaxMiddleIndex(p_EnergyVec, k);
            int bar1 = max_index1 % 22;
            int max_index2 = FindMaxMiddleIndex(p_EnergyVec, k + 1);
            int bar2 = max_index2 % 22;
            // cout << "RMS at Layer "<< k << " , " << (*p_RMSVec)[k] << endl;
            if (k % 2 == 0) {bar_info[0] = bar1; bar_info[1] = bar2; } 
            else            {bar_info[0] = bar2; bar_info[1] = bar1; }
            // std::cout << "bar_info assigned" << std::endl;
            return true;
        }
    }
    return false;
}

void FitAxisFunction(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    double bar0 = par[0];
    double cost = 0;
    for (size_t i = 0; i < g_fit_bars.size(); ++i) 
    {
        double diff = g_fit_bars[i] - bar0;
        cost += diff * diff * g_fit_energies[i] / g_fit_total_energy;
    }
    fval = cost;
}

void PrepareFitData(
    std::vector<double>* p_EnergyVec,
    int layer_start,
    int layer_end,
    std::vector<int>& bars,
    std::vector<double>& energies,
    double& total_energy)
{
    bars.clear();
    energies.clear();
    total_energy = 0.0;
    for (int layer = layer_start; layer < layer_end; layer += 2) 
    {
        int max_index = FindMaxMiddleIndex(p_EnergyVec, layer);
        int base = layer * 22;
        // cout << "layer = " << layer << " max_index = " << max_index << " max_bar = " <<  max_index % 22 << endl;
        for (int offset = -4; offset <= 4; ++offset) {
            int idx = max_index + offset;
            if (idx >= base && idx < base + 22) {
                int bar = idx % 22;
                double energy = (*p_EnergyVec)[idx];
                if (energy < 1e-2) continue;
                // cout << "layer " << layer << " bar " << bar << " , energy = " << energy << endl;
                bars.push_back(bar);
                energies.push_back(energy);
                total_energy += energy;
            }
        }
    }
}

bool Fit1DParameter(void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t),
                    double init_val, double init_err,
                    double lower_bound, double upper_bound,
                    double& fitted_val, double& fitted_err,
                    int print_level = -1) {
    TMinuit minuit(1);
    minuit.SetFCN(fcn);
    minuit.SetPrintLevel(print_level);
    minuit.SetErrorDef(1.0);  // standard chi2
    minuit.DefineParameter(0, "param", init_val, init_err, lower_bound, upper_bound);
    int status = minuit.Migrad();
    if (status != 0) {
        std::cerr << "WARNING: Fit did not converge! Status = " << status << std::endl;
        return false;
    }
    minuit.GetParameter(0, fitted_val, fitted_err);
    return true;
}

void ComputeBarEnergyInfo(const std::vector<double>* p_EnergyVec,
                          const int bar_info[2],
                          double bar_Energy_info[14],
                          double bar_Change_info[14],
                          double bar_Accumu_info[14]) 
{
    // Step 1: 计算每层中心3根条的能量和
    for (int layer = 0; layer < 14; ++layer) {
        int center_bar = (layer % 2 == 0) ? bar_info[0] : bar_info[1];
        for (int k = center_bar - 1; k <= center_bar + 1; ++k) {
            if (k >= 0 && k < 22) {
                bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
            }
        }
    }

    // Step 2: 初始化第一层
    bar_Change_info[0] = log10(bar_Energy_info[0] / 0.023);
    bar_Accumu_info[0] = bar_Energy_info[0];

    // Step 3: 累积能量和变化率
    for (int layer = 1; layer < 14; ++layer) {
        if (bar_Energy_info[layer - 1] == 0 || bar_Energy_info[layer] == 0) {
            bar_Change_info[layer - 1] = -5;
        } else {
            bar_Change_info[layer] = log10(bar_Energy_info[layer] / bar_Energy_info[layer - 1]);
        }

        bar_Accumu_info[layer]  = bar_Accumu_info[layer - 1] + bar_Energy_info[layer];
    }
}

void FillBGOHistograms(TH1D* hBGO1, TH1D* hBGO2, TH1D* hBGO3,
    const double bar_Energy_info[14],
    const double bar_Change_info[14],
    const double bar_Accumu_info[14],
    double error_scale = 0.3) 
{
    for (int layer = 0; layer < 14; ++layer) {
    hBGO1->SetBinContent(layer + 1, bar_Energy_info[layer]);hBGO1->SetBinError(layer + 1, error_scale * bar_Energy_info[layer]);
    hBGO2->SetBinContent(layer + 1, bar_Change_info[layer]);hBGO2->SetBinError(layer + 1, error_scale * bar_Change_info[layer]);
    hBGO3->SetBinContent(layer + 1, bar_Accumu_info[layer]);hBGO3->SetBinError(layer + 1, error_scale * bar_Accumu_info[layer]);
    }
}

//------------------- HIST TOOLS ----------------------------------------

void FindMaxPositiveBinSegment(TH1D* hist, double& out_sum, int& out_len, int& out_start_bin) 
{
    int nbins = hist->GetNbinsX();
    double max_sum = 0;
    int max_len = 0;
    int max_start_bin = -1;

    double curr_sum = 0;
    int curr_len = 0;
    int curr_start_bin = -1;

    for (int i = 1; i <= nbins; ++i) 
    {
        double content = hist->GetBinContent(i);
        // std::cout << "bin[" << i << "] = " << content << std::endl;

        if (content > 0) 
        {
            if (curr_len == 0) 
                curr_start_bin = i;  // 记录这一段开始的位置

            curr_sum += content;
            curr_len++;

            // if (curr_len > max_len || (curr_sum > max_sum && curr_len == max_len)) 
            if (curr_sum > max_sum || (curr_sum == max_sum && curr_len > max_len)) 
            {
                max_sum = curr_sum;
                max_len = curr_len;
                max_start_bin = curr_start_bin;  // 更新最大值对应的起始位置
            }
        } 
        else 
        {
            curr_sum = 0;
            curr_len = 0;
            curr_start_bin = -1;
        }
    }
    out_sum = max_sum;
    out_len = max_len;
    out_start_bin = max_start_bin;
}


void FindMaxValueInPositiveSegment(TH1D* hist, int start_bin, int len, double& out_max_value, int& out_max_bin)
{
    out_max_value = -1e9;  // 默认很小
    out_max_bin = -1;

    for (int i = start_bin; i < start_bin + len; ++i) 
    {
        double content = hist->GetBinContent(i);
        if (content > out_max_value) 
        {
            out_max_value = content;
            out_max_bin = i;
        }
    }
}

double MaxMinRatio(TH1D* hist, Long64_t entry) 
{
    double max_val = hist->GetMaximum();
    double min_val = hist->GetMinimum();
    if (min_val == 0) {min_val = 0.01;}
    return max_val / min_val;
}

double AccumIncreaseToPeak(TH1D* hist, int start_idx, int end_idx) 
{
    double sum = 0;
    for (int i = start_idx; i <= end_idx; ++i) {
        sum += hist->GetBinContent(i);
    }
    return sum;
}

//----------------------------------------------------------------------------------

double  Mod_Sigmoid_Percentile(double Depth2Layer, double Xmid, double Slope)
{
    double model = 1 / (1 + exp(-(Depth2Layer - Xmid) / Slope));
    return model;
}

int Inverse_Mod_sigmoid(double percentile,double Xmid, double Slope)
{
    if (percentile <= 0 || percentile >= 1) {
        std::cerr << "Percentile p must be in (0,1)" << std::endl;
        return -999;
    }
    double X =  Xmid + Slope * log(percentile / (1.0 - percentile));
    return floor(X);
}

//----------------------------------------------------------------------------------
