std::vector<int> g_fit_bars;
std::vector<double> g_fit_energies;
double g_fit_total_energy = 0;

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


double MaxMinRatio(TH1D* hist) 
{
    int nbins = hist->GetNbinsX();

    double max_val = 0;
    double min_val = 1e4;

    for (int i = 1; i <= nbins; ++i) 
    {
        double content = hist->GetBinContent(i);
        if (content > max_val) max_val = content;
        if (content < min_val && content > 0) min_val = content;
        // cout << content <<  " , min =  " << min_val <<  " ,  max " << max_val <<  endl;
    }

    if (min_val == 0) throw std::runtime_error("最小值为 0，无法计算比值");

    return max_val / min_val;
}

int FindMaxMiddleIndex(std::vector<double>* p_EnergyVec, int layer) 
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
        // cout << "layer = " << layer << " max_index = " << max_index << endl;
        for (int offset = -2; offset <= 2; ++offset) {
            int idx = max_index + offset;
            if (idx >= base && idx < base + 22) {
                int bar = idx % 22;
                double energy = (*p_EnergyVec)[idx];
                // cout << " bar " << bar << " , energy = " << energy << endl;
                bars.push_back(bar);
                energies.push_back(energy);
                total_energy += energy;
            }
        }
    }
}

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
