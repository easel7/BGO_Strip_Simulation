std::vector<int> g_fit_bars;
std::vector<double> g_fit_energies;
std::vector<double> g_fit_errors;
double g_fit_total_energy = 0;
double g_fit_E00;
double g_fit_E13;


void FindMaxPositiveSegment(const double array[], int size, double& out_sum, int& out_len, int& out_start_index) 
{
    double max_sum = 0;
    int max_len = 0;
    int max_start_index = -1;

    double curr_sum = 0;
    int curr_len = 0;
    int curr_start_index = -1;

    for (int i = 0; i < size; ++i) 
    {
        double content = array[i];
        // std::cout << "bar_Change_info[" << i << "] = " << content << std::endl;

        if (content > 0) 
        {
            if (curr_len == 0)
                curr_start_index = i ;  // 新的一段开始

            curr_sum += content;
            curr_len++;

            // if (curr_len > max_len || (curr_sum > max_sum && curr_len == max_len)) 
            if (curr_sum > max_sum || (curr_sum == max_sum && curr_len > max_len)) 

            {
                max_sum = curr_sum;
                max_len = curr_len;
                max_start_index = curr_start_index;  // 记录最大段的起点
            }
        } 
        else 
        {
            curr_sum = 0;
            curr_len = 0;
            curr_start_index = -1;
        }
    }

    out_sum = max_sum;
    out_len = max_len;
    out_start_index = max_start_index;
}


void FindMaxValueInPositiveSegment(const double array[], int start_index, int length, double& out_max_value, int& out_max_index)
{
    out_max_value = -1e9;  // 初始值很小
    out_max_index = -1;

    for (int i = start_index; i < start_index + length; ++i) 
    {
        double content = array[i];
        if (content > out_max_value) 
        {
        out_max_value = content;
        out_max_index = i;
        }
    }
}

double MaxMinRatio(const double array[], int size) 
{
    if (size <= 0) throw std::invalid_argument("数组不能为空");

    double max_val = 0;
    double min_val = 1e4;

    for (int i = 0; i < size; ++i) 
    {
        if (array[i] > max_val) max_val = array[i];
        if (array[i] < min_val && array[i] > 0) min_val = array[i];
        // cout << array[i] <<  " , min =  " << min_val <<  " ,  max " << max_val <<  endl;
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

void PrepareSigmoidData(
    const double accumu[], 
    const double error[]
)
{
    g_fit_energies.clear();
    g_fit_errors.clear();
    for (int i = 0; i < 14; ++i) 
    {
        g_fit_energies.push_back(accumu[i]);
        g_fit_errors.push_back(error[i]);

    }
}

void SigmoidFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t iflag)
{
    double chi2 = 0;

    for (size_t i = 0; i < 14; ++i) {
        double x = i + 0.5;
        double y = g_fit_energies[i];
        double err = g_fit_errors[i];
        double model = par[0] + par[4] * x + (par[1] - par[0] - par[4] * x) / (1 + exp(-(x - par[2]) / par[3]));
        if (err > 0)
            chi2 += pow((y - model) / err, 2);
    }

    fval = chi2;
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

double ComputeReducedChi2(TMinuit& minuit, void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t), int n_points, int npar) {
    double chi2 = 0;
    double params[5]; // 足够容纳参数数量（可根据需要扩展）
    double* dummy_grad = nullptr;
    int iflag = 0;

    // 获取拟合参数
    for (int i = 0; i < npar; ++i) {
        double val, err;
        minuit.GetParameter(i, val, err);
        params[i] = val;
    }

    // 调用 FCN 计算 chi²
    fcn(npar, dummy_grad, chi2, params, iflag);

    // 计算自由度
    int n_free_params = minuit.GetNumFreePars();
    int ndf = n_points - n_free_params;

    if (ndf > 0) {
        return chi2 / ndf;
    } else {
        std::cerr << "WARNING: NDF <= 0, cannot compute reduced chi2" << std::endl;
        return -1.0;
    }
}

double RunSigmoidFit(Long64_t entry, double E_L0, double Amax, int seg_peak_idx, double maxE, TMinuit*& minuit_out) {
    minuit_out = new TMinuit(5);
    minuit_out->SetFCN(SigmoidFCN);
    minuit_out->SetPrintLevel(-1);
    minuit_out->SetErrorDef(1.0);

    minuit_out->DefineParameter(0, "Ymin" , E_L0         , E_L0 * 0.1   , 0                        , maxE                    );
    minuit_out->DefineParameter(1, "Ymax" , Amax         , Amax * 0.02  , maxE                     , 1e6                      );
    minuit_out->DefineParameter(2, "Xmid" , seg_peak_idx , 0.5          , max(seg_peak_idx-3.,-1.) , min(seg_peak_idx+3.,14.) );
    minuit_out->DefineParameter(3, "Slope", 1.0          , 0.1          , 0.1                      , 10                       );
    minuit_out->DefineParameter(4, "E0"   , E_L0         , E_L0 * 0.1   , 0.1 * E_L0               , 10. * E_L0               );
    int fit_status = minuit_out->Migrad();
    double edm,  chi2, errdef;
    int nvpar, nparx, istat;
    minuit_out->mnstat(chi2, edm, errdef, nvpar, nparx, istat);
    if (fit_status != 0 && edm > 1e-3) {fit_status = minuit_out->Migrad();}
    if (fit_status != 0 && edm > 1e-2) cerr << "ERROR: Fit did not converge! Status: " << fit_status << " entry : " << entry << " edm : " << edm << endl; 
    minuit_out->mnstat(chi2, edm, errdef, nvpar, nparx, istat);
    // std::cout << "Edm = " << edm << std::endl;
    // std::cout << "chi2 = " << chi2 << std::endl;
    // std::cout << "nvpar = " << nvpar << std::endl;
    // std::cout << "nvpax = " << nparx << std::endl;
    int ndf = 14 - nvpar;
    if (ndf > 0) {
        return chi2 / ndf;
    } else {
        std::cerr << "WARNING: NDF <= 0, cannot compute reduced chi2" << std::endl;
        return -1.0;
    }
}

//------------------------------- Array Tools --------------------------------------------------//

double AccumIncreaseToPeak(const double array[], int start_idx, int end_idx) 
{
    double sum = 0;
    for (int i = start_idx; i <= end_idx; ++i) {
        sum += array[i];
    }
    return sum;
}

double FindMaxValue(const double* arr, int size) {
    return *std::max_element(arr, arr + size);
}

