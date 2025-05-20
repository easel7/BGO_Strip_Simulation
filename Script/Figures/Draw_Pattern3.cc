#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"

void Draw_Pattern3()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_L_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    int p_FH_Type;
    int p_FH_Lay;
    int p_Nhits;
    double p_FI_Dep;
    int p_FI_Lay;
    double p_Total_E;
    const char* string1;
    const char* string2 = "Deuteron_PowerLaw";
    // const char* string2 = "Proton_PowerLaw";

    auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/%s.root",string2));
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("BarEnergyVector", &p_EnergyVec);
    proton_tree->SetBranchAddress("LayerEnergyVector", &p_L_EnergyVec);
    proton_tree->SetBranchAddress("RMS"            , &p_RMSVec);
    proton_tree->SetBranchAddress("First_Had_Depth", &p_FH_Dep);
    proton_tree->SetBranchAddress("First_Had_Layer", &p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type" , &p_FH_Type);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);
    proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
    proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);
    proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);

    auto h1_p = new TH1D("h1_p", "h1_p",100,-10,0);  

    // cout  << proton_tree->GetEntries() << endl;
    Long64_t entry  = 314870;   
    int Counts = 0;
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    {
        proton_tree->GetEntry(entry);
        if (p_Nhits < 10 ) continue;
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        if (gDirectory->FindObject("hBGO1")) delete gDirectory->FindObject("hBGO1");
        if (gDirectory->FindObject("hBGO2")) delete gDirectory->FindObject("hBGO2");
        if (gDirectory->FindObject("hBGO3")) delete gDirectory->FindObject("hBGO3");
        if (gDirectory->FindObject("sigmoid")) delete gDirectory->FindObject("sigmoid");
        TH1D *hBGO1     = new TH1D("hBGO1","BGO Core Axis Energy Deposit",14,0,14); 
        TH1D *hBGO2     = new TH1D("hBGO2","Deposit Energy Change Ratio",14,0,14); 
        TH1D *hBGO3     = new TH1D("hBGO3","Accumulated Deposit Energy",14,0,14); 
        int bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double bar_Accumu_info[14] = {0};
        double bar_Accumu_error[14] = {0};
        double seg_sum      = 0;
        int    seg_len      = 0;
        int    seg_start_idx      = 0;
        double seg_peak_value      = 0;
        int    seg_peak_idx   = 0;
        int layer_start = 4;
        const double RMS_threshold = 15.0;  
        bool bar_info_assigned = AssignBarInfoFromRMS(p_RMSVec, p_EnergyVec, p_L_EnergyVec, bar_info, layer_start, RMS_threshold);

        if (!bar_info_assigned) {
            double bar_odd, bar_odd_err;
            double bar_even, bar_even_err;
            PrepareFitData(p_EnergyVec, layer_start, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            bool success_odd = Fit1DParameter(FitAxisFunction, g_fit_bars[1], 0.01, 2, 19,  bar_odd, bar_odd_err);
            if (success_odd) bar_info[0] = std::round(bar_odd);
            else std::cerr << "Failed to fit bar_odd. Check " << entry << std::endl;

            PrepareFitData(p_EnergyVec, layer_start+1, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            bool success_even = Fit1DParameter(FitAxisFunction, g_fit_bars[1], 0.01, 2, 19,  bar_even, bar_even_err);
            if (success_even) bar_info[1] = std::round(bar_even);
            else std::cerr << "Failed to fit bar_odd. Check " << entry << std::endl;
        }

        ComputeBarEnergyInfo(p_EnergyVec, bar_info, bar_Energy_info, bar_Change_info, bar_Accumu_info);
        FillBGOHistograms(hBGO1, hBGO2, hBGO3, bar_Energy_info, bar_Change_info, bar_Accumu_info,0.3);
        FindMaxPositiveBinSegment(hBGO2,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(hBGO2,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);

        TF1* sigmoid = new TF1("sigmoid", "[0]+ [4]*x + ([1]-[0] - [4]*x )/(1 + exp(-(x-[2])/[3]))", 0, 14);
        sigmoid->SetParameters(hBGO1->GetBinContent(1), 
                               hBGO3->GetBinContent(14), 
                               seg_peak_idx, 
                               1 , 
                               max(hBGO1->GetBinContent(1)*0.1,0.01)); 
        sigmoid->SetParLimits(0, 0                                     , hBGO1->GetMaximum()         );   // [0] Ymin
        sigmoid->SetParLimits(1, hBGO1->GetMaximum()                   , 1e6                         );   // [1] Ymax
        sigmoid->SetParLimits(2, max(seg_peak_idx-3, -1)               , min(seg_peak_idx + 2, 14)   );   // [2] Xmid
        sigmoid->SetParLimits(3, 0.1                                   , 10                          );   // [3] Slope，避免除0
        sigmoid->SetParLimits(4, max(hBGO1->GetBinContent(1)*0.1,0.01) , hBGO1->GetBinContent(1) * 10);   // [4] linear bias
        TFitResultPtr fitResult = hBGO3->Fit(sigmoid, "RSQ");  // R: fit range, S: return TFitResultPtr
        double percentile2 = Mod_Sigmoid_Percentile(p_FI_Dep/25.5,sigmoid->GetParameter(2),sigmoid->GetParameter(3));
        int Percent2Layer = Inverse_Mod_sigmoid(percentile2, sigmoid->GetParameter(2), sigmoid->GetParameter(3));
        h1_p->Fill(log10(percentile2));
    }
    h1_p->Draw();
}

