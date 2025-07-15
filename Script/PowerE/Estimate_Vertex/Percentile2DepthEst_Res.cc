#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 \hat{x_ine} =  - 4.789 * sigma + x_mid

void Percentile2DepthEst_Res()
{
    double Energy_Dep[8] = {8,19,43,100,320,1e3,3.2e3,1e4};
    double Energy_Bin[7] = {0};

    double Ratio[50] = {0};
    double CL95_Sur[7]  = {0};    double CL95_Sur_Err[7]  = {0};    
    double CL90_Sur[7]  = {0};    double CL90_Sur_Err[7]  = {0};    
    double CL95_Int[7]  = {0};    double CL95_Int_Err[7]  = {0};    
    double CL90_Int[7]  = {0};    double CL90_Int_Err[7]  = {0};    

    auto mean_file    = TFile::Open("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/EnergyEst_Res.root");
    auto hist_p_2     = (TH1D*)mean_file->Get("hist_p_2");
    auto Energy_Mid_p = (TH1D*)mean_file->Get("Energy_Mid_p");

    int nbins = 60;
    double xmin = 1e1;
    double xmax = 1e4;
    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    std::vector<double> bin_edges(nbins + 1);
    for (int j = 0; j <= nbins; j++) {
        bin_edges[j] = pow(10, logxmin + j * (logxmax - logxmin) / nbins);
        // cout << bin_edges[j] << endl;
    }

    auto file1 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/proton_BGO.root");
    auto hist1 = (TH1D*)file1->Get("h4");
    auto hist6 = (TH1D*)file1->Get("h1");

    auto file2 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/deuteron_BGO.root");
    auto hist2 = (TH1D*)file2->Get("h4");
    auto hist7 = (TH1D*)file2->Get("h1");
    int Nbins = hist2->GetNbinsX();
    double KN_Energy[60] = {0};
    double Deuteron[60] = {0};        double Proton[60] = {0};
    double KN_Deuteron[60] = {0};     double KN_Proton[60] = {0};
    for (int ii = 0 ; ii < Nbins ; ii++)
    {
        KN_Energy[ii]   = pow(10,hist2->GetBinCenter(ii+1)-3);
        file1->cd();Proton[ii]   = hist1->GetBinContent(ii+1); KN_Proton[ii]   = hist6->GetBinContent(ii+1)* 10.;
        file2->cd();Deuteron[ii] = hist2->GetBinContent(ii+1); KN_Deuteron[ii] = hist7->GetBinContent(ii+1)* 10.;
    }
    auto gre6 = new TGraph(Nbins,KN_Energy, KN_Proton);     gre6->SetLineColor(kRed);
    auto gre7 = new TGraph(Nbins,KN_Energy, KN_Deuteron);   gre7->SetLineColor(kBlue);
    auto line_Ntot = new TGraph();

    int p_FH_Lay; int p_FH_Type; double p_energy_res;  double p_Energy;  double p_Total_E;    int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_FH_Lay; int d_FH_Type; double d_energy_res;  double d_Energy;  double d_Total_E;    int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
    int p_FI_Lay;    double p_FI_Dep;    int p_particle;
    int d_FI_Lay;    double d_FI_Dep;    int d_particle;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("Particle"         ,&p_particle);
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
    proton_tree->SetBranchAddress("BarEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type"   ,&p_FH_Type);
    proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
    proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);
    proton_tree->SetBranchAddress("energy_res"     ,&p_energy_res);
    proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
    proton_tree->SetBranchAddress("Energy"         , &p_Energy);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);
    proton_tree->SetBranchAddress("weight"          ,&p_weight);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("Particle"         ,&d_particle);
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
    deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_FH_Lay);
    deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_FH_Type);
    deuteron_tree->SetBranchAddress("First_Ine_Depth", &d_FI_Dep);
    deuteron_tree->SetBranchAddress("First_Ine_Layer", &d_FI_Lay);
    deuteron_tree->SetBranchAddress("energy_res"     ,&d_energy_res);
    deuteron_tree->SetBranchAddress("Total_E"         ,&d_Total_E);
    deuteron_tree->SetBranchAddress("Energy"         , &d_Energy);
    deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);

    TGraphErrors* grL_p_int = new TGraphErrors(); TGraphErrors* grN_p_int = new TGraphErrors();
    TGraphErrors* grL_p_sur = new TGraphErrors(); TGraphErrors* grN_p_sur = new TGraphErrors();
    TGraphErrors* grL_d_int = new TGraphErrors(); TGraphErrors* grN_d_int = new TGraphErrors();
    TGraphErrors* grL_d_sur = new TGraphErrors(); TGraphErrors* grN_d_sur = new TGraphErrors();
    
    double Layer[14]={0};
    double Layer_Err[14]={0};
    TF1 *fitFunc_p[7][14];  
    TF1 *fitFunc_d[7][14];  
    TH1D *h1_p_Lay[7][15]; // h1_p_Lay[x][14] 用来装Passthrough的事例
    TH1D *h1_d_Lay[7][15]; // h1_d_Lay[x][14] 用来装Passthrough的事例
    TH2D *h2_p_tot[7]; 
    TH2D *h2_d_tot[7]; 
    TH1D *h1_p_int[7]; 
    TH1D *h1_d_int[7]; 
    TH1D *h1_p_inl[7]; 
    TH1D *h1_d_inl[7]; 
    TH1D *h1_p_sur[7]; 
    TH1D *h1_d_sur[7]; 
    TH1D *h1_p_lea[7]; 
    TH1D *h1_d_lea[7]; 
    TH1D *hC_p_int[7];
    TH1D *hC_p_sur[7];      
    TGraphErrors *gre_p_reso[7];
    TGraphErrors *gre_d_reso[7];
    TGraphErrors *gre_p_bias[7];
    TGraphErrors *gre_d_bias[7];

    for(int i=0 ;i <14; i++)
    {
        Layer[i] = 0.5 + i;
        Layer_Err[i] = 0.5;
    }

    for(int j=0 ;j< 7 ; j++)
    {
        h2_p_tot[j] = new TH2D(Form("h2_p_tot[%d]",j),Form("h2_p_tot[%d]",j),200,-25,375,325,-275,375);   h2_p_tot[j]->Sumw2();
        h2_d_tot[j] = new TH2D(Form("h2_d_tot[%d]",j),Form("h2_d_tot[%d]",j),200,-25,375,325,-275,375);   h2_d_tot[j]->Sumw2();
        h1_p_int[j] = new TH1D(Form("h1_p_int[%d]",j),Form("h1_p_int[%d]",j),14,0,357);                   h1_p_int[j]->Sumw2();// EST
        h1_d_int[j] = new TH1D(Form("h1_d_int[%d]",j),Form("h1_d_int[%d]",j),14,0,357);                   h1_d_int[j]->Sumw2();// EST
        h1_p_sur[j] = new TH1D(Form("h1_p_sur[%d]",j),Form("h1_p_sur[%d]",j),14,0,357);                   // EST
        h1_d_sur[j] = new TH1D(Form("h1_d_sur[%d]",j),Form("h1_d_sur[%d]",j),14,0,357);                   // EST
        hC_p_int[j] = new TH1D(Form("hC_p_int[%d]",j),Form("hC_p_int[%d]",j),14,0,357);                    hC_p_int[j]->Sumw2();
        hC_p_sur[j] = new TH1D(Form("hC_p_sur[%d]",j),Form("hC_p_sur[%d]",j),14,0,357); 

        h1_p_inl[j] = new TH1D(Form("h1_p_inl[%d]",j),Form("h1_p_inl[%d]",j),14,0,357);   h1_p_inl[j]->Sumw2();// REAL
        h1_d_inl[j] = new TH1D(Form("h1_d_inl[%d]",j),Form("h1_d_inl[%d]",j),14,0,357);   h1_d_inl[j]->Sumw2();// REAL
        h1_p_lea[j] = new TH1D(Form("h1_p_lea[%d]",j),Form("h1_p_lea[%d]",j),14,0,357);  // REAL
        h1_d_lea[j] = new TH1D(Form("h1_d_lea[%d]",j),Form("h1_d_lea[%d]",j),14,0,357);  // REAL

        gre_p_reso[j] = new TGraphErrors();
        gre_d_reso[j] = new TGraphErrors();
        gre_p_bias[j] = new TGraphErrors();
        gre_d_bias[j] = new TGraphErrors();
        for(int i =0 ; i<15 ; i++)  // Layer
        {
            h1_p_Lay[j][i] =new TH1D(Form("h1_p_Lay[%d][%d]",j,i),Form("h1_p_Lay[%d][%d]",j,i),30,-3,3); h1_p_Lay[j][i]->Sumw2(); 
            h1_d_Lay[j][i] =new TH1D(Form("h1_d_Lay[%d][%d]",j,i),Form("h1_d_Lay[%d][%d]",j,i),30,-3,3); h1_d_Lay[j][i]->Sumw2();
            if (i<14)
            {
                fitFunc_p[j][i] = new TF1(Form("fitFunc_p[%d][%d]",j, i), "gaus", -1, 1); // 注意替换范围
                fitFunc_d[j][i] = new TF1(Form("fitFunc_d[%d][%d]",j, i), "gaus", -1, 1); // 注意替换范围
            }
        }
    }

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    // for (Long64_t entry = 0; entry < 50000; entry++)
    {        
        proton_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Proton : " << entry << endl;
        int p_energy_index = int(floor((log10(p_energy_res) ) / 0.2));
        if (p_energy_index < 0 || p_energy_index > 20) continue;
        if (p_FI_Dep < 0) {p_FI_Lay = 14;}
        if (p_Nhits < 10 ) continue;
        if (p_energy_res < hist_p_2->GetXaxis()->GetXmin() || p_energy_res > hist_p_2->GetXaxis()->GetXmax()) continue;
        int BEST_FIT_MEAN_BIN  = hist_p_2->FindBin(p_energy_res);
        double BEST_FIT_MEAN   = hist_p_2->GetBinContent(BEST_FIT_MEAN_BIN);
        double sum_p = 0;
        int bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double bar_Accumu_info[14] = {0};
        if (gDirectory->FindObject("hBGO1")) delete gDirectory->FindObject("hBGO1");
        if (gDirectory->FindObject("hBGO2")) delete gDirectory->FindObject("hBGO2");
        if (gDirectory->FindObject("hBGO3")) delete gDirectory->FindObject("hBGO3");
        if (gDirectory->FindObject("sigmoid")) delete gDirectory->FindObject("sigmoid");
        auto hBGO1 = new TH1D("hBGO1","BGO Core Axis Energy Deposit",14,0,14); 
        auto hBGO2 = new TH1D("hBGO2","Deposit Energy Change Ratio",14,0,14); 
        auto hBGO3 = new TH1D("hBGO3","Accumulated Deposit Energy",14,0,14); 
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引
        int layer_start = 4;
        const double RMS_threshold = 15.0;  
        bool bar_info_assigned = AssignBarInfoFromRMS(p_RMSVec, p_EnergyVec, p_L_EnergyVec, bar_info, layer_start, RMS_threshold);
        
        if (!bar_info_assigned) {
            // cout << "No bar_info assigned, starting fit to determine cluster trajectory." << endl;
            double bar_odd, bar_odd_err;
            double bar_even, bar_even_err;
            PrepareFitData(p_EnergyVec, layer_start, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            bool success_odd = Fit1DParameter(FitAxisFunction, g_fit_bars[1], 0.01, 2, 19,  bar_odd, bar_odd_err);
            if (success_odd)              bar_info[0] = std::round(bar_odd);
            else                std::cerr << "Failed to fit bar_odd." << std::endl;

            PrepareFitData(p_EnergyVec, layer_start+1, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            bool success_even = Fit1DParameter(FitAxisFunction, g_fit_bars[1], 0.01, 2, 19,  bar_even, bar_even_err);
            if (success_even) bar_info[1] = std::round(bar_even);
            else std::cerr << "Failed to fit bar_odd." << std::endl;
        }

        ComputeBarEnergyInfo(p_EnergyVec, bar_info, bar_Energy_info, bar_Change_info, bar_Accumu_info);
        FillBGOHistograms(hBGO1, hBGO2, hBGO3, bar_Energy_info, bar_Change_info, bar_Accumu_info,0.3);
        FindMaxPositiveBinSegment(hBGO2,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(hBGO2,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);
        seg_sum_to_peak = AccumIncreaseToPeak(hBGO1,seg_start_idx,seg_peak_idx);
        seg_len_to_peak = seg_peak_idx - seg_start_idx;

        double upper_bound = std::min(seg_peak_idx+7, 14 );
        double lower_bound = std::max(seg_peak_idx-7, -1 );
        TF1 *sigmoid = new TF1("sigmoid", "[0]+ [4]*x + ([1]-[0] - [4]*x )/(1 + exp(-(x-[2])/[3]))", 0, 14);
        sigmoid->SetParameters(hBGO1->GetBinContent(1), 
                                hBGO3->GetBinContent(14), 
                                seg_peak_idx - 0.5, 
                                1 , 
                                max(hBGO1->GetBinContent(1)*0.1,0.01)); 
        sigmoid->SetParLimits(0, 0                                     , hBGO1->GetMaximum()         );   // [0] Ymin
        sigmoid->SetParLimits(1, hBGO1->GetMaximum()                   , 1e6                         );   // [1] Ymax
        sigmoid->SetParLimits(2, lower_bound                           , upper_bound                 );   // [2] Xmid
        sigmoid->SetParLimits(3, 0.01                                  , 10                          );   // [3] Slope，避免除0
        sigmoid->SetParLimits(4, max(hBGO1->GetBinContent(1)*0.1,0.01) , hBGO1->GetBinContent(1) * 10);   // [4] linear bias
        TFitResultPtr fitResult = hBGO3->Fit(sigmoid, "RSQ");  // R: fit range, S: return TFitResultPtr
        double Est_Depth = (BEST_FIT_MEAN * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5;
        // if ( Est_Depth < -40 ) cout << "Watch out " << entry << " Value " << Est_Depth <<endl;
        h2_p_tot[BEST_FIT_MEAN_BIN-1]->Fill(p_FI_Dep, Est_Depth, p_weight*1e4);   
        h1_p_Lay[BEST_FIT_MEAN_BIN-1][p_FI_Lay]->Fill( (p_FI_Dep -  Est_Depth) / p_FI_Dep , p_weight*1e4);  // bias and reso
        h1_p_int[BEST_FIT_MEAN_BIN-1]->Fill(Est_Depth, p_weight*1e4);
        h1_p_inl[BEST_FIT_MEAN_BIN-1]->Fill(p_FI_Dep, p_weight*1e4);
        hC_p_int[BEST_FIT_MEAN_BIN-1]->Fill(Est_Depth, p_weight*1e4);
    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries() ; ++entry)
    // for (Long64_t entry = 0; entry < 50000 ; ++entry)
    {
        deuteron_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Deuteron : " << entry << endl;
        int d_energy_index = int(floor((log10(d_energy_res) ) / 0.2));
        if (d_energy_index < 0 || d_energy_index > 20) continue;
        if (d_FI_Dep < 0) {d_FI_Lay = 14;}
        if (d_Nhits < 10 ) continue;
        if (d_energy_res < hist_p_2->GetXaxis()->GetXmin() || d_energy_res > hist_p_2->GetXaxis()->GetXmax()) continue;
        int BEST_FIT_MEAN_BIN = hist_p_2->FindBin(d_energy_res);
        double BEST_FIT_MEAN     = hist_p_2->GetBinContent(BEST_FIT_MEAN_BIN);
        double sum_d = 0;
        int bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double bar_Accumu_info[14] = {0};
        if (gDirectory->FindObject("hBGO1")) delete gDirectory->FindObject("hBGO1");
        if (gDirectory->FindObject("hBGO2")) delete gDirectory->FindObject("hBGO2");
        if (gDirectory->FindObject("hBGO3")) delete gDirectory->FindObject("hBGO3");
        if (gDirectory->FindObject("sigmoid")) delete gDirectory->FindObject("sigmoid");
        TH1D *hBGO1     = new TH1D("hBGO1","BGO Core Axis Energy Deposit",14,0,14); 
        TH1D *hBGO2     = new TH1D("hBGO2","Deposit Energy Change Ratio",14,0,14); 
        TH1D *hBGO3     = new TH1D("hBGO3","Accumulated Deposit Energy",14,0,14); 
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引
        int    layer_start = 4;
        const double RMS_threshold = 15.0;  
        bool bar_info_assigned = AssignBarInfoFromRMS(d_RMSVec, d_EnergyVec, d_L_EnergyVec, bar_info, layer_start, RMS_threshold);
        
        if (!bar_info_assigned) {
            // cout << "No bar_info assigned, starting fit to determine cluster trajectory." << endl;
            double bar_odd, bar_odd_err;
            double bar_even, bar_even_err;
            PrepareFitData(d_EnergyVec, layer_start, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            bool success_odd = Fit1DParameter(FitAxisFunction, g_fit_bars[1], 0.01, 2, 19,  bar_odd, bar_odd_err);
            if (success_odd)              bar_info[0] = std::round(bar_odd);
            else                std::cerr << "Failed to fit bar_odd." << std::endl;

            PrepareFitData(d_EnergyVec, layer_start+1, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            bool success_even = Fit1DParameter(FitAxisFunction, g_fit_bars[1], 0.01, 2, 19,  bar_even, bar_even_err);
            if (success_even) bar_info[1] = std::round(bar_even);
            else std::cerr << "Failed to fit bar_odd." << std::endl;
        }

        ComputeBarEnergyInfo(d_EnergyVec, bar_info, bar_Energy_info, bar_Change_info, bar_Accumu_info);
        FillBGOHistograms(hBGO1, hBGO2, hBGO3, bar_Energy_info, bar_Change_info, bar_Accumu_info,0.3);
        FindMaxPositiveBinSegment(hBGO2,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(hBGO2,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);
        seg_sum_to_peak = AccumIncreaseToPeak(hBGO1,seg_start_idx,seg_peak_idx);
        seg_len_to_peak = seg_peak_idx - seg_start_idx;

        double upper_bound = std::min(seg_peak_idx+7, 14 );
        double lower_bound = std::max(seg_peak_idx-7, -1 );
        TF1 *sigmoid = new TF1("sigmoid", "[0]+ [4]*x + ([1]-[0] -  [4]*x )/(1 + exp(-(x-[2])/[3]))", 0, 14);
        sigmoid->SetParameters(hBGO1->GetBinContent(1), 
                            hBGO3->GetBinContent(14), 
                            seg_peak_idx - 0.5, 
                            1 , 
                            max(hBGO1->GetBinContent(1)*0.1,0.01)); 
        sigmoid->SetParLimits(0, 0                                     , hBGO1->GetMaximum()         );   // [0] Ymin
        sigmoid->SetParLimits(1, hBGO1->GetMaximum()                   , 1e6                         );   // [1] Ymax
        sigmoid->SetParLimits(2, lower_bound                           , upper_bound                 );   // [2] Xmid
        sigmoid->SetParLimits(3, 0.01                                  , 10                          );   // [3] Slope，避免除0
        sigmoid->SetParLimits(4, max(hBGO1->GetBinContent(1)*0.1,0.01) , hBGO1->GetBinContent(1) * 10);   // [4] linear bias
        TFitResultPtr fitResult = hBGO3->Fit(sigmoid, "RSQ");  // R: fit range, S: return TFitResultPtr
        double Est_Depth = (BEST_FIT_MEAN * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5;
        h2_d_tot[BEST_FIT_MEAN_BIN-1]->Fill(d_FI_Dep, Est_Depth, d_weight*1e4);
        h1_d_Lay[BEST_FIT_MEAN_BIN-1][d_FI_Lay]->Fill( (d_FI_Dep -  Est_Depth) / d_FI_Dep , d_weight*1e4);
        h1_d_int[BEST_FIT_MEAN_BIN-1]->Fill(Est_Depth, d_weight*1e4);
        h1_d_inl[BEST_FIT_MEAN_BIN-1]->Fill(d_FI_Dep, d_weight*1e4);
    }

    // for(int j=6 ;j< 7 ; j++)
    for(int j=0 ;j< 7 ; j++)
    {
        CL95_Sur[j]=1e-3; int Sur_95_Tag = -1;
        CL90_Sur[j]=1e-3; int Sur_90_Tag = -1;
        CL95_Int[j]=1e-3; int Int_95_Tag = -1;
        CL90_Int[j]=1e-3; int Int_90_Tag = -1;

        auto gre_int = new TGraphErrors();  
        auto gre_sur = new TGraphErrors();
        auto gre_p_int = new TGraphErrors();  
        auto gre_p_sur = new TGraphErrors();
        auto gre_d_int = new TGraphErrors();  
        auto gre_d_sur = new TGraphErrors();
        auto chi2_int = new TGraph();
        auto chi2_sur = new TGraph();

        for(int ii = 1  ; ii<= h1_p_int[j]->GetNbinsX(); ii++)
        {
            hC_p_sur[j]->SetBinContent(ii ,(h2_p_tot[j]->Integral() - hC_p_int[j]->Integral(1,ii) ));
            h1_p_sur[j]->SetBinContent(ii, (h2_p_tot[j]->Integral() - h1_p_int[j]->Integral(1,ii) ));
            h1_d_sur[j]->SetBinContent(ii, (h2_d_tot[j]->Integral() - h1_d_int[j]->Integral(1,ii) ));
            h1_p_lea[j]->SetBinContent(ii, (h2_p_tot[j]->Integral() - h1_p_inl[j]->Integral(1,ii) ));
            h1_d_lea[j]->SetBinContent(ii, (h2_d_tot[j]->Integral() - h1_d_inl[j]->Integral(1,ii) ));
        }
        hC_p_sur[j]->Sumw2();

        auto h2_p_tot_orig = (TH2D*)h2_p_tot[j]->Clone("h2_p_tot_orig");
        auto h2_d_tot_orig = (TH2D*)h2_d_tot[j]->Clone("h2_d_tot_orig");
        auto h1_p_int_orig = (TH1D*)h1_p_int[j]->Clone("h1_p_int_orig");
        auto h1_d_int_orig = (TH1D*)h1_d_int[j]->Clone("h1_d_int_orig");
        auto h1_p_sur_orig = (TH1D*)h1_p_sur[j]->Clone("h1_p_sur_orig");
        auto h1_d_sur_orig = (TH1D*)h1_d_sur[j]->Clone("h1_d_sur_orig");

        for (int i =0; i < 50; i++) // Ratio
        // for (int i =26; i < 27; i++)
        {
            Ratio[i] = 0.01 + i *0.02;
            cout << Ratio[i] << " , " << 1-Ratio[i] <<  endl;

            // Fit for Survive
            TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x+12.75)/[1])", 80,270);fitFunc1->SetLineColor(kRed);  
            TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x+12.75)/[1])", 80,270);fitFunc2->SetLineColor(kBlue); 
            TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*exp(-(x+12.75)/[1])", 80,270);fitFunc3->SetLineColor(kBlack);

            // Fit for Interaction 
            TF1 *fitFunc5 = new TF1("fitFunc5", "25.5*[0]/[1]*exp(-x/[1])", 80,270);fitFunc5->SetLineColor(kRed);  
            TF1 *fitFunc6 = new TF1("fitFunc6", "25.5*[0]/[1]*exp(-x/[1])", 80,270);fitFunc6->SetLineColor(kBlue); 
            TF1 *fitFunc4 = new TF1("fitFunc4", "25.5*[0]/[1]*exp(-x/[1])", 80,270);fitFunc4->SetLineColor(kBlack);
            // 每次循环开始时重置成原始副本的样子
            *h1_p_int[j] = *h1_p_int_orig;
            *h2_p_tot[j] = *h2_p_tot_orig;
            *h1_d_int[j] = *h1_d_int_orig;
            *h2_d_tot[j] = *h2_d_tot_orig;
            *h1_p_sur[j] = *h1_p_sur_orig;
            *h1_d_sur[j] = *h1_d_sur_orig;

            // 然后再做缩放
            h1_p_int[j]->Scale(1 - Ratio[i]);
            h2_p_tot[j]->Scale(1 - Ratio[i]);
            h1_d_int[j]->Scale(Ratio[i]);
            h2_d_tot[j]->Scale(Ratio[i]);

            for(int ii = 1 ; ii <= h1_p_int[j]->GetNbinsX() ; ii++)
            {
                h1_p_sur[j]->SetBinContent(ii, ( h2_p_tot[j]->Integral() - h1_p_int[j]->Integral(1,ii) ) );
                h1_d_sur[j]->SetBinContent(ii, ( h2_d_tot[j]->Integral() - h1_d_int[j]->Integral(1,ii) ) );
            }

            auto h_2_int = (TH1D*)h1_p_int[j]->Clone(); h_2_int->Add(h1_d_int[j]);
            auto h_2_sur = (TH1D*)h1_p_sur[j]->Clone(); h_2_sur->Add(h1_d_sur[j]);

            TLatex latex;
            latex.SetTextSize(0.04);
            latex.SetTextFont(72);
            latex.SetTextAlign(13);  //align at top

            auto range = ROOT::Fit::DataRange();
            range.SetRange(80, 270);
            auto opt = ROOT::Fit::DataOptions();

            auto c0 = new TCanvas("c0","c0",2400,1000);
            gStyle->SetOptStat(0);
            c0->Divide(2,1);
            c0->cd(1);
            h1_p_sur[j]->SetLineColor(kRed);     h1_p_sur[j]->SetLineWidth(2);
            h1_d_sur[j]->SetLineColor(kBlue);    h1_d_sur[j]->SetLineWidth(2);
            h_2_sur->SetLineColor(kBlack);   h_2_sur->SetLineWidth(2);
            h_2_sur->SetTitle("N_{survive};Depth(mm);Counts");
            h_2_sur->GetYaxis()->SetRangeUser(0,1.2*h_2_sur->GetMaximum());
            h_2_sur->Draw("hist");
            h1_p_sur[j]->Draw("histsame");
            h1_d_sur[j]->Draw("histsame");
            fitFunc1->SetParameters(h2_p_tot[j]->Integral() , 200);
            fitFunc2->SetParameters(h2_d_tot[j]->Integral() , 200);
            fitFunc3->SetParameters((h2_p_tot[j]->Integral()+h2_d_tot[j]->Integral()), 200);
            h1_p_sur[j]->Fit(fitFunc1, "RS"); 
            h1_d_sur[j]->Fit(fitFunc2, "QSR"); 
            h_2_sur->Fit(fitFunc3, "RS"); 
            fitFunc1->Draw("same");
            fitFunc2->Draw("same");
            fitFunc3->Draw("same");

            double mixture_length      = fitFunc3->GetParameter(1);
            double mixture_length_err  = fitFunc3->GetParError(1);
            double lambda1     = fitFunc1->GetParameter(1);
            double lambda1_err = fitFunc1->GetParError(1);
            double lambda2     = fitFunc2->GetParameter(1);
            double lambda2_err = fitFunc2->GetParError(1);
            latex.DrawLatex(0,h_2_sur->GetMaximum()-1000,"Function: N_{sur} = N_{tot} #upoint e^{-x/#lambda_{x}}  ");
            latex.DrawLatex(0,h_2_sur->GetMaximum()-2000,Form("Fitted #lambda_{d}: %.2f#pm %.2f mm",mixture_length , mixture_length_err ));
            latex.DrawLatex(0,h_2_sur->GetMaximum()-4000,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda1 , lambda1_err));
            latex.DrawLatex(0,h_2_sur->GetMaximum()-5000,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda2 , lambda2_err));

            TFitResultPtr r1 = h1_p_sur[j]->Fit(fitFunc1, "RS");
            if (r1.Get() && r1->IsValid()) {
                std::cout << "p0 = " << fitFunc1->GetParameter(0)
                        << ", p1 = " << fitFunc1->GetParameter(1) << std::endl;
                fitFunc1->Draw("same");
            } else {
                std::cerr << "Fit1 failed!" << std::endl;
            }
                        
            auto data_sur = ROOT::Fit::BinData(opt, range);
            ROOT::Fit::FillData(data_sur, hC_p_sur[j]); 
            ROOT::Fit::FillData(data_sur, h_2_sur); 
            auto *f1_sur = new TF1("f1_sur","[1]*exp(-(x+12.75)/[0])",80,270);
            f1_sur->SetParameters(200, (h2_p_tot[j]->Integral()+h2_d_tot[j]->Integral()) );
            f1_sur->SetParLimits(0, 0, 500 );
            auto wf_sur = ROOT::Math::WrappedTF1(*f1_sur);
            auto fitter_sur = ROOT::Fit::Fitter() ;
            fitter_sur.SetFunction(wf_sur);
            fitter_sur.Fit(data_sur);
            auto result_sur = fitter_sur.Result();
            result_sur.Print(std::cout);
            cout << result_sur.Chi2() << endl;
            cout << fitFunc1->GetChisquare() << endl;
            cout << fitFunc3->GetChisquare() << endl;
            double Chi2_combine_sur   = result_sur.Chi2();
            double Chi2_proton_sur    = fitFunc1->GetChisquare();
            double Chi2_deuteron_sur  = fitFunc2->GetChisquare();
            double Chi2_mixture_sur   = fitFunc3->GetChisquare();

            auto c0_1 = new TCanvas("c0_1","c0_1",2400,1000);
            c0_1->cd();
            h_2_sur->Draw("hist");
            hC_p_sur[j]->Draw("histsame");
            f1_sur->Draw("same");
    
            auto lg01 = new TLegend(0.72,0.72,0.88,0.88);
            lg01->AddEntry(hC_p_sur[j],"Proton Sample","l");
            lg01->AddEntry(h_2_sur,"Mixture","l");
            lg01->AddEntry(f1_sur ,"Function share same #lambda","l");
            lg01->Draw();

            if( (Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur) ) >=2.71 && Sur_90_Tag < 0) 
            { 
                CL90_Sur[j] = Ratio[i];
                Sur_90_Tag = 1;
            }
            if( (Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur) ) >=3.84 && Sur_95_Tag < 0) 
            { 
                CL95_Sur[j] = Ratio[i];
                Sur_95_Tag = 1;
            }


            c0->cd(2);
            h1_p_int[j]->SetLineColor(kRed);     h1_p_int[j]->SetLineWidth(2);
            h1_d_int[j]->SetLineColor(kBlue);    h1_d_int[j]->SetLineWidth(2);
            h_2_int->SetLineColor(kBlack);   h_2_int->SetLineWidth(2);
            h_2_int->SetTitle("N_{interaction};Depth(mm);Counts");
            h_2_int->GetYaxis()->SetRangeUser(0,1.2*h_2_int->GetMaximum());
            h_2_int->Draw("hist");
            h1_p_int[j]->Draw("histsame");
            h1_d_int[j]->Draw("histsame");
            fitFunc5->SetParameters((h2_p_tot[j]->Integral())                     , 170);
            fitFunc6->SetParameters((h2_d_tot[j]->Integral())                     , 170);
            fitFunc4->SetParameters((h2_p_tot[j]->Integral()+h2_d_tot[j]->Integral()), 170);
            h1_p_int[j]->Fit(fitFunc5,"R");
            h1_d_int[j]->Fit(fitFunc6,"QSR");
            h_2_int->Fit(fitFunc4,"R"); 
            fitFunc4->Draw("same");
            fitFunc5->Draw("same");
            fitFunc6->Draw("same");

            double mixture_length1      = fitFunc4->GetParameter(1);
            double mixture_length1_err  = fitFunc4->GetParError(1);
            double lambda3     = fitFunc5->GetParameter(1);
            double lambda3_err = fitFunc5->GetParError(1);
            double lambda4     = fitFunc6->GetParameter(1);
            double lambda4_err = fitFunc6->GetParError(1);
            latex.DrawLatex(0,h_2_int->GetMaximum()-50,"Function: N_{int} = #frac{ N_{tot} #upoint #Delta x }{#lambda_{x} }  #upoint e^{-x/#lambda_{x} } ");
            latex.DrawLatex(0,h_2_int->GetMaximum()-200,Form("Fitted r_{d}: %.2f#pm %.2f mm",mixture_length1 , mixture_length1_err ));
            latex.DrawLatex(0,h_2_int->GetMaximum()-400,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda3 , lambda3_err));
            latex.DrawLatex(0,h_2_int->GetMaximum()-500,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda4 , lambda4_err));

            auto data_int = ROOT::Fit::BinData(opt, range);
            ROOT::Fit::FillData(data_int, hC_p_int[j]); // NDF = 13
            ROOT::Fit::FillData(data_int, h_2_int); //
            auto *f1_int = new TF1("f1_int","25.5*[1]/[0]*exp(-x/[0])",80,270);
            f1_int->SetParameters(200 , (h2_p_tot[j]->Integral()+h2_d_tot[j]->Integral()));
            f1_int->SetParLimits(0, 0, 500 );
            auto wf_int = ROOT::Math::WrappedTF1(*f1_int);
            auto fitter_int = ROOT::Fit::Fitter() ;
            fitter_int.SetFunction(wf_int);
            fitter_int.Fit(data_int);
            auto result_int = fitter_int.Result();
            result_int.Print(std::cout);
            cout << result_int.Chi2() << endl;
            cout << fitFunc4->GetChisquare() << endl;
            cout << fitFunc5->GetChisquare() << endl;
            double Chi2_combine_int   = result_int.Chi2();
            double Chi2_proton_int    = fitFunc5->GetChisquare();
            double Chi2_deuteron_int  = fitFunc6->GetChisquare();
            double Chi2_mixture_int   = fitFunc4->GetChisquare();

            if( (Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int) ) >=2.71 && Int_90_Tag < 0) 
            { 
                CL90_Int[j] = Ratio[i];
                CL95_Int[j] = Ratio[i];
                Int_90_Tag = 1;
            }
            if( (Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int) ) >=3.84 && Int_95_Tag < 0) 
            { 
                CL95_Int[j] = Ratio[i];
                Int_95_Tag = 1;
            }
            if (Ratio[i] == 0.9 && Int_95_Tag < 0 )
            {
                CL95_Int[j] = Ratio[i];
                Int_95_Tag = 1;
            }

            chi2_int->SetPoint(i,Ratio[i],Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int));
            chi2_sur->SetPoint(i,Ratio[i],Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur));

            gre_sur->SetPoint(i,Ratio[i],mixture_length);
            gre_sur->SetPointError(i,0,mixture_length_err);
            gre_int->SetPoint(i,Ratio[i],mixture_length1);
            gre_int->SetPointError(i,0,mixture_length1_err);

            gre_p_sur->SetPoint(i,Ratio[i],lambda1);
            gre_p_sur->SetPointError(i,0,lambda1_err);
            gre_p_int->SetPoint(i,Ratio[i],lambda3);
            gre_p_int->SetPointError(i,0,lambda3_err);

            gre_d_sur->SetPoint(i,Ratio[i],lambda2);
            gre_d_sur->SetPointError(i,0,lambda2_err);
            gre_d_int->SetPoint(i,Ratio[i],lambda4);
            gre_d_int->SetPointError(i,0,lambda4_err);

        }

        auto c2 = new TCanvas("c2","c2",2100,2100);
        c2->Divide(2,2);
        c2->cd(1);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gre_int->SetLineColor(kBlack);
        gre_int->SetLineWidth(2);
        gre_int->SetMarkerColor(kBlack);
        gre_int->SetMarkerStyle(20);

        gre_p_int->SetLineColor(kRed);
        gre_p_int->SetLineWidth(2);
        gre_p_int->SetMarkerColor(kRed);
        gre_p_int->SetMarkerStyle(20);

        gre_d_int->SetLineColor(kBlue);
        gre_d_int->SetLineWidth(2);
        gre_d_int->SetMarkerColor(kBlue);
        gre_d_int->SetMarkerStyle(20);

        gre_int->GetXaxis()->SetLimits(5e-4,2);
        gre_int->GetYaxis()->SetRangeUser(50,250);
        gre_int->SetTitle("Fitted From N_{int};True r_{d};#lambda (mm)");
        gre_int->Draw("AP");
        gre_p_int->Draw("PSAME");
        gre_d_int->Draw("PSAME");

        auto lg2_1 = new TLegend(0.68,0.12,0.88,0.32);
        lg2_1->AddEntry(gre_int,"Mixture","ep");
        lg2_1->AddEntry(gre_d_int,"Deuteron","ep");
        lg2_1->AddEntry(gre_p_int,"Proton","ep");
        lg2_1->Draw();

        c2->cd(2);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gre_sur->SetLineColor(kBlack);
        gre_sur->SetLineWidth(2);
        gre_sur->SetMarkerColor(kBlack);
        gre_sur->SetMarkerStyle(21);

        gre_p_sur->SetLineColor(kRed);
        gre_p_sur->SetLineWidth(2);
        gre_p_sur->SetMarkerColor(kRed);
        gre_p_sur->SetMarkerStyle(20);

        gre_d_sur->SetLineColor(kBlue);
        gre_d_sur->SetLineWidth(2);
        gre_d_sur->SetMarkerColor(kBlue);
        gre_d_sur->SetMarkerStyle(20);

        gre_sur->GetXaxis()->SetLimits(5e-4,2);
        gre_sur->GetYaxis()->SetRangeUser(50,250);
        gre_sur->SetTitle("Fitted From N_{sur};True r_{d};#lambda (mm)");
        gre_sur->Draw("AP");
        gre_p_sur->Draw("PSAME");
        gre_d_sur->Draw("PSAME");

        auto lg2_2 = new TLegend(0.68,0.12,0.88,0.32);
        lg2_2->AddEntry(gre_sur,"Mixture","ep");
        lg2_2->AddEntry(gre_d_sur,"Deuteron","ep");
        lg2_2->AddEntry(gre_p_sur,"Proton","ep");
        lg2_2->Draw();

        c2->cd(3);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        chi2_int->SetLineColor(kRed);
        chi2_int->SetLineWidth(2);
        chi2_int->GetXaxis()->SetLimits(5e-4,2);
        chi2_int->GetYaxis()->SetRangeUser(0,20);
        chi2_int->SetTitle("Fitted From N_{int};True r_{d};#Delta#chi^{2}");
        chi2_int->Draw("AL");

        c2->cd(4);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        chi2_sur->SetLineColor(kRed);
        chi2_sur->SetLineWidth(2);
        chi2_sur->GetXaxis()->SetLimits(5e-4,2);
        chi2_sur->GetYaxis()->SetRangeUser(0,20);
        chi2_sur->SetTitle("Fitted From N_{sur};True r_{d};#Delta#chi^{2}");
        chi2_sur->Draw("AL");
        // c2->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/FittedRatio_%dGeV_%dGeV.pdf",int(Energy_Dep[j]),int(Energy_Dep[j+1])));

        cout << "Survival   90 CL " << CL90_Sur[j] << " , 95 CL = " << CL95_Sur[j] << endl;
        cout << "Interaction 90 CL " << CL90_Int[j] << " , 95 CL = " << CL95_Int[j] << endl;

        // Fit for Survive
        TF1 *fitFunc1_out = new TF1("fitFunc1_out", "[0]*exp(-(x+12.75)/[1])", 80,270); fitFunc1_out->SetParameters(1e4,200); fitFunc1_out->SetLineColor(kRed);  
        TF1 *fitFunc2_out = new TF1("fitFunc2_out", "[0]*exp(-(x+12.75)/[1])", 80,270); fitFunc2_out->SetParameters(1e4,200); fitFunc2_out->SetLineColor(kBlue); 
        TF1 *fitFunc3_out = new TF1("fitFunc3_out", "[0]*exp(-(x+12.75)/[1])", 80,270); fitFunc3_out->SetParameters(1e4,200); fitFunc3_out->SetLineColor(kBlack);

        // Fit for Interaction 
        TF1 *fitFunc5_out = new TF1("fitFunc5_out", "25.5*[0]/[1]*exp(-x/[1])", 80,270); fitFunc5_out->SetParameters(1e4,200); fitFunc5_out->SetLineColor(kRed);  
        TF1 *fitFunc6_out = new TF1("fitFunc6_out", "25.5*[0]/[1]*exp(-x/[1])", 80,270); fitFunc6_out->SetParameters(1e4,170); fitFunc6_out->SetLineColor(kBlue); 
        TF1 *fitFunc4_out = new TF1("fitFunc4_out", "25.5*[0]/[1]*exp(-x/[1])", 80,270); fitFunc4_out->SetParameters(1e4,170); fitFunc4_out->SetLineColor(kBlack);

        auto c3 = new TCanvas("c3","c3",2400,1600);
        gStyle->SetOptFit(1111);
        c3->cd();
        c3->Divide(3,2);
        c3->cd(1);
        h2_p_tot_orig->SetTitle(Form("Proton Deposit %.f - %.f GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Dep[j],Energy_Dep[j+1]));
        h2_p_tot_orig->GetXaxis()->CenterTitle();
        h2_p_tot_orig->GetYaxis()->CenterTitle();
        h2_p_tot_orig->Draw("colz");

        c3->cd(2);
        gStyle->SetOptFit(1111);
        h1_p_int_orig->GetYaxis()->SetRangeUser(0,h1_p_int_orig->GetMaximum()*1.2);
        h1_p_int_orig->SetLineColor(kBlack);
        h1_p_int_orig->SetTitle("Proton N_{inteaction};Depth (mm);Counts");
        h1_p_int_orig->Draw("hist");
        h1_p_int_orig->Fit(fitFunc5_out,"QSR");  
        fitFunc5_out->Draw("same");
        h1_p_inl[j]->SetLineColor(kRed);
        h1_p_inl[j]->Draw("histsame");
        
        auto lg1 = new TLegend(0.3,0.7,0.6,0.88);
        lg1->AddEntry(h1_p_int_orig,"Estimate Depth","l");
        lg1->AddEntry(h1_p_inl[j],"Inelastic Depth","l");
        lg1->AddEntry(fitFunc5_out,"Fit on Estimated","l");
        lg1->Draw();

        c3->cd(3);
        gStyle->SetOptFit(1111);
        h1_p_sur_orig->GetYaxis()->SetRangeUser(0,h1_p_sur_orig->GetMaximum()*1.2);
        h1_p_sur_orig->SetLineColor(kBlack);
        h1_p_sur_orig->SetTitle("Proton N_{Survive};Depth (mm);Counts");
        h1_p_sur_orig->Draw("hist");
        h1_p_sur_orig->Fit(fitFunc1_out,"QSR");  
        fitFunc1_out->Draw("same");
        h1_p_lea[j]->SetLineColor(kRed);
        h1_p_lea[j]->Draw("histsame");
        lg1->Draw();

        c3->cd(4);
        h2_d_tot_orig->SetTitle(Form("Deuteron Deposit %.f - %.f GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Dep[j],Energy_Dep[j+1]));
        h2_d_tot_orig->GetXaxis()->CenterTitle();
        h2_d_tot_orig->GetYaxis()->CenterTitle();
        h2_d_tot_orig->Draw("colz");

        c3->cd(5);
        gStyle->SetOptFit(1111);
        h1_d_int_orig->GetYaxis()->SetRangeUser(0,h1_d_int_orig->GetMaximum()*1.2);
        h1_d_int_orig->SetLineColor(kBlack);
        h1_d_int_orig->SetTitle("Deuteron N_{inteaction};Depth (mm);Counts");
        h1_d_int_orig->Draw("hist");
        h1_d_int_orig->Fit(fitFunc6_out,"QSR");  
        fitFunc6_out->Draw("same");
        h1_d_inl[j]->SetLineColor(kBlue);
        h1_d_inl[j]->Draw("histsame");

        auto lg2 = new TLegend(0.3,0.7,0.6,0.88);
        lg2->AddEntry(h1_d_int_orig,"Estimate Depth","l");
        lg2->AddEntry(h1_d_inl[j],"Inelastic Depth","l");
        lg2->AddEntry(fitFunc6_out,"Fit on Estimated","l");
        lg2->Draw();

        c3->cd(6);
        gStyle->SetOptFit(1111);
        h1_d_sur_orig->GetYaxis()->SetRangeUser(0,h1_d_sur_orig->GetMaximum()*1.2);
        h1_d_sur_orig->SetLineColor(kBlack);
        h1_d_sur_orig->SetTitle("Deuteron N_{Survive};Depth (mm);Counts");
        h1_d_sur_orig->Draw("hist");
        h1_d_sur_orig->Fit(fitFunc2_out,"QSR");  
        fitFunc2_out->Draw("same");
        h1_d_lea[j]->SetLineColor(kBlue);
        h1_d_lea[j]->Draw("histsame");
        lg2->Draw();

        auto c4 = new TCanvas("c4","c4",2500,1500);
        c4->cd();
        c4->Clear();
        c4->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int k = 0; k <= 14; k++) // layer
        {
            c4->cd(k + 1);
            h1_p_Lay[j][k]->SetLineColor(kRed);   h1_p_Lay[j][k]->SetMarkerColor(kRed);  h1_p_Lay[j][k]->SetLineWidth(2);   
            h1_d_Lay[j][k]->SetLineColor(kBlue);  h1_d_Lay[j][k]->SetMarkerColor(kBlue); h1_d_Lay[j][k]->SetLineWidth(2);   
            // 获取积分，避免除以0
            double p_integral = h1_p_Lay[j][k]->Integral();
            double d_integral = h1_d_Lay[j][k]->Integral();

            // 归一化直方图
            if (p_integral > 0 && d_integral > 0) 
            {
                h1_p_Lay[j][k]->Scale(1.0 / p_integral);
                h1_d_Lay[j][k]->Scale(1.0 / d_integral);
                h1_p_Lay[j][k]->GetYaxis()->SetRangeUser(0,h1_p_Lay[j][k]->GetMaximum()*1.2);
                h1_p_Lay[j][k]->GetXaxis()->CenterTitle();
                h1_p_Lay[j][k]->GetYaxis()->CenterTitle();
                h1_p_Lay[j][k]->GetYaxis()->SetTitleSize(0.05); 
                h1_p_Lay[j][k]->GetYaxis()->SetTitleOffset(0.9); 
                h1_p_Lay[j][k]->GetXaxis()->SetTitleSize(0.05); 
                h1_p_Lay[j][k]->Draw("hist");
                h1_d_Lay[j][k]->Draw("histsame");

                if (k<14) // layer
                {
                    h1_p_Lay[j][k]->Fit(fitFunc_p[j][k],"RSQ");
                    gre_p_bias[j]->SetPoint(k,Layer[k],fitFunc_p[j][k]->GetParameter(1));
                    gre_p_bias[j]->SetPointError(k,Layer_Err[k],fitFunc_p[j][k]->GetParError(1));
                    gre_p_reso[j]->SetPoint(k,Layer[k],fitFunc_p[j][k]->GetParameter(2));
                    gre_p_reso[j]->SetPointError(k,Layer_Err[k],fitFunc_p[j][k]->GetParError(2));
                    fitFunc_p[j][k]->SetLineColor(kRed);
                    fitFunc_p[j][k]->Draw("same");

                    h1_d_Lay[j][k]->Fit(fitFunc_d[j][k],"RSQ");
                    gre_d_bias[j]->SetPoint(k,Layer[k],fitFunc_d[j][k]->GetParameter(1));
                    gre_d_bias[j]->SetPointError(k,Layer_Err[k],fitFunc_d[j][k]->GetParError(1));
                    gre_d_reso[j]->SetPoint(k,Layer[k],fitFunc_d[j][k]->GetParameter(2));
                    gre_d_reso[j]->SetPointError(k,Layer_Err[k],fitFunc_d[j][k]->GetParError(2));
                    fitFunc_d[j][k]->SetLineColor(kBlue);
                    fitFunc_d[j][k]->Draw("same");
                }
            }
            if (k<14) h1_p_Lay[j][k]->SetTitle(Form("Inelastic interaction at Layer %d;(Xine-#hat{Xest})/Xine;Normalized Count",k));
            else h1_p_Lay[j][k]->SetTitle("No Inelastic interaction, Pass through;(Xine-#hat{Xest})/Xine;Normalized Count");
        }

        gre_p_bias[j]->SetLineColor(kRed);
        gre_d_bias[j]->SetLineColor(kBlue);
        gre_p_reso[j]->SetLineColor(kRed);
        gre_d_reso[j]->SetLineColor(kBlue);
        gre_p_bias[j]->SetLineWidth(2);
        gre_d_bias[j]->SetLineWidth(2);
        gre_p_reso[j]->SetLineWidth(2);
        gre_d_reso[j]->SetLineWidth(2);

        gre_p_bias[j]->SetMarkerColor(kRed);
        gre_d_bias[j]->SetMarkerColor(kBlue);
        gre_p_reso[j]->SetMarkerColor(kRed);
        gre_d_reso[j]->SetMarkerColor(kBlue);
        gre_p_bias[j]->SetMarkerStyle(20);
        gre_d_bias[j]->SetMarkerStyle(21);
        gre_p_reso[j]->SetMarkerStyle(20);
        gre_d_reso[j]->SetMarkerStyle(21);

        auto c5 = new TCanvas("c5","c5",2500,1500);
        c5->Divide(2,1);
        c5->cd(1);
        gPad->SetGrid(0,1);
        gre_p_bias[j]->GetYaxis()->SetRangeUser(-1,1);
        gre_p_bias[j]->GetXaxis()->SetLimits(0,14);
        gre_p_bias[j]->SetTitle(";BGO Layer;Bias of (Xine-#hat{Xest})/Xine");
        gre_p_bias[j]->Draw("AP");
        gre_d_bias[j]->Draw("PSAME");
        auto lg5 = new TLegend(0.5,0.78,0.7,0.88);
        lg5->AddEntry(gre_p_bias[j],"Proton","ep");
        lg5->AddEntry(gre_d_bias[j],"Detueron","ep");
        lg5->Draw();

        c5->cd(2);
        gPad->SetGrid(0,1);
        gre_p_reso[j]->GetYaxis()->SetRangeUser(0,1);
        gre_p_reso[j]->GetXaxis()->SetLimits(0,14);
        gre_p_reso[j]->SetTitle(";BGO Layer;Resolution of (Xine-#hat{Xest})/Xine");
        gre_p_reso[j]->Draw("AP");
        gre_d_reso[j]->Draw("PSAME");
        lg5->Draw();

        grL_p_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc5_out->GetParameter(1));
        grL_p_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc5_out->GetParError(1));

        grL_p_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc1_out->GetParameter(1));
        grL_p_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc1_out->GetParError(1));

        grL_d_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc6_out->GetParameter(1));
        grL_d_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc6_out->GetParError(1));

        grL_d_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc2_out->GetParameter(1));
        grL_d_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc2_out->GetParError(1));

        grN_p_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc5_out->GetParameter(0) /  h1_p_int[j]->GetBinWidth(1));
        grN_p_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc5_out->GetParError(0) /  h1_p_int[j]->GetBinWidth(1));

        grN_p_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc1_out->GetParameter(0));
        grN_p_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc1_out->GetParError(0));

        grN_d_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc6_out->GetParameter(0) /  h1_d_int[j]->GetBinWidth(1));
        grN_d_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc6_out->GetParError(0) /  h1_d_int[j]->GetBinWidth(1));

        grN_d_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc2_out->GetParameter(0));
        grN_d_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc2_out->GetParError(0));

        line_Ntot->SetPoint(j,Energy_Mid_p->GetBinContent(j+1), h2_p_tot_orig->Integral());

        Energy_Bin[j] = Energy_Mid_p->GetBinContent(j+1);

        auto write_file = new TFile(Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_Res_%d.root",j), "RECREATE");
        write_file->cd();
        h2_p_tot[j]->Write();
        h2_d_tot[j]->Write();
        h1_p_int[j]->Write();
        h1_d_int[j]->Write();
        h1_p_inl[j]->Write();
        h1_d_inl[j]->Write();
    
        h1_p_sur[j]->Write();
        h1_d_sur[j]->Write();
        h1_p_lea[j]->Write();
        h1_d_lea[j]->Write();
        c2->Write();
        c3->Write();
        c4->Write();
        c5->Write();
        write_file->Close();

    } // J Deposit Energy
    grL_p_int->SetLineColor(kRed);  grN_p_int->SetLineColor(kRed); 
    grL_p_sur->SetLineColor(kRed);  grN_p_sur->SetLineColor(kRed); 
    grL_d_int->SetLineColor(kBlue); grN_d_int->SetLineColor(kBlue);
    grL_d_sur->SetLineColor(kBlue); grN_d_sur->SetLineColor(kBlue);

    grL_p_int->SetLineWidth(2);     grN_p_int->SetLineWidth(2);
    grL_p_sur->SetLineWidth(2);     grN_p_sur->SetLineWidth(2);
    grL_d_int->SetLineWidth(2);     grN_d_int->SetLineWidth(2);
    grL_d_sur->SetLineWidth(2);     grN_d_sur->SetLineWidth(2);

    grL_p_int->SetMarkerColor(kRed);   grN_p_int->SetMarkerColor(kRed); 
    grL_p_sur->SetMarkerColor(kRed);   grN_p_sur->SetMarkerColor(kRed); 
    grL_d_int->SetMarkerColor(kBlue);  grN_d_int->SetMarkerColor(kBlue);
    grL_d_sur->SetMarkerColor(kBlue);  grN_d_sur->SetMarkerColor(kBlue);

    grL_p_int->SetMarkerStyle(20);    grN_p_int->SetMarkerStyle(20);
    grL_p_sur->SetMarkerStyle(24);    grN_p_sur->SetMarkerStyle(24);
    grL_d_int->SetMarkerStyle(21);    grN_d_int->SetMarkerStyle(21);
    grL_d_sur->SetMarkerStyle(25);    grN_d_sur->SetMarkerStyle(25);

    auto write_file2 = new TFile("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_Energy_Res.root", "RECREATE");
    auto c6 = new TCanvas("c6","c6",2500,1200);
    c6->Divide(2,1);
    c6->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogx(1);
    grL_p_int->GetYaxis()->SetRangeUser(0,250);
    grL_p_int->GetXaxis()->SetLimits(4,2e4);
    grL_p_int->SetTitle("Fitted From N_{int};Energy(GeV);#lambda Based on Estimated Interaction Depth");
    grL_p_int->Draw("AP");
    grL_d_int->Draw("PSAME");
    gre6->Draw("LSAME");
    gre7->Draw("LSAME");

    auto lg6_1 = new TLegend(0.48,0.12,0.88,0.18);
    lg6_1->SetNColumns(2);
    lg6_1->AddEntry(grL_p_int,"Proton Fitted N_{int}" ,"pe");
    lg6_1->AddEntry(grL_d_int,"Deuteron Fitted N_{int}" ,"pe");
    lg6_1->AddEntry(gre6,"Proton GEANT4" ,"l");
    lg6_1->AddEntry(gre7,"Deuteron GEANT4" ,"l");
    lg6_1->Draw();

    c6->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogx(2);
    grL_p_sur->GetYaxis()->SetRangeUser(0,250);
    grL_p_sur->GetXaxis()->SetLimits(4,2e4);
    grL_p_sur->SetTitle("Fitted From N_{sur};Energy(GeV);#lambda Based on Estimated Interaction Depth");
    grL_p_sur->Draw("AP");
    grL_d_sur->Draw("PSAME");
    gre6->Draw("LSAME");
    gre7->Draw("LSAME");
    
    auto lg6_2 = new TLegend(0.48,0.12,0.88,0.18);
    lg6_2->SetNColumns(2);
    lg6_2->AddEntry(grL_p_sur,"Proton Fitted N_{sur}" ,"pe");
    lg6_2->AddEntry(grL_d_sur,"Deuteron Fitted N_{sur}" ,"pe");
    lg6_2->AddEntry(gre6,"Proton GEANT4" ,"l");
    lg6_2->AddEntry(gre7,"Deuteron GEANT4" ,"l");
    lg6_2->Draw();


    auto c7 = new TCanvas("c7","c7",2500,1200);
    c7->Divide(2,1);
    c7->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    double ymin_int = 8e2;
    double ymax_int = 2e8;
    grN_p_int->GetYaxis()->SetRangeUser(ymin_int, ymax_int);
    grN_p_int->GetXaxis()->SetLimits(4,2e4);
    grN_p_int->SetTitle("Fitted From N_{int};Energy(GeV);N_{0} Based on Estimated Interaction Depth");
    grN_p_int->Draw("AP");
    grN_d_int->Draw("PSAME");
    line_Ntot->Draw("Lsame");
    auto lg7_1 = new TLegend(0.48,0.72,0.88,0.88);
    lg7_1->SetNColumns(2);
    lg7_1->AddEntry(grN_p_int,"Proton Fitted N_{int}" ,"pe");
    lg7_1->AddEntry(grN_d_int,"Deuteron Fitted N_{int}" ,"pe");
    lg7_1->AddEntry(line_Ntot,"Ture Simulated N_{0}" ,"l");
    lg7_1->Draw();

    c7->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    double ymin_sur = 8e2;
    double ymax_sur = 2e8;
    grN_p_sur->GetYaxis()->SetRangeUser(ymin_sur, ymax_sur);
    grN_p_sur->GetXaxis()->SetLimits(4,2e4);
    grN_p_sur->SetTitle("Fitted From N_{sur};Energy(GeV);N_{0} Based on Estimated Interaction Depth");
    grN_p_sur->Draw("AP");
    grN_d_sur->Draw("PSAME");
    line_Ntot->Draw("Lsame");
    auto lg7_2 = new TLegend(0.48,0.72,0.88,0.88);
    lg7_2->SetNColumns(2);
    lg7_2->AddEntry(grN_p_sur,"Proton Fitted N_{sur}" ,"pe");
    lg7_2->AddEntry(grN_d_sur,"Deuteron Fitted N_{sur}" ,"pe");
    lg7_2->AddEntry(line_Ntot,"Ture Simulated N_{0}" ,"l");
    lg7_2->Draw();

    auto gre_90_sur = new TGraph(7,Energy_Bin,CL90_Sur); //,Energy_Name_Err,CL90_Sur_Err
    auto gre_95_sur = new TGraph(7,Energy_Bin,CL95_Sur); //,Energy_Name_Err,CL95_Sur_Err
    auto gre_90_int = new TGraph(7,Energy_Bin,CL90_Int); //,Energy_Name_Err,CL90_Sur_Err
    auto gre_95_int = new TGraph(7,Energy_Bin,CL95_Int); //,Energy_Name_Err,CL95_Sur_Err

    gre_90_sur->SetLineColor(kBlue);
    gre_95_sur->SetLineColor(kRed);
    gre_90_sur->SetLineWidth(2);
    gre_95_sur->SetLineWidth(2);

    gre_90_int->SetLineColor(kBlue);
    gre_95_int->SetLineColor(kRed);
    gre_90_int->SetLineWidth(2);
    gre_95_int->SetLineWidth(2);

    gre_90_int->SetName("gre_90_int");
    gre_95_int->SetName("gre_95_int");

    auto c8 = new TCanvas("c8","c8",2500,1200);
    c8->Divide(2,1);
    c8->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    gPad->SetLogx();
    gre_90_sur->GetYaxis()->SetRangeUser(9e-3,2);
    gre_90_sur->GetXaxis()->SetLimits(9e0,2e4);
    gre_90_sur->SetTitle("Sensitivity Curve from N_{sur};Energy(GeV); True r_{d}");
    gre_90_sur->Draw("ALP");
    gre_95_sur->Draw("LPSAME");
    auto lg8 = new TLegend(0.72,0.12,0.88,0.22);
    lg8->AddEntry(gre_90_sur,"90% C.L.","l");
    lg8->AddEntry(gre_95_sur,"95% C.L.","l");
    lg8->Draw();

    c8->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    gPad->SetLogx();
    gre_90_int->GetYaxis()->SetRangeUser(9e-3,2);
    gre_90_int->GetXaxis()->SetLimits(9e0,2e4);
    gre_90_int->SetTitle("Sensitivity Curve from N_{int};Energy(GeV); True r_{d}");
    gre_90_int->Draw("ALP");
    gre_95_int->Draw("LPSAME");
    lg8->Draw();


    write_file2->cd();
    gre_90_int->Write();
    gre_95_int->Write();
    c6->Write();
    c7->Write();
    c8->Write();
    write_file2->Close();
}