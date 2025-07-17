#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 \hat{x_ine} =  - 4.789 * sigma + x_mid

void Percentile2DepthEst_Dep()
{
    double Ratio[18] = {0};
    double Ratio_Edge[19];
    double Energy_Dep[8] = {8,19,43,100,320,1e3,3.2e3,1e4};

    auto mean_file    = TFile::Open("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/EnergyEst_Dep.root");
    auto hist_p_2     = (TH1D*)mean_file->Get("hist_p_2");
    auto Energy_Mid_p = (TH1D*)mean_file->Get("Energy_Mid_p");

    for (int j = 0; j < 19; j++)
    {
        if (j <= 10)       Ratio_Edge[j] = 0.005 + 0.01 * j;        // Center: 10~100 → edges: 5~105
        else               Ratio_Edge[j] = 0.15 + 0.1 * (j - 10); // Center: 2000~10000 → edges: 1500~10500    
        // cout << Ratio_Edge[j] << endl;
    }

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
    deuteron_tree->SetBranchAddress("Energy"           ,&d_Energy);
    deuteron_tree->SetBranchAddress("Nhits"            ,&d_Nhits);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);

    TGraphErrors* gre_p_int = new TGraphErrors(); TGraphErrors* grN_p_int = new TGraphErrors();
    TGraphErrors* gre_p_sur = new TGraphErrors(); TGraphErrors* grN_p_sur = new TGraphErrors();
    TGraphErrors* gre_d_int = new TGraphErrors(); TGraphErrors* grN_d_int = new TGraphErrors();
    TGraphErrors* gre_d_sur = new TGraphErrors(); TGraphErrors* grN_d_sur = new TGraphErrors();
    
    double Layer[14]={0};
    double Layer_Err[14]={0};
    TF1 *fitFunc_p[7][14];  
    TF1 *fitFunc_d[7][14];  
    TH1D *h1_p_Lay[7][15]; // h1_p_Lay[x][14] 用来装Passthrough的事例
    TH1D *h1_d_Lay[7][15]; // h1_d_Lay[x][14] 用来装Passthrough的事例
    TH2D *h2_p_int[7]; 
    TH2D *h2_d_int[7]; 
    TH1D *h1_p_int[7]; 
    TH1D *h1_d_int[7]; 
    TH1D *h1_p_inl[7]; 
    TH1D *h1_d_inl[7]; 
    TH1D *h1_p_sur[7]; 
    TH1D *h1_d_sur[7]; 
    TH1D *h1_p_lea[7]; 
    TH1D *h1_d_lea[7]; 
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
        h2_p_int[j] = new TH2D(Form("h2_p_int[%d]",j),Form("h2_p_int[%d]",j),200,-25,375,325,-275,375);
        h2_d_int[j] = new TH2D(Form("h2_d_int[%d]",j),Form("h2_d_int[%d]",j),200,-25,375,325,-275,375);
        h1_p_int[j] = new TH1D(Form("h1_p_int[%d]",j),Form("h1_p_int[%d]",j),14,0,357);  // EST
        h1_d_int[j] = new TH1D(Form("h1_d_int[%d]",j),Form("h1_d_int[%d]",j),14,0,357);  // EST
        h1_p_inl[j] = new TH1D(Form("h1_p_inl[%d]",j),Form("h1_p_inl[%d]",j),14,0,357);  // REAL
        h1_d_inl[j] = new TH1D(Form("h1_d_inl[%d]",j),Form("h1_d_inl[%d]",j),14,0,357);  // REAL
        h1_p_sur[j] = new TH1D(Form("h1_p_sur[%d]",j),Form("h1_p_sur[%d]",j),14,0,357);  // EST
        h1_d_sur[j] = new TH1D(Form("h1_d_sur[%d]",j),Form("h1_d_sur[%d]",j),14,0,357);  // EST
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
                fitFunc_p[j][i] = new TF1(Form("fitFunc_p[%d][%d]",j, i), "gaus", -2, 2); // 注意替换范围
                fitFunc_d[j][i] = new TF1(Form("fitFunc_d[%d][%d]",j, i), "gaus", -2, 2); // 注意替换范围
            }
        }
    }

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    // for (Long64_t entry = 0; entry < 50000; entry++)
    {        
        proton_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Proton : " << entry << endl;
        int p_energy_index = int(floor((log10(p_Total_E) ) / 0.2));
        if (p_energy_index < 0 || p_energy_index > 20) continue;
        if (p_FI_Dep < 0) {p_FI_Lay = 14;}
        if (p_Nhits < 10 ) continue;
        if (p_Total_E < hist_p_2->GetXaxis()->GetXmin() || p_Total_E > hist_p_2->GetXaxis()->GetXmax()) continue;
        int BEST_FIT_MEAN_BIN  = hist_p_2->FindBin(p_Total_E);
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
        h2_p_int[BEST_FIT_MEAN_BIN-1]->Fill(p_FI_Dep, Est_Depth, p_weight*1e4);   
        h1_p_Lay[BEST_FIT_MEAN_BIN-1][p_FI_Lay]->Fill( (p_FI_Dep -  Est_Depth) / p_FI_Dep , p_weight*1e4);  // bias and reso
        h1_p_int[BEST_FIT_MEAN_BIN-1]->Fill(Est_Depth, p_weight*1e4);
        h1_p_inl[BEST_FIT_MEAN_BIN-1]->Fill(p_FI_Dep, p_weight*1e4);

    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries() ; ++entry)
    // for (Long64_t entry = 0; entry < 50000 ; ++entry)
    {
        deuteron_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Deuteron : " << entry << endl;
        int d_energy_index = int(floor((log10(d_Total_E) ) / 0.2));
        if (d_energy_index < 0 || d_energy_index > 20) continue;
        if (d_FI_Dep < 0) {d_FI_Lay = 14;}
        if (d_Nhits < 10 ) continue;
        if (d_Total_E < hist_p_2->GetXaxis()->GetXmin() || d_Total_E > hist_p_2->GetXaxis()->GetXmax()) continue;
        int BEST_FIT_MEAN_BIN = hist_p_2->FindBin(d_Total_E);
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
        h2_d_int[BEST_FIT_MEAN_BIN-1]->Fill(d_FI_Dep, Est_Depth, d_weight*1e4);
        h1_d_Lay[BEST_FIT_MEAN_BIN-1][d_FI_Lay]->Fill( (d_FI_Dep -  Est_Depth) / d_FI_Dep , d_weight*1e4);
        h1_d_int[BEST_FIT_MEAN_BIN-1]->Fill(Est_Depth, d_weight*1e4);
        h1_d_inl[BEST_FIT_MEAN_BIN-1]->Fill(d_FI_Dep, d_weight*1e4);
    }

    for(int j=0 ;j< 7 ; j++)
    {
        h2_p_int[j]->Sumw2();    h1_p_int[j]->Sumw2();
        h2_d_int[j]->Sumw2();    h1_d_int[j]->Sumw2();

        // N_sur Fitting Function
        TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x+12.75)/[1])", 40,150); 
        TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x+12.75)/[1])", 40,150); 
        fitFunc1->SetParameters(h2_p_int[j]->Integral(), 200); 
        fitFunc2->SetParameters(h2_d_int[j]->Integral(), 170); 
        fitFunc1->SetLineColor(kRed); 
        fitFunc2->SetLineColor(kBlue);
        fitFunc1->SetLineStyle(2);
        fitFunc2->SetLineStyle(2);

        // N_int Fitting Function
        TF1 *fitFunc3 = new TF1("fitFunc3", "[0]/[1]* exp(-x/[1])", 40,150); 
        TF1 *fitFunc4 = new TF1("fitFunc4", "[0]/[1]* exp(-x/[1])", 40,150); 
        fitFunc3->SetParameters(h2_p_int[j]->Integral() * h1_p_int[j]->GetBinWidth(1), 200); 
        fitFunc4->SetParameters(h2_d_int[j]->Integral() * h1_d_int[j]->GetBinWidth(1), 170); 
        fitFunc3->SetLineColor(kRed); 
        fitFunc4->SetLineColor(kBlue);
        fitFunc3->SetLineStyle(2);
        fitFunc4->SetLineStyle(2);

        for(int i=1 ; i<= h1_p_int[j]->GetNbinsX(); i++)
        {
            h1_p_sur[j]->SetBinContent(i, (h2_p_int[j]->Integral() - h1_p_int[j]->Integral(1,i) ));
            h1_d_sur[j]->SetBinContent(i, (h2_d_int[j]->Integral() - h1_d_int[j]->Integral(1,i) ));
            h1_p_lea[j]->SetBinContent(i, (h2_p_int[j]->Integral() - h1_p_inl[j]->Integral(1,i) ));
            h1_d_lea[j]->SetBinContent(i, (h2_d_int[j]->Integral() - h1_d_inl[j]->Integral(1,i) ));
            // cout << h1_p_int[j]->Integral(1,i) << " , " <<  h1_p_int[j]->GetBinContent(i) << endl;
            // cout << h1_d_int->Integral(0,i) << endl;
        }
        cout << h1_p_int[j]->Integral() << " , " <<  h1_p_inl[j]->Integral() << endl;


        auto c0 = new TCanvas("c0","c0",2400,1600);
        gStyle->SetOptFit(1111);
        c0->cd();
        c0->Divide(3,2);
        c0->cd(1);
        h2_p_int[j]->SetTitle(Form("Proton Deposit %.f - %.f GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Dep[j],Energy_Dep[j+1]));
        h2_p_int[j]->GetXaxis()->CenterTitle();
        h2_p_int[j]->GetYaxis()->CenterTitle();
        h2_p_int[j]->Draw("colz");

        c0->cd(2);
        gStyle->SetOptFit(1111);
        h1_p_int[j]->GetYaxis()->SetRangeUser(0,h1_p_int[j]->GetMaximum()*1.2);
        h1_p_int[j]->SetLineColor(kBlack);
        h1_p_int[j]->SetTitle("Proton N_{inteaction};Depth (mm);Counts");
        h1_p_int[j]->Draw("hist");
        h1_p_int[j]->Fit(fitFunc3,"R");  
        fitFunc3->Draw("same");
        h1_p_inl[j]->SetLineColor(kRed);
        h1_p_inl[j]->Draw("histsame");
        
        auto lg1 = new TLegend(0.3,0.7,0.6,0.88);
        lg1->AddEntry(h1_p_int[j],"Estimate Depth","l");
        lg1->AddEntry(h1_p_inl[j],"Inelastic Depth","l");
        lg1->AddEntry(fitFunc3,"Fit on Estimated","l");
        lg1->Draw();

        c0->cd(3);
        gStyle->SetOptFit(1111);
        h1_p_sur[j]->GetYaxis()->SetRangeUser(0,h1_p_sur[j]->GetMaximum()*1.2);
        h1_p_sur[j]->SetLineColor(kBlack);
        h1_p_sur[j]->SetTitle("Proton N_{Survive};Depth (mm);Counts");
        h1_p_sur[j]->Draw("hist");
        h1_p_sur[j]->Fit(fitFunc1,"R");  
        fitFunc1->Draw("same");
        h1_p_lea[j]->SetLineColor(kRed);
        h1_p_lea[j]->Draw("histsame");
        lg1->Draw();


        c0->cd(4);
        h2_d_int[j]->SetTitle(Form("Deuteron Deposit %.f - %.f GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Dep[j],Energy_Dep[j+1]));
        h2_d_int[j]->GetXaxis()->CenterTitle();
        h2_d_int[j]->GetYaxis()->CenterTitle();
        h2_d_int[j]->Draw("colz");

        c0->cd(5);
        gStyle->SetOptFit(1111);
        h1_d_int[j]->GetYaxis()->SetRangeUser(0,h1_d_int[j]->GetMaximum()*1.2);
        h1_d_int[j]->SetLineColor(kBlack);
        h1_d_int[j]->SetTitle("Deuteron N_{inteaction};Depth (mm);Counts");
        h1_d_int[j]->Draw("hist");
        h1_d_int[j]->Fit(fitFunc4,"R");  
        fitFunc4->Draw("same");
        h1_d_inl[j]->SetLineColor(kBlue);
        h1_d_inl[j]->Draw("histsame");

        auto lg2 = new TLegend(0.3,0.7,0.6,0.88);
        lg2->AddEntry(h1_d_int[j],"Estimate Depth","l");
        lg2->AddEntry(h1_d_inl[j],"Inelastic Depth","l");
        lg2->AddEntry(fitFunc4,"Fit on Estimated","l");
        lg2->Draw();

        c0->cd(6);
        gStyle->SetOptFit(1111);
        h1_d_sur[j]->GetYaxis()->SetRangeUser(0,h1_d_sur[j]->GetMaximum()*1.2);
        h1_d_sur[j]->SetLineColor(kBlack);
        h1_d_sur[j]->SetTitle("Deuteron N_{Survive};Depth (mm);Counts");
        h1_d_sur[j]->Draw("hist");
        h1_d_sur[j]->Fit(fitFunc2,"R");  
        fitFunc2->Draw("same");
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

        gre_p_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc3->GetParameter(1));
        gre_p_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc3->GetParError(1));

        gre_p_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc1->GetParameter(1));
        gre_p_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc1->GetParError(1));

        gre_d_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc4->GetParameter(1));
        gre_d_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc4->GetParError(1));

        gre_d_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc2->GetParameter(1));
        gre_d_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc2->GetParError(1));

        grN_p_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc3->GetParameter(0) /  h1_p_int[j]->GetBinWidth(1));
        grN_p_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc3->GetParError(0) /  h1_p_int[j]->GetBinWidth(1));

        grN_p_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc1->GetParameter(0));
        grN_p_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc1->GetParError(0));

        grN_d_int->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc4->GetParameter(0) /  h1_d_int[j]->GetBinWidth(1));
        grN_d_int->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc4->GetParError(0) /  h1_d_int[j]->GetBinWidth(1));

        grN_d_sur->SetPoint(j,Energy_Mid_p->GetBinContent(j+1),fitFunc2->GetParameter(0));
        grN_d_sur->SetPointError(j,Energy_Mid_p->GetBinError(j+1),fitFunc2->GetParError(0));

        line_Ntot->SetPoint(j,Energy_Mid_p->GetBinContent(j+1), h2_p_int[j]->Integral());

        auto write_file = new TFile(Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_Dep_%d.root",j), "RECREATE");
        write_file->cd();
        h2_p_int[j]->Write();
        h2_d_int[j]->Write();
        h1_p_int[j]->Write();
        h1_d_int[j]->Write();
        h1_p_inl[j]->Write();
        h1_d_inl[j]->Write();
    
        h1_p_sur[j]->Write();
        h1_d_sur[j]->Write();
        h1_p_lea[j]->Write();
        h1_d_lea[j]->Write();
        c0->Write();
        c4->Write();
        c5->Write();
        write_file->Close();

    } // J Deposit Energy
    gre_p_int->SetLineColor(kRed);  grN_p_int->SetLineColor(kRed); 
    gre_p_sur->SetLineColor(kRed);  grN_p_sur->SetLineColor(kRed); 
    gre_d_int->SetLineColor(kBlue); grN_d_int->SetLineColor(kBlue);
    gre_d_sur->SetLineColor(kBlue); grN_d_sur->SetLineColor(kBlue);

    gre_p_int->SetLineWidth(2);     grN_p_int->SetLineWidth(2);
    gre_p_sur->SetLineWidth(2);     grN_p_sur->SetLineWidth(2);
    gre_d_int->SetLineWidth(2);     grN_d_int->SetLineWidth(2);
    gre_d_sur->SetLineWidth(2);     grN_d_sur->SetLineWidth(2);

    gre_p_int->SetMarkerColor(kRed);   grN_p_int->SetMarkerColor(kRed); 
    gre_p_sur->SetMarkerColor(kRed);   grN_p_sur->SetMarkerColor(kRed); 
    gre_d_int->SetMarkerColor(kBlue);  grN_d_int->SetMarkerColor(kBlue);
    gre_d_sur->SetMarkerColor(kBlue);  grN_d_sur->SetMarkerColor(kBlue);

    gre_p_int->SetMarkerStyle(20);    grN_p_int->SetMarkerStyle(20);
    gre_p_sur->SetMarkerStyle(24);    grN_p_sur->SetMarkerStyle(24);
    gre_d_int->SetMarkerStyle(21);    grN_d_int->SetMarkerStyle(21);
    gre_d_sur->SetMarkerStyle(25);    grN_d_sur->SetMarkerStyle(25);

    auto c6 = new TCanvas("c6","c6",2500,1200);
    c6->Divide(2,1);
    c6->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogx(1);
    gre_p_int->GetYaxis()->SetRangeUser(0,250);
    gre_p_int->GetXaxis()->SetLimits(4,2e4);
    gre_p_int->SetTitle("Fitted From N_{int};Energy(GeV);#lambda Based on Estimated Interaction Depth");
    gre_p_int->Draw("AP");
    gre_d_int->Draw("PSAME");
    gre6->Draw("LSAME");
    gre7->Draw("LSAME");

    auto lg6_1 = new TLegend(0.48,0.12,0.88,0.18);
    lg6_1->SetNColumns(2);
    lg6_1->AddEntry(gre_p_int,"Proton Fitted N_{int}" ,"pe");
    lg6_1->AddEntry(gre_d_int,"Deuteron Fitted N_{int}" ,"pe");
    lg6_1->AddEntry(gre6,"Proton GEANT4" ,"l");
    lg6_1->AddEntry(gre7,"Deuteron GEANT4" ,"l");
    lg6_1->Draw();

    c6->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogx(2);
    gre_p_sur->GetYaxis()->SetRangeUser(0,250);
    gre_p_sur->GetXaxis()->SetLimits(4,2e4);
    gre_p_sur->SetTitle("Fitted From N_{sur};Energy(GeV);#lambda Based on Estimated Interaction Depth");
    gre_p_sur->Draw("AP");
    gre_d_sur->Draw("PSAME");
    gre6->Draw("LSAME");
    gre7->Draw("LSAME");
    
    auto lg6_2 = new TLegend(0.48,0.12,0.88,0.18);
    lg6_2->SetNColumns(2);
    lg6_2->AddEntry(gre_p_sur,"Proton Fitted N_{sur}" ,"pe");
    lg6_2->AddEntry(gre_d_sur,"Deuteron Fitted N_{sur}" ,"pe");
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

    auto write_file2 = new TFile("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst.root", "RECREATE");
    write_file2->cd();
    c6->Write();
    c7->Write();
    write_file2->Close();
}