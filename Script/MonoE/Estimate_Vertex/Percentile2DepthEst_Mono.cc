#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 \hat{x_ine} =  - 4.789 * sigma + x_mid

void Percentile2DepthEst_Mono()
{
    int Energy_Name[28]={0};
    int Energy_Name_Err[28]={0};
    double Ratio[28] = {0};
    double CL95_Sur[28]  = {0};    double CL95_Sur_Err[28]  = {0};    
    double CL90_Sur[28]  = {0};    double CL90_Sur_Err[28]  = {0};    
    double CL95_Int[28]  = {0};    double CL95_Int_Err[28]  = {0};    
    double CL90_Int[28]  = {0};    double CL90_Int_Err[28]  = {0};    
    double Energy[28]= {0};    double Energy_Err[28]= {0};
    double Ratio_Edge[29];
    double Energy_Edge[29];
    auto mean_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Estimate_Vertex/DepthEst_Fitting.root");
    auto hist_p = (TH1D*)mean_file->Get("hist_p");
    for (int j = 0; j < 29; j++)
    {
        if (j <= 10)                 { Energy_Edge[j] = 5 + 10 * j;             Ratio_Edge[j] = 0.0005 + 0.001 * j;     }   // Center: 10~100 → edges: 5~105
        else if (j > 10 && j <= 19)  { Energy_Edge[j] = 150 + 100 * (j - 10);   Ratio_Edge[j] = 0.015 + 0.01 * (j - 10);      }   // Center: 10~100 → edges: 5~105  
        else                         { Energy_Edge[j] = 1500 + 1000 * (j - 19); Ratio_Edge[j] = 0.15 + 0.1 * (j - 19); } // Center: 2000~10000 → edges: 1500~10500     
    }

    auto lambda_p_sur = new TH1D("lambda_p_sur","lambda_p_sur",28,Energy_Edge);
    auto lambda_p_int = new TH1D("lambda_p_int","lambda_p_int",28,Energy_Edge);
    auto lambda_d_int = new TH1D("lambda_d_int","lambda_d_int",28,Energy_Edge);
    auto lambda_d_sur = new TH1D("lambda_d_sur","lambda_d_sur",28,Energy_Edge);

    auto N_p_sur = new TH1D("N_p_sur","N_p_sur",28,Energy_Edge);
    auto N_p_int = new TH1D("N_p_int","N_p_int",28,Energy_Edge);
    auto N_d_int = new TH1D("N_d_int","N_d_int",28,Energy_Edge);
    auto N_d_sur = new TH1D("N_d_sur","N_d_sur",28,Energy_Edge);

    for (int k =0; k < 28; k++) //. Energy
    // for (int k =1; k < 2; k++)
    // for (int k =23; k < 24; k++)
    {
        if (k < 10)      { Energy_Name[k] = (k + 1) * 10;            Energy[k] = (k + 1) * 10.;       }       // 10 ~ 100
        else if (k < 19) { Energy_Name[k] = (k - 9 + 1) * 100;       Energy[k] = (k - 9 + 1) * 100.;  }        // 200 ~ 1000
        else             { Energy_Name[k] = (k - 18 + 1) * 1000;     Energy[k] = (k - 18 + 1) * 1000.;}        // 2000 ~ 10000
        double BEST_FIT_MEAN = hist_p->GetBinContent(k+1);

        CL95_Sur[k]=1e-3; int Sur_95_Tag = -1;
        CL90_Sur[k]=1e-3; int Sur_90_Tag = -1;
        CL95_Int[k]=1e-3; int Int_95_Tag = -1;
        CL90_Int[k]=1e-3; int Int_90_Tag = -1;


        auto gre_int = new TGraphErrors();  
        auto gre_sur = new TGraphErrors();
        auto gre_p_int = new TGraphErrors();  
        auto gre_p_sur = new TGraphErrors();
        auto gre_d_int = new TGraphErrors();  
        auto gre_d_sur = new TGraphErrors();
        auto chi2_int = new TGraph();
        auto chi2_sur = new TGraph();

        cout << "Energy = " << Energy_Name[k] << " GeV ! Best Fit for Proton Estimated factor = " << BEST_FIT_MEAN << endl;
        int p_FH_Lay; int p_FH_Type; double p_Total_E;      int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
        int d_FH_Lay; int d_FH_Type; double d_Total_E;      int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
        int p_FI_Lay;    double p_FI_Dep;    int p_particle;
        int d_FI_Lay;    double d_FI_Dep;    int d_particle;

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",Energy_Name[k]));
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
        proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
        proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",Energy_Name[k]));
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
        deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
        deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);

        double Layer[14]={0};
        double Layer_Err[14]={0};

        // Depsit and Layer
        TH1D *h1_p_Lay[15]; // h1_p_Lay[14] 用来装Passthrough的事例
        TH1D *h1_d_Lay[15]; // h1_d_Lay[14] 用来装Passthrough的事例
        TF1 *fitFunc_p[14];  auto gre_p_reso = new TGraphErrors();  auto gre_p_bias = new TGraphErrors();
        TF1 *fitFunc_d[14];  auto gre_d_reso = new TGraphErrors();  auto gre_d_bias = new TGraphErrors();

        auto h2_p_tot = new TH2D("h2_p_tot","h2_p_tot",200,-25,375,325,-275,375);  h2_p_tot->Sumw2();
        auto h2_d_tot = new TH2D("h2_d_tot","h2_d_tot",200,-25,375,325,-275,375);  h2_d_tot->Sumw2();
        auto h1_p_int = new TH1D("h1_p_int","h1_p_int",14,0,357);                  h1_p_int->Sumw2(); // EST
        auto h1_d_int = new TH1D("h1_d_int","h1_d_int",14,0,357);                  h1_d_int->Sumw2(); // EST
        auto h1_p_sur = new TH1D("h1_p_sur","h1_p_sur",14,0,357);                  // h1_p_sur->Sumw2();   // EST
        auto h1_d_sur = new TH1D("h1_d_sur","h1_d_sur",14,0,357);                  // h1_d_sur->Sumw2();   // EST
        auto hC_p_int = new TH1D("hC_p_int","hC_p_int",14,0,357);                    hC_p_int->Sumw2();
        auto hC_p_sur = new TH1D("hC_p_sur","hC_p_sur",14,0,357);      

        auto h1_p_inl = new TH1D("h1_p_inl","h1_p_inl",14,0,357);                  h1_p_inl->Sumw2(); // REAL
        auto h1_d_inl = new TH1D("h1_d_inl","h1_d_inl",14,0,357);                  h1_d_inl->Sumw2(); // REAL
        auto h1_p_lea = new TH1D("h1_p_lea","h1_p_lea",14,0,357);                  // h1_p_lea->Sumw2(); // REAL
        auto h1_d_lea = new TH1D("h1_d_lea","h1_d_lea",14,0,357);                  // h1_d_lea->Sumw2(); // REAL
     
        for(int i =0 ; i<15 ; i++)  // Layer
        {
            h1_p_Lay[i] =new TH1D(Form("h1_p_Lay[%d]",i),Form("h1_p_Lay[%d]",i),30,-3,3); h1_p_Lay[i]->Sumw2(); 
            h1_d_Lay[i] =new TH1D(Form("h1_d_Lay[%d]",i),Form("h1_d_Lay[%d]",i),30,-3,3); h1_d_Lay[i]->Sumw2();
            if (i<14)
            {
                Layer[i] = 0.5 + i;
                Layer_Err[i] = 0.5;
                fitFunc_p[i] = new TF1(Form("fitFunc_p[%d]",i), "gaus", -2, 2); // 注意替换范围
                fitFunc_d[i] = new TF1(Form("fitFunc_d[%d]",i), "gaus", -2, 2); // 注意替换范围
            }
        }
        

        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
        {        
            proton_tree->GetEntry(entry);
            if (entry%1000==0) cout << " Proton : " << entry << endl;
            int p_energy_index = int(floor((log10(p_Total_E) ) / 0.2));
            if (p_energy_index < 0 || p_energy_index > 20) continue;
            if (p_FI_Dep < 0) {p_FI_Lay = 14;}
            if (p_Nhits < 10 ) continue;
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
            // if ( Est_Depth < -40 ) cout << "Watch out " << entry << " Value " << Est_Depth <<endl;
            h2_p_tot->Fill(p_FI_Dep, Est_Depth);
            h1_p_Lay[p_FI_Lay]->Fill( (p_FI_Dep -  Est_Depth) / p_FI_Dep );
            h1_p_int->Fill(Est_Depth);
            h1_p_inl->Fill(p_FI_Dep);
            hC_p_int->Fill(Est_Depth);
        }

        // cout << " Propare Deuteron" << endl; // deuteron_tree->GetEntries()
        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries() ; ++entry)
        {
            deuteron_tree->GetEntry(entry);
            if (entry%1000==0) cout << " Deuteron : " << entry << endl;
            int d_energy_index = int(floor((log10(d_Total_E) ) / 0.2));
            if (d_energy_index < 0 || d_energy_index > 20) continue;
            if (d_FI_Dep < 0) {d_FI_Lay = 14;}
            if (d_Nhits < 10 ) continue;
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
            h2_d_tot->Fill(d_FI_Dep, Est_Depth);
            h1_d_Lay[d_FI_Lay]->Fill( (d_FI_Dep -  Est_Depth) / d_FI_Dep );
            h1_d_int->Fill(Est_Depth);
            h1_d_inl->Fill(d_FI_Dep);
        }

        cout << h2_p_tot->Integral() << endl;
        cout << h2_d_tot->Integral() << endl;

        // Fit for Survive
        TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x+12.75)/[1])", 80,270); fitFunc1->SetParameters(1e4,200); fitFunc1->SetLineColor(kRed);  
        TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x+12.75)/[1])", 80,270); fitFunc2->SetParameters(1e4,200); fitFunc2->SetLineColor(kBlue); 
        TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*exp(-(x+12.75)/[1])", 80,270); fitFunc3->SetParameters(1e4,200); fitFunc3->SetLineColor(kBlack);

        // Fit for Interaction 
        TF1 *fitFunc5 = new TF1("fitFunc5", "25.5*[0]/[1]*exp(-x/[1])", 80,270); fitFunc5->SetParameters(1e4,200); fitFunc5->SetLineColor(kRed);  
        TF1 *fitFunc6 = new TF1("fitFunc6", "25.5*[0]/[1]*exp(-x/[1])", 80,270); fitFunc6->SetParameters(1e4,170); fitFunc6->SetLineColor(kBlue); 
        TF1 *fitFunc4 = new TF1("fitFunc4", "25.5*[0]/[1]*exp(-x/[1])", 80,270); fitFunc4->SetParameters(1e4,170); fitFunc4->SetLineColor(kBlack);

        for(int ii = 1 ; ii <= hC_p_int->GetNbinsX() ; ii++)
        {
            hC_p_sur->SetBinContent(ii ,(h2_p_tot->Integral() - hC_p_int->Integral(1,ii) ));
            h1_p_sur->SetBinContent(ii, (h2_p_tot->Integral() - h1_p_int->Integral(1,ii) ));
            h1_d_sur->SetBinContent(ii, (h2_d_tot->Integral() - h1_d_int->Integral(1,ii) ));
            h1_p_lea->SetBinContent(ii, (h2_p_tot->Integral() - h1_p_inl->Integral(1,ii) ));
            h1_d_lea->SetBinContent(ii, (h2_d_tot->Integral() - h1_d_inl->Integral(1,ii) ));
                cout << ii << " , " << hC_p_sur->GetBinContent(ii) << endl;
        }
        hC_p_sur->Sumw2();

        auto h2_p_tot_orig = (TH2D*)h2_p_tot->Clone("h2_p_tot_orig");
        auto h2_d_tot_orig = (TH2D*)h2_d_tot->Clone("h2_d_tot_orig");

        auto h1_p_int_orig = (TH1D*)h1_p_int->Clone("h1_p_int_orig");
        auto h1_d_int_orig = (TH1D*)h1_d_int->Clone("h1_d_int_orig");

        auto h1_p_sur_orig = (TH1D*)h1_p_sur->Clone("h1_p_sur_orig");
        auto h1_d_sur_orig = (TH1D*)h1_d_sur->Clone("h1_d_sur_orig");


        for (int i =0; i < 27; i++)
        // for (int i =18; i < 19; i++)
        {
            if      (i < 10)          Ratio[i] = (i + 1) * 0.001;     
            else if (i < 19)          Ratio[i] = (i - 9 + 1) * 0.01;
            else                      Ratio[i] = (i - 18 + 1) * 0.1;  
            cout << Ratio[i] << " , " << 1-Ratio[i] <<  endl;
            // 每次循环开始时重置成原始副本的样子
            *h1_p_int = *h1_p_int_orig;
            *h2_p_tot = *h2_p_tot_orig;
            *h1_d_int = *h1_d_int_orig;
            *h2_d_tot = *h2_d_tot_orig;
            *h1_p_sur = *h1_p_sur_orig;
            *h1_d_sur = *h1_d_sur_orig;

            // 然后再做缩放
            h1_p_int->Scale(1 - Ratio[i]);
            h2_p_tot->Scale(1 - Ratio[i]);
            h1_d_int->Scale(Ratio[i]);
            h2_d_tot->Scale(Ratio[i]);

            for(int ii = 1 ; ii <= h1_p_int->GetNbinsX() ; ii++)
            {
                h1_p_sur->SetBinContent(ii, ( h2_p_tot->Integral() - h1_p_int->Integral(1,ii) ) );
                h1_d_sur->SetBinContent(ii, ( h2_d_tot->Integral() - h1_d_int->Integral(1,ii) ) );
            }

            auto h_2_int = (TH1D*)h1_p_int->Clone(); h_2_int->Add(h1_d_int);
            auto h_2_sur = (TH1D*)h1_p_sur->Clone(); h_2_sur->Add(h1_d_sur);

            // h_2_int->Sumw2();    h1_p_int->Sumw2();
            // h_2_sur->Sumw2();    h1_d_sur->Sumw2();

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
            h1_p_sur->SetLineColor(kRed);     h1_p_sur->SetLineWidth(2);
            h1_d_sur->SetLineColor(kBlue);    h1_d_sur->SetLineWidth(2);
            h_2_sur->SetLineColor(kBlack);   h_2_sur->SetLineWidth(2);
            h_2_sur->SetTitle("N_{survive};Depth(mm);Counts");
            h_2_sur->GetYaxis()->SetRangeUser(0,1.2*h_2_sur->GetMaximum());
            h_2_sur->Draw("hist");
            h1_p_sur->Draw("histsame");
            h1_d_sur->Draw("histsame");
            fitFunc1->SetParameters(h2_p_tot->Integral()                       , 200);
            fitFunc2->SetParameters(h2_d_tot->Integral()                       , 200);
            fitFunc3->SetParameters((h2_p_tot->Integral()+h2_d_tot->Integral()), 200);
            h1_p_sur->Fit(fitFunc1, "QSR"); 
            h1_d_sur->Fit(fitFunc2, "QSR"); 
            h_2_sur->Fit(fitFunc3, "QSR"); 
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
            
            auto data_sur = ROOT::Fit::BinData(opt, range);
            ROOT::Fit::FillData(data_sur, hC_p_sur); 
            ROOT::Fit::FillData(data_sur, h_2_sur); 
            auto *f1_sur = new TF1("f1_sur","[1]*exp(-(x+12.75)/[0])",80,270);
            f1_sur->SetParameters(200, (h2_p_tot->Integral()+h2_d_tot->Integral()) );
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
            hC_p_sur->Draw("histsame");
            f1_sur->Draw("same");
    
            auto lg01 = new TLegend(0.72,0.72,0.88,0.88);
            lg01->AddEntry(hC_p_sur,"Proton Sample","l");
            lg01->AddEntry(h_2_sur,"Mixture","l");
            lg01->AddEntry(f1_sur ,"Function share same #lambda","l");
            lg01->Draw();

            if( (Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur) ) >=2.71 && Sur_90_Tag < 0) 
            { 
                CL90_Sur[k] = Ratio[i];
                Sur_90_Tag = 1;
            }
            if( (Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur) ) >=3.84 && Sur_95_Tag < 0) 
            { 
                CL95_Sur[k] = Ratio[i];
                Sur_95_Tag = 1;
            }


            c0->cd(2);
            h1_p_int->SetLineColor(kRed);     h1_p_int->SetLineWidth(2);
            h1_d_int->SetLineColor(kBlue);    h1_d_int->SetLineWidth(2);
            h_2_int->SetLineColor(kBlack);   h_2_int->SetLineWidth(2);
            h_2_int->SetTitle("N_{interaction};Depth(mm);Counts");
            h_2_int->GetYaxis()->SetRangeUser(0,1.2*h_2_int->GetMaximum());
            h_2_int->Draw("hist");
            h1_p_int->Draw("histsame");
            h1_d_int->Draw("histsame");
            fitFunc5->SetParameters((h2_p_tot->Integral())                     , 170);
            fitFunc6->SetParameters((h2_d_tot->Integral())                     , 170);
            fitFunc4->SetParameters((h2_p_tot->Integral()+h2_d_tot->Integral()), 170);
            h1_p_int->Fit(fitFunc5,"QSR");
            h1_d_int->Fit(fitFunc6,"QSR");
            h_2_int->Fit(fitFunc4,"QSR"); 
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
            ROOT::Fit::FillData(data_int, hC_p_int); // NDF = 13
            ROOT::Fit::FillData(data_int, h_2_int); //
            auto *f1_int = new TF1("f1_int","25.5*[1]/[0]*exp(-x/[0])",80,270);
            f1_int->SetParameters(200 , (h2_p_tot->Integral()+h2_d_tot->Integral()));
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
                CL90_Int[k] = Ratio[i];
                CL95_Int[k] = Ratio[i];
                Int_90_Tag = 1;
            }
            if( (Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int) ) >=3.84 && Int_95_Tag < 0) 
            { 
                CL95_Int[k] = Ratio[i];
                Int_95_Tag = 1;
            }
            if (Ratio[i] == 0.9 && Int_95_Tag < 0 )
            {
                CL95_Int[k] = Ratio[i];
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
        c2->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Estimate_Vertex/FittedRatio_%dGeV.pdf",Energy_Name[k]));

        cout << "Survival   90 CL " << CL90_Sur[k] << " , 95 CL = " << CL95_Sur[k] << endl;
        cout << "Interaction 90 CL " << CL90_Int[k] << " , 95 CL = " << CL95_Int[k] << endl;

        auto c3 = new TCanvas("c3","c3",2400,1600);
        gStyle->SetOptFit(1111);
        c3->cd();
        c3->Divide(3,2);
        c3->cd(1);
        h2_p_tot_orig->SetTitle(Form("Incident Proton %d GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Name[k]));
        h2_p_tot_orig->GetXaxis()->CenterTitle();
        h2_p_tot_orig->GetYaxis()->CenterTitle();
        h2_p_tot_orig->Draw("colz");

        c3->cd(2);
        gStyle->SetOptFit(1111);
        h1_p_int_orig->GetYaxis()->SetRangeUser(0,h1_p_int_orig->GetMaximum()*1.2);
        h1_p_int_orig->SetLineColor(kBlack);
        h1_p_int_orig->SetTitle("Proton N_{inteaction};Depth (mm);Counts");
        h1_p_int_orig->Draw("hist");
        fitFunc5->SetParameters(1e4, 200);
        h1_p_int_orig->Fit(fitFunc5,"QSR");  
        fitFunc5->Draw("same");
        h1_p_inl->SetLineColor(kRed);
        h1_p_inl->Draw("histsame");
        
        auto lg1 = new TLegend(0.3,0.7,0.6,0.88);
        lg1->AddEntry(h1_p_int_orig,"Estimate Depth","l");
        lg1->AddEntry(h1_p_inl,"Inelastic Depth","l");
        lg1->AddEntry(fitFunc5,"Fit on Estimated","l");
        lg1->Draw();

        c3->cd(3);
        gStyle->SetOptFit(1111);
        h1_p_sur_orig->GetYaxis()->SetRangeUser(0,h1_p_sur_orig->GetMaximum()*1.2);
        h1_p_sur_orig->SetLineColor(kBlack);
        h1_p_sur_orig->SetTitle("Proton N_{Survive};Depth (mm);Counts");
        h1_p_sur_orig->Draw("hist");
        fitFunc1->SetParameters(1e4,200);
        h1_p_sur_orig->Fit(fitFunc1,"QSR"); 
        fitFunc1->Draw("same");
        h1_p_lea->SetLineColor(kRed);
        h1_p_lea->Draw("histsame");
        lg1->Draw();

        c3->cd(4);
        h2_d_tot_orig->SetTitle(Form("Incident Deuteron %d GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Name[k]));
        h2_d_tot_orig->GetXaxis()->CenterTitle();
        h2_d_tot_orig->GetYaxis()->CenterTitle();
        h2_d_tot_orig->Draw("colz");

        c3->cd(5);
        gStyle->SetOptFit(1111);
        h1_d_int_orig->GetYaxis()->SetRangeUser(0,h1_d_int_orig->GetMaximum()*1.2);
        h1_d_int_orig->SetLineColor(kBlack);
        h1_d_int_orig->SetTitle("Deuteron N_{inteaction};Depth (mm);Counts");
        h1_d_int_orig->Draw("hist");
        fitFunc6->SetParameters(1e4, 200);
        h1_d_int_orig->Fit(fitFunc6,"QSR");  
        fitFunc6->Draw("same");
        h1_d_inl->SetLineColor(kBlue);
        h1_d_inl->Draw("histsame");

        auto lg2 = new TLegend(0.3,0.7,0.6,0.88);
        lg2->AddEntry(h1_d_int_orig,"Estimate Depth","l");
        lg2->AddEntry(h1_d_inl,"Inelastic Depth","l");
        lg2->AddEntry(fitFunc6,"Fit on Estimated","l");
        lg2->Draw();

        c3->cd(6);
        gStyle->SetOptFit(1111);
        h1_d_sur_orig->GetYaxis()->SetRangeUser(0,h1_d_sur_orig->GetMaximum()*1.2);
        h1_d_sur_orig->SetLineColor(kBlack);
        h1_d_sur_orig->SetTitle("Deuteron N_{Survive};Depth (mm);Counts");
        h1_d_sur_orig->Draw("hist");
        fitFunc2->SetParameters(1e4,170);
        h1_d_sur_orig->Fit(fitFunc2,"QSR");  
        fitFunc2->Draw("same");
        h1_d_lea->SetLineColor(kBlue);
        h1_d_lea->Draw("histsame");
        lg2->Draw();

        lambda_p_int->SetBinContent(k+1, fitFunc5->GetParameter(1)); 
        lambda_p_sur->SetBinContent(k+1, fitFunc1->GetParameter(1)); 
        lambda_d_int->SetBinContent(k+1, fitFunc6->GetParameter(1)); 
        lambda_d_sur->SetBinContent(k+1, fitFunc2->GetParameter(1)); 

        lambda_p_int->SetBinError(k+1, fitFunc5->GetParError(1));
        lambda_p_sur->SetBinError(k+1, fitFunc1->GetParError(1));
        lambda_d_int->SetBinError(k+1, fitFunc6->GetParError(1));
        lambda_d_sur->SetBinError(k+1, fitFunc2->GetParError(1));

        N_p_int->SetBinContent(k+1, fitFunc5->GetParameter(0) / h1_p_int->GetBinWidth(1));
        N_p_sur->SetBinContent(k+1, fitFunc1->GetParameter(0)  );
        N_d_int->SetBinContent(k+1, fitFunc6->GetParameter(0) / h1_d_int->GetBinWidth(1));
        N_d_sur->SetBinContent(k+1, fitFunc2->GetParameter(0)  );

        N_p_int->SetBinError(k+1, fitFunc5->GetParError(0) / h1_p_int->GetBinWidth(1));
        N_p_sur->SetBinError(k+1, fitFunc1->GetParError(0));
        N_d_int->SetBinError(k+1, fitFunc6->GetParError(0) / h1_d_int->GetBinWidth(1));
        N_d_sur->SetBinError(k+1, fitFunc2->GetParError(0));

        // Xine = Xmid - 4.879 * Slope
        auto c4 = new TCanvas("c4","c4",2500,1500);
        c4->cd();
        c4->Clear();
        c4->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j <= 14; j++) // layer
        {
            c4->cd(j + 1);
            h1_p_Lay[j]->SetLineColor(kRed);   h1_p_Lay[j]->SetMarkerColor(kRed);  h1_p_Lay[j]->SetLineWidth(2);   
            h1_d_Lay[j]->SetLineColor(kBlue);  h1_d_Lay[j]->SetMarkerColor(kBlue); h1_d_Lay[j]->SetLineWidth(2);   
            h1_p_Lay[j]->Scale(1.0/h1_p_Lay[j]->Integral()); 
            h1_d_Lay[j]->Scale(1.0/h1_d_Lay[j]->Integral()); 
            h1_p_Lay[j]->GetYaxis()->SetRangeUser(0,h1_p_Lay[j]->GetMaximum()*1.2);
            h1_p_Lay[j]->GetXaxis()->CenterTitle();
            h1_p_Lay[j]->GetYaxis()->CenterTitle();
            h1_p_Lay[j]->Draw("hist");
            h1_d_Lay[j]->Draw("histsame");
            if (j<14)
            {
                h1_p_Lay[j]->Fit(fitFunc_p[j],"RSQ");
                h1_d_Lay[j]->Fit(fitFunc_d[j],"RSQ");
                gre_p_bias->SetPoint(j,Layer[j],fitFunc_p[j]->GetParameter(1));
                gre_d_bias->SetPoint(j,Layer[j],fitFunc_d[j]->GetParameter(1));
                gre_p_bias->SetPointError(j,Layer_Err[j],fitFunc_p[j]->GetParError(1));
                gre_d_bias->SetPointError(j,Layer_Err[j],fitFunc_d[j]->GetParError(1));
                gre_p_reso->SetPoint(j,Layer[j],fitFunc_p[j]->GetParameter(2));
                gre_d_reso->SetPoint(j,Layer[j],fitFunc_d[j]->GetParameter(2));
                gre_p_reso->SetPointError(j,Layer_Err[j],fitFunc_p[j]->GetParError(2));
                gre_d_reso->SetPointError(j,Layer_Err[j],fitFunc_d[j]->GetParError(2));
                h1_p_Lay[j]->SetTitle(Form("Inelastic interaction at Layer %d;(Xine-#hat{Xest})/Xine;Normalized Count",j));
                fitFunc_p[j]->SetLineColor(kRed);
                fitFunc_d[j]->SetLineColor(kBlue);
                fitFunc_p[j]->Draw("same");
                fitFunc_d[j]->Draw("same");
            }
            else h1_p_Lay[j]->SetTitle("No Inelastic interaction, Pass through;(Xine-#hat{Xest})/Xine;Normalized Count");
            h1_p_Lay[j]->GetYaxis()->SetTitleSize(0.05); 
            h1_p_Lay[j]->GetYaxis()->SetTitleOffset(0.9); 
            h1_p_Lay[j]->GetXaxis()->SetTitleSize(0.05); 
        }


        gre_p_bias->SetLineColor(kRed);
        gre_d_bias->SetLineColor(kBlue);
        gre_p_reso->SetLineColor(kRed);
        gre_d_reso->SetLineColor(kBlue);
        gre_p_bias->SetLineWidth(2);
        gre_d_bias->SetLineWidth(2);
        gre_p_reso->SetLineWidth(2);
        gre_d_reso->SetLineWidth(2);

        gre_p_bias->SetMarkerColor(kRed);
        gre_d_bias->SetMarkerColor(kBlue);
        gre_p_reso->SetMarkerColor(kRed);
        gre_d_reso->SetMarkerColor(kBlue);
        gre_p_bias->SetMarkerStyle(20);
        gre_d_bias->SetMarkerStyle(21);
        gre_p_reso->SetMarkerStyle(20);
        gre_d_reso->SetMarkerStyle(21);

        auto c5 = new TCanvas("c5","c5",2500,1500);
        c5->Divide(2,1);
        c5->cd(1);
        gPad->SetGrid(0,1);
        gre_p_bias->GetYaxis()->SetRangeUser(-1,1);
        gre_p_bias->GetXaxis()->SetLimits(0,14);
        gre_p_bias->SetTitle(";BGO Layer;Bias of (Xine-#hat{Xest})/Xine");
        gre_p_bias->Draw("AP");
        gre_d_bias->Draw("PSAME");
        auto lg5 = new TLegend(0.5,0.78,0.7,0.88);
        lg5->AddEntry(gre_p_bias,"Proton","ep");
        lg5->AddEntry(gre_d_bias,"Detueron","ep");
        lg5->Draw();

        c5->cd(2);
        gPad->SetGrid(0,1);
        gre_p_reso->GetYaxis()->SetRangeUser(0,1);
        gre_p_reso->GetXaxis()->SetLimits(0,14);
        gre_p_reso->SetTitle(";BGO Layer;Resolution of (Xine-#hat{Xest})/Xine");
        gre_p_reso->Draw("AP");
        gre_d_reso->Draw("PSAME");
        lg5->Draw();

        auto write_file = new TFile(Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Estimate_Vertex/DepthEst_%dGeV.root",Energy_Name[k]), "RECREATE");
        write_file->cd();
        h2_p_tot->Write();
        h2_d_tot->Write();

        h1_p_int->Write();
        h1_d_int->Write();
        h1_p_inl->Write();
        h1_d_inl->Write();

        h1_p_sur->Write();
        h1_d_sur->Write();
        h1_p_lea->Write();
        h1_d_lea->Write();
        c2->Write();
        c3->Write();
        c4->Write();
        c5->Write();
        write_file->Close();
    } //. Energy k

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


    TGraphErrors* grL_p_int = new TGraphErrors(lambda_p_int);     TGraphErrors* grN_p_int = new TGraphErrors(N_p_int);
    TGraphErrors* grL_p_sur = new TGraphErrors(lambda_p_sur);     TGraphErrors* grN_p_sur = new TGraphErrors(N_p_sur);
    TGraphErrors* grL_d_int = new TGraphErrors(lambda_d_int);     TGraphErrors* grN_d_int = new TGraphErrors(N_d_int);
    TGraphErrors* grL_d_sur = new TGraphErrors(lambda_d_sur);     TGraphErrors* grN_d_sur = new TGraphErrors(N_d_sur);

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


    auto c6 = new TCanvas("c6","c6",2500,1200);
    c6->Divide(2,1);
    c6->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogx(1);
    grL_p_int->GetYaxis()->SetRangeUser(100,250);
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
    grL_p_sur->GetYaxis()->SetRangeUser(100,250);
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
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    grN_p_int->GetYaxis()->SetRangeUser(8e3,2e4);
    grN_p_int->GetXaxis()->SetLimits(4,2e4);
    grN_p_int->SetTitle("Fitted From N_{int};Energy(GeV);N_{0} Based on Estimated Interaction Depth");
    grN_p_int->Draw("AP");
    grN_d_int->Draw("PSAME");
    auto line_Ntot = new TLine(4,1e4,2e4,1e4);
    line_Ntot->Draw("same");
    auto lg7_1 = new TLegend(0.48,0.72,0.88,0.88);
    lg7_1->SetNColumns(2);
    lg7_1->AddEntry(grN_p_int,"Proton Fitted N_{int}" ,"pe");
    lg7_1->AddEntry(grN_d_int,"Deuteron Fitted N_{int}" ,"pe");
    lg7_1->AddEntry(line_Ntot,"Ture Simulated N_{0}" ,"l");
    lg7_1->Draw();

    c7->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    grN_p_sur->GetYaxis()->SetRangeUser(8e3,2e4);
    grN_p_sur->GetXaxis()->SetLimits(4,2e4);
    grN_p_sur->SetTitle("Fitted From N_{sur};Energy(GeV);N_{0} Based on Estimated Interaction Depth");
    grN_p_sur->Draw("AP");
    grN_d_sur->Draw("PSAME");
    line_Ntot->Draw("same");
    auto lg7_2 = new TLegend(0.48,0.72,0.88,0.88);
    lg7_2->SetNColumns(2);
    lg7_2->AddEntry(grN_p_sur,"Proton Fitted N_{sur}" ,"pe");
    lg7_2->AddEntry(grN_d_sur,"Deuteron Fitted N_{sur}" ,"pe");
    lg7_2->AddEntry(line_Ntot,"Ture Simulated N_{0}" ,"l");
    lg7_2->Draw();

    auto gre_90_sur = new TGraph(28,Energy,CL90_Sur); //,Energy_Name_Err,CL90_Sur_Err
    auto gre_95_sur = new TGraph(28,Energy,CL95_Sur); //,Energy_Name_Err,CL95_Sur_Err

    auto gre_90_int = new TGraph(28,Energy,CL90_Int); //,Energy_Name_Err,CL90_Sur_Err
    auto gre_95_int = new TGraph(28,Energy,CL95_Int); //,Energy_Name_Err,CL95_Sur_Err

    gre_90_sur->SetLineColor(kBlue);
    gre_95_sur->SetLineColor(kRed);
    gre_90_sur->SetLineWidth(2);
    gre_95_sur->SetLineWidth(2);

    gre_90_int->SetLineColor(kBlue);
    gre_95_int->SetLineColor(kRed);
    gre_90_int->SetLineWidth(2);
    gre_95_int->SetLineWidth(2);

    auto c8 = new TCanvas("c8","c8",2500,1200);
    c8->Divide(2,1);
    c8->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    gPad->SetLogx();
    gre_90_sur->GetYaxis()->SetRangeUser(9e-4,2);
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
    gre_90_int->GetYaxis()->SetRangeUser(9e-4,2);
    gre_90_int->GetXaxis()->SetLimits(9e0,2e4);
    gre_90_int->SetTitle("Sensitivity Curve from N_{int};Energy(GeV); True r_{d}");
    gre_90_int->Draw("ALP");
    gre_95_int->Draw("LPSAME");
    lg8->Draw();

    auto write_file2 = new TFile("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Estimate_Vertex/DepthEst_Energy.root", "RECREATE");
    write_file2->cd();
    c6->Write();
    c7->Write();
    c8->Write();
    write_file2->Close();
}