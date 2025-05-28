#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 \hat{x_ine} =  - 4.789 * sigma + x_mid

void Percentile2DepthEst()
{
    int Energy_Name[28]={0};     
    auto mean_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_Fitting.root");
    auto hist_p = (TH1D*)mean_file->Get("hist_p");
    for (int k =18; k < 19; k++)
    {
        if (k < 10)      Energy_Name[k] = (k + 1) * 10;               // 10 ~ 100
        else if (k < 19) Energy_Name[k] = (k - 9 + 1) * 100;           // 200 ~ 1000
        else             Energy_Name[k] = (k - 18 + 1) * 1000;         // 2000 ~ 10000
        double BEST_FIT_MEAN = hist_p->GetBinContent(k+1);
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

        auto h2_p_int = new TH2D("h2_p_int","h2_p_int",200,-25,375,325,-275,375);  
        auto h2_d_int = new TH2D("h2_d_int","h2_d_int",200,-25,375,325,-275,375);  
        auto h1_p_int = new TH1D("h1_p_int","h1_p_int",14,0,357);  
        auto h1_d_int = new TH1D("h1_d_int","h1_d_int",14,0,357);  
        auto h1_p_inl = new TH1D("h1_p_inl","h1_p_inl",14,0,357);  
        auto h1_d_inl = new TH1D("h1_d_inl","h1_d_inl",14,0,357);  
        auto h1_p_sur = new TH1D("h1_p_sur","h1_p_sur",14,0,357);  
        auto h1_d_sur = new TH1D("h1_d_sur","h1_d_sur",14,0,357);  
        auto h1_p_lea = new TH1D("h1_p_lea","h1_p_lea",14,0,357);  
        auto h1_d_lea = new TH1D("h1_d_lea","h1_d_lea",14,0,357);  

        for(int i =0 ; i<15 ; i++)  // Layer
        {
            h1_p_Lay[i] =new TH1D(Form("h1_p_Lay[%d]",i),Form("h1_p_Lay[%d]",i),60,-3,3);  
            h1_d_Lay[i] =new TH1D(Form("h1_d_Lay[%d]",i),Form("h1_d_Lay[%d]",i),60,-3,3); 
            if (i<14)
            {
                Layer[i] = 0.5 + i;
                Layer_Err[i] = 0.5;
                fitFunc_p[i] = new TF1(Form("fitFunc_p[%d]",i), "gaus", -1, 1); // 注意替换范围
                fitFunc_d[i] = new TF1(Form("fitFunc_d[%d]",i), "gaus", -1, 1); // 注意替换范围
            }
        }
        

        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
        {        
            proton_tree->GetEntry(entry);
            if (entry%1000==0) cout << " Proton : " << entry << endl;
            int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
            // if (p_energy_index < 0 || p_energy_index > 14) continue;
            if (p_FI_Dep < 0) {p_FI_Lay = 14;}
            // if (p_Nhits < 10 ) continue;
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
            h2_p_int->Fill(p_FI_Dep, Est_Depth);
            h1_p_Lay[p_FI_Lay]->Fill( (p_FI_Dep -  Est_Depth) / p_FI_Dep );
            h1_p_int->Fill(Est_Depth);
            h1_p_inl->Fill(p_FI_Dep);
        }

        // cout << " Propare Deuteron" << endl; // deuteron_tree->GetEntries()
        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries() ; ++entry)
        {
            deuteron_tree->GetEntry(entry);
            if (entry%1000==0) cout << " Deuteron : " << entry << endl;
            int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
            // if (d_energy_index < 0 || d_energy_index > 14) continue;
            if (d_FI_Dep < 0) {d_FI_Lay = 14;}
            // if (d_Nhits < 10 ) continue;
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
            h2_d_int->Fill(d_FI_Dep, Est_Depth);
            h1_d_Lay[d_FI_Lay]->Fill( (d_FI_Dep -  Est_Depth) / d_FI_Dep );
            h1_d_int->Fill(Est_Depth);
            h1_d_inl->Fill(d_FI_Dep);
        }

        h2_p_int->Sumw2();    h1_p_int->Sumw2();
        h2_d_int->Sumw2();    h1_d_int->Sumw2();

        // N_sur Fitting Function
        TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x)/[1])", 80,250); 
        TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x)/[1])", 80,250); 
        fitFunc1->SetParameters(1e4, 200); 
        fitFunc2->SetParameters(1e4, 170); 
        fitFunc1->SetLineColor(kRed); 
        fitFunc2->SetLineColor(kBlue);
        fitFunc1->SetLineStyle(2);
        fitFunc2->SetLineStyle(2);
        // fitFunc1->FixParameter(0, 1e4);
        // fitFunc2->FixParameter(0, 1e4);

        // N_int Fitting Function
        TF1 *fitFunc3 = new TF1("fitFunc3", "[0]/[1]* exp(-x/[1])", 80,250); 
        TF1 *fitFunc4 = new TF1("fitFunc4", "[0]/[1]* exp(-x/[1])", 80,250); 
        fitFunc3->SetParameters(1e4 * h1_p_int->GetBinWidth(1), 200); 
        fitFunc4->SetParameters(1e4 * h1_d_int->GetBinWidth(1), 170); 
        cout << h1_p_int->GetBinWidth(1) << endl;
        fitFunc3->SetLineColor(kRed); 
        fitFunc4->SetLineColor(kBlue);
        fitFunc3->SetLineStyle(2);
        fitFunc4->SetLineStyle(2);
        // fitFunc3->FixParameter(0, 1e4 * h1_p_int->GetBinWidth(1) );
        // fitFunc4->FixParameter(0, 1e4 * h1_d_int->GetBinWidth(1) );

        cout << h2_p_int->Integral() << endl;
        cout << h2_d_int->Integral() << endl;
        for(int i=1 ; i<= h1_p_int->GetNbinsX(); i++)
        {
            h1_p_sur->SetBinContent(i, (h2_p_int->Integral() - h1_p_int->Integral(1,i) ));
            h1_d_sur->SetBinContent(i, (h2_d_int->Integral() - h1_d_int->Integral(1,i) ));
            h1_p_lea->SetBinContent(i, (h2_p_int->Integral() - h1_p_inl->Integral(1,i) ));
            h1_d_lea->SetBinContent(i, (h2_d_int->Integral() - h1_d_inl->Integral(1,i) ));
            cout << h1_p_int->Integral(1,i) << " , " <<  h1_p_int->GetBinContent(i) << endl;

            // cout << h1_d_int->Integral(0,i) << endl;
        }

        auto c3 = new TCanvas("c3","c3",2400,1600);
        gStyle->SetOptFit(1111);
        c3->cd();
        c3->Divide(3,2);
        c3->cd(1);
        h2_p_int->SetTitle(Form("Incident Proton %d GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Name[k]));
        h2_p_int->GetXaxis()->CenterTitle();
        h2_p_int->GetYaxis()->CenterTitle();
        h2_p_int->Draw("colz");

        c3->cd(2);
        gStyle->SetOptFit(1111);
        h1_p_int->GetYaxis()->SetRangeUser(0,h1_p_int->GetMaximum()*1.2);
        h1_p_int->SetLineColor(kBlack);
        h1_p_int->SetTitle("Proton N_{inteaction};Depth (mm);Counts");
        h1_p_int->Draw("hist");
        h1_p_int->Fit(fitFunc3,"R");  
        fitFunc3->Draw("same");
        h1_p_inl->SetLineColor(kRed);
        h1_p_inl->Draw("histsame");
        
        auto lg1 = new TLegend(0.3,0.7,0.6,0.88);
        lg1->AddEntry(h1_p_int,"Estimate Depth","l");
        lg1->AddEntry(h1_p_inl,"Inelastic Depth","l");
        lg1->AddEntry(fitFunc3,"Fit on Estimated","l");
        lg1->Draw();

        c3->cd(3);
        gStyle->SetOptFit(1111);
        h1_p_sur->GetYaxis()->SetRangeUser(0,h1_p_sur->GetMaximum()*1.2);
        h1_p_sur->SetLineColor(kBlack);
        h1_p_sur->SetTitle("Proton N_{Survive};Depth (mm);Counts");
        h1_p_sur->Draw("hist");
        h1_p_sur->Fit(fitFunc1,"R");  
        fitFunc1->Draw("same");
        h1_p_lea->SetLineColor(kRed);
        h1_p_lea->Draw("histsame");

        lg1->Draw();


        c3->cd(4);
        h2_d_int->SetTitle(Form("Incident Deuteron %d GeV;True Inelastic Depth (mm);Estimate Depth (mm)", Energy_Name[k]));
        h2_d_int->GetXaxis()->CenterTitle();
        h2_d_int->GetYaxis()->CenterTitle();
        h2_d_int->Draw("colz");

        c3->cd(5);
        gStyle->SetOptFit(1111);
        h1_d_int->GetYaxis()->SetRangeUser(0,h1_d_int->GetMaximum()*1.2);
        h1_d_int->SetLineColor(kBlack);
        h1_d_int->SetTitle("Deuteron N_{inteaction};Depth (mm);Counts");
        h1_d_int->Draw("hist");
        h1_d_int->Fit(fitFunc4,"R");  
        fitFunc4->Draw("same");
        h1_d_inl->SetLineColor(kBlue);
        h1_d_inl->Draw("histsame");

        auto lg2 = new TLegend(0.3,0.7,0.6,0.88);
        lg2->AddEntry(h1_d_int,"Estimate Depth","l");
        lg2->AddEntry(h1_d_inl,"Inelastic Depth","l");
        lg2->AddEntry(fitFunc4,"Fit on Estimated","l");
        lg2->Draw();

        c3->cd(6);
        gStyle->SetOptFit(1111);
        h1_d_sur->GetYaxis()->SetRangeUser(0,h1_d_sur->GetMaximum()*1.2);
        h1_d_sur->SetLineColor(kBlack);
        h1_d_sur->SetTitle("Deuteron N_{Survive};Depth (mm);Counts");
        h1_d_sur->Draw("hist");
        h1_d_sur->Fit(fitFunc2,"R");  
        fitFunc2->Draw("same");
        h1_d_lea->SetLineColor(kBlue);
        h1_d_lea->Draw("histsame");

        lg2->Draw();

        // Xine = Xmid - 4.879 * Slope
        auto c4 = new TCanvas("c4","c4",2500,1500);
        c4->cd();
        c4->Clear();
        c4->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j <= 14; j++) // layer
        {
            c4->cd(j + 1);
            h1_p_Lay[j]->SetLineColor(kRed);   h1_p_Lay[j]->SetMarkerColor(kRed);  h1_p_Lay[j]->SetLineWidth(2);   h1_p_Lay[j]->Sumw2();
            h1_d_Lay[j]->SetLineColor(kBlue);  h1_d_Lay[j]->SetMarkerColor(kBlue); h1_d_Lay[j]->SetLineWidth(2);   h1_d_Lay[j]->Sumw2();
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
        gre_p_bias->GetXaxis()->SetLimits(0,14);
        gre_p_bias->SetTitle(";BGO Layer;Bias of (Xine-#hat{Xest})/Xine");
        gre_p_bias->Draw("AP");
        gre_d_bias->Draw("PSAME");

        c5->cd(2);
        gPad->SetGrid(0,1);
        gre_p_reso->GetXaxis()->SetLimits(0,14);
        gre_p_reso->SetTitle(";BGO Layer;Resolution of (Xine-#hat{Xest})/Xine");
        gre_p_reso->Draw("AP");
        gre_d_reso->Draw("PSAME");

        auto write_file = new TFile(Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_%dGeV.root",Energy_Name[k]), "RECREATE");
        write_file->cd();
        h2_p_int->Write();
        h2_d_int->Write();

        h1_p_int->Write();
        h1_d_int->Write();
        h1_p_inl->Write();
        h1_d_inl->Write();

        h1_p_sur->Write();
        h1_d_sur->Write();
        h1_p_lea->Write();
        h1_d_lea->Write();
        c3->Write();
        c4->Write();
        c5->Write();
        write_file->Close();
    }

}