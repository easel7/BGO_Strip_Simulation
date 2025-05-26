#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 \hat{x_ine} =  - 4.789 * sigma + x_mid

void Percentile2DepthEst()
{
    int Energy_Name[29]={0};                 
    for (int i =18; i < 29; i++)
    {
        if (i < 10)      Energy_Name[i] = (i + 1) * 10;               // 10 ~ 100
        else if (i < 19) Energy_Name[i] = (i - 9 + 1) * 100;           // 200 ~ 1000
        else             Energy_Name[i] = (i - 18 + 1) * 1000;         // 2000 ~ 10000

        cout << "Energy = " << Energy_Name[i] << " GeV !" << endl;
        int p_FH_Lay; int p_FH_Type; double p_Total_E;      int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
        int d_FH_Lay; int d_FH_Type; double d_Total_E;      int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
        int p_FI_Lay;    double p_FI_Dep;    int p_particle;
        int d_FI_Lay;    double d_FI_Dep;    int d_particle;

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",Energy_Name[i]));
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
        // proton_tree->SetBranchAddress("weight"          ,&p_weight);
        proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",Energy_Name[i]));
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
        // deuteron_tree->SetBranchAddress("weight"           ,&d_weight);
        deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);

        double Energy[15]={0};
        double Energy_LL[15]={0};      
        double Energy_UL[15]={0};

        double Layer[14]={0};
        double Layer_Err[14]={0};

        // Depsit and Layer
        TH1D *h1_p[15][14];     TH1D *h1_p_inter[15];  TH1D *h1_p_Lay[14];
        TH1D *h1_d[15][14];     TH1D *h1_d_inter[15];  TH1D *h1_d_Lay[14];

        auto h2_p_int = new TH2D("h2_p_int","h2_p_int",200,-40,360,200,-40,360);  
        auto h2_d_int = new TH2D("h2_d_int","h2_d_int",200,-40,360,200,-40,360);  

        auto h1_p_int = new TH1D("h1_p_int","h1_p_int",18,0,360);  
        auto h1_d_int = new TH1D("h1_d_int","h1_d_int",18,0,360);  

        auto h1_p_inl = new TH1D("h1_p_inl","h1_p_inl",18,0,360);  
        auto h1_d_inl = new TH1D("h1_d_inl","h1_d_inl",18,0,360);  

        auto h1_p_sur = new TH1D("h1_p_sur","h1_p_sur",18,0,360);  
        auto h1_d_sur = new TH1D("h1_d_sur","h1_d_sur",18,0,360);  

        auto h1_p_lea = new TH1D("h1_p_lea","h1_p_lea",18,0,360);  
        auto h1_d_lea = new TH1D("h1_d_lea","h1_d_lea",18,0,360);  

        for(int i =0 ; i<15 ; i++)  // Deposit Energy Bin
        {
            Energy[i]    = 1.1 + 0.2 * i;
            Energy_LL[i] = 1.0 + 0.2 * i;
            Energy_UL[i] = 1.2 + 0.2 * i;

            h1_p_inter[i] =new TH1D(Form("h1_p_inter[%d]",i),Form("h1_p_inter[%d]",i),120,-10,5);  
            h1_d_inter[i] =new TH1D(Form("h1_d_inter[%d]",i),Form("h1_d_inter[%d]",i),120,-10,5);  
            for( int j= 0; j<14 ;j++)
            {
                h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),120,-10,5);  
                h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),120,-10,5);  
                Layer[j] = 0.5 + j;
                Layer_Err[j] = 0.5;
            }
            if (i<14)
            {
                h1_p_Lay[i] =new TH1D(Form("h1_p_Lay[%d]",i),Form("h1_p_Lay[%d]",i),60,-3,3);  
                h1_d_Lay[i] =new TH1D(Form("h1_d_Lay[%d]",i),Form("h1_d_Lay[%d]",i),60,-3,3);  
            }
        }
    
        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
        {        
            proton_tree->GetEntry(entry);
            int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
            // if(p_energy_index < 0 || p_energy_index > 14) continue;
            if(p_FI_Dep < 0) continue;
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

            // h1_p[p_energy_index][p_FI_Lay]->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
            // h1_p_inter[p_energy_index]->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3)) ;
            h2_p_int->Fill(p_FI_Dep, (-4.789 * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5);
            h1_p_Lay[p_FI_Lay]->Fill( (p_FI_Dep -  (-4.789 * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5) / p_FI_Dep );
            h1_p_int->Fill((-4.789 * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5);
            h1_p_inl->Fill(p_FI_Dep);
        }

        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
        {
            deuteron_tree->GetEntry(entry);
            int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
            // if(d_energy_index < 0 || d_energy_index > 14) continue;
            if(d_FI_Dep < 0) continue;
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

            // h1_d[d_energy_index][d_FI_Lay]->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
            // h1_d_inter[d_energy_index]->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3)) ;
            h2_d_int->Fill(d_FI_Dep, (-4.789 * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5);
            h1_d_Lay[d_FI_Lay]->Fill( (d_FI_Dep -  (-4.789 * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5) / d_FI_Dep );
            h1_d_int->Fill((-4.789 * sigmoid->GetParameter(3) + sigmoid->GetParameter(2)) * 25.5);
            h1_d_inl->Fill(d_FI_Dep);
        }

        // for (int i = 9; i < 10; i++) // Deposit Energy Bin
        // {
        //     auto c1 = new TCanvas("c1","c1",2500,1500);
        //     c1->Clear();
        //     c1->Divide(5,3);
        //     gStyle->SetOptStat(0);

        //     double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};     
        //     double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0};  
        //     double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0};  
        //     for (int j = 0; j < 14; j++) // layer
        //     {
        //         h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
        //         h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
                
        //         c1->cd(j + 1);

        //         // h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral());
        //         // h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral());
        //         // h1_p[i][j]->GetYaxis()->SetRangeUser(0,h1_p[i][j]->GetMaximum()*1.2);

        //         h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV]Interaction happened in L%d;(Xine-Xmid)/Slope;Normalized Count", pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j ));
        //         h1_p[i][j]->Draw("hist");
        //         h1_d[i][j]->Draw("histsame");

        //     }
        //     c1->cd(15);
        //     TLatex *tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));tex->SetNDC();tex->Draw(); 
        //     auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        //     legend1->AddEntry(h1_p[i][0], "Proton", "l");
        //     legend1->AddEntry(h1_d[i][0], "Deuteron", "l");     
        //     legend1->Draw();       
        //     // c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Longti_PowerE/PDF/Longti_EnergyVec_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        //     auto c2 = new TCanvas("c2","c2",1000,1000);
        //     c2->cd();
        //     h1_p_inter[i]->SetLineColor(kRed);   h1_p_inter[i]->SetMarkerColor(kRed);  h1_p_inter[i]->SetLineWidth(2);   h1_p_inter[i]->Sumw2();
        //     h1_d_inter[i]->SetLineColor(kBlue);  h1_d_inter[i]->SetMarkerColor(kBlue); h1_d_inter[i]->SetLineWidth(2);   h1_d_inter[i]->Sumw2();
        //     // h1_p_inter[i]->Scale(1.0/h1_p_inter[i]->Integral()); 
        //     // h1_d_inter[i]->Scale(1.0/h1_d_inter[i]->Integral()); 
        //     // h1_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
        //     h1_p_inter[i]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] Stack Multi Layer;(Xine-Xmid)/Slope;Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        //     h1_p_inter[i]->Draw("hist");
        //     h1_d_inter[i]->Draw("histsame");
        // }

        h2_p_int->Sumw2();    h1_p_int->Sumw2();
        h2_d_int->Sumw2();    h1_d_int->Sumw2();
        TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-x/[1])", 60, 260); 
        TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 60, 260); 
        fitFunc1->SetParameters(1e4, 200); 
        fitFunc2->SetParameters(1e4, 170); 
        fitFunc1->SetLineColor(kRed); 
        fitFunc2->SetLineColor(kBlue);
        // fitFunc1->FixParameter(0, h2_p_int->Integral());
        // fitFunc2->FixParameter(0, h2_d_int->Integral());

        TF1 *fitFunc3 = new TF1("fitFunc3", "[0]/[1]*exp(-x/[1])", 60, 260); 
        TF1 *fitFunc4 = new TF1("fitFunc4", "[0]/[1]*exp(-x/[1])", 60, 260); 
        fitFunc3->SetParameters(1e3, 200); 
        fitFunc4->SetParameters(1e3, 170); 
        fitFunc3->SetLineColor(kRed); 
        fitFunc4->SetLineColor(kBlue);
        // fitFunc3->FixParameter(0, h2_p_int->Integral() / h2_p_int->GetBinWidth(1));
        // fitFunc4->FixParameter(0, h2_d_int->Integral() / h2_d_int->GetBinWidth(1));
        // fitFunc3->FixParameter(2, h2_p_int->GetBinWidth(1));
        // fitFunc4->FixParameter(2, h2_d_int->GetBinWidth(1));

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
        c3->cd();
        c3->Divide(3,2);
        c3->cd(1);
        h2_p_int->SetTitle("Proton Stack Multi Layer;True Inelastic Depth;Depth (mm)");
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
        lg1->AddEntry(h1_p_int,"Estimated","l");
        lg1->AddEntry(h1_p_inl,"Inelastic","l");
        lg1->Draw();

        c3->cd(3);
        gStyle->SetOptFit(1111);
        h1_p_sur->GetYaxis()->SetRangeUser(0,h1_p_sur->GetMaximum()*1.2);
        h1_p_sur->SetLineColor(kBlack);
        h1_p_sur->SetTitle("Proton N_{Survive};Depth (mm);Counts");
        h1_p_sur->Draw("hist");
        h1_p_lea->SetLineColor(kRed);
        h1_p_lea->Draw("histsame");
        h1_p_lea->Fit(fitFunc1,"R");  
        fitFunc1->Draw("same");
        lg1->Draw();


        c3->cd(4);
        h2_d_int->SetTitle("Deuteron Stack Multi Layer;True Inelastic Depth;Depth (mm);");
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
        lg2->AddEntry(h1_d_int,"Estimated","l");
        lg2->AddEntry(h1_d_inl,"Inelastic","l");
        lg2->Draw();

        c3->cd(6);
        gStyle->SetOptFit(1111);
        h1_d_sur->GetYaxis()->SetRangeUser(0,h1_d_sur->GetMaximum()*1.2);
        h1_d_sur->SetLineColor(kBlack);
        h1_d_sur->SetTitle("Deuteron N_{Survive};Depth (mm);Counts");
        h1_d_sur->Draw("hist");
        h1_d_lea->SetLineColor(kBlue);
        h1_d_lea->Draw("histsame");
        h1_d_lea->Fit(fitFunc2,"R");  
        fitFunc2->Draw("same");
        lg2->Draw();

    // Xine = Xmid - 4.879 * Slope
        auto c4 = new TCanvas("c4","c4",1000,1000);
        c4->cd();
        c4->Clear();
        c4->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j < 14; j++) // layer
        {
            c4->cd(j + 1);
            h1_p_Lay[j]->SetLineColor(kRed);   h1_p_Lay[j]->SetMarkerColor(kRed);  h1_p_Lay[j]->SetLineWidth(2);   h1_p_Lay[j]->Sumw2();
            h1_d_Lay[j]->SetLineColor(kBlue);  h1_d_Lay[j]->SetMarkerColor(kBlue); h1_d_Lay[j]->SetLineWidth(2);   h1_d_Lay[j]->Sumw2();
            h1_p_Lay[j]->Scale(1.0/h1_p_Lay[j]->Integral()); 
            h1_d_Lay[j]->Scale(1.0/h1_d_Lay[j]->Integral()); 
            h1_p_Lay[j]->GetYaxis()->SetRangeUser(0,h1_p_Lay[j]->GetMaximum()*1.2);
            h1_p_Lay[j]->SetTitle(" Stack Multi Layer;(Xine-#hat{Xest})/Xine;Normalized Count");
            h1_p_Lay[j]->Draw("hist");
            h1_d_Lay[j]->Draw("histsame");
        }


        auto write_file = new TFile(Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_%dGeV.root",Energy_Name[i]), "RECREATE");
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
        write_file->Close();
    }

}