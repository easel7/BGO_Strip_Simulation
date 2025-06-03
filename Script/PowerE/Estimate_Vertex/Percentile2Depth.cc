#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 (x_ine - x_mid)

void Percentile2Depth()
{
    int Energy_Name[28]={0};
    int Energy_Name_Err[28]={0};
    auto gre_p = new TGraphErrors();
    auto gre_d = new TGraphErrors();
    if (gDirectory->FindObject("hist_p")) delete gDirectory->FindObject("hist_p");
    if (gDirectory->FindObject("hist_d")) delete gDirectory->FindObject("hist_d");
    auto hist_p = new TH1D("hist_p","hist_p",28,0,28);
    auto hist_d = new TH1D("hist_d","hist_d",28,0,28);

    // for (int k =0; k < 28; k++)
    for (int k =18; k < 19; k++)
    {
        if (k < 10)      {Energy_Name[k] = (k + 1) * 10;        }       // 10 ~ 100
        else if (k < 19) {Energy_Name[k] = (k - 9 + 1) * 100;   }        // 200 ~ 1000
        else             {Energy_Name[k] = (k - 18 + 1) * 1000; }        // 2000 ~ 10000
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
        // proton_tree->SetBranchAddress("weight"          ,&p_weight);
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

        auto h1_p_int = new TH1D("h1_p_int","h1_p_int",120,-10,5);  
        auto h1_d_int = new TH1D("h1_d_int","h1_d_int",120,-10,5);  

        auto h1_p_int_Sel = new TH1D("h1_p_int_Sel","h1_p_int_Sel",120,-10,5);  // Selected Ranged [3 - 12]
        auto h1_d_int_Sel = new TH1D("h1_d_int_Sel","h1_d_int_Sel",120,-10,5);  // Selected Ranged [3 - 12]


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
                h1_p_Lay[i] =new TH1D(Form("h1_p_Lay[%d]",i),Form("h1_p_Lay[%d]",i),120,-10,5);  
                h1_d_Lay[i] =new TH1D(Form("h1_d_Lay[%d]",i),Form("h1_d_Lay[%d]",i),120,-10,5);  
            }
        }
    
        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
        {        
            proton_tree->GetEntry(entry);
            if (entry%1000==0) cout << " Proton : " << entry << endl;
            int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
            // if(p_energy_index < 0 || p_energy_index > 14) continue;
            if(p_FI_Dep < 0) continue;
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

            // h1_p[p_energy_index][p_FI_Lay]->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
            // h1_p_inter[p_energy_index]->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3)) ;
            if(p_FI_Lay>= 0 && p_FI_Lay <= 4 ) {h1_p_int_Sel->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));}
            h1_p_int->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
            h1_p_Lay[p_FI_Lay]->Fill((p_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
        }

        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
        {
            deuteron_tree->GetEntry(entry);
            if (entry%1000==0) cout << " Deuteron : " << entry << endl;
            int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
            // if(d_energy_index < 0 || d_energy_index > 14) continue;
            if(d_FI_Dep < 0) continue;
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

            // h1_d[d_energy_index][d_FI_Lay]->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
            // h1_d_inter[d_energy_index]->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3)) ;
            if(d_FI_Lay>= 0 && d_FI_Lay <= 4 ) {h1_d_int_Sel->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));}
            h1_d_int->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
            h1_d_Lay[d_FI_Lay]->Fill((d_FI_Dep/25.5 - sigmoid->GetParameter(2) )/sigmoid->GetParameter(3));
        }

        for (int i = 9; i < 10; i++) // Deposit Energy Bin
        {
            auto c1 = new TCanvas("c1","c1",2500,1500);
            c1->Clear();
            c1->Divide(5,3);
            gStyle->SetOptStat(0);

            double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};     
            double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0};  
            double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0};  
            for (int j = 0; j < 14; j++) // layer
            {
                h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
                h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
                
                c1->cd(j + 1);

                // h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral());
                // h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral());
                // h1_p[i][j]->GetYaxis()->SetRangeUser(0,h1_p[i][j]->GetMaximum()*1.2);

                h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV]Interaction happened in L%d;(Xine-Xmid)/Slope;Normalized Count", pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j ));
                h1_p[i][j]->Draw("hist");
                h1_d[i][j]->Draw("histsame");

            }
            c1->cd(15);
            TLatex *tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));tex->SetNDC();tex->Draw(); 
            auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
            legend1->AddEntry(h1_p[i][0], "Proton", "l");
            legend1->AddEntry(h1_d[i][0], "Deuteron", "l");     
            legend1->Draw();       
            // c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Longti_PowerE/PDF/Longti_EnergyVec_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

            auto c2 = new TCanvas("c2","c2",1000,1000);
            c2->cd();
            h1_p_inter[i]->SetLineColor(kRed);   h1_p_inter[i]->SetMarkerColor(kRed);  h1_p_inter[i]->SetLineWidth(2);   h1_p_inter[i]->Sumw2();
            h1_d_inter[i]->SetLineColor(kBlue);  h1_d_inter[i]->SetMarkerColor(kBlue); h1_d_inter[i]->SetLineWidth(2);   h1_d_inter[i]->Sumw2();
            // h1_p_inter[i]->Scale(1.0/h1_p_inter[i]->Integral()); 
            // h1_d_inter[i]->Scale(1.0/h1_d_inter[i]->Integral()); 
            // h1_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
            h1_p_inter[i]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] Stack Multi Layer;(Xine-Xmid)/Slope;Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
            h1_p_inter[i]->Draw("hist");
            h1_d_inter[i]->Draw("histsame");

        }

        auto c3 = new TCanvas("c3","c3",1000,1000);
        c3->cd();
        gStyle->SetOptFit(1111);
        auto fitFunc1 = new TF1("fitFunc1", "gaus", -10, 5); // 注意替换范围
        auto fitFunc2 = new TF1("fitFunc2", "gaus", -10, 5); // 注意替换范围

        fitFunc1->SetLineColor(kRed);
        fitFunc1->SetLineWidth(2);
        fitFunc2->SetLineColor(kBlue);
        fitFunc2->SetLineWidth(2);

        h1_p_int_Sel->SetLineColor(kRed);   h1_p_int_Sel->SetMarkerColor(kRed);  h1_p_int_Sel->SetLineWidth(2);   h1_p_int_Sel->Sumw2();
        h1_d_int_Sel->SetLineColor(kBlue);  h1_d_int_Sel->SetMarkerColor(kBlue); h1_d_int_Sel->SetLineWidth(2);   h1_d_int_Sel->Sumw2();
        h1_p_int_Sel->Scale(1.0/h1_p_int_Sel->Integral()); 
        h1_d_int_Sel->Scale(1.0/h1_d_int_Sel->Integral()); 
        h1_p_int_Sel->GetYaxis()->SetRangeUser(0,h1_p_int_Sel->GetMaximum()*1.2);
        h1_p_int_Sel->SetTitle(" Stack Multi Layer;(Xine-Xmid)/Slope;Normalized Count");

        h1_p_int_Sel->Fit(fitFunc1,"R");  
        h1_d_int_Sel->Fit(fitFunc2,"R");  
        h1_p_int_Sel->Draw("hist");
        h1_d_int_Sel->Draw("histsame");
        fitFunc1->Draw("same");
        fitFunc2->Draw("same");
        gre_p->AddPointError(Energy_Name[k],fitFunc1->GetParameter(1),Energy_Name_Err[k],fitFunc1->GetParError(1));
        gre_d->AddPointError(Energy_Name[k],fitFunc2->GetParameter(1),Energy_Name_Err[k],fitFunc2->GetParError(1));
        hist_p->SetBinContent(k+1, fitFunc1->GetParameter(1)); hist_p->SetBinError(k+1, fitFunc1->GetParError(1));
        hist_d->SetBinContent(k+1, fitFunc2->GetParameter(1)); hist_d->SetBinError(k+1, fitFunc2->GetParError(1));
        hist_p->SetTitle("Proton;Energy Bin;(Xine-Xmid)/Slope");
        hist_d->SetTitle("Deuteron;Energy Bin;(Xine-Xmid)/Slope");

        auto c4 = new TCanvas("c4","c4",1000,1000);
        c4->cd();
        c4->Clear();
        c4->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int i = 0; i < 14; i++) // layer
        {
            c4->cd(i + 1);
            h1_p_Lay[i]->SetLineColor(kRed);   h1_p_Lay[i]->SetMarkerColor(kRed);  h1_p_Lay[i]->SetLineWidth(2);   h1_p_Lay[i]->Sumw2();
            h1_d_Lay[i]->SetLineColor(kBlue);  h1_d_Lay[i]->SetMarkerColor(kBlue); h1_d_Lay[i]->SetLineWidth(2);   h1_d_Lay[i]->Sumw2();
            h1_p_Lay[i]->Scale(1.0/h1_p_Lay[i]->Integral()); 
            h1_d_Lay[i]->Scale(1.0/h1_d_Lay[i]->Integral()); 
            h1_p_Lay[i]->GetYaxis()->SetRangeUser(0,h1_p_Lay[i]->GetMaximum()*1.2);
            h1_p_Lay[i]->SetTitle(" Stack Multi Layer;(Xine-Xmid)/Slope;Normalized Count");
            h1_p_Lay[i]->Draw("hist");
            h1_d_Lay[i]->Draw("histsame");
        }

    }

    auto write_file = new TFile("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_Fitting.root", "RECREATE");
    write_file->cd();

    auto c0 = new TCanvas("c0","c0",1600,1200);
    c0->cd();
    gPad->SetLogx();
    gre_p->SetTitle(";Energy(GeV);(Xine-Xmid)/Slope");
    gre_p->SetLineColor(kRed);
    gre_p->SetLineWidth(2);
    gre_p->SetMarkerColor(kRed);
    gre_p->SetMarkerStyle(20);

    gre_d->SetTitle(";Energy(GeV);(Xine-Xmid)/Slope");
    gre_d->SetLineColor(kBlue);
    gre_d->SetLineWidth(2);
    gre_d->SetMarkerColor(kBlue);
    gre_d->SetMarkerStyle(21);

    gre_p->Draw("ALP");
    gre_d->Draw("LPSAME");

    auto lg0 = new TLegend(0.7,0.7,0.88,0.88);
    lg0->AddEntry(gre_p,"Proton","lep");
    lg0->AddEntry(gre_d,"Deuteron","lep");
    lg0->Draw();
    hist_p->Write();
    hist_d->Write();
    gre_p->Write();
    gre_d->Write();
    c0->Write();
    write_file->Close();
}