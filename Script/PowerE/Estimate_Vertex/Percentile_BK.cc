#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti.hh"

void Percentile()
{
    int p_FH_Lay; int p_FH_Type; double p_Total_E;      int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_FH_Lay; int d_FH_Type; double d_Total_E;      int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
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
    proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
    proton_tree->SetBranchAddress("weight"          ,&p_weight);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);


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
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);
    deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);

    double Energy[15]={0};
    double Energy_LL[15]={0};      
    double Energy_UL[15]={0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[15][14];
    TH1D *h1_d[15][14];

    TH1D *h1_p_inter[15];
    TH1D *h1_d_inter[15];


    for(int i =0 ; i<15 ; i++)  // Deposit Energy Bin
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;

        h1_p_inter[i] =new TH1D(Form("h1_p_inter[%d]",i),Form("h1_p_inter[%d]",i),100,-10,0);  
        h1_d_inter[i] =new TH1D(Form("h1_d_inter[%d]",i),Form("h1_d_inter[%d]",i),100,-10,0);  
        for( int j= 0; j<14 ;j++)
        {
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),100,-10,0);  
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),100,-10,0);  
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
        }
    }
  
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    {        
        proton_tree->GetEntry(entry);
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        if(p_FI_Dep < 0) continue;
        if (p_Nhits < 10 ) continue;
        double sum_p = 0;
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double bar_Accumu_info[14] = {0};
        double bar_Accumu_error[14] = {0};
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引
        int layer_start = 4;
        const double RMS_threshold = 15.0;  // 自定义阈值，越小越“直”，你可以调整
        bool bar_info_assigned = false;     // 标志变量，判断是否已赋值
        for (int k = layer_start; k <= 12; k ++) 
        {  // 每次两层作为一个窗口
            if((*p_RMSVec)[k]<=RMS_threshold && (*p_RMSVec)[k+1]<=RMS_threshold)
            {
                int max_index1 = FindMaxMiddleIndex(p_EnergyVec, k);
                int bar1 = max_index1 % 22;
                int max_index2 = FindMaxMiddleIndex(p_EnergyVec, k+1);
                int bar2 = max_index2 % 22;
                if (k % 2 == 0) {
                    bar_info[0] = bar1; // odd
                    bar_info[1] = bar2; // even
                } else {
                    bar_info[0] = bar2; // odd
                    bar_info[1] = bar1; // even
                }
                // cout << "Directly determined bar_info: " 
                // << "bar1 = " << bar_info[0] << ", bar2 = " << bar_info[1] << endl;
                bar_info_assigned = true;
                break;  // 一旦赋值，跳出循环
            }
        }
        
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

        for(int layer = 0 ; layer<14 ; layer++)
        {   
            int center_bar = (layer % 2 == 0) ? bar_info[0] : bar_info[1];  // select center bar
            for (int k = center_bar - 1; k <= center_bar + 1; k++)
            {
                bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
            }
        }
        bar_Change_info[0] = log10(bar_Energy_info[0] / 0.023);
        bar_Accumu_info[0] = bar_Energy_info[0];
        bar_Accumu_error[0] = 0.3 * bar_Accumu_info[0];
        for(int layer = 1 ; layer<14 ; layer++)
        {
            if( bar_Energy_info[layer-1] == 0 || bar_Energy_info[layer] == 0) { bar_Change_info[layer-1] = -5 ;  }//  cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
            else {bar_Change_info[layer] = log10(bar_Energy_info[layer]/bar_Energy_info[layer-1]); }// cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
            bar_Accumu_info[layer]  += bar_Accumu_info[layer-1] + bar_Energy_info[layer];
            bar_Accumu_error[layer] = 0.3 * bar_Accumu_info[layer];
        }

        FindMaxPositiveSegment(bar_Change_info,14,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(bar_Change_info,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);
        rate_max_min = MaxMinRatio(bar_Energy_info,14);
        seg_sum_to_peak = AccumIncreaseToPeak(bar_Change_info,seg_start_idx,seg_peak_idx);
        seg_len_to_peak = seg_peak_idx - seg_start_idx;

        double E_L0 = bar_Energy_info[0];
        double maxE = FindMaxValue(bar_Energy_info, 14);
        double Amax = bar_Accumu_info[13];
        PrepareSigmoidData(bar_Accumu_info,bar_Accumu_error);
        TMinuit* myMinuit = nullptr;
        double reducedChi2 = RunSigmoidFit(entry, E_L0, Amax, seg_peak_idx, maxE, myMinuit);

        double Ymin, Ymin_err, Ymax, Ymax_err;
        double Slope, Slope_err, Xmid, Xmid_err;
        double E0, E0_err;
        myMinuit->GetParameter(0, Ymin, Ymin_err);
        myMinuit->GetParameter(1, Ymax, Ymax_err);
        myMinuit->GetParameter(2, Xmid, Xmid_err);
        myMinuit->GetParameter(3, Slope, Slope_err);
        myMinuit->GetParameter(4, E0, E0_err);
        double percentile = Mod_Sigmoid_Percentile(p_FI_Dep/25.5,Xmid,Slope);

        h1_p[p_energy_index][p_FI_Lay]->Fill(log10(percentile));
        h1_p_inter[p_energy_index]->Fill(log10(percentile));
        delete myMinuit;
    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry);
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 14) continue;
        if(d_FI_Dep < 0) continue;
        if (d_Nhits < 10 ) continue;
        double sum_d = 0;
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double bar_Accumu_info[14] = {0};
        double bar_Accumu_error[14] = {0};
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引
        int layer_start = 4;
        const double RMS_threshold = 15.0;  // 自定义阈值，越小越“直”，你可以调整
        bool bar_info_assigned = false;     // 标志变量，判断是否已赋值
        for (int k = layer_start; k <= 12; k ++) 
        {  // 每次两层作为一个窗口
            if((*d_RMSVec)[k]<=RMS_threshold && (*d_RMSVec)[k+1]<=RMS_threshold)
            {
                int max_index1 = FindMaxMiddleIndex(d_EnergyVec, k);
                int bar1 = max_index1 % 22;
                int max_index2 = FindMaxMiddleIndex(d_EnergyVec, k+1);
                int bar2 = max_index2 % 22;
                if (k % 2 == 0) {
                    bar_info[0] = bar1; // odd
                    bar_info[1] = bar2; // even
                } else {
                    bar_info[0] = bar2; // odd
                    bar_info[1] = bar1; // even
                }
                // cout << "Directly determined bar_info: " 
                // << "bar1 = " << bar_info[0] << ", bar2 = " << bar_info[1] << endl;
                bar_info_assigned = true;
                break;  // 一旦赋值，跳出循环
            }
        }
        
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

        for(int layer = 0 ; layer<14 ; layer++)
        {   
            int center_bar = (layer % 2 == 0) ? bar_info[0] : bar_info[1];  // select center bar
            for (int k = center_bar - 1; k <= center_bar + 1; k++)
            {
                bar_Energy_info[layer] += (*d_EnergyVec)[layer * 22 + k];
            }
        }
        bar_Change_info[0] = log10(bar_Energy_info[0] / 0.023);
        bar_Accumu_info[0] = bar_Energy_info[0];
        bar_Accumu_error[0] = 0.3 * bar_Accumu_info[0];
        for(int layer = 1 ; layer<14 ; layer++)
        {
            if( bar_Energy_info[layer-1] == 0 || bar_Energy_info[layer] == 0) { bar_Change_info[layer-1] = -5 ;  }
            else {bar_Change_info[layer] = log10(bar_Energy_info[layer]/bar_Energy_info[layer-1]); }
            bar_Accumu_info[layer]  += bar_Accumu_info[layer-1] + bar_Energy_info[layer];
            bar_Accumu_error[layer] = 0.3 * bar_Accumu_info[layer];
        }

        FindMaxPositiveSegment(bar_Change_info,14,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(bar_Change_info,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);
        rate_max_min = MaxMinRatio(bar_Energy_info,14);
        seg_sum_to_peak = AccumIncreaseToPeak(bar_Change_info,seg_start_idx,seg_peak_idx);
        seg_len_to_peak = seg_peak_idx - seg_start_idx;

        double E_L0 = bar_Energy_info[0];
        double maxE = FindMaxValue(bar_Energy_info, 14);
        double Amax = bar_Accumu_info[13];
        PrepareSigmoidData(bar_Accumu_info,bar_Accumu_error);
        TMinuit* myMinuit = nullptr;
        double reducedChi2 = RunSigmoidFit(entry, E_L0, Amax, seg_peak_idx, maxE, myMinuit);

        double Ymin, Ymin_err, Ymax, Ymax_err;
        double Slope, Slope_err, Xmid, Xmid_err;
        double E0, E0_err;
        myMinuit->GetParameter(0, Ymin, Ymin_err);
        myMinuit->GetParameter(1, Ymax, Ymax_err);
        myMinuit->GetParameter(2, Xmid, Xmid_err);
        myMinuit->GetParameter(3, Slope, Slope_err);
        myMinuit->GetParameter(4, E0, E0_err);
        double percentile = Mod_Sigmoid_Percentile(d_FI_Dep/25.5,Xmid,Slope);

        h1_d[d_energy_index][d_FI_Lay]->Fill(log10(percentile));
        h1_d_inter[d_energy_index]->Fill(log10(percentile)) ;
        delete myMinuit;
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

            h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV]Interaction happened in L%d;log10(Inelastic Percentile);Normalized Count", pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j ));
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
        h1_p_inter[i]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] Stack Multi Layer;log10(Inelastic Percentile);Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        h1_p_inter[i]->Draw("hist");
        h1_d_inter[i]->Draw("histsame");

    }
}