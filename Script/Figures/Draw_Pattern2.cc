#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti.hh"

void Draw_Pattern2()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    double p_Total_E;
    int p_FH_Type;
    int p_FH_Lay;
    int p_Nhits;    
    int p_FH_Second;
    double p_FI_Dep;
    int p_FI_Lay;

    const char* string1;
    const char* string2 = "Proton_1000GeV";
    // const char* string2 = "Proton_10000GeV";
    // const char* string2 = "Deuteron_1000GeV";
    // const char* string2 = "Deuteron_10000GeV";
    auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/%s.root",string2));
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("BarEnergyVector", &p_EnergyVec);
    proton_tree->SetBranchAddress("RMS"            , &p_RMSVec);
    proton_tree->SetBranchAddress("First_Had_Depth", &p_FH_Dep);
    proton_tree->SetBranchAddress("First_Had_Layer", &p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type" , &p_FH_Type);    
    proton_tree->SetBranchAddress("First_Had_Second", &p_FH_Second);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);
    proton_tree->SetBranchAddress("Total_E"        , &p_Total_E);
    proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
    proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);
    int nbins = 30;
    double xmin = 1e-2;
    double xmax = 1e1;

    double logxmin = log10(xmin);
    double logxmax = log10(xmax);

    std::vector<double> bin_edges(nbins + 1);
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * (logxmax - logxmin) / nbins);
    }

    auto h_max_min0 = new TH1D("h_max_min0","h_max_min0",60,-1,5); // Until Inelastic
    auto h_max_min1 = new TH1D("h_max_min1","h_max_min1",60,-1,5); // Inelastic
    auto h_max_min2 = new TH1D("h_max_min2","h_max_min2",60,-1,5); // Elastic
    auto h_max_min3 = new TH1D("h_max_min3","h_max_min3",60,-1,5); // Pass

    auto h_change_0 = new TH1D("h_change_0","h_change_0",nbins, bin_edges.data()); // Until Inelastic
    auto h_change_1 = new TH1D("h_change_1","h_change_1",nbins, bin_edges.data()); // Inelastic
    auto h_change_2 = new TH1D("h_change_2","h_change_2",nbins, bin_edges.data()); // Elastic
    auto h_change_3 = new TH1D("h_change_3","h_change_3",nbins, bin_edges.data()); // Pass

    auto h_contin_0 = new TH1D("h_contin_0","h_contin_0",14,0,14); // Until Inelastic
    auto h_contin_1 = new TH1D("h_contin_1","h_contin_1",14,0,14); // Inelastic
    auto h_contin_2 = new TH1D("h_contin_2","h_contin_2",14,0,14); // Elastic
    auto h_contin_3 = new TH1D("h_contin_3","h_contin_3",14,0,14); // Pass

    auto g_sum_len0 = new TGraph();  int point_counter   = 0; int pc   = 0;
    auto g_sum_len1 = new TGraph();  int point_counter_i = 0; int pc_i = 0;
    auto g_sum_len2 = new TGraph();  int point_counter_e = 0; int pc_e = 0;
    auto g_sum_len3 = new TGraph();  int point_counter_p = 0; int pc_p = 0;

    auto h_peak_had0 = new TH2I("h_peak_had0","h_peak_had0",14,0,14,16,-2,14);
    auto h_peak_had1 = new TH2I("h_peak_had1","h_peak_had1",14,0,14,16,-2,14);
    auto h_peak_had2 = new TH2I("h_peak_had2","h_peak_had2",14,0,14,16,-2,14);
    auto h_peak_had3 = new TH2I("h_peak_had3","h_peak_had3",14,0,14,16,-2,14);

    auto h_start_had0 = new TH2I("h_start_had0","h_start_had0",14,0,14,16,-2,14);
    auto h_start_had1 = new TH2I("h_start_had1","h_start_had1",14,0,14,16,-2,14);
    auto h_start_had2 = new TH2I("h_start_had2","h_start_had2",14,0,14,16,-2,14);
    auto h_start_had3 = new TH2I("h_start_had3","h_start_had3",14,0,14,16,-2,14);

    auto h_len_had0 = new TH2I("h_len_had0","h_len_had0",14,0,14,16,-2,14);
    auto h_len_had1 = new TH2I("h_len_had1","h_len_had1",14,0,14,16,-2,14);
    auto h_len_had2 = new TH2I("h_len_had2","h_len_had2",14,0,14,16,-2,14);
    auto h_len_had3 = new TH2I("h_len_had3","h_len_had3",14,0,14,16,-2,14);

    auto h_Lay_MM1  = new TH2D("h_Lay_MM1","h_Lay_MM1",60,-1,5,16,-2,14);
    auto h_Lay_MM2  = new TH2D("h_Lay_MM2","h_Lay_MM2",60,-1,5,16,-2,14);
    auto h_Lay_MM3  = new TH2D("h_Lay_MM3","h_Lay_MM3",60,-1,5,16,-2,14);

    auto h_CR_MM1  = new TH2D("h_CR_MM1","h_CR_MM1",60,-1,5,60,-1,5);
    auto h_CR_MM2  = new TH2D("h_CR_MM2","h_CR_MM2",60,-1,5,60,-1,5);
    auto h_CR_MM3  = new TH2D("h_CR_MM3","h_CR_MM3",60,-1,5,60,-1,5);

    auto h_peak_val1 = new TH1D("h_peak_val1","h_peak_val1",50,0,5); // Inelastic
    auto h_peak_val2 = new TH1D("h_peak_val2","h_peak_val2",50,0,5); // Elastic
    auto h_peak_val3 = new TH1D("h_peak_val3","h_peak_val3",50,0,5); // Pass

    auto h_sp_rate1 = new TH1D("h_sp_rate1","h_sp_rate1",60,-1,5); // Inelastic
    auto h_sp_rate2 = new TH1D("h_sp_rate2","h_sp_rate2",60,-1,5); // Elastic
    auto h_sp_rate3 = new TH1D("h_sp_rate3","h_sp_rate3",60,-1,5); // Pass

    auto h_sp_bin1 = new TH1D("h_sp_bin1","h_sp_bin1",14,0,14); // Inelastic
    auto h_sp_bin2 = new TH1D("h_sp_bin2","h_sp_bin2",14,0,14); // Elastic
    auto h_sp_bin3 = new TH1D("h_sp_bin3","h_sp_bin3",14,0,14); // Pass

    auto h_int      = new TH1D("h_int","h_int",14,0,14);
    auto h_sur      = new TH1D("h_sur","h_sur",14,0,14);

    auto h_peak_Ine = new TH2I("h_peak_Ine","h_peak_Ine",14,0,14,14,0,14);
    auto h_Ine_Reso = new TH2I("h_Ine_Reso","h_Ine_Reso",14,0,14,14,0,14);

    // cout  << proton_tree->GetEntries() << endl;
    Long64_t entry  = 162;   
    // for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    // for (Long64_t entry = 0; entry < 100; entry++)
    {        
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
        
        proton_tree->GetEntry(entry);
        // if (p_Nhits < 10 ) continue;
        int layer_start = 4;
        const double RMS_threshold = 15.0;  // 自定义阈值，越小越“直”，你可以调整
        bool bar_info_assigned = false;  // 标志变量，判断是否已赋值
        // cout << " Layer " << layer_start << " , RMS = " <<  (*p_RMSVec)[layer_start] << " , Next Layer RMS = " << (*p_RMSVec)[layer_start+1] << endl;
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
            PrepareFitData(p_EnergyVec, layer_start, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            TMinuit minuit0(1);
            minuit0.SetFCN(FitAxisFunction);
            minuit0.SetPrintLevel(-1);
            minuit0.DefineParameter(0, "bar_odd", g_fit_bars[1], 0.01, 2, 19);
            minuit0.Migrad();
            double bar_odd, bar_odd_err;
            minuit0.GetParameter(0, bar_odd, bar_odd_err);
            // cout << " bar_odd = " << bar_odd << endl;
            bar_info[0] = std::round(bar_odd);

            PrepareFitData(p_EnergyVec, layer_start+1, 14, g_fit_bars, g_fit_energies, g_fit_total_energy);
            TMinuit minuit1(1);
            minuit1.SetFCN(FitAxisFunction);
            minuit1.SetPrintLevel(-1);
            minuit1.DefineParameter(0, "bar_even", g_fit_bars[1], 0.04, 2, 19);
            minuit1.Migrad();
            double bar_even, bar_even_err;
            minuit1.GetParameter(0, bar_even, bar_even_err);
            // cout << " bar_even = " << bar_even << endl;
            bar_info[1] = std::round(bar_even);
        }

        for(int layer = 0 ; layer<14 ; layer++)
        {   
            int center_bar = (layer % 2 == 0) ? bar_info[0] : bar_info[1];  // select center bar
            for (int k = center_bar - 1; k <= center_bar + 1; k++)
            {
                bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
            }
            if(layer==0)
            {
                bar_Change_info[0] = log10(bar_Energy_info[0] / 0.023);
            }
            else // (layer>0) 
            {
                if( bar_Energy_info[layer-1] == 0 || bar_Energy_info[layer] == 0) { bar_Change_info[layer-1] = -5 ;  }//  cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
                else {bar_Change_info[layer] = log10(bar_Energy_info[layer]/bar_Energy_info[layer-1]); }// cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
            }
            bar_Accumu_info[layer]  += bar_Accumu_info[layer-1] + bar_Energy_info[layer];
            bar_Accumu_error[layer] = 0.3 * bar_Accumu_info[layer];
        }

        PrepareSigmoidData(bar_Accumu_info,bar_Accumu_error);
        TMinuit minuit(4);
        minuit.SetFCN(SigmoidFCN);
        minuit.SetPrintLevel(-1); // 静默输出
        minuit.SetErrorDef(1.0);  // Δχ² = 1 规则
        minuit.DefineParameter(0, "Ymin", bar_Accumu_info[0], 1, 0, bar_Accumu_info[13]); // initVal, initErr, LowerL, UpperL
        minuit.DefineParameter(1, "Ymax", bar_Accumu_info[13], 1, 0, bar_Accumu_info[13]*1.5); 
        minuit.DefineParameter(2, "Xmid", seg_peak_idx, 0.1, 1, 14); // 拐点
        minuit.DefineParameter(3, "Slope", 1.0, 0.1, 0.1, 5); // 斜率
        minuit.FixParameter(0);
        minuit.FixParameter(1);
        minuit.Migrad();
        double Ymin, Ymin_err, Ymax, Ymax_err;
        double Slope, Slope_err, Xmid, Xmid_err;
        minuit.GetParameter(0, Ymin, Ymin_err);
        minuit.GetParameter(1, Ymax, Ymax_err);
        minuit.GetParameter(2, Xmid, Xmid_err);
        minuit.GetParameter(3, Slope, Slope_err);
        cout << "Ymin: " << Ymin << " ± " << Ymin_err    << endl;
        cout << "Ymax: " << Ymax << " ± " << Ymax_err    << endl;
        cout << "Xmid: " << Xmid << " ± " << Xmid_err    << endl;
        cout << "Slope: " << Slope << " ± " << Slope_err << endl;
        double percentile = Mod_Sigmoid_Percentile(p_FI_Lay,Xmid,Slope);
        cout << "perenctile: " << percentile << endl; 

        FindMaxPositiveSegment(bar_Change_info,14,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(bar_Change_info,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);
        rate_max_min = MaxMinRatio(bar_Energy_info,14);
        seg_sum_to_peak = AccumIncreaseToPeak(bar_Change_info,seg_start_idx,seg_peak_idx);
        seg_len_to_peak = seg_peak_idx - seg_start_idx;
        // cout << " begin to increase bin = " << seg_start_idx << endl;
        // cout << "Max Positive Bin Length = " << seg_len << endl;
        // cout << "Max Increase Rate = " <<  seg_peak_value << endl;
        // cout << "Max Increase Rate Bin = " << seg_peak_idx << endl;
        g_sum_len0->SetPoint(point_counter++,seg_sum,seg_len); 
        
        
        // if (rate_max_min >  1e2 && seg_len > 5 ) 
        {h_int->Fill(seg_peak_idx);}
        // if (rate_max_min < 100 && rate_max_min>60  && p_FH_Type == 1 && p_FH_Lay>11) { cout << entry << " , " <<  rate_max_min <<  endl; }
        if (rate_max_min >  1e2 )
        {}
        h_peak_Ine->Fill(seg_peak_idx,p_FI_Lay);
        h_Ine_Reso->Fill(p_FI_Lay,abs(seg_peak_idx-p_FI_Lay));

        
        if(p_FH_Type == 1)       {  string1 = "Inelastic"; 
            h_max_min1->Fill(log10(rate_max_min)); 
            h_change_1->Fill(seg_sum); 
            h_contin_1->Fill(seg_len); 
            g_sum_len1->SetPoint(point_counter_i++,seg_sum,seg_len+0.1); 
            h_peak_had1->Fill(seg_peak_idx,p_FH_Lay); 
            h_start_had1->Fill(seg_start_idx,p_FH_Lay); 
            h_len_had1->Fill(seg_len,p_FH_Lay); 
            h_Lay_MM1->Fill(log10(rate_max_min),p_FH_Lay);
            h_CR_MM1->Fill(log10(rate_max_min),seg_sum);
            h_sp_rate1->Fill(seg_sum_to_peak);
            h_sp_bin1->Fill(seg_len_to_peak);
            h_peak_val1->Fill(seg_peak_value);
        } 
        else if (p_FH_Type == 2) {  string1 = "Elastic";   
            h_max_min2->Fill(log10(rate_max_min)); 
            h_change_2->Fill(seg_sum); 
            h_contin_2->Fill(seg_len); 
            g_sum_len2->SetPoint(point_counter_e++,seg_sum,seg_len); 
            h_peak_had2->Fill(seg_peak_idx,p_FH_Lay);
            h_start_had2->Fill(seg_start_idx,p_FH_Lay); 
            h_len_had2->Fill(seg_len,p_FH_Lay); 

            h_Lay_MM2->Fill(log10(rate_max_min),p_FH_Lay);
            h_CR_MM2->Fill(log10(rate_max_min),seg_sum);

            h_sp_rate2->Fill(seg_sum_to_peak);
            h_sp_bin2->Fill(seg_len_to_peak);
            h_peak_val2->Fill(seg_peak_value);
            if(p_FI_Dep>0)
            {   
                h_max_min0->Fill(log10(rate_max_min)); 
                h_change_0->Fill(seg_sum); 
                h_contin_0->Fill(seg_len); 
            }
        }
        else                     {  string1 = "Pass";      
            h_max_min3->Fill(log10(rate_max_min)); 
            h_change_3->Fill(seg_sum); 
            h_contin_3->Fill(seg_len); 
            g_sum_len3->SetPoint(point_counter_p++,seg_sum,seg_len-0.1);
            h_peak_had3->Fill(seg_peak_idx,p_FH_Lay);
            h_start_had3->Fill(seg_start_idx,p_FH_Lay); 
            h_len_had3->Fill(seg_len,p_FH_Lay); 

            h_Lay_MM3->Fill(log10(rate_max_min),p_FH_Lay);
            h_CR_MM3->Fill(log10(rate_max_min),seg_sum);

            h_sp_rate3->Fill(seg_sum_to_peak);
            h_sp_bin3->Fill(seg_len_to_peak);
            h_peak_val3->Fill(seg_peak_value);
        }
    }
    
    auto c0    = new TCanvas("c0","c0",1200,1200);
    c0->Divide(2,2);
    c0->cd(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_max_min1->GetXaxis()->SetLimits(-1, 5);        // X 轴范围
    h_max_min1->SetTitle(";log10(Emax/Emin);Counts");
    h_max_min1->SetLineColor(kRed);
    h_max_min2->SetLineColor(kBlue);
    h_max_min3->SetLineColor(kOrange-3);
    h_max_min0->SetLineColor(kBlack);
    h_max_min1->Draw("hist");
    h_max_min2->Draw("histsame");
    h_max_min3->Draw("histsame");
    h_max_min0->Draw("histsame");


    auto legend0 = new TLegend(0.12,0.7,0.32,0.88);
    legend0->AddEntry(h_max_min1, "FH Inlastic","l");
    legend0->AddEntry(h_max_min2, "FH Elastic","l");
    legend0->AddEntry(h_max_min3, "Pass through","l");
    legend0->AddEntry(h_max_min0, "Until Inelastic","l");
    legend0->Draw();

    c0->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    // h_change_1->GetXaxis()->SetLimits(1e-2, 3e3);
    h_change_1->SetTitle(";#sum log10(Change Rate);Counts");
    h_change_1->SetLineColor(kRed);
    h_change_2->SetLineColor(kBlue);
    h_change_3->SetLineColor(kOrange-3);
    h_change_0->SetLineColor(kBlack);
    h_change_1->Draw("hist");
    h_change_2->Draw("histsame");
    h_change_3->Draw("histsame");
    h_change_0->Draw("histsame");

    legend0->Draw();

    c0->cd(3);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    // h_contin_1->GetYaxis()->SetRangeUser(1e0, 3e3); 
    h_contin_1->SetTitle(";# Continues Positive Bins;Counts");
    h_contin_1->SetLineColor(kRed);
    h_contin_2->SetLineColor(kBlue);
    h_contin_3->SetLineColor(kOrange-3);
    h_contin_0->SetLineColor(kBlack);
    h_contin_1->Draw("hist");
    h_contin_2->Draw("histsame");
    h_contin_3->Draw("histsame");
    h_contin_0->Draw("histsame");
    legend0->Draw();

    c0->cd(4);
    gPad->SetLogx();
    g_sum_len0->SetTitle(";#sum log10(Change Rate);# Continues Positive Bins");
    g_sum_len0->GetXaxis()->SetLimits(1e-2, 100);        // X 轴范围
    g_sum_len0->GetYaxis()->SetRangeUser(0, 14);    // Y 轴范围
    g_sum_len0->SetMarkerStyle(20);  
    g_sum_len0->SetMarkerColorAlpha(kBlack, 0);  
    g_sum_len0->SetMarkerSize(0);

    g_sum_len1->SetMarkerStyle(20);  
    g_sum_len1->SetMarkerColorAlpha(kRed, 0.1); 
    g_sum_len1->SetMarkerSize(0.8);

    g_sum_len2->SetMarkerStyle(21);  
    g_sum_len2->SetMarkerColorAlpha(kBlue, 0.1);
    g_sum_len2->SetMarkerSize(0.8);

    g_sum_len3->SetMarkerStyle(22);  
    g_sum_len3->SetMarkerColorAlpha(kOrange-3, 0.1);
    g_sum_len3->SetMarkerSize(0.8);

    g_sum_len0->Draw("AP"); 
    g_sum_len1->Draw("PSAME");
    g_sum_len2->Draw("PSAME");
    g_sum_len3->Draw("PSAME");

    auto legend1 = new TLegend(0.7,0.7,0.88,0.88);
    legend1->AddEntry(g_sum_len1, "FH Inelastic","p");
    legend1->AddEntry(g_sum_len2, "FH Elastic","p");
    legend1->AddEntry(g_sum_len3, "Pass through","p");
    legend1->Draw();
    std::cout << "Number of valid points: " << g_sum_len0->GetN() << std::endl;
    c0->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/Figures/%s_FIG.pdf",string2));
    
    //////////////////////////////////////////////

    double minVal_1 = 1 ;       double minVal_2 = 1 ;
    double maxVal_1 = 1e4;      double maxVal_2 = 2e2; 
    h_peak_had1->SetMinimum(minVal_1);h_start_had1->SetMinimum(minVal_1); h_Lay_MM1->SetMinimum(minVal_2);
    h_peak_had1->SetMaximum(maxVal_1);h_start_had1->SetMaximum(maxVal_1); h_Lay_MM1->SetMaximum(maxVal_2);
    h_peak_had2->SetMinimum(minVal_1);h_start_had2->SetMinimum(minVal_1); h_Lay_MM2->SetMinimum(minVal_2);
    h_peak_had2->SetMaximum(maxVal_1);h_start_had2->SetMaximum(maxVal_1); h_Lay_MM2->SetMaximum(maxVal_2);
    h_peak_had3->SetMinimum(minVal_1);h_start_had3->SetMinimum(minVal_1); h_Lay_MM3->SetMinimum(minVal_2);
    h_peak_had3->SetMaximum(maxVal_1);h_start_had3->SetMaximum(maxVal_1); h_Lay_MM3->SetMaximum(maxVal_2);

    auto c1 = new TCanvas("c1","c1",1800,1200);
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogz();
    h_peak_had1->SetTitle("Inelastic;Bin of Maximum Change Ratio; First Hadronic Layer");
    h_peak_had1->Draw("colz");
    double cov1 = h_peak_had1->GetCorrelationFactor();
    cout << " cov1 = " << cov1 << endl;

    c1->cd(2);
    gPad->SetLogz();
    h_peak_had2->SetTitle("Elastic;Bin of Maximum Change Ratio; First Hadronic Layer");
    h_peak_had2->Draw("colz");


    c1->cd(3);
    gPad->SetLogz();
    h_peak_had3->SetTitle("Pass through;Bin of Maximum Change Ratio; First Hadronic Layer");
    h_peak_had3->Draw("colz");

    c1->cd(4);
    gPad->SetLogz();
    h_peak_Ine->SetMinimum(minVal_1);
    h_peak_Ine->SetMaximum(maxVal_1);
    h_peak_Ine->SetTitle(";Bin of Maximum Change Ratio; First Inelastic Layer");
    h_peak_Ine->Draw("colz");
    double cov0 = h_peak_Ine->GetCorrelationFactor();
    cout << " cov0 = " << cov0 << endl;

    c1->cd(5);
    gPad->SetLogz();
    h_Ine_Reso->SetMinimum(minVal_1);
    h_Ine_Reso->SetMaximum(maxVal_1);
    h_Ine_Reso->SetTitle(";First Inelastic Layer; Bin of Maximum Change Ratio - First Inelastic Layer");
    h_Ine_Reso->Draw("colz");
    TH1D *h1_p_reso[14];
    for (int j = 0; j < 14; j++) // layer
    {
        h1_p_reso[j] = h_Ine_Reso->ProjectionY(Form("h1_p_reso[%d]",j), j+1, j+1, "");
        cout <<h1_p_reso[j]->Integral()<< endl;
        h1_p_reso[j]->Scale(1/h1_p_reso[j]->Integral());
        h1_p_reso[j]->SetTitle(Form("Inelastic in Layer %d;|Bin of Maximum Change Ratio - First Inelastic Layer|; Normalized Count",j));
    }
    ///////////////////////////////

    double Proton_Acc0[14]={0};   
    double Proton_Acc1[14]={0};   
    double Proton_Acc2[14]={0};   

    double Layer[14]={0};   


    auto c1_2 = new TCanvas("c1_2","c1_2",2000,1200);
    c1_2->Divide(5,3);
    for (int j = 0; j < 14; j++) // layer
    {
        c1_2->cd(j+1);
        gPad->SetLogy();
        h1_p_reso[j]->Draw("hist");

        Proton_Acc0[j]   = h1_p_reso[j]->Integral(1,1); 
        Proton_Acc1[j]   = h1_p_reso[j]->Integral(2,2); 
        Proton_Acc2[j]   = h1_p_reso[j]->Integral(3,14); 

        cout << Proton_Acc0[j] << " , " << Proton_Acc1[j] << " , " << Proton_Acc2[j] << endl;
        Layer[j] = 0.5 + j;
    }
    c1_2->cd(15);
    auto gre_p = new TGraph(14,Layer,Proton_Acc0);
    auto gre_d = new TGraph(14,Layer,Proton_Acc1);
    auto gre_c = new TGraph(14,Layer,Proton_Acc2);


    gre_p->SetTitle(Form("Incident %s ; BGO Layer; Percentile",string2));
    gre_p->SetMarkerStyle(22);
    gre_p->GetXaxis()->SetLimits(0,14);
    gre_p->GetYaxis()->SetRangeUser(0,1);
    gre_p->SetMarkerStyle(20);
    gre_p->SetMarkerColor(kRed);
    gre_p->SetLineColor(kRed);
    gre_d->SetMarkerStyle(21);
    gre_d->SetMarkerColor(kBlue);
    gre_d->SetLineColor(kBlue);
    gre_c->SetMarkerStyle(22);
    gre_c->SetMarkerColor(kOrange-3);
    gre_c->SetLineColor(kOrange-3);
    gre_p->Draw("ALP");
    gre_d->Draw("LPSAME");
    gre_c->Draw("LPSAME");

    auto legend1_2 = new TLegend(0.2,0.7,0.7,0.88);
    legend1_2->AddEntry(gre_p, "#Delta Layer = 0","p");
    legend1_2->AddEntry(gre_d, "#Delta Layer = 1","p");
    legend1_2->AddEntry(gre_c, "#Delta Layer >= 2","p");
    legend1_2->Draw();

    ///////////////////////////////

    for(int ii = 1 ; ii <= 14 ; ii++)
    {
        h_sur->SetBinContent(ii, ( 1e4 - h_int->Integral(0,ii) ) );
    }

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top
    TF1 *fitFunc0 = new TF1("fitFunc0", "[0]*exp(-x/[1])", 0, 5); fitFunc0->SetParameters(100, 10); fitFunc0->SetLineColor(kBlue);

    auto c2 = new TCanvas("c2","c2",1200,600);
    c2->Divide(2,1);
    c2->cd(1);
    h_int->SetTitle(";Layer;N_{interact} (Bin of Maximum Change Ratio)");
    h_int->Draw("");
    c2->cd(2);
    gPad->SetLogy();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    h_sur->GetYaxis()->SetRangeUser(1e2,2e4);
    h_sur->SetTitle(";Layer;N_{survive}");
    h_sur->Draw("");
    h_sur->Fit(fitFunc0, "R"); // 进行拟合
    double constant2   = fitFunc0->GetParameter(0);
    double lambda2     = fitFunc0->GetParameter(1);
    double lambda2_err = fitFunc0->GetParError(1);
    double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    double hi_section = 1 / (lambda2*25) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hi_section_err = hi_section * lambda2_err/lambda2; // barn
    latex.DrawLatex(0,pow(10,3.3),"Fitting Function: N_{leave} =N_{total} *exp(-x/#lambda)");
    latex.DrawLatex(0,pow(10,3.0),Form("Deuteron Fitting #lambda: %.2f mm",lambda2*25));

    //////////////////////////////////////////////

    auto c3 = new TCanvas("c3","c3",1200,1200);
    c3->Divide(2,2);
    c3->cd(1);
    gStyle->SetOptStat(0);
    h_Lay_MM1->SetTitle("Inelastic;log10(Emax/Emin);First Hadronic Layer");
    h_Lay_MM1->Draw("colz");

    c3->cd(2);
    gStyle->SetOptStat(0);
    h_Lay_MM2->SetTitle("Elastic;log10(Emax/Emin);First Hadronic Layer");
    h_Lay_MM2->Draw("colz");

    c3->cd(3);
    gStyle->SetOptStat(0);
    h_CR_MM3->SetTitle("Pass through;log10(Emax/Emin);First Hadronic Layer");
    h_CR_MM3->Draw("colz");

    ////////////////////////////////////////

    auto c4 = new TCanvas("c4","c4",600,600);
    c4->Divide(2,2);
    c4->cd(1);
    gStyle->SetOptStat(0);
    h_CR_MM1->SetTitle("Inelastic;log10(Emax/Emin);#sum log10(Change Rate)");
    h_CR_MM1->Draw("colz");

    c4->cd(2);
    gStyle->SetOptStat(0);
    h_CR_MM2->SetTitle("Elastic;log10(Emax/Emin);#sum log10(Change Rate)");
    h_CR_MM2->Draw("colz");

    c4->cd(3);
    gStyle->SetOptStat(0);
    h_CR_MM3->SetTitle("Pass through;log10(Emax/Emin);#sum log10(Change Rate)");
    h_CR_MM3->Draw("colz");


    //////////////////////////////////////////

    auto c5 = new TCanvas("c5","c5",1200,1200);
    c5->Divide(2,2);
    c5->cd(1);
    gPad->SetLogz();
    h_start_had1->SetTitle("Inelastic;Bin where continuous energy increase begins ; First Hadronic Layer");
    h_start_had1->Draw("colz");
    double cov2 = h_start_had1->GetCorrelationFactor();
    cout << " cov2 = " << cov2 << endl;

    c5->cd(2);
    gPad->SetLogz();
    h_start_had2->SetTitle("Elastic;Bin where continuous energy increase begins ; First Hadronic Layer");
    h_start_had2->Draw("colz");

    c5->cd(3);
    gPad->SetLogz();
    h_start_had3->SetTitle("Pass through;Bin where continuous energy increase begins ; First Hadronic Layer");
    h_start_had3->Draw("colz");

    //////////////////////////////////////////

    auto c6 = new TCanvas("c5","c5",1200,1200);
    c6->Divide(2,2);
    c6->cd(1);
    gPad->SetLogz();
    h_len_had1->SetTitle("Inelastic;# Continues Positive Bins; First Hadronic Layer");
    h_len_had1->Draw("colz");
 
    c6->cd(2);
    gPad->SetLogz();
    h_len_had2->SetTitle("Elastic;# Continues Positive Bins; First Hadronic Layer");
    h_len_had2->Draw("colz");

    c6->cd(3);
    gPad->SetLogz();
    h_len_had3->SetTitle("Pass through;# Continues Positive Bins; First Hadronic Layer");
    h_len_had3->Draw("colz");


    auto c7 = new TCanvas("c7","c7",1200,1200);
    c7->Divide(2,2);
    c7->cd(1);
    h_sp_rate1->SetTitle(";Energy increase from start to peak;Counts");
    h_sp_rate1->SetLineColor(kRed);
    h_sp_rate2->SetLineColor(kBlue);
    h_sp_rate3->SetLineColor(kOrange-3);
    h_sp_rate1->Draw("hist");
    h_sp_rate2->Draw("histsame");
    h_sp_rate3->Draw("histsame");
    legend0->Draw();

    c7->cd(2);
    h_sp_bin1->SetTitle(";Bin length from start to peak;Counts");
    h_sp_bin1->SetLineColor(kRed);
    h_sp_bin2->SetLineColor(kBlue);
    h_sp_bin3->SetLineColor(kOrange-3);
    h_sp_bin1->Draw("hist");
    h_sp_bin2->Draw("histsame");
    h_sp_bin3->Draw("histsame");
    legend0->Draw();

    c7->cd(3);
    h_peak_val1->SetTitle(";Maximun Change Ratio;Counts");
    h_peak_val1->SetLineColor(kRed);
    h_peak_val2->SetLineColor(kBlue);
    h_peak_val3->SetLineColor(kOrange-3);
    h_peak_val1->Draw("hist");
    h_peak_val2->Draw("histsame");
    h_peak_val3->Draw("histsame");
    legend0->Draw();




}