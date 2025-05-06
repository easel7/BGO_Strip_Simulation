std::vector<int> g_fit_bars;
std::vector<double> g_fit_energies;
double g_fit_total_energy = 0;

int FindMaxMiddleIndex(std::vector<double>* p_EnergyVec, int layer);

void FindMaxPositiveSegment(const double array[], int size, double& out_sum, int& out_len, int& out_start_index);

void FindMaxValueInPositiveSegment(const double array[], int start_index, int length, double& out_max_value, int& out_max_index);

double MaxMinRatio(const double array[], int size);

void PrepareFitData(
    std::vector<double>* p_EnergyVec,
    int layer_start,
    int layer_end,
    std::vector<int>& bars,
    std::vector<double>& energies,
    double& total_energy);

void FitAxisFunction(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

double AccumIncreaseToPeak(const double array[], int start_idx, int end_idx);

void Draw_Pattern2_PowerLaw()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    double p_Total_E;
    int p_FH_Type;
    int p_FH_Lay;
    int p_Nhits;
    int p_FH_Second;


    const char* string1;
    // const char* string2 = "Proton_PowerLaw";
    const char* string2 = "Deuteron_PowerLaw";

    auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Weight/%s.root",string2));
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("BarEnergyVector", &p_EnergyVec);
    proton_tree->SetBranchAddress("RMS"            , &p_RMSVec);
    proton_tree->SetBranchAddress("First_Had_Depth", &p_FH_Dep);
    proton_tree->SetBranchAddress("First_Had_Layer", &p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type" , &p_FH_Type);
    proton_tree->SetBranchAddress("First_Had_Second", &p_FH_Second);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);
    proton_tree->SetBranchAddress("Total_E"        , &p_Total_E);


    int nbins = 30;
    double xmin = 1e-2;
    double xmax = 1e1;

    double logxmin = log10(xmin);
    double logxmax = log10(xmax);

    std::vector<double> bin_edges(nbins + 1);
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * (logxmax - logxmin) / nbins);
    }

    auto h_max_min1 = new TH1D("h_max_min1","h_max_min1",60,-1,5); // Inelastic
    auto h_max_min2 = new TH1D("h_max_min2","h_max_min2",60,-1,5); // Elastic
    auto h_max_min3 = new TH1D("h_max_min3","h_max_min3",60,-1,5); // Pass

    auto h_change_1 = new TH1D("h_change_1","h_change_1",60,-1,5); // Inelastic
    auto h_change_2 = new TH1D("h_change_2","h_change_2",60,-1,5); // Elastic
    auto h_change_3 = new TH1D("h_change_3","h_change_3",60,-1,5); // Pass

    auto h_contin_1 = new TH1D("h_contin_1","h_contin_1",14,0,14); // Inelastic
    auto h_contin_2 = new TH1D("h_contin_2","h_contin_2",14,0,14); // Elastic
    auto h_contin_3 = new TH1D("h_contin_3","h_contin_3",14,0,14); // Pass

    auto h_layer_1 = new TH1D("h_layer_1","h_layer_1",14,0,14); // Inelastic
    auto h_layer_2 = new TH1D("h_layer_2","h_layer_2",14,0,14); // Elastic
    auto h_layer_3 = new TH1D("h_layer_3","h_layer_3",14,0,14); // Pass

    auto h_sec_1 = new TH1D("h_sec_1","h_sec_1",20,0,200); // Inelastic
    auto h_sec_2 = new TH1D("h_sec_2","h_sec_2",20,0,200); // Elastic
    auto h_sec_3 = new TH1D("h_sec_3","h_sec_3",20,0,200); // Pass

    auto h_RMS_1 = new TH1D("h_RMS_1","h_RMS_1",20,0,200); // Inelastic
    auto h_RMS_2 = new TH1D("h_RMS_2","h_RMS_2",20,0,200); // Elastic
    auto h_RMS_3 = new TH1D("h_RMS_3","h_RMS_3",20,0,200); // Pass

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
    auto h_sur        = new TH1D("h_sur","h_sur",14,0,14);

    auto h_peak_rate1 = new TH1D("h_peak_rate1","h_peak_rate1",50,-5,0); // Inelastic
    auto h_peak_rate2 = new TH1D("h_peak_rate2","h_peak_rate2",50,-5,0); // Elastic
    auto h_peak_rate3 = new TH1D("h_peak_rate3","h_peak_rate3",50,-5,0); // Pass

    int through =0; 
    int el = 0;
  
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    {        
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引

        proton_tree->GetEntry(entry);
        if (p_Nhits < 10 ) continue;
        if (log10(p_Total_E) < 2.8 || log10(p_Total_E)> 3.0 ) continue;
        
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
            minuit1.DefineParameter(0, "bar_even", g_fit_bars[1], 0.01, 2, 19);
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
        }
        bar_Change_info[0] = log10(bar_Energy_info[0] / 0.023);
        for(int layer = 1 ; layer<14 ; layer++)
        {
            if( bar_Energy_info[layer-1] == 0 || bar_Energy_info[layer] == 0) { bar_Change_info[layer-1] = -5 ;  }//  cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
            else {bar_Change_info[layer] = log10(bar_Energy_info[layer]/bar_Energy_info[layer-1]); }// cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
        }

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
        h_peak_had0->Fill(seg_peak_idx,p_FH_Lay);
        
        if (rate_max_min >  1e2 && seg_len > 5 ) 
        {h_int->Fill(seg_peak_idx); }
        // if (rate_max_min < 100 && rate_max_min>60  && p_FH_Type == 1 && p_FH_Lay>11) { cout << entry << " , " <<  rate_max_min <<  endl; }
        
        if ((*p_RMSVec)[13] < 15 &&  (*p_RMSVec)[12] < 15 )
        { through ++ ;}
        // cout << entry << " , " << p_Total_E << endl;
        
        
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
            h_layer_1->Fill(p_FH_Lay);
            h_sec_1->Fill(p_FH_Second);
            h_peak_rate1->Fill(seg_sum_to_peak - seg_sum);

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
            el ++ ;
            h_layer_2->Fill(p_FH_Lay);
            h_sec_2->Fill(p_FH_Second);

            h_peak_rate2->Fill(seg_sum_to_peak - seg_sum);


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
            h_layer_3->Fill(p_FH_Lay);
            h_sec_3->Fill(p_FH_Second);

            h_peak_rate2->Fill(seg_sum_to_peak - seg_sum);

        }
        
    }

    cout << "through " << through << endl;
    cout << "through el " << el << endl;

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
    h_max_min1->Draw("hist");
    h_max_min2->Draw("histsame");
    h_max_min3->Draw("histsame");

    auto legend0 = new TLegend(0.12,0.7,0.32,0.88);
    legend0->AddEntry(h_max_min1, "FH Inlastic","l");
    legend0->AddEntry(h_max_min2, "FH Elastic","l");
    legend0->AddEntry(h_max_min3, "Pass through","l");
    legend0->Draw();

    c0->cd(2);
    gStyle->SetOptStat(0);
    // gPad->SetLogx(1);
    gPad->SetLogy(1);
    // h_change_1->GetXaxis()->SetLimits(1e-2, 3e3);
    h_change_1->SetTitle(";#sum log10(Change Rate);Counts");
    h_change_1->SetLineColor(kRed);
    h_change_2->SetLineColor(kBlue);
    h_change_3->SetLineColor(kOrange-3);
    h_change_1->Draw("hist");
    h_change_2->Draw("histsame");
    h_change_3->Draw("histsame");

    legend0->Draw();

    c0->cd(3);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    // h_contin_1->GetYaxis()->SetRangeUser(1e0, 3e3); 
    h_contin_1->SetTitle(";# Continues Positive Bins;Counts");
    h_contin_1->SetLineColor(kRed);
    h_contin_2->SetLineColor(kBlue);
    h_contin_3->SetLineColor(kOrange-3);
    h_contin_1->Draw("hist");
    h_contin_2->Draw("histsame");
    h_contin_3->Draw("histsame");

    legend0->Draw();

    c0->cd(4);
    // gPad->SetLogx();
    g_sum_len0->SetTitle(";#sum log10(Change Rate);# Continues Positive Bins");
    g_sum_len0->GetXaxis()->SetLimits(1, 5);        // X 轴范围
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

    auto c1 = new TCanvas("c1","c1",1200,1200);
    c1->Divide(2,2);
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

    ///////////////////////////////

    for(int ii = 1 ; ii <= 14 ; ii++)
    {
        h_sur->SetBinContent(ii, ( point_counter - h_int->Integral(0,ii) ) );
    }

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top
    TF1 *fitFunc0 = new TF1("fitFunc0", "[0]*exp(-x/[1])", 1, 5); fitFunc0->SetParameters(100, 10); fitFunc0->SetLineColor(kBlue);

    auto c2 = new TCanvas("c2","c2",1200,600);
    c2->Divide(2,1);
    c2->cd(1);
    h_layer_1->SetTitle(";Layer;N_{interact}");
    h_layer_1->SetLineColor(kRed);
    h_layer_2->SetLineColor(kBlue);
    h_layer_3->SetLineColor(kOrange-3);
    h_layer_1->Draw("hist");
    h_layer_2->Draw("histsame");
    h_layer_3->Draw("histsame");

    h_int->SetLineColor(kBlack);
    h_int->SetLineStyle(2);
    h_int->Draw("histsame");

    auto legend2 = new TLegend(0.4,0.65,0.88,0.88);
    legend2->AddEntry(h_layer_1, "FH Inelastic","l");
    legend2->AddEntry(h_layer_2, "FH Elastic","l");
    legend2->AddEntry(h_layer_3, "Pass through","l");
    legend2->AddEntry(h_int, "Bin of Maximum Change Ratio","l");
    legend2->Draw();

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
    latex.DrawLatex(0,pow(10,3.0),Form("Fitting #lambda: %.2f mm",lambda2*25));

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
    gPad->SetLogy();
    h_sp_rate1->SetTitle(";log10(Energy increase from start to peak);Counts");
    h_sp_rate1->SetLineColor(kRed);
    h_sp_rate2->SetLineColor(kBlue);
    h_sp_rate3->SetLineColor(kOrange-3);
    h_sp_rate1->Draw("hist");
    h_sp_rate2->Draw("histsame");
    h_sp_rate3->Draw("histsame");
    legend0->Draw();

    c7->cd(2);
    gPad->SetLogy();
    h_sp_bin1->SetTitle(";Bin length from start to peak;Counts");
    h_sp_bin1->SetLineColor(kRed);
    h_sp_bin2->SetLineColor(kBlue);
    h_sp_bin3->SetLineColor(kOrange-3);
    h_sp_bin1->Draw("hist");
    h_sp_bin2->Draw("histsame");
    h_sp_bin3->Draw("histsame");
    legend0->Draw();

    c7->cd(3);
    gPad->SetLogy();
    h_peak_val1->SetTitle(";Peak Change Ratio;Counts");
    h_peak_val1->SetLineColor(kRed);
    h_peak_val2->SetLineColor(kBlue);
    h_peak_val3->SetLineColor(kOrange-3);
    h_peak_val1->Draw("hist");
    h_peak_val2->Draw("histsame");
    h_peak_val3->Draw("histsame");
    legend0->Draw();


    c7->cd(4);
    h_peak_rate1->SetTitle(";(Start-Peak)/(Start - End);Counts");
    h_peak_rate1->SetLineColor(kRed);
    h_peak_rate2->SetLineColor(kBlue);
    h_peak_rate3->SetLineColor(kOrange-3);
    h_peak_rate1->Draw("hist");
    h_peak_rate2->Draw("histsame");
    h_peak_rate3->Draw("histsame");
    legend0->Draw();

    auto c8 = new TCanvas("c8","c8",1200,1200);
    c8->Divide(2,2);
    c8->cd(1);
    gPad->SetLogy();
    h_sec_1->SetTitle(";Second;Counts");
    h_sec_1->SetLineColor(kRed);
    h_sec_2->SetLineColor(kBlue);
    h_sec_3->SetLineColor(kOrange-3);
    h_sec_1->Draw("hist");
    h_sec_2->Draw("histsame");
    h_sec_3->Draw("histsame");
    legend0->Draw();
}

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

double AccumIncreaseToPeak(const double array[], int start_idx, int end_idx) 
{
    double sum = 0;
    for (int i = start_idx; i <= end_idx; ++i) {
        sum += array[i];
    }
    return sum;
}
