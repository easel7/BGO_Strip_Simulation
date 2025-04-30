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

void Draw_Pattern2()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    int p_FH_Type;
    int p_FH_Lay;
    int p_Nhits;

    const char* string1;
    const char* string2 = "Proton_1000GeV";
    // const char* string2 = "Proton_10000GeV";
    // const char* string2 = "Deuteron_1000GeV";
    // const char* string2 = "Deuteron_10000GeV";

    auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/%s.root",string2));
    // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("BarEnergyVector", &p_EnergyVec);
    proton_tree->SetBranchAddress("RMS"            , &p_RMSVec);
    proton_tree->SetBranchAddress("First_Had_Depth", &p_FH_Dep);
    proton_tree->SetBranchAddress("First_Had_Layer", &p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type" , &p_FH_Type);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);

    int nbins = 30;
    double xmin = 1e-2;
    double xmax = 1e1;

    double logxmin = log10(xmin);
    double logxmax = log10(xmax);

    std::vector<double> bin_edges(nbins + 1);
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * (logxmax - logxmin) / nbins);
    }

    auto h_max_min0 = new TH1D("h_max_min0","h_max_min0",60,-1,5); 
    auto h_max_min1 = new TH1D("h_max_min1","h_max_min1",60,-1,5); // Inelastic
    auto h_max_min2 = new TH1D("h_max_min2","h_max_min2",60,-1,5); // Elastic
    auto h_max_min3 = new TH1D("h_max_min3","h_max_min3",60,-1,5); // Pass

    auto h_contin_0 = new TH1D("h_contin_0","h_contin_0",14,0,14); 
    auto h_contin_1 = new TH1D("h_contin_1","h_contin_1",14,0,14); // Inelastic
    auto h_contin_2 = new TH1D("h_contin_2","h_contin_2",14,0,14); // Elastic
    auto h_contin_3 = new TH1D("h_contin_3","h_contin_3",14,0,14); // Pass

    auto h_change_0 = new TH1D("h_change_0","h_change_0",60,-1,5); 
    auto h_change_1 = new TH1D("h_change_1","h_change_1",60,-1,5); // Inelastic
    auto h_change_2 = new TH1D("h_change_2","h_change_2",60,-1,5); // Elastic
    auto h_change_3 = new TH1D("h_change_3","h_change_3",60,-1,5); // Pass


    auto g_sum_len0 = new TGraph();      auto h_poi_had0 = new TH2I("h_poi_had0","h_poi_had0",14,0,14,16,-2,14);
    auto g_sum_len1 = new TGraph();      auto h_poi_had1 = new TH2I("h_poi_had1","h_poi_had1",14,0,14,16,-2,14);
    auto g_sum_len2 = new TGraph();      auto h_poi_had2 = new TH2I("h_poi_had2","h_poi_had2",14,0,14,16,-2,14);
    auto g_sum_len3 = new TGraph();      auto h_poi_had3 = new TH2I("h_poi_had3","h_poi_had3",14,0,14,16,-2,14);

    int point_counter   = 0; int pc   = 0;
    int point_counter_i = 0; int pc_i = 0;
    int point_counter_e = 0; int pc_e = 0;
    int point_counter_p = 0; int pc_p = 0;

    auto h1         = new TH2D("h1","h1",60,-1,5,16,-2,14);
    auto h2         = new TH1D("h2","h2",14,0,14);
    auto hC2        = new TH1D("hC2","hC2",14,0,14);
    auto h3         = new TH2D("h3","h3",60,-1,5,60,-1,5);
    auto h4         = new TH1D("h4","h4",40,0,2);

    // cout  << proton_tree->GetEntries() << endl;
    // Long64_t entry  = 99;   
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    // for (Long64_t entry = 0; entry < 100; entry++)
    {        
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double rate_max_min  = 0;
        double rate_sum      = 0;
        int    rate_len      = 0;
        int    rate_poi      = 0;
        double max_rate      = 0;
        int max_rate_index   = 0;

        proton_tree->GetEntry(entry);
        if (p_Nhits < 10 ) continue;
        // cout << entry << endl;
        // cout << " Depth " << p_FH_Dep << " layer " <<  p_FH_Lay << " Type =" << p_FH_Type << endl;
        
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
            if(layer==0)
            {
                bar_Change_info[0] = log10(bar_Energy_info[0] / 0.02);
            }
            else // (layer>0) 
            {
                if( bar_Energy_info[layer-1] == 0) { bar_Change_info[layer-1] = -5 ;  }//  cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
                else if( bar_Energy_info[layer] == 0) { bar_Change_info[layer-1] = -4 ; }//  cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
                else {bar_Change_info[layer] = log10(bar_Energy_info[layer]/bar_Energy_info[layer-1]); }// cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
            }
        }

        FindMaxPositiveSegment(bar_Change_info,14,rate_sum,rate_len,rate_poi);
        FindMaxValueInPositiveSegment(bar_Change_info,rate_poi,rate_len,max_rate,max_rate_index);
        // cout << " begin to increase bin = " << rate_poi << endl;
        // cout << "Max Positive Bin Length = " << rate_len << endl;
        // cout << "Max Increase Rate = " <<  max_rate << endl;
        // cout << "Max Increase Rate Bin = " << max_rate_index << endl;
        rate_max_min = MaxMinRatio(bar_Energy_info,14);    
        h_max_min0->Fill(log10(rate_max_min)); 
        h_change_0->Fill(log10(rate_sum)); 
        h_contin_0->Fill(rate_len); 
        g_sum_len0->SetPoint(point_counter++,rate_sum,rate_len); 
        h_poi_had0->Fill(max_rate_index,p_FH_Lay);
        
        
        if (rate_max_min >  1e2 && rate_len > 5 ) 
        {h2->Fill(max_rate_index); }
        // if (rate_max_min < 100 && rate_max_min>60  && p_FH_Type == 1 && p_FH_Lay>11) { cout << entry << " , " <<  rate_max_min <<  endl; }
        if (rate_max_min >  1e2 )
        {}
        
        
        if(p_FH_Type == 1)       {  string1 = "Inelastic"; 
            h_max_min1->Fill(log10(rate_max_min)); 
            h_change_1->Fill(log10(rate_sum)); 
            h_contin_1->Fill(rate_len); 
            g_sum_len1->SetPoint(point_counter_i++,rate_sum,rate_len+0.1); 
            h_poi_had1->Fill(max_rate_index,p_FH_Lay); 
            h4->Fill(max_rate);
            h3->Fill(log10(rate_max_min),rate_sum);
        } 
        else if (p_FH_Type == 2) {  string1 = "Elastic";   
            h_max_min2->Fill(log10(rate_max_min)); 
            h_change_2->Fill(log10(rate_sum)); 
            h_contin_2->Fill(rate_len); 
            g_sum_len2->SetPoint(point_counter_e++,rate_sum,rate_len); 
            h_poi_had2->Fill(max_rate_index,p_FH_Lay);
        }
        else                     {  string1 = "Pass";      
            h_max_min3->Fill(log10(rate_max_min)); 
            h_change_3->Fill(log10(rate_sum)); 
            h_contin_3->Fill(rate_len); 
            g_sum_len3->SetPoint(point_counter_p++,rate_sum,rate_len-0.1);
            h_poi_had3->Fill(max_rate_index,p_FH_Lay);

            h1->Fill(log10(rate_max_min),p_FH_Lay);
        }
        
    }
    auto c2    = new TCanvas("c2","c2",2100,2100);
    c2->Divide(3,3);
    c2->cd(1);
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

    auto legend1 = new TLegend(0.7,0.7,0.88,0.88);
    legend1->AddEntry(h_max_min1, "FH Inlastic","l");
    legend1->AddEntry(h_max_min2, "FH Elastic","l");
    legend1->AddEntry(h_max_min3, "Pass through","l");
    legend1->Draw();

    c2->cd(2);
    gStyle->SetOptStat(0);
    // gPad->SetLogx(1);
    gPad->SetLogy(1);
    h_change_1->GetXaxis()->SetLimits(1e-2, 1e2);
    h_change_1->SetTitle(";#sum log10(Change Rate);Counts");
    h_change_1->SetLineColor(kRed);
    h_change_2->SetLineColor(kBlue);
    h_change_3->SetLineColor(kOrange-3);
    // h_change_1->Draw("hist");
    // h_change_2->Draw("histsame");
    // h_change_3->Draw("histsame");

    legend1->Draw();

    c2->cd(3);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    h_contin_1->GetYaxis()->SetRangeUser(1e0, 3e3); 
    h_contin_1->SetTitle(";# Continues Positive Bins;Counts");
    h_contin_1->SetLineColor(kRed);
    h_contin_2->SetLineColor(kBlue);
    h_contin_3->SetLineColor(kOrange-3);
    h_contin_1->Draw("hist");
    h_contin_2->Draw("histsame");
    h_contin_3->Draw("histsame");

    legend1->Draw();


    c2->cd(4);
    gStyle->SetOptStat(0);
    // h1->SetTitle("Inelastic;log10(Emax/Emin);Bin of Maximum Change Ratio");
    h1->SetTitle("Pass throught;log10(Emax/Emin);First Hadronic Layer");
    h1->Draw("colz");


    c2->cd(5);
    gPad->SetLogx();
    g_sum_len0->SetTitle(";#sum log10(Change Rate);# Continues Positive Bins");
    g_sum_len0->GetXaxis()->SetLimits(1e-2, 100);        // X 轴范围
    g_sum_len0->GetYaxis()->SetRangeUser(0, 14);    // Y 轴范围
    g_sum_len0->SetMarkerStyle(20);  
    g_sum_len0->SetMarkerColorAlpha(kBlack, 0);  // 0.0 = fully transparent, 1.0 = fully opaque
    g_sum_len0->SetMarkerSize(0);

    g_sum_len1->SetMarkerStyle(20);  
    g_sum_len1->SetMarkerColorAlpha(kRed, 0.1);  // 0.0 = fully transparent, 1.0 = fully opaque
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

    auto legend2 = new TLegend(0.7,0.7,0.88,0.88);
    legend2->AddEntry(g_sum_len1, "FH Inelastic","p");
    legend2->AddEntry(g_sum_len2, "FH Elastic","p");
    legend2->AddEntry(g_sum_len3, "Pass through","p");
    legend2->Draw();


    c2->cd(6);
    // gPad->SetLogy();
    gStyle->SetOptStat(0);
    h3->SetTitle("Elastic;log10(Emax/Emin);#sum log10(Change Rate)");
    h3->Draw("colz");
    
    c2->cd(7);
    gPad->SetLogz();
    h_poi_had1->SetTitle("Inelastic;Bin of Maximum Change Ratio; First Hadronic Layer");
    h_poi_had1->Draw("colz");
    double cov1 = h_poi_had1->GetCorrelationFactor();
    cout << " cov1 = " << cov1 << endl;

    c2->cd(8);
    gPad->SetLogz();
    h_poi_had2->SetMinimum(h_poi_had2->GetMinimum(1e-10)); // 取非零最小值
    h_poi_had2->SetMaximum(h_poi_had2->GetMaximum());
    h_poi_had2->SetTitle("Elastic;Bin of Maximum Change Ratio; First Hadronic Layer");
    h_poi_had2->Draw("colz");
    double cov2 = h_poi_had2->GetCorrelationFactor();
    cout << " cov2 = " << cov2 << endl;

    c2->cd(9);
    gPad->SetLogz();
    h_poi_had3->SetMinimum(h_poi_had3->GetMinimum(1e-10)); // 取非零最小值
    h_poi_had3->SetMaximum(h_poi_had3->GetMaximum());
    h_poi_had3->SetTitle("Pass through;Bin of Maximum Change Ratio; First Hadronic Layer");
    h_poi_had3->Draw("colz");

    c2->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/Figures/%s_FIG.pdf",string2));
    std::cout << "Number of valid points: " << g_sum_len0->GetN() << std::endl;

    for(int ii = 1 ; ii <= 14 ; ii++)
    {
        hC2->SetBinContent(ii, ( 1e4 - h2->Integral(0,ii) ) );
    }

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 0, 6); fitFunc2->SetParameters(100, 10); fitFunc2->SetLineColor(kBlue);

    auto c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1);
    c1->cd(1);
    h2->SetTitle(";Layer;N_{interact}");
    h2->Draw("");
    c1->cd(2);
    gPad->SetLogy();
    // gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    hC2->GetYaxis()->SetRangeUser(1e2,2e4);
    hC2->SetTitle(";Layer;N_{survive}");
    hC2->Draw("");
    hC2->Fit(fitFunc2, "R"); // 进行拟合
    double constant2   = fitFunc2->GetParameter(0);
    double lambda2     = fitFunc2->GetParameter(1);
    double lambda2_err = fitFunc2->GetParError(1);
    double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    double hi_section = 1 / (lambda2*25) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hi_section_err = hi_section * lambda2_err/lambda2; // barn
    latex.DrawLatex(0,pow(10,3.3),"Fitting Function: N_{leave} =N_{total} *exp(-x/#lambda)");
    latex.DrawLatex(0,pow(10,3.0),Form("Deuteron Fitting #lambda: %.2f mm",lambda2*25));

    auto c3 = new TCanvas("c3","c3",600,600);
    h4->Draw();
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
