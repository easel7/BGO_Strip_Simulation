std::vector<int> g_fit_bars;
std::vector<double> g_fit_energies;
double g_fit_total_energy = 0;

void FindMaxPositiveBinSegment(TH1D* hist, double& out_sum, int& out_len, int& out_start_bin);

void FindMaxValueInPositiveSegment(TH1D* hist, int start_bin, int len, double& out_max_value, int& out_max_bin);

void FitAxisFunction(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

void PrepareFitData(
    std::vector<double>* p_EnergyVec,
    int layer_start,
    int layer_end,
    std::vector<int>& bars,
    std::vector<double>& energies,
    double& total_energy);

int FindMaxMiddleIndex(std::vector<double>* p_EnergyVec, int layer);

double MaxMinRatio(TH1D* hist);

void Draw_Pattern()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    int p_FH_Type;
    int p_FH_Lay;
    int p_Nhits;

    const char* string1;
    const char* string2 = "Proton_1000GeV";
    // const char* string2 = "Deuteron_1000GeV";

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
    auto h_max_min0 = new TH1D("h_max_min0","h_max_min0",60,-1,5); 
    auto h_max_min1 = new TH1D("h_max_min1","h_max_min1",60,-1,5); // Inelastic
    auto h_max_min2 = new TH1D("h_max_min2","h_max_min2",60,-1,5); // Elastic
    auto h_max_min3 = new TH1D("h_max_min3","h_max_min3",60,-1,5); // Pass

    auto g_sum_len0 = new TGraph();
    auto g_sum_len1 = new TGraph();
    auto g_sum_len2 = new TGraph();
    auto g_sum_len3 = new TGraph();

    int point_counter   = 0;
    int point_counter_i = 0;
    int point_counter_e = 0;
    int point_counter_p = 0;

    // cout  << proton_tree->GetEntries() << endl;
    Long64_t entry  = 5668;   
    // for (Long64_t entry = 0; entry < 100; entry++)
    {
        auto c1    = new TCanvas("c1","c1",1400,1400);
        auto hXZ   = new TH2D("hXZ","BGO X-Z Plane",22,-11,11,14,0,14);
        auto hYZ   = new TH2D("hYZ","BGO Y-Z Plane",22,-11,11,14,0,14);
        auto hBGO1 = new TH1D("hBGO1","BGO Core Axis Energy Deposit",14,0,14); 
        auto hBGO2 = new TH1D("hBGO2","Deposit Energy Change Ratio",14,0,14); 
        
        auto *g_core0 = new TGraph();          auto *g_ch_0 = new TGraph();
        auto *g_core1 = new TGraph();          auto *g_ch_1 = new TGraph();
        auto *g_core2 = new TGraph();          auto *g_ch_2 = new TGraph();
        auto *g_core3 = new TGraph();          auto *g_ch_3 = new TGraph();

        auto    *box0 = new TBox();
        auto    *box1 = new TBox();
        
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double rate_sum      = 0;
        int    rate_len      = 0;
        int    rate_poi      = 0;
        double max_rate      = 0;
        int max_rate_index   = 0;

        proton_tree->GetEntry(entry);
        // if (p_Nhits < 10 ) continue;
        // cout << " RMS at 0 layer = " << p_FH_Lay << endl;

        for(size_t i = 0; i < p_EnergyVec->size(); i++) 
        {
            int layer = i / 22; 
            int bar = i % 22;
            double energy = 0;
            if ((*p_EnergyVec)[i] < 1e-2)  {energy = -5;              }
            else                           {energy = log10((*p_EnergyVec)[i]);}

            // cout << "ECEnergyVector[" << i << "] = " << (*p_EnergyVec)[i] << " , layer " <<  layer  << " , bar " <<  bar << " Calculated Energy " << energy << " , "<< layer % 2 << endl;
            if(layer % 2 == 0)
            {
                hXZ->SetBinContent(bar+1, 14-layer, energy);
                hYZ->SetBinContent(bar+1, 14-layer, -5);                
            }
            else
            {
                hYZ->SetBinContent(bar+1, 14-layer, energy);
                hXZ->SetBinContent(bar+1, 14-layer, -5);
            }
        }

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
                    bar_info[0] = bar1; // 奇数层
                    bar_info[1] = bar2; // 偶数层
                } else {
                    bar_info[0] = bar2; // 奇数层
                    bar_info[1] = bar1; // 偶数层
                }
                // 输出确定的 bar_info
                // cout << "Directly determined bar_info: " 
                // << "bar1 = " << bar_info[0] << ", bar2 = " << bar_info[1] << endl;

                // 设置标志为 true，表示已经赋值
                bar_info_assigned = true;
                break;  // 一旦赋值，跳出循环
            }
        }
        
        if (!bar_info_assigned) {
            cout << "No bar_info assigned, starting fit to determine cluster trajectory." << endl;
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


        box0->SetX1(bar_info[0] - 9);
        box0->SetY1(0);
        box0->SetX2(bar_info[0] - 12);
        box0->SetY2(14 - layer_start );

        box1->SetX1(bar_info[1] - 9);
        box1->SetY1(0);
        box1->SetX2(bar_info[1] - 12);
        box1->SetY2(14 - layer_start);

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
                if( bar_Energy_info[layer-1] == 0) { bar_Change_info[layer-1] = -5 ;   cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
                else if( bar_Energy_info[layer] == 0) { bar_Change_info[layer-1] = -4 ;   cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
                else {bar_Change_info[layer] = log10(bar_Energy_info[layer]/bar_Energy_info[layer-1]);  cout << "entry = " << entry << " , layer "<< layer-1 << " , rate " << bar_Change_info[layer-1] << endl;}
            }
            hBGO1->SetBinContent(layer+1, bar_Energy_info[layer]);
            hBGO2->SetBinContent(layer+1, bar_Change_info[layer]);
        }
        FindMaxPositiveBinSegment(hBGO2,rate_sum,rate_len,rate_poi);
        FindMaxValueInPositiveSegment(hBGO2,rate_poi,rate_len,max_rate,max_rate_index);

        cout << "Begin to increase bin = " << rate_poi << endl;
        cout << "Max Positive Bin Length = " << rate_len << endl;
        cout << "Max Increase Rate = " <<  max_rate << endl;
        cout << "Max Increase Rate Bin = " << max_rate_index << endl;


        auto *line1 = new TLine(-11,14 - p_FH_Lay,11,14 - p_FH_Lay);
        auto *line2 = new TLine(p_FH_Lay,hBGO1->GetMinimum(),p_FH_Lay,hBGO1->GetMaximum());
        auto *line3 = new TLine(p_FH_Lay,hBGO2->GetMinimum(),p_FH_Lay,hBGO2->GetMaximum());
        auto *line4 = new TLine(0,0,14,0);
        g_core0->SetPoint(0, bar_info[0] - 10.5, 15 - rate_poi );                       g_ch_0->SetPoint(0, bar_info[0] - 10.5, 15 - max_rate_index );              
        g_core1->SetPoint(0, bar_info[1] - 10.5, 15 - rate_poi );                       g_ch_1->SetPoint(0, bar_info[1] - 10.5, 15 - max_rate_index );              
        g_core2->SetPoint(0, rate_poi - 0.5, hBGO1->GetBinContent(rate_poi) );          g_ch_2->SetPoint(0, max_rate_index - 0.5, hBGO1->GetBinContent(max_rate_index) );
        g_core3->SetPoint(0, rate_poi - 0.5, hBGO2->GetBinContent(rate_poi) );          g_ch_3->SetPoint(0, max_rate_index - 0.5, hBGO2->GetBinContent(max_rate_index) );
        line4->SetLineColor(kBlack);
        line4->SetLineWidth(2);
        line4->SetLineStyle(9);
        h_max_min0->Fill(log10(MaxMinRatio(hBGO1)));
        g_sum_len0->SetPoint(point_counter++,rate_sum,rate_len);
        if(p_FH_Type == 1)       {  string1 = "Inelastic";  g_core0->SetMarkerColor(kRed);    g_core1->SetMarkerColor(kRed);    g_core2->SetMarkerColor(kRed);    g_core3->SetMarkerColor(kRed);      g_ch_0->SetMarkerColor(kRed);      g_ch_1->SetMarkerColor(kRed);      g_ch_2->SetMarkerColor(kRed);      g_ch_3->SetMarkerColor(kRed);  box0->SetLineColor(kRed);     box1->SetLineColor(kRed);       line1->SetLineColor(kRed);     line2->SetLineColor(kRed);     line3->SetLineColor(kRed);     h_max_min1->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));  g_sum_len1->SetPoint(point_counter_i++,rate_sum,rate_len+0.1);}
        else if (p_FH_Type == 2) {  string1 = "Elastic";    g_core0->SetMarkerColor(kMagenta);g_core1->SetMarkerColor(kMagenta);g_core2->SetMarkerColor(kMagenta);g_core3->SetMarkerColor(kMagenta);  g_ch_0->SetMarkerColor(kMagenta);  g_ch_1->SetMarkerColor(kMagenta);  g_ch_2->SetMarkerColor(kMagenta);  g_ch_3->SetMarkerColor(kRed);  box0->SetLineColor(kMagenta); box1->SetLineColor(kMagenta);   line1->SetLineColor(kMagenta); line2->SetLineColor(kMagenta); line3->SetLineColor(kMagenta); h_max_min2->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));  g_sum_len2->SetPoint(point_counter_e++,rate_sum,rate_len);    }
        else                     {  string1 = "Pass";       g_core0->SetMarkerColor(kGray);   g_core1->SetMarkerColor(kGray);   g_core2->SetMarkerColor(kGray);   g_core3->SetMarkerColor(kGray);     g_ch_0->SetMarkerColor(kGray);     g_ch_1->SetMarkerColor(kGray);     g_ch_2->SetMarkerColor(kGray);     g_ch_3->SetMarkerColor(kRed);  box0->SetLineColor(kGray);    box1->SetLineColor(kGray);      line1->SetLineColor(kGray);    line2->SetLineColor(kGray);    line3->SetLineColor(kGray);    h_max_min3->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));  g_sum_len3->SetPoint(point_counter_p++,rate_sum,rate_len-0.1);}
        
        box0->SetFillStyle(0);        box0->SetLineWidth(2);
        box1->SetFillStyle(0);        box1->SetLineWidth(2);
        g_core0->SetMarkerStyle(29);  g_core0->SetMarkerSize(2);    g_ch_0->SetMarkerStyle(22);   g_ch_0->SetMarkerSize(2);
        g_core1->SetMarkerStyle(29);  g_core1->SetMarkerSize(2);    g_ch_1->SetMarkerStyle(22);   g_ch_1->SetMarkerSize(2);
        g_core2->SetMarkerStyle(29);  g_core2->SetMarkerSize(2);    g_ch_2->SetMarkerStyle(22);   g_ch_2->SetMarkerSize(2);
        g_core3->SetMarkerStyle(29);  g_core3->SetMarkerSize(2);    g_ch_3->SetMarkerStyle(22);   g_ch_3->SetMarkerSize(2);
        line1->SetLineWidth(2);       line4->SetLineColor(kBlack);
        line2->SetLineWidth(2);       line4->SetLineWidth(2);
        line3->SetLineWidth(2);       line4->SetLineStyle(9);

        for (int layer = 0; layer < 14; layer++) {
            hXZ->SetBinContent(1, 14-layer, -5);
            hXZ->SetBinContent(22, 14-layer, -5);
            hYZ->SetBinContent(1, 14-layer, -5);
            hYZ->SetBinContent(22, 14-layer, -5);
        }
        double minVal = -2;  // 设置最小值
        double maxVal = 2;  // 设置最大值

        hXZ->SetMinimum(minVal); 
        hXZ->SetMaximum(maxVal); 
        hYZ->SetMinimum(minVal); 
        hYZ->SetMaximum(maxVal); 

        c1->Clear();
        TLatex latex;
        latex.SetTextSize(0.03);
        latex.SetTextColor(kBlack);
        c1->Divide(2,2);
        c1->cd(1);
        gStyle->SetOptStat(0);
        hXZ->Draw("colz");
        line1->Draw("same");
        g_core0->Draw("psame"); 
        g_ch_0->Draw("psame"); 

        box0->Draw("psame");

        
        gPad->Update();
        auto *palette1 = (TPaletteAxis*) hXZ->GetListOfFunctions()->FindObject("palette");
        if(palette1)
        {
            palette1->SetX1NDC(0.90);
            palette1->SetX2NDC(0.93);
            gPad->Modified();
            gPad->Update();
        }
        latex.DrawLatexNDC(0.82, 0.92, "log_{10}(E/GeV)");

        c1->cd(2);
        gStyle->SetOptStat(0);
        hYZ->Draw("colz");
        line1->Draw("same");
        g_core1->Draw("psame"); 
        g_ch_1->Draw("psame"); 

        box1->Draw("psame");

        gPad->Update();
        auto  *palette2 = (TPaletteAxis*) hYZ->GetListOfFunctions()->FindObject("palette");
        if(palette2)
        {
            palette2->SetX1NDC(0.90);
            palette2->SetX2NDC(0.93);
            gPad->Modified();
            gPad->Update();
        }
        latex.DrawLatexNDC(0.82, 0.92, "log_{10}(E/GeV)");


        c1->cd(3);
        gPad->SetLogy();
        hBGO1->SetLineColor(kBlack);
        hBGO1->SetLineWidth(2);
        hBGO1->SetTitle("BGO Core Axis Energy Deposit;BGO Layer; 3 Bars Energy / GeV");
        hBGO1->GetYaxis()->CenterTitle();
        hBGO1->GetYaxis()->SetTitleOffset(0.97);
        hBGO1->GetXaxis()->CenterTitle();
        hBGO1->Draw("hist");
        g_core2->Draw("psame");
        g_ch_2->Draw("psame"); 

        if(p_FH_Lay%2 != -1 ) {line2->Draw("");}

        c1->cd(4);
        hBGO2->SetLineColor(kBlack);
        hBGO2->SetLineWidth(2);
        hBGO2->SetTitle("BGO Core Axis Energy Deposit Change Rate;BGO Layer; log10(^{}E^{}_{i}/^{}E^{}_{i-1})");
        hBGO2->GetYaxis()->CenterTitle();
        hBGO2->GetYaxis()->SetRangeUser(-0.5,2.5);
        hBGO2->GetYaxis()->SetTitleOffset(1.05);
        hBGO2->GetXaxis()->CenterTitle();
        hBGO2->Draw("hist");
        line4->Draw();
        g_core3->Draw("psame");
        g_ch_3->Draw("psame"); 

        if(p_FH_Lay%2 != -1 ) {line3->Draw("");}
        if (entry < 100) { c1->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/Figures/%s/%s/%lld.pdf",string1,string2,entry)); }
    }
    auto c2    = new TCanvas("c2","c2",1400,700);
    c2->Divide(2,1);
    c2->cd(1);
    h_max_min1->SetTitle(";log10(Emax/Emin);Counts");

    h_max_min1->SetLineColor(kRed);
    h_max_min2->SetLineColor(kBlue);
    h_max_min3->SetLineColor(kOrange-3);

    h_max_min1->Draw("hist");
    h_max_min2->Draw("histsame");
    h_max_min3->Draw("histsame");

    c2->cd(2);
    gPad->SetLogx();
    g_sum_len0->SetTitle("; #sum log10(Change Rate);# Continues Positive Bins");
    g_sum_len0->GetXaxis()->SetLimits(1e-2, 100);        // X 轴范围
    g_sum_len0->GetYaxis()->SetRangeUser(0, 14);    // Y 轴范围
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

    std::cout << "Number of valid points: " << g_sum_len0->GetN() << std::endl;
}

void FindMaxPositiveBinSegment(TH1D* hist, double& out_sum, int& out_len, int& out_start_bin) 
{
    int nbins = hist->GetNbinsX();
    double max_sum = 0;
    int max_len = 0;
    int max_start_bin = -1;

    double curr_sum = 0;
    int curr_len = 0;
    int curr_start_bin = -1;

    for (int i = 1; i <= nbins; ++i) 
    {
        double content = hist->GetBinContent(i);
        // std::cout << "bin[" << i << "] = " << content << std::endl;

        if (content > 0) 
        {
            if (curr_len == 0) 
                curr_start_bin = i;  // 记录这一段开始的位置

            curr_sum += content;
            curr_len++;

            if (curr_sum > max_sum || (curr_sum == max_sum && curr_len > max_len)) 
            {
                max_sum = curr_sum;
                max_len = curr_len;
                max_start_bin = curr_start_bin;  // 更新最大值对应的起始位置
            }
        } 
        else 
        {
            curr_sum = 0;
            curr_len = 0;
            curr_start_bin = -1;
        }
    }

    out_sum = max_sum;
    out_len = max_len;
    out_start_bin = max_start_bin;
}


void FindMaxValueInPositiveSegment(TH1D* hist, int start_bin, int len, double& out_max_value, int& out_max_bin)
{
    out_max_value = -1e9;  // 默认很小
    out_max_bin = -1;

    for (int i = start_bin; i < start_bin + len; ++i) 
    {
        double content = hist->GetBinContent(i);
        if (content > out_max_value) 
        {
        out_max_value = content;
        out_max_bin = i;
        }
    }
}


double MaxMinRatio(TH1D* hist) 
{
    int nbins = hist->GetNbinsX();

    double max_val = 0;
    double min_val = 1e4;

    for (int i = 1; i <= nbins; ++i) 
    {
        double content = hist->GetBinContent(i);
        if (content > max_val) max_val = content;
        if (content < min_val && content > 0) min_val = content;
        // cout << content <<  " , min =  " << min_val <<  " ,  max " << max_val <<  endl;
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
