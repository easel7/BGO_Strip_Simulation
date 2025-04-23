void FindMaxPositiveBinSegment(TH1D* hist, double& out_sum, int& out_len) {
    int nbins = hist->GetNbinsX();
    double max_sum = 0;
    int max_len = 0;
    double curr_sum = 0;
    int curr_len = 0;

    for (int i = 1; i <= nbins; ++i) 
    {
        double content = hist->GetBinContent(i);
        if (content > 0) 
        {
            curr_sum += content;
            curr_len++;
            if (curr_sum > max_sum || (curr_sum == max_sum && curr_len > max_len)) 
            {
                max_sum = curr_sum;
                max_len = curr_len;
            }
        } else {
            curr_sum = 0;
            curr_len = 0;
        }
    }
    out_sum = max_sum;
    out_len = max_len;
}

void Draw_Pattern()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    int p_FH_Type;
    int p_FH_Lay;

    const char* string1;
    const char* string2 = "Proton_1000GeV";

    auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/%s.root",string2));
    // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("BarEnergyVector", &p_EnergyVec);
    proton_tree->SetBranchAddress("RMS"            , &p_RMSVec);
    proton_tree->SetBranchAddress("First_Had_Depth", &p_FH_Dep);
    proton_tree->SetBranchAddress("First_Had_Layer", &p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type" , &p_FH_Type);

    auto h_max_min0 = new TH1D("h_max_min0","h_max_min0",50,0,5); 
    auto h_max_min1 = new TH1D("h_max_min1","h_max_min1",50,0,5); // Inelastic
    auto h_max_min2 = new TH1D("h_max_min2","h_max_min2",50,0,5); // Elastic
    auto h_max_min3 = new TH1D("h_max_min3","h_max_min3",50,0,5); // Pass

    auto g_sum_len0 = new TGraph();
    auto g_sum_len1 = new TGraph();
    auto g_sum_len2 = new TGraph();
    auto g_sum_len3 = new TGraph();

    int point_counter   = 0;
    int point_counter_i = 0;
    int point_counter_e = 0;
    int point_counter_p = 0;


    // cout  << proton_tree->GetEntries() << endl;
    // Long64_t entry  = 9;   
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    // for (Long64_t entry = 0; entry < 100; entry++)
    {
        auto c1    = new TCanvas("c1","c1",1400,1400);
        auto hXZ   = new TH2D("hXZ","BGO X-Z Plane",22,-11,11,14,0,14);
        auto hYZ   = new TH2D("hYZ","BGO Y-Z Plane",22,-11,11,14,0,14);
        auto hBGO1 = new TH1D("hBGO1","BGO Core Axis Energy Deposit",14,0,14); 
        auto hBGO2 = new TH1D("hBGO2","Deposit Energy Change Ratio",13,1,14); 
        auto   *line1 = new TLine(-11,14 - p_FH_Lay,11,14 - p_FH_Lay);
        
        auto *g_core0 = new TGraph();
        auto *g_core1 = new TGraph();
        auto    *box0 = new TBox();
        auto    *box1 = new TBox();
        
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[13] = {0};
        double rate_sum      = 0;
        int    rate_len      = 0;

        proton_tree->GetEntry(entry);
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

        // cout << " Depth " << p_FH_Dep << " layer " <<  p_FH_Lay << " Type =" << p_FH_Type << endl;
        for(int layer = 10 ; layer<12 ; layer++)
        {
            if(layer == 10) // %2 == 0
            {
                auto p_start = p_EnergyVec->begin() + layer * 22;  auto p_end = p_EnergyVec->begin() + (layer+1) * 22;
                auto max_iter = max_element(p_start, p_end);
                double max_val = *max_iter;
                int max_index = distance(p_EnergyVec->begin(), max_iter);
                // cout << "Max value in [" << layer * 22 << "," << (layer+1) * 22 << ") = " << max_val  << " at index " << max_index << endl;
                // cout << "Core point at bar = " << max_index % 22 << ", depth = " << 14. - p_FH_Dep / 25.5 << endl;
                bar_info[0] = max_index % 22;
                g_core0->SetPoint(0, bar_info[0] - 10.5, 14 - 10 );
                box0->SetX1(bar_info[0] - 9);
                box0->SetY1(0);
                box0->SetX2(bar_info[0] - 12);
                box0->SetY2(14 - 10 );
            }
            if(layer == 11)
            {
                auto p_start = p_EnergyVec->begin() + layer * 22;  auto p_end = p_EnergyVec->begin() + (layer+1) * 22;
                auto max_iter = max_element(p_start, p_end);
                double max_val = *max_iter;
                int max_index = distance(p_EnergyVec->begin(), max_iter);
                // cout << "Max value in [" << layer * 22 << "," << (layer+1) * 22 << ") = " << max_val  << " at index " << max_index << endl;
                // cout << "Core point at bar = " << max_index % 22 << endl;
                bar_info[1] = max_index % 22;
                g_core1->SetPoint(0, bar_info[1] - 10.5, 14 - 10 );
                box1->SetX1(bar_info[1] - 9);
                box1->SetY1(0);
                box1->SetX2(bar_info[1] - 12);
                box1->SetY2(14 - 10);
            }
        }

        for(int layer = 0 ; layer<14 ; layer++)
        {   
            if (layer %2 == 0)
            {
                for (int k = bar_info[0] - 1 ; k <= bar_info[0] + 1 ; k++) 
                {
                    bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                    // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                }
                if(layer>0) 
                {
                    bar_Change_info[layer-1] = bar_Energy_info[layer]/bar_Energy_info[layer-1];
                }
            }
            else
            {
                for (int k = bar_info[1] - 1 ; k <= bar_info[1] + 1 ; k++) 
                {
                    bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                    // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                }
                bar_Change_info[layer-1] = bar_Energy_info[layer]/bar_Energy_info[layer-1];
            } 
            hBGO1->SetBinContent(layer+1, bar_Energy_info[layer]);
            // cout << layer << " , " << bar_Energy_info[layer] << endl;
            if(layer>0) 
            {
                hBGO2->SetBinContent(layer, log10(bar_Change_info[layer-1]));
            }
            if( bar_Energy_info[layer] == 0) 
            {
                hBGO2->SetBinContent(layer, 0);
            }
        }

        FindMaxPositiveBinSegment(hBGO2,rate_sum,rate_len);
        h_max_min0->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));        g_sum_len0->SetPoint(point_counter++,rate_sum,rate_len);
        if(p_FH_Type == 1)       {  string1 = "Inelastic";  g_core0->SetMarkerColor(kRed);     g_core1->SetMarkerColor(kRed);     box0->SetLineColor(kRed);     box1->SetLineColor(kRed);       line1->SetLineColor(kRed);      h_max_min1->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));        g_sum_len1->SetPoint(point_counter_i++,rate_sum,rate_len);}
        else if (p_FH_Type == 2) {  string1 = "Elastic";    g_core0->SetMarkerColor(kMagenta); g_core1->SetMarkerColor(kMagenta); box0->SetLineColor(kMagenta); box1->SetLineColor(kMagenta);   line1->SetLineColor(kMagenta);  h_max_min2->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));        g_sum_len2->SetPoint(point_counter_e++,rate_sum,rate_len);}
        else                     {  string1 = "Pass";       g_core0->SetMarkerColor(kGray);    g_core1->SetMarkerColor(kGray);    box0->SetLineColor(kGray);    box1->SetLineColor(kGray);      line1->SetLineColor(kGray);     h_max_min3->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));        g_sum_len3->SetPoint(point_counter_p++,rate_sum,rate_len);}
        
        // g_core0->SetMarkerStyle(29);   box0->SetFillStyle(0);    
        // g_core0->SetMarkerSize(2);     box0->SetLineWidth(2);
        // g_core1->SetMarkerStyle(29);   box1->SetFillStyle(0);    
        // g_core1->SetMarkerSize(2);     box1->SetLineWidth(2);
        // line1->SetLineWidth(2);

        // for (int layer = 0; layer < 14; layer++) {
        //     hXZ->SetBinContent(1, 14-layer, -5);
        //     hXZ->SetBinContent(22, 14-layer, -5);
        //     hYZ->SetBinContent(1, 14-layer, -5);
        //     hYZ->SetBinContent(22, 14-layer, -5);
        // }
        // double minVal = -2;  // 设置最小值
        // double maxVal = 2;  // 设置最大值

        // hXZ->SetMinimum(minVal); 
        // hXZ->SetMaximum(maxVal); 
        // hYZ->SetMinimum(minVal); 
        // hYZ->SetMaximum(maxVal); 

        // c1->Clear();
        // TLatex latex;
        // latex.SetTextSize(0.03);
        // latex.SetTextColor(kBlack);
        // c1->Divide(2,2);
        // c1->cd(1);
        // gStyle->SetOptStat(0);
        // hXZ->Draw("colz");
        // line1->Draw("same");
        // g_core0->Draw("psame"); 
        // box0->Draw("psame");

        
        // gPad->Update();
        // auto *palette1 = (TPaletteAxis*) hXZ->GetListOfFunctions()->FindObject("palette");
        // if(palette1)
        // {
        //     palette1->SetX1NDC(0.90);
        //     palette1->SetX2NDC(0.93);
        //     gPad->Modified();
        //     gPad->Update();
        // }
        // latex.DrawLatexNDC(0.82, 0.92, "log_{10}(E/GeV)");


        // c1->cd(2);
        // gStyle->SetOptStat(0);
        // hYZ->Draw("colz");
        // line1->Draw("same");
        // g_core1->Draw("psame"); 
        // box1->Draw("psame");

        // gPad->Update();
        // auto  *palette2 = (TPaletteAxis*) hYZ->GetListOfFunctions()->FindObject("palette");
        // if(palette2)
        // {
        //     palette2->SetX1NDC(0.90);
        //     palette2->SetX2NDC(0.93);
        //     gPad->Modified();
        //     gPad->Update();
        // }
        // latex.DrawLatexNDC(0.82, 0.92, "log_{10}(E/GeV)");


        // c1->cd(3);
        // gPad->SetLogy();
        // hBGO1->SetLineColor(kBlack);
        // hBGO1->SetLineWidth(2);
        // hBGO1->SetTitle("BGO Core Axis Energy Deposit;BGO Layer; 3 Bars Energy / GeV");
        // hBGO1->GetYaxis()->CenterTitle();
        // hBGO1->GetYaxis()->SetTitleOffset(0.97);
        // hBGO1->GetXaxis()->CenterTitle();
        // hBGO1->Draw("hist");
        // auto *line2 = new TLine(p_FH_Lay,hBGO1->GetMinimum(),p_FH_Lay,hBGO1->GetMaximum());
        // line2->SetLineColor(kRed);
        // line2->SetLineWidth(2);
        // if(p_FH_Lay%2 != -1 ) {line2->Draw("");}

        // c1->cd(4);
        // hBGO2->SetLineColor(kBlack);
        // hBGO2->SetLineWidth(2);
        // hBGO2->SetTitle("BGO Core Axis Energy Deposit Change Rate;BGO Layer; log10(^{}E^{}_{i}/^{}E^{}_{i-1})");
        // hBGO2->GetYaxis()->CenterTitle();
        // hBGO2->GetYaxis()->SetRangeUser(-1.5,1.5);
        // hBGO2->GetYaxis()->SetTitleOffset(1.05);
        // hBGO2->GetXaxis()->CenterTitle();
        // hBGO2->Draw("hist");
        // auto *line3 = new TLine(p_FH_Lay,hBGO2->GetMinimum(),p_FH_Lay,hBGO2->GetMaximum());
        // line3->SetLineColor(kRed);
        // line3->SetLineWidth(2);
        // if(p_FH_Lay%2 != -1 ) {line3->Draw("");}
        // c1->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/Figures/%s/%s/%lld.pdf",string1,string2,entry));
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
    g_sum_len1->GetXaxis()->SetLimits(0, 5);  // 根据你的数据范围改
    g_sum_len1->GetYaxis()->SetRangeUser(0, 10); 
    g_sum_len1->SetTitle("; #sum log10(Change Rate);# Continues Positive Bins");
    g_sum_len1->SetMarkerStyle(20);  
    g_sum_len1->SetMarkerColorAlpha(kRed, 0.1);  // 0.0 = fully transparent, 1.0 = fully opaque
    g_sum_len1->SetMarkerSize(0.8);

    g_sum_len2->SetMarkerStyle(21);  
    g_sum_len2->SetMarkerColorAlpha(kBlue, 0.1);
    g_sum_len2->SetMarkerSize(0.8);

    g_sum_len3->SetMarkerStyle(22);  
    g_sum_len3->SetMarkerColorAlpha(kOrange-3, 0.1);
    g_sum_len3->SetMarkerSize(0.8);

    g_sum_len1->Draw("AP");
    g_sum_len2->Draw("PSAME");
    g_sum_len3->Draw("PSAME");

    std::cout << "Number of valid points: " << g_sum_len0->GetN() << std::endl;


}
