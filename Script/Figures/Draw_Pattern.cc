#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"

void Draw_Pattern()
{
    vector<double>* p_EnergyVec = nullptr;
    vector<double>* p_RMSVec = nullptr;
    double p_FH_Dep;
    int p_FH_Type;
    int p_FH_Lay;
    int p_Nhits;

    double p_FI_Dep;
    int p_FI_Lay;

    const char* string1;
    const char* string2 = "Proton_1000GeV";

    auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/%s.root",string2));
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("BarEnergyVector", &p_EnergyVec);
    proton_tree->SetBranchAddress("RMS"            , &p_RMSVec);
    proton_tree->SetBranchAddress("First_Had_Depth", &p_FH_Dep);
    proton_tree->SetBranchAddress("First_Had_Layer", &p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type" , &p_FH_Type);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);
    proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
    proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);

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
    Long64_t entry  = 162;   
    // for (Long64_t entry = 0; entry < 100; entry++)
    {
        auto c1    = new TCanvas("c1","c1",2100,1400);
        auto hXZ   = new TH2D("hXZ","BGO X-Z Plane",22,-11,11,14,0,14);
        auto hYZ   = new TH2D("hYZ","BGO Y-Z Plane",22,-11,11,14,0,14);
        auto hBGO1 = new TH1D("hBGO1","BGO Core Axis Energy Deposit",14,0,14); 
        auto hBGO2 = new TH1D("hBGO2","Deposit Energy Change Ratio",14,0,14); 
        auto hBGO3 = new TH1D("hBGO3","Accumulated Deposit Energy",14,0,14); 
        auto  hResidual = new TH1D("hResidual","Accumulated Deposit Energy - Fit",14,0,14); 

        auto *g_core0 = new TGraph();          auto *g_ch_0 = new TGraph();            auto *g_vert0 = new TGraph();
        auto *g_core1 = new TGraph();          auto *g_ch_1 = new TGraph();            auto *g_vert1 = new TGraph();
        auto *g_core2 = new TGraph();          auto *g_ch_2 = new TGraph();            auto *g_vert2 = new TGraph();
        auto *g_core3 = new TGraph();          auto *g_ch_3 = new TGraph();            auto *g_vert3 = new TGraph();

        auto    *box0 = new TBox();
        auto    *box1 = new TBox();
        
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        double bar_Change_info[14] = {0};
        double bar_Accumu_info[14] = {0};
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
        const double RMS_threshold = 15.0;  
        bool bar_info_assigned = false;  
        // cout << " Layer " << layer_start << " , RMS = " <<  (*p_RMSVec)[layer_start] << " , Next Layer RMS = " << (*p_RMSVec)[layer_start+1] << endl;
        for (int k = layer_start; k <= 12; k ++) // two layer as a search window
        {  
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
                break;  
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
            cout << "RMS " << layer << " , " << (*p_RMSVec)[layer] << endl;
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
            bar_Accumu_info[layer] += bar_Accumu_info[layer-1] + bar_Energy_info[layer];
            hBGO3->SetBinContent(layer+1, bar_Accumu_info[layer]);
            hBGO3->SetBinError(layer+1, 0.3 * bar_Accumu_info[layer]);
            cout << "Info " << bar_Accumu_info[layer] << " ± " << 0.3 * bar_Accumu_info[layer] << endl;

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

        double vert_layer = p_FI_Dep / 25.5; 
        cout << "vert_layer = " << vert_layer << endl;
        g_vert0->SetPoint(0, bar_info[0] - 10.5, 15 - vert_layer);
        g_vert1->SetPoint(0, bar_info[1] - 10.5, 15 - vert_layer);   
        g_vert2->SetPoint(0, vert_layer, hBGO1->GetBinContent(ceil(vert_layer)) );
        g_vert3->SetPoint(0, vert_layer, hBGO2->GetBinContent(ceil(vert_layer)) );

        h_max_min0->Fill(log10(MaxMinRatio(hBGO1)));
        g_sum_len0->SetPoint(point_counter++,rate_sum,rate_len);
        if(p_FH_Type == 1)       {  h_max_min1->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));  g_sum_len1->SetPoint(point_counter_i++,rate_sum,rate_len+0.1);}
        else if (p_FH_Type == 2) {  h_max_min2->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));  g_sum_len2->SetPoint(point_counter_e++,rate_sum,rate_len);    }
        else                     {  h_max_min3->Fill(log10(hBGO1->GetMaximum()/hBGO1->GetMinimum()));  g_sum_len3->SetPoint(point_counter_p++,rate_sum,rate_len-0.1);}

        if(p_FH_Type == 1)       {  string1 = "Inelastic";  g_core0->SetMarkerColor(kRed);    g_core1->SetMarkerColor(kRed);    g_core2->SetMarkerColor(kRed);    g_core3->SetMarkerColor(kRed);      g_ch_0->SetMarkerColor(kRed);      g_ch_1->SetMarkerColor(kRed);      g_ch_2->SetMarkerColor(kRed);      g_ch_3->SetMarkerColor(kRed);     box0->SetLineColor(kRed);     box1->SetLineColor(kRed);       line1->SetLineColor(kRed);     line2->SetLineColor(kRed);     line3->SetLineColor(kRed);     }
        else if (p_FH_Type == 2) {  string1 = "Elastic";    g_core0->SetMarkerColor(kMagenta);g_core1->SetMarkerColor(kMagenta);g_core2->SetMarkerColor(kMagenta);g_core3->SetMarkerColor(kMagenta);  g_ch_0->SetMarkerColor(kMagenta);  g_ch_1->SetMarkerColor(kMagenta);  g_ch_2->SetMarkerColor(kMagenta);  g_ch_3->SetMarkerColor(kMagenta); box0->SetLineColor(kMagenta); box1->SetLineColor(kMagenta);   line1->SetLineColor(kMagenta); line2->SetLineColor(kMagenta); line3->SetLineColor(kMagenta); }
        else                     {  string1 = "Pass";       g_core0->SetMarkerColor(kGray);   g_core1->SetMarkerColor(kGray);   g_core2->SetMarkerColor(kGray);   g_core3->SetMarkerColor(kGray);     g_ch_0->SetMarkerColor(kGray);     g_ch_1->SetMarkerColor(kGray);     g_ch_2->SetMarkerColor(kGray);     g_ch_3->SetMarkerColor(kGray);    box0->SetLineColor(kGray);    box1->SetLineColor(kGray);      line1->SetLineColor(kGray);    line2->SetLineColor(kGray);    line3->SetLineColor(kGray);    }
        
        box0->SetFillStyle(0);        box0->SetLineWidth(2);
        box1->SetFillStyle(0);        box1->SetLineWidth(2);
        g_core0->SetMarkerStyle(29);  g_core0->SetMarkerSize(2);    g_ch_0->SetMarkerStyle(22);   g_ch_0->SetMarkerSize(2); g_vert0->SetMarkerStyle(3);   g_vert0->SetMarkerColor(3);
        g_core1->SetMarkerStyle(29);  g_core1->SetMarkerSize(2);    g_ch_1->SetMarkerStyle(22);   g_ch_1->SetMarkerSize(2); g_vert1->SetMarkerStyle(3);   g_vert1->SetMarkerColor(3);
        g_core2->SetMarkerStyle(29);  g_core2->SetMarkerSize(2);    g_ch_2->SetMarkerStyle(22);   g_ch_2->SetMarkerSize(2); g_vert2->SetMarkerStyle(3);   g_vert2->SetMarkerColor(3);
        g_core3->SetMarkerStyle(29);  g_core3->SetMarkerSize(2);    g_ch_3->SetMarkerStyle(22);   g_ch_3->SetMarkerSize(2); g_vert3->SetMarkerStyle(3);   g_vert3->SetMarkerColor(3);

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
        c1->Divide(3,2);
        c1->cd(1);
        gStyle->SetOptStat(0);
        hXZ->Draw("colz");
        line1->Draw("same");
        g_core0->Draw("psame"); 
        g_ch_0->Draw("psame"); 
        g_vert0->Draw("psame"); 
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
        g_vert1->Draw("psame"); 

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


        c1->cd(4);
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
        g_vert2->Draw("psame"); 
        if(p_FH_Lay%2 != -1 ) {line2->Draw("");}

        c1->cd(5);
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
        g_vert3->Draw("psame"); 
        if(p_FH_Lay%2 != -1 ) {line3->Draw("");}


        c1->cd(6);
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->SetGrid(1, 1);      // Vertical ,Horizontal grid
        pad1->SetLogy(1);
        pad1->SetLogx(0);
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        hBGO3->SetLineColor(kBlack);
        hBGO3->SetLineWidth(2);
        hBGO3->SetTitle("Accumulated Energy Deposit;BGO Layer; Energy / GeV");
        hBGO3->GetYaxis()->CenterTitle();
        hBGO3->GetYaxis()->SetTitleOffset(0.97);
        hBGO3->GetXaxis()->CenterTitle();
        hBGO3->Draw("hist");
        
        TF1 *sigmoid = new TF1("sigmoid", "[0]+ 0.02*x + ([1]-[0] - 0.02*x )/(1 + exp(-(x-[2])/[3]))", 0, 14);
        cout << max_rate_index << endl;
        sigmoid->SetParameters(hBGO3->GetMinimum(), hBGO3->GetMaximum(), max_rate_index - 1, 1); 
        sigmoid->FixParameter(0,hBGO3->GetMinimum()); 
        sigmoid->FixParameter(1,hBGO3->GetMaximum()); 
        hBGO3->Fit(sigmoid, "R");
        sigmoid->Draw("same");
        double percentile = (sigmoid->Eval(p_FH_Lay) - hBGO3->GetMinimum() - 0.02 * p_FH_Lay) / (hBGO3->GetMaximum() - hBGO3->GetMinimum() - 0.02 * p_FH_Lay);
        cout << "Percentile = " << percentile << endl;

        double percentile2 = Mod_Sigmoid_Percentile(p_FI_Lay,sigmoid->GetParameter(2),sigmoid->GetParameter(3));
        cout << "perenctile2: " << percentile2 << endl; 
        latex.DrawLatexNDC(0.12, 0.82, Form("Ine Vertex Sigmoid Percentile:%.2f %%", percentile*100));

        c1->cd(6);
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.2);
        pad2->SetGrid(1, 1);      // Vertical ,Horizontal grid
        pad2->SetLogy(0);
        pad2->SetLogx(0);
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad
        hResidual->SetTitle(";BGO Layer; (Data - Fit)/Err");
        hResidual->GetXaxis()->SetLabelSize(0.06);
        hResidual->GetXaxis()->SetTitleSize(0.1);
        hResidual->GetYaxis()->SetLabelSize(0.06);
        hResidual->GetYaxis()->SetTitleOffset(0.5);
        hResidual->GetYaxis()->SetTitleSize(0.07);

        for (int i = 1; i <= 14; ++i) {
            double x =     hBGO3->GetBinCenter(i);
            double yData = hBGO3->GetBinContent(i);
            double yFit  = sigmoid->Eval(x);
            double yErr = hBGO3->GetBinError(i);
            if (yErr > 0) hResidual->SetBinContent(i, (yData - yFit) / yErr);
        }
        // hResidual->GetYaxis()->SetRangeUser(-1,1);
        hResidual->SetLineColor(kBlack);
        hResidual->SetLineWidth(2);
        hResidual->Draw("hist");
        auto line0 = new TLine(0,0,14,0);
        line0->SetLineColor(kRed);
        line0->SetLineWidth(2);
        line0->SetLineStyle(2);
        line0->Draw("same");


        if (entry < 100) { c1->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/Figures/%s/%s/%lld.pdf",string1,string2,entry)); }
    }
}

