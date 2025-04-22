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

    // cout  << proton_tree->GetEntries() << endl;
    // Long64_t entry  = 484;   
    // for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    for (Long64_t entry = 0; entry < 1000; entry++)
    {
        auto c1    = new TCanvas("c1","c1",1800,600);
        auto hXZ   = new TH2D("hXZ","BGO X-Z Plane",22,-11,11,14,0,14);
        auto hYZ   = new TH2D("hYZ","BGO Y-Z Plane",22,-11,11,14,0,14);
        auto hBGO  = new TH1D("hBGO","BGO Core Axis Energy Deposit",14,0,14); 
        
        auto *g_core1 = new TGraph();
        auto *g_core2 = new TGraph();
        auto    *box1 = new TBox();
        auto    *box2 = new TBox();
        double bar_info[2] = {0};
        double bar_Energy_info[14] = {0};
        proton_tree->GetEntry(entry);
        // cout << " RMS at 0 layer = " << (*p_RMSVec)[0] << endl;

        for(size_t i = 0; i < p_EnergyVec->size(); i++) 
        {
            int layer = i / 22; 
            int bar = i % 22;
            double energy = 0;
            if ((*p_EnergyVec)[i] < 1e-2)  {energy = -5;              }
            else                         {energy = log10((*p_EnergyVec)[i]);}

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

        for(int layer = 0 ; layer<14 ; layer++)
        {
            if (p_FH_Lay>-1)
            {
                if(layer == p_FH_Lay)
                {
                    auto p_start = p_EnergyVec->begin() + layer * 22;  auto p_end = p_EnergyVec->begin() + (layer+1) * 22;
                    auto max_iter = max_element(p_start, p_end);
                    double max_val = *max_iter;
                    int max_index = distance(p_EnergyVec->begin(), max_iter);
                    // cout << "Max value in [" << layer * 22 << "," << (layer+1) * 22 << ") = " << max_val  << " at index " << max_index << " Depth " << p_FH_Dep << " layer " <<  p_FH_Lay << " Type =" << p_FH_Type << endl;
                    // cout << "Core point at bar = " << max_index % 22 << ", depth = " << 14. - p_FH_Dep / 25.5 << endl;
                    bar_info[0] = max_index % 22;
                    g_core1->SetPoint(0, bar_info[0] - 10.5, 14. - p_FH_Dep / 25.5);
                    box1->SetX1(bar_info[0] - 9);
                    box1->SetY1(0);
                    box1->SetX2(bar_info[0] - 12);
                    box1->SetY2(14 - p_FH_Lay);
                }
                if(layer == p_FH_Lay + 1)
                {
                    auto p_start = p_EnergyVec->begin() + layer * 22;  auto p_end = p_EnergyVec->begin() + (layer+1) * 22;
                    auto max_iter = max_element(p_start, p_end);
                    double max_val = *max_iter;
                    int max_index = distance(p_EnergyVec->begin(), max_iter);
                    // cout << "Max value in [" << layer * 22 << "," << (layer+1) * 22 << ") = " << max_val  << " at index " << max_index << " Depth " << endl;
                    // cout << "Core point at bar = " << max_index % 22 << endl;
                    bar_info[1] = max_index % 22;
                    g_core2->SetPoint(0, bar_info[1] - 10.5, 14. - p_FH_Dep / 25.5);
                    box2->SetX1(bar_info[1] - 9);
                    box2->SetY1(0);
                    box2->SetX2(bar_info[1] - 12);
                    box2->SetY2(14 - p_FH_Lay);
                }
            }
            else // Pass
            {
                auto p_end0 = p_EnergyVec->begin();  
                auto p_end1 = p_EnergyVec->begin() + 22;
                auto p_end2 = p_EnergyVec->begin() + 44;
                auto max_iter1 = max_element(p_end0, p_end1);
                double max_val1 = *max_iter1;
                auto max_iter2 = max_element(p_end1, p_end2);
                double max_val2 = *max_iter2;
                int max_index1 = distance(p_EnergyVec->begin(), max_iter1);
                int max_index2 = distance(p_EnergyVec->begin(), max_iter2);
                bar_info[0] = max_index1 % 22;
                bar_info[1] = max_index2 % 22;
                // cout << "Max value in [" << layer * 22 << "," << (layer+1) * 22 << ") = " << max_val  << " at index " << max_index << " Depth " << endl;
                // cout << "Core point at bar = " << max_index % 22 << endl;
            }
           
        }

        for(int layer = 0 ; layer<14 ; layer++)
        {
            if(p_FH_Lay%2 == 0 ) 
            {
                if (layer %2 == 0)
                {
                    for (int k = bar_info[0] - 1 ; k <= bar_info[0] + 1 ; k++) 
                    {
                        bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                        // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                    }
                    // cout << endl;
                }
                else
                {
                    for (int k = bar_info[1] - 1 ; k <= bar_info[1] + 1 ; k++) 
                    {
                        bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                        // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                    }
                    // cout << endl;
                } 
            }
            if (p_FH_Lay%2 == 1) 
            {
                if (layer %2 == 0)
                {
                    for (int k = bar_info[1] - 1 ; k <= bar_info[1] + 1 ; k++) 
                    {
                        bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                        // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                    }
                    // cout << endl;
                }
                else
                {
                    for (int k = bar_info[0] - 1 ; k <= bar_info[0] + 1 ; k++) 
                    {
                        bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                        // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                    }
                    // cout << endl;
                }
            }
            if (p_FH_Lay%2 == -1)// Pass
            {
                if (layer %2 == 0)
                {
                    for (int k = bar_info[0] - 1 ; k <= bar_info[0] + 1 ; k++) 
                    {
                        bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                        // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                    }
                    // cout << endl;
                }
                else
                {
                    for (int k = bar_info[1] - 1 ; k <= bar_info[1] + 1 ; k++) 
                    {
                        bar_Energy_info[layer] += (*p_EnergyVec)[layer * 22 + k];
                        // cout << "Layer " << layer << " , bar " << k << " , Energy " << (*p_EnergyVec)[layer * 22 + k] << " , ";
                    }
                    // cout << endl;
                } 
            }
            hBGO->SetBinContent(layer+1, bar_Energy_info[layer]);
        }

        if(p_FH_Type == 1)       {  string1 = "Inelastic";  g_core1->SetMarkerColor(kRed);     g_core2->SetMarkerColor(kRed);     box1->SetLineColor(kRed);     box2->SetLineColor(kRed);}
        else if (p_FH_Type == 2) {  string1 = "Elastic";  g_core1->SetMarkerColor(kMagenta); g_core2->SetMarkerColor(kMagenta); box1->SetLineColor(kMagenta); box2->SetLineColor(kMagenta);}
        else                     {  string1 = "Pass";  g_core1->SetMarkerColor(kGray);    g_core2->SetMarkerColor(kGray);    box1->SetLineColor(kGray);    box2->SetLineColor(kGray);}
        
        g_core1->SetMarkerStyle(29);         box1->SetFillStyle(0);    
        g_core1->SetMarkerSize(2);           box1->SetLineWidth(2);
        g_core2->SetMarkerStyle(29);         box2->SetFillStyle(0);    
        g_core2->SetMarkerSize(2);           box2->SetLineWidth(2);


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
        c1->Divide(3,1);
        c1->cd(1);
        gStyle->SetOptStat(0);
        hXZ->Draw("colz");

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

        if(p_FH_Lay%2 ==0 ) 
        { 
            c1->cd(1); g_core1->Draw("psame"); box1->Draw("psame");
            c1->cd(2); g_core2->Draw("psame"); box2->Draw("psame");
            
        }
        if(p_FH_Lay%2 ==1 ) 
        {
            c1->cd(1); g_core2->Draw("psame"); box2->Draw("psame");
            c1->cd(2); g_core1->Draw("psame"); box1->Draw("psame");
        }

        c1->cd(3);
        gPad->SetLogy();
        hBGO->SetLineColor(kBlack);
        hBGO->SetTitle("BGO Core Axis Energy Deposit;BGO Layer; 3 Bars Energy / GeV");
        hBGO->GetYaxis()->CenterTitle();
        hBGO->GetYaxis()->SetTitleOffset(0.97);
        hBGO->GetXaxis()->CenterTitle();
        hBGO->Draw("hist");
        auto *line = new TLine(p_FH_Lay,hBGO->GetMinimum(),p_FH_Lay,hBGO->GetMaximum());
        line->SetLineColor(kRed);
        if(p_FH_Lay%2 != -1 ) {line->Draw("");}
        c1->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/Figures/%s/%s/%lld.pdf",string1,string2,entry));

    }

}
