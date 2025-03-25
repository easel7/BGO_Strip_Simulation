void Longti_RMS_Single()
{
    auto c1 = new TCanvas("c1","c1",1200,900);
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_100GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    std::vector<double>* RMSVec = nullptr;
    std::vector<double>* EnergyVec = nullptr;
    std::vector<double>* Efrac = nullptr;

    int First_Had_Layer;
    proton_tree->SetBranchAddress("RMS", &RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector", &EnergyVec);
    proton_tree->SetBranchAddress("Efrac", &Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer", &First_Had_Layer);

    cout  << proton_tree->GetEntries() << endl;
    auto h1_p = new TH1D("h1_p","h1_p",50,-5,0);  

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);
        if ((*EnergyVec)[0] > 0.23 && (*EnergyVec)[1] > 0.23 && (*EnergyVec)[2] > 0.23 && (*EnergyVec)[0] > 0.046) // HET
        // if ((*EnergyVec)[0] > 0.0092 && (*EnergyVec)[1] > 0.0092 ) // UBT
        {
            h1_p->Fill(log10((*Efrac)[0]));
        }
    }
    h1_p->Sumw2(); 
    h1_p->Scale(1.0/h1_p->Integral());
    h1_p->SetLineColor(kRed);     h1_p->SetMarkerColor(kRed);     h1_p->SetLineWidth(2);h1_p->GetYaxis()->SetRangeUser(0,0.25);h1_p->SetTitle("100 GeV incident EdepRatio Distrubution in L0;og_{10}(Energy Deposit in Layer/ Total Deposit);Normalized Count");
    h1_p->Draw();
    double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
    double p_values[3];  
    h1_p->GetQuantiles(3, p_values, quantiles);

    TLine *l_p[3];
    for (int ii = 0 ;ii< 3 ; ii++)
    {
        l_p[ii] = new TLine(p_values[ii],0,p_values[ii],h1_p->GetBinContent(h1_p->FindBin(p_values[ii])));l_p[ii]->SetLineColor(kRed);     l_p[ii]->SetLineWidth(2);l_p[ii]->Draw();
    }
    auto legend1 = new TLegend(0.12, 0.68, 0.28, 0.88);
    legend1->AddEntry(h1_p, "Proton", "el");
    legend1->Draw();       
    // Proton_Ratio = p_values[1];     Proton_Ratio_LL  = p_values[1] - p_values[0];   Proton_Ratio_UL= p_values[2] - p_values[1];      
    // cout << Proton_Ratio  << " , " << Proton_Ratio_LL<< " , " << Proton_Ratio_UL << endl; 
    // cout << Layer  << " , " << Layer_Err << " , " << Layer_Err << endl;         
}