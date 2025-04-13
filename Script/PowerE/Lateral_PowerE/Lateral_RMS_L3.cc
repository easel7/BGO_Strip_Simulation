void Lateral_RMS_L3()
{
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");

    auto proton_tree = (TTree*)proton_file->Get("B4");
    std::vector<double>* RMSVec = nullptr;
    std::vector<double>* EnergyVec = nullptr;
    int First_Had_Layer;
    double Total_E;
    double Weight;
    proton_tree->SetBranchAddress("RMS", &RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector", &EnergyVec);
    proton_tree->SetBranchAddress("First_Had_Layer", &First_Had_Layer);
    proton_tree->SetBranchAddress("Total_E", &Total_E);
    proton_tree->SetBranchAddress("weight", &Weight);

    cout  << proton_tree->GetEntries() << endl;
    auto c1    = new TCanvas("c1","c1",1800,1200);

    auto h1_p = new TH1D("h1_p","BGO Layer 1st ;RMS;Counts",80,0,160);
    auto h2_p = new TH1D("h2_p","BGO Layer 5th ;RMS;Counts",80,0,160);
    auto h3_p = new TH1D("h3_p","BGO Layer 10th ;RMS;Counts",80,0,160);
    auto h4_p = new TH1D("h4_p","BGO Layer 14th ;RMS;Counts",80,0,160);

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);
        if ((*EnergyVec)[0] > 0.23 && (*EnergyVec)[1] > 0.23 && (*EnergyVec)[2] > 0.23 && (*EnergyVec)[0] > 0.046 && Total_E > 1E3 && Total_E<1580) //
        {
            h1_p->Fill((*RMSVec)[0],Weight);
            h2_p->Fill((*RMSVec)[4],Weight);
            h3_p->Fill((*RMSVec)[9],Weight);
            h4_p->Fill((*RMSVec)[13],Weight);
        }
    }   
    h1_p->Sumw2();
    h2_p->Sumw2();
    h3_p->Sumw2();
    h4_p->Sumw2();
    h1_p->Scale(1/h1_p->Integral());
    h1_p->SetLineWidth(2);
    h2_p->Scale(1/h2_p->Integral());
    h2_p->SetLineWidth(2);
    h3_p->Scale(1/h3_p->Integral());
    h3_p->SetLineWidth(2);
    h4_p->Scale(1/h4_p->Integral());
    h4_p->SetLineWidth(2);

    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    // gPad->SetLogy();
    h1_p->Draw("hist");

    c1->cd(2);
    h2_p->Draw("hist");
    
    c1->cd(3);
    h3_p->Draw("hist");

    
    c1->cd(4);
    h4_p->Draw("hist");
  

}