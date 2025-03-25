void Edep_TrackLength_Fit()
{
    std::vector<double>* EnergyVec = nullptr;
    std::vector<double>* LengthVec = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_100GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("LayerEnergyVector", &EnergyVec);
    proton_tree->SetBranchAddress("LayerLengthVector", &LengthVec);

    int Counter = 0;
    auto c1 = new TCanvas("c1","c1",1200,1200);
    auto L0 = new TH2D("L0","L0",50,0,25,50,0,25);
    auto L1 = new TProfile("L1","L1",50,0,25,0,25,"s");

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);
        if ((*EnergyVec)[0] > 0.23 && (*EnergyVec)[1] > 0.23 && (*EnergyVec)[2] > 0.23 && (*EnergyVec)[0] > 0.046)
        {
            L0->Fill((*LengthVec)[3],(*EnergyVec)[3]);
            L1->Fill((*LengthVec)[3],(*EnergyVec)[3]);
            Counter++;
        } // HET
    }
    cout << Counter << endl;
    
    c1->Clear();

    gStyle->SetOptStat(1111);
    L0->SetTitle("100 GeV Proton Edep vs TrackLength in L6;TrackLength(m);Energy(GeV)");
    L0->Draw("colz");
    L1->SetLineColor(kRed);
    L1->SetLineWidth(2);
    L1->Draw("profsame");
    auto fitFunc = new TF1("fitFunc", "pol1", 5, 15);
    L1->Fit(fitFunc, "", "", 5, 15);

    // Get slope parameter from fit
    double slope = fitFunc->GetParameter(1);
    double slope_err = fitFunc->GetParError(1);

    // Add slope text to the canvas
    TLatex *tex = new TLatex();
    tex->SetNDC();  // Use normalized coordinates (0-1)
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.3, 0.7, Form("Fitted Slope = %.2f #pm %.2f", slope,slope_err));


}