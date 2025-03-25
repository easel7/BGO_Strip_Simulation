void Lateral_RMS()
{
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    std::vector<double>* RMSVec = nullptr;
    std::vector<double>* EnergyVec = nullptr;
    int First_Had_Layer;
    proton_tree->SetBranchAddress("RMS", &RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector", &EnergyVec);
    proton_tree->SetBranchAddress("First_Had_Layer", &First_Had_Layer);


    cout  << proton_tree->GetEntries() << endl;
    auto c1    = new TCanvas("c1","c1",1200,1200);

    auto h1_p = new TH1D("h1_p","BGO Layer 0 ;RMS;Counts",75,0,150);
    auto h2_p = new TH1D("h2_p","BGO Layer 4 ;RMS;Counts",75,0,150);
    auto h3_p = new TH1D("h3_p","BGO Layer 9 ;RMS;Counts",75,0,150);
    auto h4_p = new TH1D("h4_p","BGO Layer 13;RMS;Counts",75,0,150);

    auto h1_F = new TH1D("h1_F","BGO Layer 0 ;RMS;Counts",75,0,150);
    auto h1_N = new TH1D("h1_N","BGO Layer 0 ;RMS;Counts",75,0,150);


    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);
        if ((*EnergyVec)[0] > 0.23 && (*EnergyVec)[1] > 0.23 && (*EnergyVec)[2] > 0.23 && (*EnergyVec)[0] > 0.046) // )
        {
            h1_p->Fill((*RMSVec)[0]);
            h2_p->Fill((*RMSVec)[4]);
            h3_p->Fill((*RMSVec)[9]);
            h4_p->Fill((*RMSVec)[13]);

            if( First_Had_Layer == 0)
            {
                h1_F->Fill((*RMSVec)[0]);
            }
            else           
            {
                h1_N->Fill((*RMSVec)[0]);
            }
        }
    }   
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    h1_p->Draw();
    h1_F->SetLineColor(kRed);
    h1_N->SetLineColor(kBlue);
    h1_F->Draw("same");
    h1_N->Draw("same");
    c1->cd(2);
    h2_p->Draw();
    c1->cd(3);
    h3_p->Draw();
    c1->cd(4);
    h4_p->Draw();

}