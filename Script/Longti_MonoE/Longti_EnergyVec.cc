void Longti_EnergyVec()
{
    double Energy[19]={0};
    for (int i = 0; i < 1; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        std::vector<double>* e_EnergyVec = nullptr;    double e_Total_E; double e_Energy;
        auto c1 = new TCanvas("c1","c1",2500,1500);
        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
        electron_tree->SetBranchAddress("Total_E",&e_Total_E);
        electron_tree->SetBranchAddress("Energy",&e_Energy);
        TH1D *h1_e[14];

        for (int j = 0; j < 14; j++)
        {
            h1_e[j] = new TH1D(Form("h1_e[%d]",j),Form("h1_e[%d]",j),100,0,10);
        }

        for (Long64_t entry = 0; entry < electron_tree->GetEntries(); ++entry)
        {
            electron_tree->GetEntry(entry);
            {
                for (int j = 0; j < 14; j++) // Sub Fig
                {
                    h1_e[j]->Fill((*e_EnergyVec)[j]);
                }
            }
        }
        c1->Clear();
        c1->Divide(5,3);
        // gStyle->SetOptStat(0);
        for (int j = 0; j < 15; j++) // Sub Fig
        {
            c1->cd(j+1);
            if(j<=13)
            {
                h1_e[j]->Sumw2(); h1_e[j]->SetLineColor(kRed);     h1_e[j]->SetMarkerColor(kOrange-3);     h1_e[j]->SetLineWidth(2); // h1_e[j]->Scale(1.0/h1_e[j]->Integral());
                h1_e[j]->Draw();
            }
            else // (j==14)    
            {
                TLatex *tex = new TLatex(0.1,0.9,Form("IncidentParicle %dGeV",int(Energy[i])));tex->SetNDC();tex->Draw(); 
                auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
                legend1->AddEntry(h1_e[0], "Electrπ", "l");
                legend1->Draw();       
            }
        }
    }
    
}