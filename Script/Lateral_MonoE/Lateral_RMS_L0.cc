void Lateral_RMS_L0()
{
    double Energy[19]={0};
    for (int i = 0; i < 1; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        auto c1 = new TCanvas("c1","c1",2500,1500);
        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        std::vector<double>* RMSVec = nullptr;
        std::vector<double>* EnergyVec = nullptr;
        int First_Had_Layer;
        proton_tree->SetBranchAddress("RMS", &RMSVec);
        proton_tree->SetBranchAddress("LayerEnergyVector", &EnergyVec);
        proton_tree->SetBranchAddress("First_Had_Layer", &First_Had_Layer);
        cout  << proton_tree->GetEntries() << endl;
        TH1D *h1_p[14];  TF1  *fitFunc_p[14];

        double Proton_Ratio[14]={0};   
        double Proton_Ratio_LL[14]={0};
        double Proton_Ratio_UL[14]={0};
        double Layer[14]={0};
        double Layer_Err[14]={0};
        for (int j = 0; j < 14; j++) // Sub Fig
        {
            h1_p[j] = new TH1D(Form("h1_p[%d]",j), Form("h1_p[%d]",j),75,0,150); 
        }
        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
        {
            proton_tree->GetEntry(entry);
            if ((*EnergyVec)[0] > 0.23 && (*EnergyVec)[1] > 0.23 && (*EnergyVec)[2] > 0.23 && (*EnergyVec)[0] > 0.046) //
            {
                for (int j = 0; j < 14; j++) // Sub Fig
                {
                    h1_p[j]->Fill((*RMSVec)[j]);
                    double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
                    double p_values[3];  h1_p[j]->GetQuantiles(3, p_values, quantiles);
                    Proton_Ratio[j] = p_values[1];     Proton_Ratio_LL[j]  = p_values[1] - p_values[0];   Proton_Ratio_UL[j]= p_values[2] - p_values[1];                
                    Layer[j] = 0.5 + j;
                    Layer_Err[j] = 0.5;
                    // cout << Proton_Ratio[j]  << " , " << Proton_Ratio_LL[j]<< " , " << Proton_Ratio_UL[j] << endl; 
                    // cout << Layer[j]  << " , " << Layer_Err[j] << " , " << Layer_Err[j] << endl; 
                }
            }
        }   
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);

        for (int j = 0; j < 15; j++) // Sub Fig
        {
            c1->cd(j+1);
            if(j<=13)
            {
                h1_p[j]->Draw();
            }
            else // (j==14)    
            {
                TLatex *tex = new TLatex(0.1,0.9,Form("IncidentParicle %dGeV",int(Energy[i])));tex->SetNDC();tex->Draw(); 
                auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
                legend1->AddEntry(h1_p[0], "Proton", "l");
                legend1->Draw();       
            }
        }

    }

}