void Lateral_RMS2()
{
    int p_First_Had_Layer;  double p_Total_E;  std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_First_Had_Layer;  double d_Total_E;  std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr; double d_weight;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);
    proton_tree->SetBranchAddress("weight"           ,&p_weight);


    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);

    
    double Energy[15] = {0};
    double Energy_LL[15] = {0};
    double Energy_UL[15] = {0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[15][14];  TH1D *hC_p[15][14];       TH2D *h2_p[15];      TGraph *gr_proton[15];          int point_counter_p = 0;
    TH1D *h1_d[15][14];  TH1D *hC_d[15][14];       TH2D *h2_d[15];      TGraph *gr_deuteron[15];        int point_counter_d = 0;



    for(int i =0 ; i<15 ; i++)
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;
        gr_proton[i] = new TGraph();
        gr_deuteron[i] = new TGraph();
        h2_p[i] = new TH2D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),14,0,14,75,0,150);
        h2_d[i] = new TH2D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),14,0,14,75,0,150);
        for( int j= 0; j<14 ;j++)
        {
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),75,0,150); 
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),75,0,150); 
            hC_p[i][j] = new TH1D(Form("hC_p[%d][%d]",i,j), Form("hC_p[%d][%d]",i,j),75,0,150);  
            hC_d[i][j] = new TH1D(Form("hC_d[%d][%d]",i,j), Form("hC_d[%d][%d]",i,j),75,0,150);  
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
        }
    }
    
    // proton_file->cd();
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 10; ++entry)
    {
        proton_tree->GetEntry(entry);   
        double sum_p = 0;
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        // if( (*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 && (*p_RMSVec)[2]<45 && (*p_RMSVec)[3]<45 )  // 
        if ((*p_EnergyVec)[0] > 0.23 && (*p_EnergyVec)[1] > 0.23 && (*p_EnergyVec)[2] > 0.23 && (*p_EnergyVec)[0] > 0.046)
        {
            for (int k = 0; k < p_EnergyVec->size(); k ++)
            {
                // h1_p[p_energy_index][k]->Fill((*p_RMSVec)[k]);
                h1_p[p_energy_index][k]->Fill((*p_RMSVec)[k],p_weight);
            }
        }
    }

    deuteron_file->cd();
    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        deuteron_tree->GetEntry(entry); 
        double sum_d = 0;
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 14) continue;
        // cout << d_Total_E << "," <<  d_energy_index << endl;
        // cout << d_EnergyVec->size() << endl;
        // if((*d_RMSVec)[0]>15 && (*d_RMSVec)[1]>15 && (*d_RMSVec)[2]<45 && (*d_RMSVec)[3]<45) // 
        if ((*d_EnergyVec)[0] > 0.23 && (*d_EnergyVec)[1] > 0.23 && (*d_EnergyVec)[2] > 0.23 && (*d_EnergyVec)[0] > 0.046)
        {
            for (int k = 0; k < d_EnergyVec->size(); k ++)
            {
                int index = int(k / 22);  // Get the Layer
                // h1_d[d_energy_index][k]->Fill((*d_RMSVec)[k]);
                h1_d[d_energy_index][k]->Fill((*d_RMSVec)[k],d_weight);
            }
            // cout << endl;
        }
        // cout << "Next " << endl;
    }

    for(int i = 0 ; i<15 ; i++) // Energy
    {
        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
    
        double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};     double SUM_Proton_Ratio[14]={0};     double SUM_Deuteron_Ratio[14]={0};    
        double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0};  double SUM_Proton_Ratio_LL[14]={0};  double SUM_Deuteron_Ratio_LL[14]={0}; 
        double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0};  double SUM_Proton_Ratio_UL[14]={0};  double SUM_Deuteron_Ratio_UL[14]={0}; 
        for( int j= 0; j<14 ;j++) // layer
        {
            h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
            h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
        
            
            c1->cd(j+1);
            h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral()); 
            h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral()); 

            h1_p[i][j]->GetYaxis()->SetRangeUser(0,0.2);
            h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] in L%d;log10(Rm) = log10(Max Energy Deposit bar/ Total Deposit);Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j));
            h1_p[i][j]->Draw();
            h1_p[i][j]->Draw("hist");
            h1_d[i][j]->Draw("histsame");

            double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
            double p_values[3];      h1_p[i][j]->GetQuantiles(3, p_values, quantiles);      Proton_Ratio[j]       = p_values[1];     Proton_Ratio_LL[j] = p_values[1] - p_values[0];   Proton_Ratio_UL[j]   = p_values[2] - p_values[1];
            double d_values[3];      h1_d[i][j]->GetQuantiles(3, d_values, quantiles);      Deuteron_Ratio[j]     = d_values[1];   Deuteron_Ratio_LL[j] = d_values[1] - d_values[0];   Deuteron_Ratio_UL[j] = d_values[2] - d_values[1];
        }
        auto tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));

        c1->cd(15);
        tex->SetNDC();
        tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[i][0], "Proton", "l");
        legend1->AddEntry(h1_d[i][0], "Deuteron", "l");
        legend1->Draw();       
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Lateral_PowerE/PDF/Lateral_RMS_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));


        auto c3 = new TCanvas("c3","c3",1000,1000);
        auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_Ratio  ,Layer_Err,Layer_Err,Proton_Ratio_LL  ,Proton_Ratio_UL);
        auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_Ratio,Layer_Err,Layer_Err,Deuteron_Ratio_LL,Deuteron_Ratio_UL);
        gre_p->SetMarkerStyle(20);gre_p->SetMarkerColor(kRed);     gre_p->SetLineColor(kRed);     gre_p->SetLineWidth(2);
        gre_d->SetMarkerStyle(21);gre_d->SetMarkerColor(kBlue);    gre_d->SetLineColor(kBlue);    gre_d->SetLineWidth(2);  
        gre_d->SetTitle(Form("Deposit Energy [%.2fGeV, %.2fGeV]; BGO Layer; RMS (mm)",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));

        gre_d->GetYaxis()->SetRangeUser(0,150);
        gre_d->GetXaxis()->SetLimits(0,14);
        gre_d->Draw("AP");
        gre_p->Draw("PSAME");

        auto legend3 = new TLegend(0.70, 0.80, 0.88, 0.88);
        legend3->SetNColumns(2);
        legend3->AddEntry(gre_p, "Proton", "ep");
        legend3->AddEntry(gre_d, "Deuteron", "ep");
        legend3->Draw();
        c3->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Lateral_PowerE/PDF/Lateral_RMS_BGOLayer_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

    }
}
