void Rm2_Interaction()
{
    int p_First_Had_Layer; int p_First_Had_Type; double p_Total_E;  std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;
    int d_First_Had_Layer; int d_First_Had_Type; double d_Total_E;  std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("BarEnergyVector"  ,&p_EnergyVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("First_Had_Type"  ,&p_First_Had_Type);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_First_Had_Type);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    
    double Energy[15] = {0};
    double Energy_LL[15] = {0};
    double Energy_UL[15] = {0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[15][14];  TH1D *hC_p[15][14];       TH2D *h2_p[15];      TGraph *gr_proton[15];          int point_counter_p = 0;
    TH1D *h1_d[15][14];  TH1D *hC_d[15][14];       TH2D *h2_d[15];      TGraph *gr_deuteron[15];        int point_counter_d = 0;

    TH1D *h1_p_inter[15];     TH1D *hC_p_inter[15];
    TH1D *h1_d_inter[15];     TH1D *hC_d_inter[15];

    for(int i =0 ; i<15 ; i++)
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;
        gr_proton[i] = new TGraph();
        gr_deuteron[i] = new TGraph();
        h2_p[i] = new TH2D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),14,0,14,50,-5,0);
        h2_d[i] = new TH2D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),14,0,14,50,-5,0);

        h1_p_inter[i] =new TH1D(Form("h1_p_inter[%d]",i),Form("h1_p_inter[%d]",i), 50,-5,0);  
        h1_d_inter[i] =new TH1D(Form("h1_d_inter[%d]",i),Form("h1_d_inter[%d]",i), 50,-5,0);  
        hC_p_inter[i] =new TH1D(Form("hC_p_inter[%d]",i),Form("hC_p_inter[%d]",i), 50,-5,0);  
        hC_d_inter[i] =new TH1D(Form("hC_d_inter[%d]",i),Form("hC_d_inter[%d]",i), 50,-5,0);  
        for( int j= 0; j<14 ;j++)
        {
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),50,-5,0); 
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),50,-5,0); 
            hC_p[i][j] = new TH1D(Form("hC_p[%d][%d]",i,j), Form("hC_p[%d][%d]",i,j),50,-5,0); 
            hC_d[i][j] = new TH1D(Form("hC_d[%d][%d]",i,j), Form("hC_d[%d][%d]",i,j),50,-5,0); 
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
            
        }
    }
    
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);   
        double sum_p = 0;
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        if(p_First_Had_Type != 1) continue;
        {
            for (size_t k = 0; k <= p_First_Had_Layer * 22 ; k += 22)
            {
                int index = int(k / 22);  // Get the Layer
                auto p_start = p_EnergyVec->begin() + k;  auto p_end = (k + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  
                double p_maxVal = *std::max_element(p_start, p_end); 
                sum_p += p_maxVal;
                if (index == p_First_Had_Layer)
                {
                    h1_p[p_energy_index][index]->Fill(log10(p_maxVal/p_Total_E)); 
                    hC_p[p_energy_index][index]->Fill(log10(sum_p / p_Total_E));
                    h1_p_inter[p_energy_index]->Fill(log10(p_maxVal/p_Total_E)); 
                    hC_p_inter[p_energy_index]->Fill(log10(sum_p/p_Total_E)); 
                }
            }
        }
    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry); 
        double sum_d = 0;
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 14) continue;
        if(d_First_Had_Type != 1) continue;
        {
            // if (d_energy_index==10) cout << "deuteron entry = " << entry << endl;
            for (size_t k = 0; k <= d_First_Had_Layer * 22; k += 22)
            {
                int index = int(k / 22);  // Get the Layer
                auto d_start = d_EnergyVec->begin() + k;  auto d_end = (k + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end(); 
                double d_maxVal = *std::max_element(d_start, d_end); 
                sum_d += d_maxVal;
                if (index == d_First_Had_Layer)
                {
                    h1_d[d_energy_index][index]->Fill(log10(d_maxVal/d_Total_E)); 
                    hC_d[d_energy_index][index]->Fill(log10(sum_d / d_Total_E));
                    h1_d_inter[d_energy_index]->Fill(log10(d_maxVal/d_Total_E));
                    hC_d_inter[d_energy_index]->Fill(log10(sum_d /d_Total_E)); 
                }
            }
        }
    }

    for(int i = 9 ; i<10 ; i++) // Energy
    {
        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
    
        auto c2 = new TCanvas("c2","c2",2500,1500);
        c2->Clear();
        c2->Divide(5,3);
        gStyle->SetOptStat(0);
        double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};     double SUM_Proton_Ratio[14]={0};     double SUM_Deuteron_Ratio[14]={0};    
        double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0};  double SUM_Proton_Ratio_LL[14]={0};  double SUM_Deuteron_Ratio_LL[14]={0}; 
        double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0};  double SUM_Proton_Ratio_UL[14]={0};  double SUM_Deuteron_Ratio_UL[14]={0}; 
        for( int j= 0; j<14 ;j++) // layer
        {
            h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
            h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
            hC_p[i][j]->SetLineColor(kRed);   hC_p[i][j]->SetMarkerColor(kRed);  hC_p[i][j]->SetLineWidth(2);   hC_p[i][j]->Sumw2(); 
            hC_d[i][j]->SetLineColor(kBlue);  hC_d[i][j]->SetMarkerColor(kBlue); hC_d[i][j]->SetLineWidth(2);   hC_d[i][j]->Sumw2(); 
            
            c1->cd(j+1);
            h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral()); 
            h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral()); 

            h1_p[i][j]->GetYaxis()->SetRangeUser(0,h1_p[i][j]->GetMaximum()*1.2);
            h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] in L%d;log10(Rm) = log10(Max Energy Deposit bar/ Total Deposit);Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j));
            h1_p[i][j]->Draw("hist");
            h1_d[i][j]->Draw("histsame");
        
            c2->cd(j + 1);
            hC_p[i][j]->Scale(1.0/hC_p[i][j]->Integral());
            hC_d[i][j]->Scale(1.0/hC_d[i][j]->Integral());
            hC_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] in L%d;#sum_{0}^{%d} (dE/dx) / Edep",pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j,j));
            hC_p[i][j]->GetYaxis()->SetRangeUser(0,hC_p[i][j]->GetMaximum()*1.2);
            hC_p[i][j]->Draw("hist");
            hC_d[i][j]->Draw("histsame");
        }
        auto tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));

        c1->cd(15);
        tex->SetNDC();
        tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[i][0], "Proton", "l");
        legend1->AddEntry(h1_d[i][0], "Deuteron", "l");
        legend1->Draw();       
        // c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Rm_PowerE/PDF/Rm_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        c2->cd(15);
        tex->Draw();
        legend1->Draw();       
        // c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Rm_PowerE/CDF/Rm_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        auto c3 = new TCanvas("c3","c3",1000,1000);
        c3->cd();
        h1_p_inter[i]->SetLineColor(kRed);   h1_p_inter[i]->SetMarkerColor(kRed);  h1_p_inter[i]->SetLineWidth(2);   h1_p_inter[i]->Sumw2();
        h1_d_inter[i]->SetLineColor(kBlue);  h1_d_inter[i]->SetMarkerColor(kBlue); h1_d_inter[i]->SetLineWidth(2);   h1_d_inter[i]->Sumw2();
        h1_p_inter[i]->Scale(1.0/h1_p_inter[i]->Integral()); 
        h1_d_inter[i]->Scale(1.0/h1_d_inter[i]->Integral()); 
        h1_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
        h1_p_inter[i]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] Stack Multi Layer;log10(Rm) = log10(Max Energy Deposit bar/ Total Deposit);Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        h1_p_inter[i]->Draw("hist");
        h1_d_inter[i]->Draw("histsame");

        auto c4 = new TCanvas("c4","c4",1000,1000);
        c4->cd();
        hC_p_inter[i]->SetLineColor(kRed);   hC_p_inter[i]->SetMarkerColor(kRed);  hC_p_inter[i]->SetLineWidth(2);   hC_p_inter[i]->Sumw2();
        hC_d_inter[i]->SetLineColor(kBlue);  hC_d_inter[i]->SetMarkerColor(kBlue); hC_d_inter[i]->SetLineWidth(2);   hC_d_inter[i]->Sumw2();
        hC_p_inter[i]->Scale(1.0/hC_p_inter[i]->Integral()); 
        hC_d_inter[i]->Scale(1.0/hC_d_inter[i]->Integral()); 
        hC_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
        hC_p_inter[i]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] Stack Multi Layer;#sum_{0}^{Interaction Layer} (dE/dx) / Edep;Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        hC_p_inter[i]->Draw("hist");
        hC_d_inter[i]->Draw("histsame");
    }
}
