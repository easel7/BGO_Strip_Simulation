void Rm2()
{
    int p_First_Had_Layer;  double p_Total_E;  std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr;
    int d_First_Had_Layer;  double d_Total_E;  std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr;
    int e_First_Had_Layer;  double e_Total_E;  std::vector<double>* e_RMSVec = nullptr;    std::vector<double>* e_EnergyVec = nullptr;    std::vector<double>* e_Efrac = nullptr;
    int h_First_Had_Layer;  double h_Total_E;  std::vector<double>* h_RMSVec = nullptr;    std::vector<double>* h_EnergyVec = nullptr;    std::vector<double>* h_Efrac = nullptr;
    int H_First_Had_Layer;  double H_Total_E;  std::vector<double>* H_RMSVec = nullptr;    std::vector<double>* H_EnergyVec = nullptr;    std::vector<double>* H_Efrac = nullptr;
    int c_First_Had_Layer;  double c_Total_E;  std::vector<double>* c_RMSVec = nullptr;    std::vector<double>* c_EnergyVec = nullptr;    std::vector<double>* c_Efrac = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    
    double Proton_Ratio[15]={0};     double Deuteron_Ratio[15]={0};    double Energy[15]    = {0};
    double Proton_Ratio_LL[15]={0};  double Deuteron_Ratio_LL[15]={0}; double Energy_LL[15] = {0};
    double Proton_Ratio_UL[15]={0};  double Deuteron_Ratio_UL[15]={0}; double Energy_UL[15] = {0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[15][14];
    TH1D *h1_d[15][14];

    double p_maxVal[15][14]={0};
    double d_maxVal[15][14]={0};

    auto c1 = new TCanvas("c1","c1",500,500);
    // auto c1 = new TCanvas("c1","c1",1200,900);
    c1->Clear();
    c1->Divide(5,3);
    gStyle->SetOptStat(0);

    // auto c2 = new TCanvas("c2","c2",2500,1500);
    // c2->Clear();
    // c2->Divide(5,3);
    // gStyle->SetOptStat(0);

    for(int i =0 ; i<15 ; i++)
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.2 + 0.2 * i;
        Energy_UL[i] = 1.0 + 0.2 * i;
        for( int j= 0; j<14 ;j++)
        {
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),50,-5,0); 
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),50,-5,0); 
        }
        if (i < 14)
        {
            Layer[i]      =  i + 0.5;
            Layer_Err[i]  =  0.5; 
        }

    }
    
    // for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    for (Long64_t entry = 0; entry < 1; ++entry)
    {
        proton_tree->GetEntry(entry);   
        deuteron_tree->GetEntry(entry); 
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        cout << p_Total_E << "," <<  p_energy_index << endl;
        cout << d_Total_E << "," <<  d_energy_index << endl;

        for (int k = 0; k < 14; k++ ) // layer
        {  
            int i = 22*k;
            auto p_start = p_EnergyVec->begin() + i;  auto p_end = (i + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  

            if(p_energy_index >= 0 && p_energy_index <= 14 && (*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 && (*p_RMSVec)[0]<40 && (*p_RMSVec)[1]<40) 
            {    

                p_maxVal[p_energy_index][k] = *std::max_element(p_start, p_end); 
                h1_p[p_energy_index][k]->Fill(log10(p_maxVal[p_energy_index][k]/p_Total_E)); 
            }
            auto d_start = d_EnergyVec->begin() + i;  auto d_end = (i + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end(); 
            // for (auto it = d_start; it != d_end; ++it) {
            //     cout << *it << " ";
            // }
            if(d_energy_index >= 0 && d_energy_index <= 14 && (*d_RMSVec)[0]>15 && (*d_RMSVec)[1]>15 && (*d_RMSVec)[0]<40 && (*d_RMSVec)[1]<40) 
            { 
                // d_maxVal[d_energy_index][k] = *std::max_element(d_start, d_end); 
                // h1_d[d_energy_index][k]->Fill(log10(d_maxVal[d_energy_index][k]/d_Total_E)); 
            }
        }
    }

    // for(int i =0 ; i<15 ; i++) // Energy
    // {
    //     for( int j= 0; j<14 ;j++) // layer
    //     {
    //         c1->cd(j+1);
    //         h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
    //         h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
    //         h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral()); 
    //         h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral()); 

    //         h1_p[i][j]->GetYaxis()->SetRangeUser(0,0.6);
    //         h1_p[i][j]->SetTitle(Form("%.f GeV incident in L%d;Rm = Max Energy Deposit bar/ Total Deposit;Normalized Count",Energy[i],j));
    //         h1_p[i][j]->Draw();
    //         h1_p[i][j]->Draw("hist");
    //         h1_d[i][j]->Draw("histsame");
    //     }
    // }
    
}

//TLatex *tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));tex->SetNDC();tex->Draw(); 
// auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
// legend1->AddEntry(h1_p[i][0], "Proton", "l");
// legend1->AddEntry(h1_d[i][0], "Deuteron", "l");
// legend1->Draw();       
// // c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Longti_PowerE/Longti_Efrac_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));


// auto c2 = new TCanvas("c2","c2",2500,1500);
// // c2->cd();
// auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_Ratio  ,Layer_Err,Layer_Err,Proton_Ratio_LL  ,Proton_Ratio_UL);
// auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_Ratio,Layer_Err,Layer_Err,Deuteron_Ratio_LL,Deuteron_Ratio_UL);


// gre_p->SetTitle(Form("Deposit Energy 10^{%.1f} - 10^{%.1f} GeV ; BGO Layer; log10(Deposit Energy/GeV)",Energy_LL[i],Energy_UL[i]));
// gre_p->SetMarkerStyle(20);gre_p->SetMarkerColor(kRed);     gre_p->SetLineColor(kRed);     gre_p->SetLineWidth(2);
// gre_d->SetMarkerStyle(21);gre_d->SetMarkerColor(kBlue);    gre_d->SetLineColor(kBlue);    gre_d->SetLineWidth(2);

// gre_d->Draw("AP");
// gre_p->Draw("PSAME");

// auto legend2 = new TLegend(0.42, 0.12, 0.58, 0.32);
// legend2->SetNColumns(2);
// legend2->AddEntry(gre_p, "Proton", "ep");
// legend2->AddEntry(gre_d, "Deuteron", "ep");
// legend2->Draw();

// c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Longti_PowerE/Longti_Efrac_BGOLayer_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));