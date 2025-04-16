void Lateral_RMS()
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
    deuteron_tree_tree->SetBranchAddress("Total_E"     ,&d_Total_E);

    auto electron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Electron_PowerLaw.root");
    auto electron_tree = (TTree*)electron_file->Get("B4");
    electron_tree->SetBranchAddress("RMS"              ,&e_RMSVec);
    electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
    electron_tree->SetBranchAddress("Efrac"            ,&e_Efrac);
    electron_tree->SetBranchAddress("First_Had_Layer"  ,&e_First_Had_Layer);
    electron_tree->SetBranchAddress("Total_E"          ,&e_Total_E);

    auto helium4_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium4_PowerLaw.root");
    auto helium4_tree = (TTree*)helium4_file->Get("B4");
    helium4_tree->SetBranchAddress("RMS"              ,&h_RMSVec);
    helium4_tree->SetBranchAddress("LayerEnergyVector",&h_EnergyVec);
    helium4_tree->SetBranchAddress("Efrac"            ,&h_Efrac);
    helium4_tree->SetBranchAddress("First_Had_Layer"  ,&h_First_Had_Layer);
    helium4_tree->SetBranchAddress("Total_E"          ,&h_Total_E);

    auto helium3_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium3_PowerLaw.root");
    auto helium3_tree = (TTree*)helium3_file->Get("B4");
    helium3_tree->SetBranchAddress("RMS"              ,&H_RMSVec);
    helium3_tree->SetBranchAddress("LayerEnergyVector",&H_EnergyVec);
    helium3_tree->SetBranchAddress("Efrac"            ,&H_Efrac);
    helium3_tree->SetBranchAddress("First_Had_Layer"  ,&H_First_Had_Layer);
    helium3_tree->SetBranchAddress("Total_E"          ,&H_Total_E);

    auto carbon_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Carbon_PowerLaw.root");
    auto carbon_tree = (TTree*)carbon_file->Get("B4");
    carbon_tree->SetBranchAddress("RMS"              ,&c_RMSVec);
    carbon_tree->SetBranchAddress("LayerEnergyVector",&c_EnergyVec);
    carbon_tree->SetBranchAddress("Efrac"            ,&c_Efrac);
    carbon_tree->SetBranchAddress("First_Had_Layer"  ,&c_First_Had_Layer);
    carbon_tree->SetBranchAddress("Total_E"          ,&c_Total_E);
    
    double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};    double Electron_Ratio[14]={0};    double Helium4_Ratio[14]={0};     double Helium3_Ratio[14]={0};     double Carbon_Ratio[14]={0};        double Energy[13]={0};
    double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0}; double Electron_Ratio_LL[14]={0}; double Helium4_Ratio_LL[14]={0};  double Helium3_Ratio_LL[14]={0};  double Carbon_Ratio_LL[14]={0};     double Energy_LL[13]={0};      
    double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0}; double Electron_Ratio_UL[14]={0}; double Helium4_Ratio_UL[14]={0};  double Helium3_Ratio_UL[14]={0};  double Carbon_Ratio_UL[14]={0};     double Energy_UL[13]={0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[13][14];
    TH1D *h1_d[13][14];
    TH1D *h1_e[13][14];
    TH1D *h1_h[13][14];
    TH1D *h1_H[13][14];
    TH1D *h1_c[13][14];

    for (int i = 0; i < 1; i++) // Deposit Energy Bin
    {
        Energy[i]    =  0.2*i+1.1;
        Energy_LL[i] =  0.2*i+1.0;
        Energy_UL[i] =  0.2*i+1.2;

        double Layer[14]={0};
        double Layer_Err[14]={0};

        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j < 14; j++) // layer
        {
            c1->cd(j + 1);
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),75,0,150);  
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),75,0,150);  
            h1_e[i][j] = new TH1D(Form("h1_e[%d][%d]",i,j), Form("h1_e[%d][%d]",i,j),75,0,150);  
            h1_h[i][j] = new TH1D(Form("h1_h[%d][%d]",i,j), Form("h1_h[%d][%d]",i,j),75,0,150);  
            h1_H[i][j] = new TH1D(Form("h1_H[%d][%d]",i,j), Form("h1_H[%d][%d]",i,j),75,0,150);  
            h1_c[i][j] = new TH1D(Form("h1_c[%d][%d]",i,j), Form("h1_c[%d][%d]",i,j),75,0,150);  

            for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
            {
                proton_tree->GetEntry(entry);   if ( log10(p_Total_E) > 0.2*k && log10(p_Total_E) <= 0.2*(k+1) ) h1_p[i][j]->Fill(log10((*p_Efrac)[j]));
                deuteron_tree->GetEntry(entry); if ( log10(d_Total_E) > 0.2*k && log10(d_Total_E) <= 0.2*(k+1) ) h1_d[i][j]->Fill(log10((*d_Efrac)[j]));
                electron_tree->GetEntry(entry); if ( log10(e_Total_E) > 0.2*k && log10(e_Total_E) <= 0.2*(k+1) ) h1_e[i][j]->Fill(log10((*e_Efrac)[j]));
                helium4_tree->GetEntry(entry);  if ( log10(h_Total_E) > 0.2*k && log10(h_Total_E) <= 0.2*(k+1) ) h1_h[i][j]->Fill(log10((*h_Efrac)[j]));
                helium3_tree->GetEntry(entry);  if ( log10(H_Total_E) > 0.2*k && log10(H_Total_E) <= 0.2*(k+1) ) h1_H[i][j]->Fill(log10((*H_Efrac)[j]));
                carbon_tree->GetEntry(entry);   if ( log10(c_Total_E) > 0.2*k && log10(c_Total_E) <= 0.2*(k+1) ) h1_c[i][j]->Fill(log10((*c_Efrac)[j]));
            }

            h1_p[i][j]->Sumw2(); h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral());h1_p[i][j]->SetLineColor(kRed);     h1_p[i][j]->SetMarkerColor(kRed);     h1_p[i][j]->SetLineWidth(2);
            h1_d[i][j]->Sumw2(); h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral());h1_d[i][j]->SetLineColor(kBlue);    h1_d[i][j]->SetMarkerColor(kBlue);    h1_d[i][j]->SetLineWidth(2);
            h1_e[i][j]->Sumw2(); h1_e[i][j]->Scale(1.0/h1_e[i][j]->Integral());h1_e[i][j]->SetLineColor(kOrange-3);h1_e[i][j]->SetMarkerColor(kOrange-3);h1_e[i][j]->SetLineWidth(2);
            h1_h[i][j]->Sumw2(); h1_h[i][j]->Scale(1.0/h1_h[i][j]->Integral());h1_h[i][j]->SetLineColor(kGreen-3); h1_h[i][j]->SetMarkerColor(kGreen-3); h1_h[i][j]->SetLineWidth(2);
            h1_H[i][j]->Sumw2(); h1_H[i][j]->Scale(1.0/h1_H[i][j]->Integral());h1_H[i][j]->SetLineColor(kGreen-3); h1_H[i][j]->SetMarkerColor(kGreen-3); h1_H[i][j]->SetLineWidth(2);
            h1_c[i][j]->Sumw2(); h1_c[i][j]->Scale(1.0/h1_c[i][j]->Integral());h1_c[i][j]->SetLineColor(kMagenta); h1_c[i][j]->SetMarkerColor(kMagenta); h1_c[i][j]->SetLineWidth(2);

            h1_e[i][j]->GetYaxis()->SetRangeUser(0,0.25);h1_e[j]->SetTitle(Form("10^%.f - %.f GeV Deposit EdepRatio Distrubution in L%d;RMS Distribution;Normalized Count", 1+0.2*k ,1.2+0.2*k ,j ));
            h1_e[i][j]->Draw();
            h1_H[i][j]->Draw("same");
            h1_c[i][j]->Draw("same");
            h1_h[i][j]->Draw("same");
            h1_p[i][j]->Draw("same");
            h1_d[i][j]->Draw("same");

            double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
            double p_values[3];  h1_p[i][j]->GetQuantiles(3, p_values, quantiles);
            double d_values[3];  h1_d[i][j]->GetQuantiles(3, d_values, quantiles);
            double e_values[3];  h1_e[i][j]->GetQuantiles(3, e_values, quantiles);
            double h_values[3];  h1_h[i][j]->GetQuantiles(3, h_values, quantiles);
            double H_values[3];  h1_H[i][j]->GetQuantiles(3, H_values, quantiles);
            double c_values[3];  h1_c[i][j]->GetQuantiles(3, c_values, quantiles);

            Proton_Ratio[j] = p_values[1];     Proton_Ratio_LL[j]  = p_values[1] - p_values[0];   Proton_Ratio_UL[j]= p_values[2] - p_values[1];
            Deuteron_Ratio[j] = d_values[1];   Deuteron_Ratio_LL[j]= d_values[1] - d_values[0];   Deuteron_Ratio_UL[j]= d_values[2] - d_values[1];
            Electron_Ratio[j] = e_values[1];   Electron_Ratio_LL[j]= e_values[1] - e_values[0];   Electron_Ratio_UL[j]= e_values[2] - e_values[1];
            Helium4_Ratio[j] = h_values[1];    Helium4_Ratio_LL[j] = h_values[1] - h_values[0];   Helium4_Ratio_UL[j]= h_values[2] - h_values[1];
            Helium3_Ratio[j] = H_values[1];    Helium3_Ratio_LL[j] = H_values[1] - H_values[0];   Helium3_Ratio_UL[j]= H_values[2] - H_values[1];
            Carbon_Ratio[j] = c_values[1];     Carbon_Ratio_LL[j] = c_values[1] - c_values[0];    Carbon_Ratio_UL[j]= c_values[2] - c_values[1];

            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
        }
        c1->cd(15);
        TLatex *tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));tex->SetNDC();tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[i][0], "Proton", "l");
        legend1->AddEntry(h1_d[i][0], "Deuteron", "l");
        legend1->AddEntry(h1_e[i][0], "Electron", "l");
        legend1->AddEntry(h1_h[i][0], "Helium4", "l");         
        legend1->AddEntry(h1_H[i][0], "Helium3", "l");         
        legend1->AddEntry(h1_c[i][0], "Carbon", "l");         
        legend1->Draw();       
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Lateral_PowerE/Lateral_RMS_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));


        auto c2 = new TCanvas("c2","c2",2500,1500);
        // c2->cd();
        auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_Ratio  ,Layer_Err,Layer_Err,Proton_Ratio_LL  ,Proton_Ratio_UL);
        auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_Ratio,Layer_Err,Layer_Err,Deuteron_Ratio_LL,Deuteron_Ratio_UL);
        auto gre_e = new TGraphAsymmErrors(14,Layer,Electron_Ratio,Layer_Err,Layer_Err,Electron_Ratio_LL,Electron_Ratio_UL);
        auto gre_h = new TGraphAsymmErrors(14,Layer,Helium4_Ratio ,Layer_Err,Layer_Err,Helium4_Ratio_LL ,Helium4_Ratio_UL);
        auto gre_H = new TGraphAsymmErrors(14,Layer,Helium3_Ratio ,Layer_Err,Layer_Err,Helium3_Ratio_LL ,Helium3_Ratio_UL);
        auto gre_c = new TGraphAsymmErrors(14,Layer,Carbon_Ratio ,Layer_Err,Layer_Err,Carbon_Ratio_LL ,Carbon_Ratio_UL);

        gre_e->SetTitle(Form("RMS Distribution 10^{%.1f} - 10^{%.1f} GeV ; BGO Layer; RMS Distribution",Energy_LL[i],Energy_UL[i]));
        gre_e->SetMarkerStyle(22);
        
        gre_e->SetMarkerColor(kOrange-3);
        gre_e->SetLineColor(kOrange-3);
        gre_p->SetMarkerStyle(20);
        gre_p->SetMarkerColor(kRed);
        gre_p->SetLineColor(kRed);
        gre_d->SetMarkerStyle(21);
        gre_d->SetMarkerColor(kBlue);
        gre_d->SetLineColor(kBlue);
        gre_h->SetMarkerStyle(23);
        gre_h->SetMarkerColor(kGreen-3);
        gre_h->SetLineColor(kGreen-3);
        gre_H->SetMarkerColor(kGreen-3);
        gre_H->SetLineColor(kGreen-3);
        gre_c->SetMarkerColor(kMagenta);
        gre_c->SetLineColor(kMagenta);

        gre_p->SetLineWidth(2);
        gre_d->SetLineWidth(2);
        gre_e->SetLineWidth(2);
        gre_h->SetLineWidth(2);
        gre_H->SetLineWidth(2);gre_H->SetLineStyle(2);
        gre_c->SetLineWidth(2);
        
        gre_e->Draw("AP");
    
        gre_h->Draw("PSAME");
        gre_H->Draw("PSAME");
        gre_c->Draw("PSAME");
        gre_d->Draw("PSAME");
        gre_p->Draw("PSAME");

        auto legend2 = new TLegend(0.42, 0.12, 0.58, 0.32);
        legend2->SetNColumns(2);
        legend2->AddEntry(gre_p, "HET Proton", "ep");
        legend2->AddEntry(gre_d, "HET Deuteron", "ep");
        legend2->AddEntry(gre_e, "HET Electron", "ep");
        legend2->AddEntry(gre_h, "HET Helium4", "ep");
        legend2->AddEntry(gre_H, "HET Helium3", "ep");
        legend2->AddEntry(gre_c, "HET Carbon", "ep");
        legend2->Draw();

        c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Lateral_PowerE/Lateral_RMS_BGOLayer_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));
    }
    
}