void Fval()
{
    double Energy[19]={0};
    for (int i = 0; i < 19; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        int p_First_Had_Layer;    std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Fval = nullptr;
        int d_First_Had_Layer;    std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Fval = nullptr;
        int e_First_Had_Layer;    std::vector<double>* e_RMSVec = nullptr;    std::vector<double>* e_EnergyVec = nullptr;    std::vector<double>* e_Fval = nullptr;
        int h_First_Had_Layer;    std::vector<double>* h_RMSVec = nullptr;    std::vector<double>* h_EnergyVec = nullptr;    std::vector<double>* h_Fval = nullptr;
        int H_First_Had_Layer;    std::vector<double>* H_RMSVec = nullptr;    std::vector<double>* H_EnergyVec = nullptr;    std::vector<double>* H_Fval = nullptr;
        int c_First_Had_Layer;    std::vector<double>* c_RMSVec = nullptr;    std::vector<double>* c_EnergyVec = nullptr;    std::vector<double>* c_Fval = nullptr;

        TH1D *h1_p[14];  TF1  *fitFunc_p[14];
        TH1D *h1_d[14];  TF1  *fitFunc_d[14];
        TH1D *h1_e[14];  TF1  *fitFunc_e[14];
        TH1D *h1_h[14];  TF1  *fitFunc_h[14];
        TH1D *h1_H[14];  TF1  *fitFunc_H[14];
        TH1D *h1_c[14];  TF1  *fitFunc_c[14];

        double Proton_RMS[14]={0};     double Deuteron_RMS[14]={0};    double Electron_RMS[14]={0};    double Helium4_RMS[14]={0};     double Helium3_RMS[14]={0};     double Carbon_RMS[14]={0};    
        double Proton_RMS_LL[14]={0};  double Deuteron_RMS_LL[14]={0}; double Electron_RMS_LL[14]={0}; double Helium4_RMS_LL[14]={0};  double Helium3_RMS_LL[14]={0};  double Carbon_RMS_LL[14]={0};       
        double Proton_RMS_UL[14]={0};  double Deuteron_RMS_UL[14]={0}; double Electron_RMS_UL[14]={0}; double Helium4_RMS_UL[14]={0};  double Helium3_RMS_UL[14]={0};  double Carbon_RMS_UL[14]={0};       

        double Layer[14]={0};
        double Layer_Err[14]={0};

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
        proton_tree->SetBranchAddress("Fval"             ,&p_Fval);
        proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
        deuteron_tree->SetBranchAddress("Fval"             ,&d_Fval);
        deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        electron_tree->SetBranchAddress("RMS"              ,&e_RMSVec);
        electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
        electron_tree->SetBranchAddress("Fval"             ,&e_Fval);
        electron_tree->SetBranchAddress("First_Had_Layer"  ,&e_First_Had_Layer);
        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        helium4_tree->SetBranchAddress("RMS"              ,&h_RMSVec);
        helium4_tree->SetBranchAddress("LayerEnergyVector",&h_EnergyVec);
        helium4_tree->SetBranchAddress("Fval"             ,&h_Fval);
        helium4_tree->SetBranchAddress("First_Had_Layer"  ,&h_First_Had_Layer);
        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        helium3_tree->SetBranchAddress("RMS"              ,&H_RMSVec);
        helium3_tree->SetBranchAddress("LayerEnergyVector",&H_EnergyVec);
        helium3_tree->SetBranchAddress("Fval"             ,&H_Fval);
        helium3_tree->SetBranchAddress("First_Had_Layer"  ,&H_First_Had_Layer);
        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        carbon_tree->SetBranchAddress("RMS"              ,&c_RMSVec);
        carbon_tree->SetBranchAddress("LayerEnergyVector",&c_EnergyVec);
        carbon_tree->SetBranchAddress("Fval"             ,&c_Fval);
        carbon_tree->SetBranchAddress("First_Had_Layer"  ,&c_First_Had_Layer);

        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j < 14; j++) // Sub Fig
        {
            c1->cd(j + 1);
            h1_p[j] = new TH1D(Form("h1_p[%d]",j), Form("h1_p[%d]",j),100,-3,2);  
            h1_d[j] = new TH1D(Form("h1_d[%d]",j), Form("h1_d[%d]",j),100,-3,2);  
            h1_e[j] = new TH1D(Form("h1_e[%d]",j), Form("h1_e[%d]",j),100,-3,2);  
            h1_h[j] = new TH1D(Form("h1_h[%d]",j), Form("h1_h[%d]",j),100,-3,2);  
            h1_H[j] = new TH1D(Form("h1_H[%d]",j), Form("h1_H[%d]",j),100,-3,2);  
            h1_c[j] = new TH1D(Form("h1_c[%d]",j), Form("h1_c[%d]",j),100,-3,2);  

            for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
            {
                proton_tree->GetEntry(entry);  /*if ((*p_EnergyVec)[0] > 0.23 && (*p_EnergyVec)[1] > 0.23 && (*p_EnergyVec)[2] > 0.23 && (*p_EnergyVec)[0] > 0.046)*/  h1_p[j]->Fill(log10((*p_Fval)[j]));
                deuteron_tree->GetEntry(entry);/*if ((*d_EnergyVec)[0] > 0.23 && (*d_EnergyVec)[1] > 0.23 && (*d_EnergyVec)[2] > 0.23 && (*d_EnergyVec)[0] > 0.046)*/  h1_d[j]->Fill(log10((*d_Fval)[j]));
                electron_tree->GetEntry(entry);/*if ((*e_EnergyVec)[0] > 0.23 && (*e_EnergyVec)[1] > 0.23 && (*e_EnergyVec)[2] > 0.23 && (*e_EnergyVec)[0] > 0.046)*/  h1_e[j]->Fill(log10((*e_Fval)[j]));
                helium4_tree->GetEntry(entry); /*if ((*h_EnergyVec)[0] > 0.23 && (*h_EnergyVec)[1] > 0.23 && (*h_EnergyVec)[2] > 0.23 && (*h_EnergyVec)[0] > 0.046)*/  h1_h[j]->Fill(log10((*h_Fval)[j]));
                helium3_tree->GetEntry(entry); /*if ((*H_EnergyVec)[0] > 0.23 && (*H_EnergyVec)[1] > 0.23 && (*H_EnergyVec)[2] > 0.23 && (*H_EnergyVec)[0] > 0.046)*/  h1_H[j]->Fill(log10((*H_Fval)[j]));
                carbon_tree->GetEntry(entry);  /*if ((*c_EnergyVec)[0] > 0.23 && (*c_EnergyVec)[1] > 0.23 && (*c_EnergyVec)[2] > 0.23 && (*c_EnergyVec)[0] > 0.046)*/  h1_c[j]->Fill(log10((*c_Fval)[j]));
            }
            h1_p[j]->Sumw2(); h1_p[j]->Scale(1.0/h1_p[j]->Integral());h1_p[j]->SetLineColor(kRed);     h1_p[j]->SetMarkerColor(kRed);     h1_p[j]->SetLineWidth(2);
            h1_d[j]->Sumw2(); h1_d[j]->Scale(1.0/h1_d[j]->Integral());h1_d[j]->SetLineColor(kBlue);    h1_d[j]->SetMarkerColor(kBlue);    h1_d[j]->SetLineWidth(2);
            h1_e[j]->Sumw2(); h1_e[j]->Scale(1.0/h1_e[j]->Integral());h1_e[j]->SetLineColor(kOrange-3);h1_e[j]->SetMarkerColor(kOrange-3);h1_e[j]->SetLineWidth(2);
            h1_h[j]->Sumw2(); h1_h[j]->Scale(1.0/h1_h[j]->Integral());h1_h[j]->SetLineColor(kGreen-3); h1_h[j]->SetMarkerColor(kGreen-3); h1_h[j]->SetLineWidth(2);
            h1_H[j]->Sumw2(); h1_H[j]->Scale(1.0/h1_H[j]->Integral());h1_H[j]->SetLineColor(kGreen-3); h1_H[j]->SetMarkerColor(kGreen-3); h1_H[j]->SetLineWidth(2);
            h1_c[j]->Sumw2(); h1_c[j]->Scale(1.0/h1_c[j]->Integral());h1_c[j]->SetLineColor(kMagenta); h1_c[j]->SetMarkerColor(kMagenta); h1_c[j]->SetLineWidth(2);

            h1_e[j]->GetYaxis()->SetRangeUser(0,0.25);h1_e[j]->SetTitle(Form("Lateral RMS Distribution in L%d;RMS (mm);Normalized Count",j));
            h1_e[j]->Draw();
            h1_H[j]->Draw("same");
            h1_c[j]->Draw("same");
            h1_h[j]->Draw("same");
            h1_p[j]->Draw("same");
            h1_d[j]->Draw("same");

            double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
            double p_values[3];  h1_p[j]->GetQuantiles(3, p_values, quantiles);
            double d_values[3];  h1_d[j]->GetQuantiles(3, d_values, quantiles);
            double e_values[3];  h1_e[j]->GetQuantiles(3, e_values, quantiles);
            double h_values[3];  h1_h[j]->GetQuantiles(3, h_values, quantiles);
            double H_values[3];  h1_H[j]->GetQuantiles(3, H_values, quantiles);
            double c_values[3];  h1_c[j]->GetQuantiles(3, c_values, quantiles);

            Proton_RMS[j] = p_values[1];     Proton_RMS_LL[j]  = p_values[1] - p_values[0];   Proton_RMS_UL[j]= p_values[2] - p_values[1];
            Deuteron_RMS[j] = d_values[1];   Deuteron_RMS_LL[j]= d_values[1] - d_values[0];   Deuteron_RMS_UL[j]= d_values[2] - d_values[1];
            Electron_RMS[j] = e_values[1];   Electron_RMS_LL[j]= e_values[1] - e_values[0];   Electron_RMS_UL[j]= e_values[2] - e_values[1];
            Helium4_RMS[j] = h_values[1];    Helium4_RMS_LL[j] = h_values[1] - h_values[0];   Helium4_RMS_UL[j]= h_values[2] - h_values[1];
            Helium3_RMS[j] = H_values[1];    Helium3_RMS_LL[j] = H_values[1] - H_values[0];   Helium3_RMS_UL[j]= H_values[2] - H_values[1];
            Carbon_RMS[j] = c_values[1];     Carbon_RMS_LL[j] = c_values[1] - c_values[0];    Carbon_RMS_UL[j]= c_values[2] - c_values[1];

            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;

        }   

        c1->cd(15);
        TLatex *tex = new TLatex(0.1,0.9,Form("IncidentParicle %dGeV",int(Energy[i])));tex->SetNDC();tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[0], "Proton", "l");
        legend1->AddEntry(h1_d[0], "Deuteron", "l");
        legend1->AddEntry(h1_e[0], "Electron", "l");
        legend1->AddEntry(h1_h[0], "Helium4", "l");         
        legend1->AddEntry(h1_H[0], "Helium3", "l");         
        legend1->AddEntry(h1_c[0], "Carbon", "l");         
        legend1->Draw();       
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Fval_MonoE/Fval_%dGeV.pdf",int(Energy[i])));

        auto c2 = new TCanvas("c2","c2",1500,1000);
        c2->cd();
        auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_RMS  ,Layer_Err,Layer_Err,Proton_RMS_LL  ,Proton_RMS_UL);
        auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_RMS,Layer_Err,Layer_Err,Deuteron_RMS_LL,Deuteron_RMS_UL);
        auto gre_e = new TGraphAsymmErrors(14,Layer,Electron_RMS,Layer_Err,Layer_Err,Electron_RMS_LL,Electron_RMS_UL);
        auto gre_h = new TGraphAsymmErrors(14,Layer,Helium4_RMS ,Layer_Err,Layer_Err,Helium4_RMS_LL ,Helium4_RMS_UL);
        auto gre_H = new TGraphAsymmErrors(14,Layer,Helium3_RMS ,Layer_Err,Layer_Err,Helium3_RMS_LL ,Helium3_RMS_UL);
        auto gre_c = new TGraphAsymmErrors(14,Layer,Carbon_RMS ,Layer_Err,Layer_Err,Carbon_RMS_LL ,Carbon_RMS_UL);

        gre_e->SetTitle(Form("Incident Energy %d GeV ; BGO Layer; RMS (mm)",int(Energy[i])));
        gre_e->SetMarkerStyle(22);
        gre_e->GetYaxis()->SetRangeUser(0,150);
        
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

        auto legend2 = new TLegend(0.72, 0.68, 0.88, 0.88);
        legend2->SetNColumns(2);
        legend2->AddEntry(gre_p, "HET Proton", "ep");
        legend2->AddEntry(gre_d, "HET Deuteron", "ep");
        legend2->AddEntry(gre_e, "HET Electron", "ep");
        legend2->AddEntry(gre_h, "HET Helium4", "ep");
        legend2->AddEntry(gre_H, "HET Helium3", "ep");
        legend2->AddEntry(gre_c, "HET Carbon", "ep");
        legend2->Draw();

        c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Fval_MonoE/Fval_BGOLayer_%dGeV.pdf",int(Energy[i])));
    }
}