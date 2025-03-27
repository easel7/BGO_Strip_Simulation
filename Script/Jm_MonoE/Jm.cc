void Jm()
{
    double Energy[19]={0};
    for (int i = 0; i < 19; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        int p_First_Had_Layer; double p_E_total;   std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;
        int d_First_Had_Layer; double d_E_total;   std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;
        int e_First_Had_Layer; double e_E_total;   std::vector<double>* e_RMSVec = nullptr;    std::vector<double>* e_EnergyVec = nullptr;    std::vector<double>* e_L_EnergyVec = nullptr;
        int h_First_Had_Layer; double h_E_total;   std::vector<double>* h_RMSVec = nullptr;    std::vector<double>* h_EnergyVec = nullptr;    std::vector<double>* h_L_EnergyVec = nullptr;
        int H_First_Had_Layer; double H_E_total;   std::vector<double>* H_RMSVec = nullptr;    std::vector<double>* H_EnergyVec = nullptr;    std::vector<double>* H_L_EnergyVec = nullptr;
        int c_First_Had_Layer; double c_E_total;   std::vector<double>* c_RMSVec = nullptr;    std::vector<double>* c_EnergyVec = nullptr;    std::vector<double>* c_L_EnergyVec = nullptr;

        TH1D *h1_p[14];  TF1  *fitFunc_p[14]; double p_maxVal[14]={0};
        TH1D *h1_d[14];  TF1  *fitFunc_d[14]; double d_maxVal[14]={0};
        TH1D *h1_e[14];  TF1  *fitFunc_e[14]; double e_maxVal[14]={0};
        TH1D *h1_h[14];  TF1  *fitFunc_h[14]; double h_maxVal[14]={0};
        TH1D *h1_H[14];  TF1  *fitFunc_H[14]; double H_maxVal[14]={0};
        TH1D *h1_c[14];  TF1  *fitFunc_c[14]; double c_maxVal[14]={0};

        double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};    double Electron_Ratio[14]={0};    double Helium4_Ratio[14]={0};     double Helium3_Ratio[14]={0};     double Carbon_Ratio[14]={0};    
        double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0}; double Electron_Ratio_LL[14]={0}; double Helium4_Ratio_LL[14]={0};  double Helium3_Ratio_LL[14]={0};  double Carbon_Ratio_LL[14]={0};       
        double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0}; double Electron_Ratio_UL[14]={0}; double Helium4_Ratio_UL[14]={0};  double Helium3_Ratio_UL[14]={0};  double Carbon_Ratio_UL[14]={0};       

        double Layer[14]={0};
        double Layer_Err[14]={0};

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
        proton_tree->SetBranchAddress("BarEnergyVector",&p_EnergyVec);
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
        proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
        proton_tree->SetBranchAddress("Total_E"          ,&p_E_total);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
        deuteron_tree->SetBranchAddress("BarEnergyVector",&d_EnergyVec);
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
        deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
        deuteron_tree->SetBranchAddress("Total_E"          ,&d_E_total);

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        electron_tree->SetBranchAddress("RMS"              ,&e_RMSVec);
        electron_tree->SetBranchAddress("BarEnergyVector",&e_EnergyVec);
        electron_tree->SetBranchAddress("LayerEnergyVector",&e_L_EnergyVec);
        electron_tree->SetBranchAddress("First_Had_Layer"  ,&e_First_Had_Layer);
        electron_tree->SetBranchAddress("Total_E"          ,&e_E_total);

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        helium4_tree->SetBranchAddress("RMS"              ,&h_RMSVec);
        helium4_tree->SetBranchAddress("BarEnergyVector",&h_EnergyVec);
        helium4_tree->SetBranchAddress("LayerEnergyVector",&h_L_EnergyVec);
        helium4_tree->SetBranchAddress("First_Had_Layer"  ,&h_First_Had_Layer);
        helium4_tree->SetBranchAddress("Total_E"          ,&h_E_total);
        
        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        helium3_tree->SetBranchAddress("RMS"              ,&H_RMSVec);
        helium3_tree->SetBranchAddress("BarEnergyVector",&H_EnergyVec);
        helium3_tree->SetBranchAddress("LayerEnergyVector",&H_L_EnergyVec);
        helium3_tree->SetBranchAddress("First_Had_Layer"  ,&H_First_Had_Layer);
        helium3_tree->SetBranchAddress("Total_E"          ,&H_E_total);

        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        carbon_tree->SetBranchAddress("RMS"              ,&c_RMSVec);
        carbon_tree->SetBranchAddress("BarEnergyVector",&c_EnergyVec);
        carbon_tree->SetBranchAddress("LayerEnergyVector",&c_L_EnergyVec);
        carbon_tree->SetBranchAddress("First_Had_Layer"  ,&c_First_Had_Layer);
        carbon_tree->SetBranchAddress("Total_E"          ,&c_E_total);

        auto c1 = new TCanvas("c1","c1",2500,1500);
        // auto c1 = new TCanvas("c1","c1",1200,900);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j < 14; j++)
        {
            h1_p[j] = new TH1D(Form("h1_p[%d]",j), Form("h1_p[%d]",j),40,-1,0);  
            h1_d[j] = new TH1D(Form("h1_d[%d]",j), Form("h1_d[%d]",j),40,-1,0);  
            h1_e[j] = new TH1D(Form("h1_e[%d]",j), Form("h1_e[%d]",j),40,-1,0);  
            h1_h[j] = new TH1D(Form("h1_h[%d]",j), Form("h1_h[%d]",j),40,-1,0);  
            h1_H[j] = new TH1D(Form("h1_H[%d]",j), Form("h1_H[%d]",j),40,-1,0);  
            h1_c[j] = new TH1D(Form("h1_c[%d]",j), Form("h1_c[%d]",j),40,-1,0);  
        }

        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
        // for (Long64_t entry = 0; entry < 1; ++entry)
        {
            proton_tree->GetEntry(entry);   
            deuteron_tree->GetEntry(entry); 
            electron_tree->GetEntry(entry); 
            helium4_tree->GetEntry(entry);  
            helium3_tree->GetEntry(entry);  
            carbon_tree->GetEntry(entry);   
            // cout << p_E_total << endl;
            // cout << p_EnergyVec->size() << endl;
            for (size_t k = 0; k < p_EnergyVec->size(); k += 22)
            {
                auto p_start = p_EnergyVec->begin() + k;  auto p_end = (k + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal[int(k/22)] = *std::max_element(p_start, p_end);  h1_p[int(k/22)]->Fill(log10(p_maxVal[int(k/22)]/(*p_L_EnergyVec)[int(k/22)]));
                auto d_start = d_EnergyVec->begin() + k;  auto d_end = (k + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end();  d_maxVal[int(k/22)] = *std::max_element(d_start, d_end);  h1_d[int(k/22)]->Fill(log10(d_maxVal[int(k/22)]/(*d_L_EnergyVec)[int(k/22)]));
                auto e_start = e_EnergyVec->begin() + k;  auto e_end = (k + 22 < e_EnergyVec->size() ) ? e_start + 22 : e_EnergyVec->end();  e_maxVal[int(k/22)] = *std::max_element(e_start, e_end);  h1_e[int(k/22)]->Fill(log10(e_maxVal[int(k/22)]/(*e_L_EnergyVec)[int(k/22)]));
                auto h_start = h_EnergyVec->begin() + k;  auto h_end = (k + 22 < h_EnergyVec->size() ) ? h_start + 22 : h_EnergyVec->end();  h_maxVal[int(k/22)] = *std::max_element(h_start, h_end);  h1_h[int(k/22)]->Fill(log10(h_maxVal[int(k/22)]/(*h_L_EnergyVec)[int(k/22)]));
                auto H_start = H_EnergyVec->begin() + k;  auto H_end = (k + 22 < H_EnergyVec->size() ) ? H_start + 22 : H_EnergyVec->end();  H_maxVal[int(k/22)] = *std::max_element(H_start, H_end);  h1_H[int(k/22)]->Fill(log10(H_maxVal[int(k/22)]/(*H_L_EnergyVec)[int(k/22)]));
                auto c_start = c_EnergyVec->begin() + k;  auto c_end = (k + 22 < c_EnergyVec->size() ) ? c_start + 22 : c_EnergyVec->end();  c_maxVal[int(k/22)] = *std::max_element(c_start, c_end);  h1_c[int(k/22)]->Fill(log10(c_maxVal[int(k/22)]/(*c_L_EnergyVec)[int(k/22)]));
                // cout << " bar " << k << " Layer " << int(k/22) <<endl; 
                // std::copy(p_start, p_end, std::ostream_iterator<double>(std::cout, ", "));
                // cout << " max = "<< p_maxVal[int(k/22)] <<  endl;
            }

        }

        for (int j = 0; j < 14; j++)
        {
            c1->cd(j + 1);
            h1_p[j]->Sumw2(); h1_p[j]->Scale(1.0/h1_p[j]->Integral()); h1_p[j]->SetLineColor(kRed);     h1_p[j]->SetMarkerColor(kRed);     h1_p[j]->SetLineWidth(2);
            h1_d[j]->Sumw2(); h1_d[j]->Scale(1.0/h1_d[j]->Integral()); h1_d[j]->SetLineColor(kBlue);    h1_d[j]->SetMarkerColor(kBlue);    h1_d[j]->SetLineWidth(2);
            h1_e[j]->Sumw2(); h1_e[j]->Scale(1.0/h1_e[j]->Integral()); h1_e[j]->SetLineColor(kOrange-3);h1_e[j]->SetMarkerColor(kOrange-3);h1_e[j]->SetLineWidth(2);
            h1_h[j]->Sumw2(); h1_h[j]->Scale(1.0/h1_h[j]->Integral()); h1_h[j]->SetLineColor(kGreen-3); h1_h[j]->SetMarkerColor(kGreen-3); h1_h[j]->SetLineWidth(2);
            h1_H[j]->Sumw2(); h1_H[j]->Scale(1.0/h1_H[j]->Integral()); h1_H[j]->SetLineColor(kGreen-3); h1_H[j]->SetMarkerColor(kGreen-3); h1_H[j]->SetLineWidth(2);
            h1_c[j]->Sumw2(); h1_c[j]->Scale(1.0/h1_c[j]->Integral()); h1_c[j]->SetLineColor(kMagenta); h1_c[j]->SetMarkerColor(kMagenta); h1_c[j]->SetLineWidth(2);
            h1_e[j]->GetYaxis()->SetRangeUser(0,0.6);
            h1_e[j]->SetTitle(Form("%.f GeV incident in L%d;log_{10}(Jm) = log_{10}(Max Energy Deposit bar in L%d/ Energy Deposit in L%d);Normalized Count",Energy[i],j,j,j));
            h1_e[j]->Draw();
            h1_H[j]->Draw("same");
            h1_c[j]->Draw("same");
            h1_h[j]->Draw("same");
            h1_p[j]->Draw("histsame");
            h1_d[j]->Draw("histsame");

            double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
            double p_values[3];  h1_p[j]->GetQuantiles(3, p_values, quantiles);
            double d_values[3];  h1_d[j]->GetQuantiles(3, d_values, quantiles);
            double e_values[3];  h1_e[j]->GetQuantiles(3, e_values, quantiles);
            double h_values[3];  h1_h[j]->GetQuantiles(3, h_values, quantiles);
            double H_values[3];  h1_H[j]->GetQuantiles(3, H_values, quantiles);
            double c_values[3];  h1_c[j]->GetQuantiles(3, c_values, quantiles);

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
        TLatex *tex = new TLatex(0.1,0.9,Form("IncidentParicle %dGeV",int(Energy[i])));tex->SetNDC();tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[0], "Proton", "l");
        legend1->AddEntry(h1_d[0], "Deuteron", "l");
        legend1->AddEntry(h1_e[0], "Electron", "l");
        legend1->AddEntry(h1_h[0], "Helium4", "l");         
        legend1->AddEntry(h1_H[0], "Helium3", "l");         
        legend1->AddEntry(h1_c[0], "Carbon", "l");         
        legend1->Draw();       
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Jm_MonoE/Jm_%dGeV.pdf",int(Energy[i])));

        auto c2 = new TCanvas("c2","c2",1000,1000);
        // c2->cd();
        auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_Ratio  ,Layer_Err,Layer_Err,Proton_Ratio_LL  ,Proton_Ratio_UL);
        auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_Ratio,Layer_Err,Layer_Err,Deuteron_Ratio_LL,Deuteron_Ratio_UL);
        auto gre_e = new TGraphAsymmErrors(14,Layer,Electron_Ratio,Layer_Err,Layer_Err,Electron_Ratio_LL,Electron_Ratio_UL);
        auto gre_h = new TGraphAsymmErrors(14,Layer,Helium4_Ratio ,Layer_Err,Layer_Err,Helium4_Ratio_LL ,Helium4_Ratio_UL);
        auto gre_H = new TGraphAsymmErrors(14,Layer,Helium3_Ratio ,Layer_Err,Layer_Err,Helium3_Ratio_LL ,Helium3_Ratio_UL);
        auto gre_c = new TGraphAsymmErrors(14,Layer,Carbon_Ratio ,Layer_Err,Layer_Err,Carbon_Ratio_LL ,Carbon_Ratio_UL);

        gre_e->SetTitle(Form("Incident Energy %d GeV ; BGO Layer (Rad Length); log10(Jm) = log10( (dE/dx) / Edep)",int(Energy[i])));
        gre_e->SetMarkerStyle(22);
        gre_e->GetYaxis()->SetRangeUser(-1,0);
        
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

        c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/Jm_MonoE/Jm_BGOLayer_%dGeV.pdf",int(Energy[i])));
    }
    
}