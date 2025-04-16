void Sum_Rm_L1()
{
    int p_First_Had_Layer;   
    int p_First_Had_Type; double p_E_total;    
    std::vector<double>* p_RMSVec = nullptr;    
    std::vector<double>* p_EnergyVec = nullptr;    
    std::vector<double>* p_L_EnergyVec = nullptr;
    std::vector<double>* p_Efrac = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
    // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_1000GeV.root");

    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("BarEnergyVector"  ,&p_EnergyVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("First_Had_Type"   ,&p_First_Had_Type);
    proton_tree->SetBranchAddress("Total_E"          ,&p_E_total);
  
    cout  << proton_tree->GetEntries() << endl;
    TH1D *h1_p[14];  
    TH1D *h1_d[14];  
    TH1D *h1_e[14];  
    TH1D *h1_E[14];  
    TH1D *h1_h[14];  
    TH1D *h1_H[14];  
    TH1D *h1_c[14];  


    for (int k = 0; k < 14; k++)
    {
        h1_p[k] = new TH1D(Form("h1_p[%d]",k),Form("h1_p[%d]",k),100,0,1);  
        h1_d[k] = new TH1D(Form("h1_d[%d]",k),Form("h1_d[%d]",k),100,0,1);  
        h1_e[k] = new TH1D(Form("h1_e[%d]",k),Form("h1_e[%d]",k),100,0,1);  
        h1_E[k] = new TH1D(Form("h1_E[%d]",k),Form("h1_E[%d]",k),100,0,1);  
        h1_h[k] = new TH1D(Form("h1_h[%d]",k),Form("h1_h[%d]",k),100,0,1);  
        h1_H[k] = new TH1D(Form("h1_H[%d]",k),Form("h1_H[%d]",k),100,0,1);  
        h1_c[k] = new TH1D(Form("h1_c[%d]",k),Form("h1_c[%d]",k),100,0,1);  
    }

    int GL = 3; // Goal Layer
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        proton_tree->GetEntry(entry);    
        double sum_p = 0;
        double p_maxVal[14] = {0};
        if((*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 && (*p_RMSVec)[2]<45 && (*p_RMSVec)[3]<45  ) //  
        {
            for (size_t j = 0; j < p_EnergyVec->size(); j += 22)
            {
                int k = int(j / 22);  // Get the Layer
                auto p_start = p_EnergyVec->begin() + j;  auto p_end = (j + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal[k] = *std::max_element(p_start, p_end); 
                sum_p       += p_maxVal[k];
                if(k>0 && k<13)
                {
                                                                                                                            h1_p[k]->Fill((sum_p / p_E_total));
                    if(p_First_Had_Layer==k)                                                                              { h1_d[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer==-1)                                                                             { h1_h[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=floor((k-1)/2) && p_First_Had_Layer<=k-1             && p_First_Had_Type == 1)  { h1_e[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=0              && p_First_Had_Layer<floor((k-1)/2)   && p_First_Had_Type == 1)  { h1_H[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer<=k-1 && p_First_Had_Type == 2)                                                   { h1_c[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=k+1)                                                                            { h1_E[k]->Fill((sum_p / p_E_total)); }
                }
                else if (k==0)
                {
                                                                                                                            h1_p[k]->Fill((sum_p / p_E_total));
                    if(p_First_Had_Layer==k)                                                                              { h1_d[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer==-1)                                                                             { h1_h[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=k+1 && p_First_Had_Type == 1)                                                   { h1_e[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=k+1 && p_First_Had_Type == 2)                                                   { h1_c[k]->Fill((sum_p / p_E_total)); }
                }
                else if (k==13)
                {
                                                                                                                            h1_p[k]->Fill((sum_p / p_E_total));
                    if(p_First_Had_Layer==k)                                                                              { h1_d[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer==-1)                                                                             { h1_h[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=floor((k-1)/2) && p_First_Had_Layer<=k-1             && p_First_Had_Type == 1)  { h1_e[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer>=0              && p_First_Had_Layer<floor((k-1)/2)   && p_First_Had_Type == 1)  { h1_H[k]->Fill((sum_p / p_E_total)); }
                    if(p_First_Had_Layer<=k-1 && p_First_Had_Type == 2)                                                   { h1_c[k]->Fill((sum_p / p_E_total)); }
                }
            }   
        }
    }

    for(int k = 0; k<14 ; k++)
    {
        h1_p[k]->Sumw2(); /*h1_p[k]->Scale(1.0/h1_p[k]->Integral());*/ h1_p[k]->SetLineColor(kRed);     h1_p[k]->SetMarkerColor(kRed);     h1_p[k]->SetLineWidth(2);
        h1_d[k]->Sumw2(); /*h1_d[k]->Scale(1.0/h1_d[k]->Integral());*/ h1_d[k]->SetLineColor(kBlue);    h1_d[k]->SetMarkerColor(kBlue);    h1_d[k]->SetLineWidth(2);
        h1_e[k]->Sumw2(); /*h1_e[k]->Scale(1.0/h1_e[k]->Integral());*/ h1_e[k]->SetLineColor(kOrange-3);h1_e[k]->SetMarkerColor(kOrange-3);h1_e[k]->SetLineWidth(2);
        h1_h[k]->Sumw2(); /*h1_h[k]->Scale(1.0/h1_h[k]->Integral());*/ h1_h[k]->SetLineColor(kGreen-3); h1_h[k]->SetMarkerColor(kGreen-3); h1_h[k]->SetLineWidth(2);
        h1_c[k]->Sumw2(); /*h1_c[k]->Scale(1.0/h1_c[k]->Integral());*/ h1_c[k]->SetLineColor(kMagenta); h1_c[k]->SetMarkerColor(kMagenta); h1_c[k]->SetLineWidth(2);
        h1_H[k]->Sumw2(); /*h1_H[k]->Scale(1.0/h1_H[k]->Integral());*/ h1_H[k]->SetLineColor(kBlack);   h1_H[k]->SetMarkerColor(kBlack);   h1_H[k]->SetLineWidth(2);
        h1_E[k]->Sumw2(); /*h1_H[k]->Scale(1.0/h1_H[k]->Integral());*/ h1_E[k]->SetLineColor(kCyan);    h1_E[k]->SetMarkerColor(kCyan);    h1_E[k]->SetLineWidth(2);
    }
   
    // h1_p->GetYaxis()->SetRangeUser(0,0.5);h1_p->SetTitle("1000 GeV incident; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L0/ Total Deposit);Normalized Count");
    h1_p[GL]->GetYaxis()->SetRangeUser(0,200);
    h1_p[GL]->SetTitle(Form("1000 GeV incident in L%d; #sum Rm = #sum_{0}^{L%d} Max Energy Deposit bar/ Total Deposit;Count",GL,GL));
    
    auto c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    c1->Clear();
    h1_p[GL]->Draw("hist");h1_p[GL]->SetStats(kFALSE); 
    h1_d[GL]->Draw("same");h1_d[GL]->SetStats(kFALSE); 
    h1_e[GL]->Draw("same");h1_e[GL]->SetStats(kFALSE); 
    h1_h[GL]->Draw("same");h1_h[GL]->SetStats(kFALSE); 
    h1_c[GL]->Draw("same");h1_c[GL]->SetStats(kFALSE); 

    auto legend1 = new TLegend(0.42, 0.58, 0.88, 0.88);

    if(GL>0 && GL<13)
    {
        h1_H[GL]->Draw("same");h1_H[GL]->SetStats(kFALSE); 
        h1_E[GL]->Draw("same");h1_c[GL]->SetStats(kFALSE); 
        legend1->AddEntry(h1_p[GL], "All", "el");
        legend1->AddEntry(h1_d[GL], Form("Had Interaction at L%d",GL), "el");
        legend1->AddEntry(h1_e[GL], Form("Inelastic Interaction at shallower Layer [L%d,L%d]",int(floor((GL-1)/2)),GL-1), "el");
        legend1->AddEntry(h1_H[GL], Form("Inelastic Interaction at shallower Layer [L%d,L%d]",0,int(floor((GL-1)/2))), "el");  
        legend1->AddEntry(h1_c[GL], "Elastic Interaction at shallower Layer", "el");  
        legend1->AddEntry(h1_E[GL], "Had Interaction at deeper Layer", "el");  
        legend1->AddEntry(h1_h[GL], "Pass Through", "el");       
    }
    else if (GL==0)
    {
        legend1->AddEntry(h1_p[GL], "All", "el");
        legend1->AddEntry(h1_d[GL], Form("Had Interaction at L%d",GL), "el");
        legend1->AddEntry(h1_e[GL], "Inelastic Interaction at deeper Layer", "el");
        legend1->AddEntry(h1_c[GL], "Elastic Interaction at deeper Layer", "el");  
        legend1->AddEntry(h1_h[GL], "Pass Through", "el");    
    }
    else if (GL==13)
    {
        h1_H[GL]->Draw("same");h1_H[GL]->SetStats(kFALSE); 
        legend1->AddEntry(h1_p[GL], "All", "el");
        legend1->AddEntry(h1_d[GL], Form("Had Interaction at L%d",GL), "el");
        legend1->AddEntry(h1_e[GL], Form("Inelastic Interaction at shallower Layer [L%d,L%d]",int(floor((GL-1)/2)),GL-1), "el");
        legend1->AddEntry(h1_H[GL], Form("Inelastic Interaction at shallower Layer [L%d,L%d]",0,int(floor((GL-1)/2))), "el");  
        legend1->AddEntry(h1_c[GL], "Elastic Interaction at shallower Layer", "el");  
        legend1->AddEntry(h1_h[GL], "Pass Through", "el");       
    }

    legend1->Draw();       
}