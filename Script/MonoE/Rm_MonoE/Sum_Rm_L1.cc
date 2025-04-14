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
    auto h1_p = new TH1D("h1_p","h1_p",50,-5,0);  
    auto h1_d = new TH1D("h1_d","h1_d",50,-5,0);  
    auto h1_e = new TH1D("h1_e","h1_e",50,-5,0);  
    auto h1_E = new TH1D("h1_E","h1_E",50,-5,0);  
    auto h1_h = new TH1D("h1_h","h1_h",50,-5,0);  
    auto h1_H = new TH1D("h1_H","h1_H",50,-5,0);  
    auto h1_c = new TH1D("h1_c","h1_c",50,-5,0);  


    for (int j = 0; j < 14; j++)
    {
    }

    int k = 3; // Goal Layer
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        proton_tree->GetEntry(entry);    
        double sum_p = 0;
        double sum_p_total = 0;
        double p_maxVal[14] = {0};
        for (size_t j = 0; j < p_EnergyVec->size(); j += 22)
        {
            int k = int(j / 22);  // Get the Layer
            auto p_start = p_EnergyVec->begin() + j;  auto p_end = (j + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal[k] = *std::max_element(p_start, p_end); 
            if((*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15) //  && (*p_RMSVec)[0]<40 && (*p_RMSVec)[1]<40)
            {
                sum_p       += p_maxVal[k];
                sum_p_total += (*p_L_EnergyVec)[k];

                if(k>0 && k<13)
                {
                    h1_p->Fill(log10(sum_p / sum_p_total));
                    if(p_First_Had_Layer==k)                                                                              { h1_d->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer==-1)                                                                             { h1_h->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=floor((k-1)/2) && p_First_Had_Layer<=k-1              && p_First_Had_Type == 1) { h1_e->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=0              && p_First_Had_Layer<floor((k-1)/2)   && p_First_Had_Type == 1)  { h1_H->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer<=k-1 && p_First_Had_Type == 2)                                                   { h1_c->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=k+1)                                                                            { h1_E->Fill(log10(sum_p / sum_p_total)); }
                    
                }
                else if (k==0)
                {
                    h1_p->Fill(log10(sum_p / sum_p_total));
                    if(p_First_Had_Layer==k)                                                                              { h1_d->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer==-1)                                                                             { h1_h->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=k+1 && p_First_Had_Type == 1)                                                   { h1_e->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=k+1 && p_First_Had_Type == 2)                                                   { h1_E->Fill(log10(sum_p / sum_p_total)); }
                }
                else if (k==13)
                {
                    h1_p->Fill(log10(sum_p / sum_p_total));
                    if(p_First_Had_Layer==k)                                                                              { h1_d->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer==-1)                                                                             { h1_h->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=floor((k-1)/2) && p_First_Had_Layer<=k-1             && p_First_Had_Type == 1)  { h1_e->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer>=0              && p_First_Had_Layer<floor((k-1)/2)   && p_First_Had_Type == 1)  { h1_H->Fill(log10(sum_p / sum_p_total)); }
                    if(p_First_Had_Layer<=k-1 && p_First_Had_Type == 2)                                                   { h1_c->Fill(log10(sum_p / sum_p_total)); }
                }
                
            }   
        }
    }


    h1_p->Sumw2(); /*h1_p->Scale(1.0/h1_p->Integral());*/ h1_p->SetLineColor(kRed);     h1_p->SetMarkerColor(kRed);     h1_p->SetLineWidth(2);
    h1_d->Sumw2(); /*h1_d->Scale(1.0/h1_d->Integral());*/ h1_d->SetLineColor(kBlue);    h1_d->SetMarkerColor(kBlue);    h1_d->SetLineWidth(2);
    h1_e->Sumw2(); /*h1_e->Scale(1.0/h1_e->Integral());*/ h1_e->SetLineColor(kOrange-3);h1_e->SetMarkerColor(kOrange-3);h1_e->SetLineWidth(2);
    h1_h->Sumw2(); /*h1_h->Scale(1.0/h1_h->Integral());*/ h1_h->SetLineColor(kGreen-3); h1_h->SetMarkerColor(kGreen-3); h1_h->SetLineWidth(2);
    h1_c->Sumw2(); /*h1_c->Scale(1.0/h1_c->Integral());*/ h1_c->SetLineColor(kMagenta); h1_c->SetMarkerColor(kMagenta); h1_c->SetLineWidth(2);
    h1_H->Sumw2(); /*h1_H->Scale(1.0/h1_H->Integral());*/ h1_H->SetLineColor(kBlack);   h1_H->SetMarkerColor(kBlack);   h1_H->SetLineWidth(2);
    h1_E->Sumw2(); /*h1_H->Scale(1.0/h1_H->Integral());*/ h1_E->SetLineColor(kCyan);    h1_E->SetMarkerColor(kCyan);     h1_E->SetLineWidth(2);

   
    // h1_p->GetYaxis()->SetRangeUser(0,0.5);h1_p->SetTitle("1000 GeV incident; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L0/ Total Deposit);Normalized Count");
    h1_p->GetYaxis()->SetRangeUser(0,300);h1_p->SetTitle(Form("1000 GeV incident; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L%d/ Total Deposit);Count",k));
    
    auto c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    c1->Clear();
    h1_p->Draw("hist");h1_p->SetStats(kFALSE); 
    h1_d->Draw("same");h1_d->SetStats(kFALSE); 
    h1_e->Draw("same");h1_e->SetStats(kFALSE); 
    h1_h->Draw("same");h1_h->SetStats(kFALSE); 
    h1_c->Draw("same");h1_c->SetStats(kFALSE); 

    auto legend1 = new TLegend(0.12, 0.58, 0.68, 0.88);

    if(k>0 && k<13)
    {
        h1_H->Draw("same");h1_H->SetStats(kFALSE); 
        h1_E->Draw("same");h1_c->SetStats(kFALSE); 

        legend1->AddEntry(h1_p, "All", "el");
        legend1->AddEntry(h1_d, Form("Had Interaction at L%d",k), "el");
        legend1->AddEntry(h1_e, Form("Inelastic Interaction at shallower Layer [L%d,L%d]",int(floor((k-1)/2)),k-1), "el");
        legend1->AddEntry(h1_H, Form("Inelastic Interaction at shallower Layer [L%d,L%d]",0,int(floor((k-1)/2))), "el");  
        legend1->AddEntry(h1_c, "Elastic Interaction at shallower Layer", "el");  
        legend1->AddEntry(h1_E, "Had Interaction at deeper Layer", "el");  
        legend1->AddEntry(h1_h, "Pass Through", "el");       
    }
    else if (k==0)
    {
        legend1->AddEntry(h1_p, "All", "el");
        legend1->AddEntry(h1_d, Form("Had Interaction at L%d",k), "el");
        legend1->AddEntry(h1_e, "Inelastic Interaction at deeper Layer", "el");
        legend1->AddEntry(h1_c, "Elastic Interaction at deeper Layer", "el");  
        legend1->AddEntry(h1_h, "Pass Through", "el");    
    }
    else if (k==13)
    {
        h1_H->Draw("same");h1_H->SetStats(kFALSE); 
        legend1->AddEntry(h1_p, "All", "el");
        legend1->AddEntry(h1_d, Form("Had Interaction at L%d",k), "el");
        legend1->AddEntry(h1_e, Form("Inelastic Interaction at shallower Layer [L%d,L%d]",int(floor((k-1)/2)),k-1), "el");
        legend1->AddEntry(h1_H, Form("Inelastic Interaction at shallower Layer [L%d,L%d]",0,int(floor((k-1)/2))), "el");  
        legend1->AddEntry(h1_c, "Elastic Interaction at shallower Layer", "el");  
        legend1->AddEntry(h1_h, "Pass Through", "el");       
    }


    legend1->Draw();       
}