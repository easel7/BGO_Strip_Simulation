void Jm_Single()
{
    double p_E_total;    std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;
    double d_E_total;    std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;
    double e_E_total;    std::vector<double>* e_RMSVec = nullptr;    std::vector<double>* e_EnergyVec = nullptr;    std::vector<double>* e_L_EnergyVec = nullptr;
    double h_E_total;    std::vector<double>* h_RMSVec = nullptr;    std::vector<double>* h_EnergyVec = nullptr;    std::vector<double>* h_L_EnergyVec = nullptr;
    double H_E_total;    std::vector<double>* H_RMSVec = nullptr;    std::vector<double>* H_EnergyVec = nullptr;    std::vector<double>* H_L_EnergyVec = nullptr;
    double c_E_total;    std::vector<double>* c_RMSVec = nullptr;    std::vector<double>* c_EnergyVec = nullptr;    std::vector<double>* c_L_EnergyVec = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_10GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("BarEnergyVector"  ,&p_EnergyVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
    proton_tree->SetBranchAddress("Total_E"          ,&p_E_total);
    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_10GeV.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_E_total);
    auto electron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Electron_10GeV.root");
    auto electron_tree = (TTree*)electron_file->Get("B4");
    electron_tree->SetBranchAddress("RMS"              ,&e_RMSVec);
    electron_tree->SetBranchAddress("BarEnergyVector"  ,&e_EnergyVec);
    electron_tree->SetBranchAddress("LayerEnergyVector",&e_L_EnergyVec);
    electron_tree->SetBranchAddress("Total_E"          ,&e_E_total);
    auto helium4_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium4_10GeV.root");
    auto helium4_tree = (TTree*)helium4_file->Get("B4");
    helium4_tree->SetBranchAddress("RMS"              ,&h_RMSVec);
    helium4_tree->SetBranchAddress("BarEnergyVector"  ,&h_EnergyVec);
    helium4_tree->SetBranchAddress("LayerEnergyVector",&h_L_EnergyVec);
    helium4_tree->SetBranchAddress("Total_E"          ,&h_E_total);
    auto helium3_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium3_10GeV.root");
    auto helium3_tree = (TTree*)helium3_file->Get("B4");
    helium3_tree->SetBranchAddress("RMS"              ,&H_RMSVec);
    helium3_tree->SetBranchAddress("BarEnergyVector"  ,&H_EnergyVec);
    helium3_tree->SetBranchAddress("LayerEnergyVector",&H_L_EnergyVec);
    helium3_tree->SetBranchAddress("Total_E"          ,&H_E_total);
    auto carbon_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Carbon_10GeV.root");
    auto carbon_tree = (TTree*)carbon_file->Get("B4");
    carbon_tree->SetBranchAddress("RMS"              ,&c_RMSVec);
    carbon_tree->SetBranchAddress("BarEnergyVector"  ,&c_EnergyVec);
    carbon_tree->SetBranchAddress("LayerEnergyVector",&c_L_EnergyVec);
    carbon_tree->SetBranchAddress("Total_E"          ,&c_E_total);

    cout  << proton_tree->GetEntries() << endl;
    auto h1_p = new TH1D("h1_p","h1_p",40,-4,0);  
    auto h1_d = new TH1D("h1_d","h1_d",40,-4,0);  
    auto h1_e = new TH1D("h1_e","h1_e",40,-4,0);  
    auto h1_h = new TH1D("h1_h","h1_h",40,-4,0);  
    auto h1_H = new TH1D("h1_H","h1_H",40,-4,0);  
    auto h1_c = new TH1D("h1_c","h1_c",40,-4,0);  

    double p_maxVal;
    double d_maxVal;
    double e_maxVal;
    double h_maxVal;
    double H_maxVal;
    double c_maxVal;

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        proton_tree->GetEntry(entry);    
        deuteron_tree->GetEntry(entry); 
        electron_tree->GetEntry(entry); 
        helium4_tree->GetEntry(entry);  
        helium3_tree->GetEntry(entry);  
        carbon_tree->GetEntry(entry);   
        if (!p_EnergyVec || p_EnergyVec->empty()) cout << "branch fail"<< endl;
        if (!d_EnergyVec || d_EnergyVec->empty()) cout << "branch fail"<< endl;
        if (!e_EnergyVec || e_EnergyVec->empty()) cout << "branch fail"<< endl;
        if (!h_EnergyVec || h_EnergyVec->empty()) cout << "branch fail"<< endl;
        if (!H_EnergyVec || H_EnergyVec->empty()) cout << "branch fail"<< endl;
        if (!c_EnergyVec || c_EnergyVec->empty()) cout << "branch fail"<< endl;
        cout << p_E_total << endl;
        // for (size_t i = 0; i < 22; i += 22)
        for (size_t i = 0; i < p_EnergyVec->size(); i += 22)
        {
            auto p_start = p_EnergyVec->begin() + i;  auto p_end = (i + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal = *std::max_element(p_start, p_end);
            auto d_start = d_EnergyVec->begin() + i;  auto d_end = (i + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end();  d_maxVal = *std::max_element(d_start, d_end);
            auto e_start = e_EnergyVec->begin() + i;  auto e_end = (i + 22 < e_EnergyVec->size() ) ? e_start + 22 : e_EnergyVec->end();  e_maxVal = *std::max_element(e_start, e_end);
            auto h_start = h_EnergyVec->begin() + i;  auto h_end = (i + 22 < h_EnergyVec->size() ) ? h_start + 22 : h_EnergyVec->end();  h_maxVal = *std::max_element(h_start, h_end);
            auto H_start = H_EnergyVec->begin() + i;  auto H_end = (i + 22 < H_EnergyVec->size() ) ? H_start + 22 : H_EnergyVec->end();  H_maxVal = *std::max_element(H_start, H_end);
            auto c_start = c_EnergyVec->begin() + i;  auto c_end = (i + 22 < c_EnergyVec->size() ) ? c_start + 22 : c_EnergyVec->end();  c_maxVal = *std::max_element(c_start, c_end);
            
            // cout << " bar " << i << " Layer " << int(i/22) <<endl; 
            // std::copy(p_start, p_end, std::ostream_iterator<double>(std::cout, " "));
            // cout << " max = "<< p_maxVal <<  endl;
        }
        h1_p->Fill(log10(p_maxVal/(*p_L_EnergyVec)[0]));
        h1_d->Fill(log10(d_maxVal/(*p_L_EnergyVec)[0]));
        h1_e->Fill(log10(e_maxVal/(*p_L_EnergyVec)[0]));
        h1_h->Fill(log10(h_maxVal/(*p_L_EnergyVec)[0]));
        h1_H->Fill(log10(H_maxVal/(*p_L_EnergyVec)[0]));
        h1_c->Fill(log10(c_maxVal/(*p_L_EnergyVec)[0]));
    }

    
    h1_p->Sumw2(); h1_p->Scale(1.0/h1_p->Integral()); h1_p->SetLineColor(kRed);     h1_p->SetMarkerColor(kRed);     h1_p->SetLineWidth(2);h1_p->GetYaxis()->SetRangeUser(0,0.4);h1_p->SetTitle("100 GeV incident EdepRatio Distrubution in L0; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L0/ Energy Deposit in L0);Normalized Count");
    h1_d->Sumw2(); h1_d->Scale(1.0/h1_d->Integral()); h1_d->SetLineColor(kBlue);    h1_d->SetMarkerColor(kBlue);    h1_d->SetLineWidth(2);
    h1_e->Sumw2(); h1_e->Scale(1.0/h1_e->Integral()); h1_e->SetLineColor(kOrange-3);h1_e->SetMarkerColor(kOrange-3);h1_e->SetLineWidth(2);
    h1_h->Sumw2(); h1_h->Scale(1.0/h1_h->Integral()); h1_h->SetLineColor(kGreen-3); h1_h->SetMarkerColor(kGreen-3); h1_h->SetLineWidth(2);
    h1_H->Sumw2(); h1_H->Scale(1.0/h1_H->Integral()); h1_H->SetLineColor(kGreen-3); h1_H->SetMarkerColor(kGreen-3); h1_H->SetLineWidth(2);
    h1_c->Sumw2(); h1_c->Scale(1.0/h1_c->Integral()); h1_c->SetLineColor(kMagenta); h1_c->SetMarkerColor(kGreen-3); h1_c->SetLineWidth(2);
    
    auto c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    c1->Clear();
    h1_p->Draw("hist");      h1_p->SetStats(kTRUE); 
    h1_d->Draw("histsame");h1_d->SetStats(kFALSE); 
    h1_e->Draw("same");h1_e->SetStats(kFALSE); 
    h1_h->Draw("same");h1_h->SetStats(kFALSE); 
    // h1_H->Draw("same");h1_H->SetStats(kFALSE); 
    h1_c->Draw("same");h1_c->SetStats(kFALSE); 

    double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
    double p_values[3];      h1_p->GetQuantiles(3, p_values, quantiles);    TLine *l_p[3];
    double d_values[3];      h1_d->GetQuantiles(3, d_values, quantiles);    TLine *l_d[3];
    double e_values[3];      h1_e->GetQuantiles(3, e_values, quantiles);    TLine *l_e[3];
    double h_values[3];      h1_h->GetQuantiles(3, h_values, quantiles);    TLine *l_h[3];
    double H_values[3];      h1_H->GetQuantiles(3, H_values, quantiles);    TLine *l_H[3];
    double c_values[3];      h1_c->GetQuantiles(3, c_values, quantiles);    TLine *l_c[3];

    for (int ii = 0 ;ii< 3 ; ii++)
    {
        l_p[ii] = new TLine(p_values[ii],0,p_values[ii],h1_p->GetBinContent(h1_p->FindBin(p_values[ii])));l_p[ii]->SetLineColor(kRed);     l_p[ii]->SetLineWidth(2);l_p[ii]->Draw();
        l_d[ii] = new TLine(d_values[ii],0,d_values[ii],h1_d->GetBinContent(h1_d->FindBin(d_values[ii])));l_d[ii]->SetLineColor(kBlue);    l_d[ii]->SetLineWidth(2);l_d[ii]->Draw();
        l_e[ii] = new TLine(e_values[ii],0,e_values[ii],h1_e->GetBinContent(h1_e->FindBin(e_values[ii])));l_e[ii]->SetLineColor(kOrange-3);l_e[ii]->SetLineWidth(2);l_e[ii]->Draw();
        l_h[ii] = new TLine(h_values[ii],0,h_values[ii],h1_h->GetBinContent(h1_h->FindBin(h_values[ii])));l_h[ii]->SetLineColor(kGreen-3); l_h[ii]->SetLineWidth(2);l_h[ii]->Draw();
        l_H[ii] = new TLine(H_values[ii],0,H_values[ii],h1_H->GetBinContent(h1_H->FindBin(H_values[ii])));l_H[ii]->SetLineColor(kGreen-3); l_H[ii]->SetLineWidth(2);l_H[ii]->Draw();l_H[ii]->SetLineStyle(2);
        l_c[ii] = new TLine(c_values[ii],0,c_values[ii],h1_c->GetBinContent(h1_c->FindBin(c_values[ii])));l_c[ii]->SetLineColor(kMagenta); l_c[ii]->SetLineWidth(2);l_c[ii]->Draw();
    }

    auto legend1 = new TLegend(0.12, 0.68, 0.28, 0.88);
    legend1->AddEntry(h1_p, "Proton", "el");
    legend1->AddEntry(h1_d, "Deuteron", "el");
    legend1->AddEntry(h1_e, "Electron", "el");
    legend1->AddEntry(h1_h, "Helium4", "el");         
    // legend1->AddEntry(h1_H, "Helium3", "el");         
    legend1->AddEntry(h1_c, "Carbon", "el");         
    legend1->Draw();       
}