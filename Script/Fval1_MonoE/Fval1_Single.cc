void Fval1_Single()
{
    int p_First_Had_Layer;    std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Fval = nullptr;    std::vector<double>* p_Efrac = nullptr;
    int d_First_Had_Layer;    std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Fval = nullptr;    std::vector<double>* d_Efrac = nullptr;
    int e_First_Had_Layer;    std::vector<double>* e_RMSVec = nullptr;    std::vector<double>* e_EnergyVec = nullptr;    std::vector<double>* e_Fval = nullptr;    std::vector<double>* e_Efrac = nullptr;
    int h_First_Had_Layer;    std::vector<double>* h_RMSVec = nullptr;    std::vector<double>* h_EnergyVec = nullptr;    std::vector<double>* h_Fval = nullptr;    std::vector<double>* h_Efrac = nullptr;
    int H_First_Had_Layer;    std::vector<double>* H_RMSVec = nullptr;    std::vector<double>* H_EnergyVec = nullptr;    std::vector<double>* H_Fval = nullptr;    std::vector<double>* H_Efrac = nullptr;
    int c_First_Had_Layer;    std::vector<double>* c_RMSVec = nullptr;    std::vector<double>* c_EnergyVec = nullptr;    std::vector<double>* c_Fval = nullptr;    std::vector<double>* c_Efrac = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_100GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Fval"             ,&p_Fval);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_100GeV.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Fval"             ,&d_Fval);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    auto electron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Electron_100GeV.root");
    auto electron_tree = (TTree*)electron_file->Get("B4");
    electron_tree->SetBranchAddress("RMS"              ,&e_RMSVec);
    electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
    electron_tree->SetBranchAddress("Fval"             ,&e_Fval);
    electron_tree->SetBranchAddress("First_Had_Layer"  ,&e_First_Had_Layer);
    auto helium4_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium4_100GeV.root");
    auto helium4_tree = (TTree*)helium4_file->Get("B4");
    helium4_tree->SetBranchAddress("RMS"              ,&h_RMSVec);
    helium4_tree->SetBranchAddress("LayerEnergyVector",&h_EnergyVec);
    helium4_tree->SetBranchAddress("Fval"             ,&h_Fval);
    helium4_tree->SetBranchAddress("First_Had_Layer"  ,&h_First_Had_Layer);
    auto helium3_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium3_100GeV.root");
    auto helium3_tree = (TTree*)helium3_file->Get("B4");
    helium3_tree->SetBranchAddress("RMS"              ,&H_RMSVec);
    helium3_tree->SetBranchAddress("LayerEnergyVector",&H_EnergyVec);
    helium3_tree->SetBranchAddress("Fval"             ,&H_Fval);
    helium3_tree->SetBranchAddress("First_Had_Layer"  ,&H_First_Had_Layer);
    auto carbon_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Carbon_100GeV.root");
    auto carbon_tree = (TTree*)carbon_file->Get("B4");
    carbon_tree->SetBranchAddress("RMS"              ,&c_RMSVec);
    carbon_tree->SetBranchAddress("LayerEnergyVector",&c_EnergyVec);
    carbon_tree->SetBranchAddress("Fval"             ,&c_Fval);
    carbon_tree->SetBranchAddress("First_Had_Layer"  ,&c_First_Had_Layer);

    cout  << proton_tree->GetEntries() << endl;
    auto h1_p = new TH1D("h1_p","h1_p",100,0,10);  
    auto h1_d = new TH1D("h1_d","h1_d",100,0,10);  
    auto h1_e = new TH1D("h1_e","h1_e",100,0,10);  
    auto h1_h = new TH1D("h1_h","h1_h",100,0,10);  
    auto h1_H = new TH1D("h1_H","h1_H",100,0,10);  
    auto h1_c = new TH1D("h1_c","h1_c",100,0,10);  

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);    h1_p->Fill(log10(pow((*p_RMSVec)[4],3)));  // (*p_Fval)[4]*
        deuteron_tree->GetEntry(entry);  h1_d->Fill(log10(pow((*d_RMSVec)[4],3)));  // (*d_Fval)[4]*
        electron_tree->GetEntry(entry);  h1_e->Fill(log10(pow((*e_RMSVec)[4],3)));  // (*e_Fval)[4]*
        helium4_tree->GetEntry(entry);   h1_h->Fill(log10(pow((*h_RMSVec)[4],3)));  // (*h_Fval)[4]*
        helium3_tree->GetEntry(entry);   h1_H->Fill(log10(pow((*H_RMSVec)[4],3)));  // (*H_Fval)[4]*
        carbon_tree->GetEntry(entry);    h1_c->Fill(log10(pow((*c_RMSVec)[4],3)));  // (*c_Fval)[4]*
    }
    h1_p->Sumw2(); h1_p->Scale(1.0/h1_p->Integral());h1_p->SetLineColor(kRed);     h1_p->SetMarkerColor(kRed);     h1_p->SetLineWidth(2);h1_p->GetYaxis()->SetRangeUser(0,0.25);h1_p->SetTitle("100 GeV incident Particle; Fval ;Normalized Count");
    h1_d->Sumw2(); h1_d->Scale(1.0/h1_d->Integral());h1_d->SetLineColor(kBlue);    h1_d->SetMarkerColor(kBlue);    h1_d->SetLineWidth(2);
    h1_e->Sumw2(); h1_e->Scale(1.0/h1_e->Integral());h1_e->SetLineColor(kOrange-3);h1_e->SetMarkerColor(kOrange-3);h1_e->SetLineWidth(2);
    h1_h->Sumw2(); h1_h->Scale(1.0/h1_h->Integral());h1_h->SetLineColor(kGreen-3); h1_h->SetMarkerColor(kGreen-3); h1_h->SetLineWidth(2);
    h1_H->Sumw2(); h1_H->Scale(1.0/h1_H->Integral());h1_H->SetLineColor(kGreen-3); h1_H->SetMarkerColor(kGreen-3); h1_H->SetLineWidth(2);
    h1_c->Sumw2(); h1_c->Scale(1.0/h1_c->Integral());h1_c->SetLineColor(kMagenta); h1_c->SetMarkerColor(kGreen-3); h1_c->SetLineWidth(2);
    auto c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    c1->Clear();
    h1_p->Draw();
    h1_d->Draw("same");
    h1_e->Draw("same");
    h1_h->Draw("same");
    h1_H->Draw("same");
    h1_c->Draw("same");

    double quantiles[4] = {0.16, 0.50, 0.84};  // Percentiles
    double p_values[4];      h1_p->GetQuantiles(3, p_values, quantiles);        TLine *l_p[4];
    double d_values[4];      h1_d->GetQuantiles(3, d_values, quantiles);        TLine *l_d[4];
    double e_values[4];      h1_e->GetQuantiles(3, e_values, quantiles);        TLine *l_e[4];
    double h_values[4];      h1_h->GetQuantiles(3, h_values, quantiles);        TLine *l_h[4];
    double H_values[4];      h1_H->GetQuantiles(3, H_values, quantiles);        TLine *l_H[4];
    double c_values[4];      h1_c->GetQuantiles(3, c_values, quantiles);        TLine *l_c[4];

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
    legend1->AddEntry(h1_H, "Helium3", "el");         
    legend1->AddEntry(h1_c, "Carbon", "el");         

    legend1->Draw();       
   
}