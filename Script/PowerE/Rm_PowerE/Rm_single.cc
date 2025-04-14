void Rm_Single()
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

    cout  << proton_tree->GetEntries() << endl;
    auto h1_p = new TH1D("h1_p","h1_p",50,-5,0);      double p_maxVal;
    auto h1_d = new TH1D("h1_d","h1_d",50,-5,0);      double d_maxVal;

    int k = 0; // Layer
    int j = 10; // Energy Range log10
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1000; ++entry)
    {
        proton_tree->GetEntry(entry);   
        deuteron_tree->GetEntry(entry); 
        for (size_t i = 22*k; i < 22*(k+1); i += 22) // layer
        {  
            auto p_start = p_EnergyVec->begin() + i;  auto p_end = (i + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal = *std::max_element(p_start, p_end); 
            if(p_Total_E > pow(10,1+0.2*j) && p_Total_E <= pow(10,1.2+0.2*j) /* && p_First_Had_Layer==0 */) 
            { 
                // cout << i << " , " << p_Total_E << " , " << p_maxVal << " , " <<  log10(p_maxVal/p_Total_E) <<  endl;
                h1_p->Fill(log10(p_maxVal/p_Total_E)); 
            }
            auto d_start = d_EnergyVec->begin() + i;  auto d_end = (i + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end();  d_maxVal = *std::max_element(d_start, d_end); 
            if(d_Total_E > pow(10,1+0.2*j) && d_Total_E <= pow(10,1.2+0.2*j) /* && d_First_Had_Layer==0)*/) 
            { 
                h1_d->Fill(log10(d_maxVal/d_Total_E)); 
            }
        }
    }
    h1_p->Sumw2(); h1_p->Scale(1.0/h1_p->Integral()); h1_p->SetLineColor(kRed);     h1_p->SetMarkerColor(kRed);     h1_p->SetLineWidth(2);
    h1_d->Sumw2(); h1_d->Scale(1.0/h1_d->Integral()); h1_d->SetLineColor(kBlue);    h1_d->SetMarkerColor(kBlue);    h1_d->SetLineWidth(2);
    h1_p->GetYaxis()->SetRangeUser(0,0.5);h1_p->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV]; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L%d/ Total Deposit);Count", 1+0.2*j ,1.2+0.2*j ,k ));

    auto c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    c1->Clear();
    h1_p->Draw("hist");      
    h1_d->Draw("same");
    h1_p->SetStats(kTRUE); 
    h1_d->SetStats(kFALSE); 
    double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
    double p_values[3];      h1_p->GetQuantiles(3, p_values, quantiles);    TLine *l_p[3];
    double d_values[3];      h1_d->GetQuantiles(3, d_values, quantiles);    TLine *l_d[3];

    auto legend1 = new TLegend(0.12, 0.68, 0.28, 0.88);
    legend1->AddEntry(h1_p, "Proton", "el");
    legend1->AddEntry(h1_d, "Deuteron", "el");
    legend1->Draw();       
}