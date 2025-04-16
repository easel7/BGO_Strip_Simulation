void Rm_Single()
{
    int p_First_Had_Layer; double p_E_total;    std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr;
    int d_First_Had_Layer; double d_E_total;    std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr;
    int e_First_Had_Layer; double e_E_total;    std::vector<double>* e_RMSVec = nullptr;    std::vector<double>* e_EnergyVec = nullptr;    std::vector<double>* e_Efrac = nullptr;
    int h_First_Had_Layer; double h_E_total;    std::vector<double>* h_RMSVec = nullptr;    std::vector<double>* h_EnergyVec = nullptr;    std::vector<double>* h_Efrac = nullptr;
    int H_First_Had_Layer; double H_E_total;    std::vector<double>* H_RMSVec = nullptr;    std::vector<double>* H_EnergyVec = nullptr;    std::vector<double>* H_Efrac = nullptr;
    int c_First_Had_Layer; double c_E_total;    std::vector<double>* c_RMSVec = nullptr;    std::vector<double>* c_EnergyVec = nullptr;    std::vector<double>* c_Efrac = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("BarEnergyVector"  ,&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("Total_E"          ,&p_E_total);
    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_1000GeV.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_E_total);
    auto electron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Electron_1000GeV.root");
    auto electron_tree = (TTree*)electron_file->Get("B4");
    electron_tree->SetBranchAddress("RMS"              ,&e_RMSVec);
    electron_tree->SetBranchAddress("BarEnergyVector"  ,&e_EnergyVec);
    electron_tree->SetBranchAddress("Efrac"            ,&e_Efrac);
    electron_tree->SetBranchAddress("First_Had_Layer"  ,&e_First_Had_Layer);
    electron_tree->SetBranchAddress("Total_E"          ,&e_E_total);
    auto helium4_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium4_1000GeV.root");
    auto helium4_tree = (TTree*)helium4_file->Get("B4");
    helium4_tree->SetBranchAddress("RMS"              ,&h_RMSVec);
    helium4_tree->SetBranchAddress("BarEnergyVector"  ,&h_EnergyVec);
    helium4_tree->SetBranchAddress("Efrac"            ,&h_Efrac);
    helium4_tree->SetBranchAddress("First_Had_Layer"  ,&h_First_Had_Layer);
    helium4_tree->SetBranchAddress("Total_E"          ,&h_E_total);
    auto helium3_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Helium3_1000GeV.root");
    auto helium3_tree = (TTree*)helium3_file->Get("B4");
    helium3_tree->SetBranchAddress("RMS"              ,&H_RMSVec);
    helium3_tree->SetBranchAddress("BarEnergyVector"  ,&H_EnergyVec);
    helium3_tree->SetBranchAddress("Efrac"            ,&H_Efrac);
    helium3_tree->SetBranchAddress("First_Had_Layer"  ,&H_First_Had_Layer);
    helium3_tree->SetBranchAddress("Total_E"          ,&H_E_total);
    auto carbon_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Carbon_1000GeV.root");
    auto carbon_tree = (TTree*)carbon_file->Get("B4");
    carbon_tree->SetBranchAddress("RMS"              ,&c_RMSVec);
    carbon_tree->SetBranchAddress("BarEnergyVector"  ,&c_EnergyVec);
    carbon_tree->SetBranchAddress("Efrac"            ,&c_Efrac);
    carbon_tree->SetBranchAddress("First_Had_Layer"  ,&c_First_Had_Layer);
    carbon_tree->SetBranchAddress("Total_E"          ,&c_E_total);

    cout  << proton_tree->GetEntries() << endl;
    auto h1_p = new TH1D("h1_p","h1_p",50,-5,0);      double p_maxVal;
    auto h1_d = new TH1D("h1_d","h1_d",50,-5,0);      double d_maxVal;
    auto h1_e = new TH1D("h1_e","h1_e",50,-5,0);      double e_maxVal;
    auto h1_h = new TH1D("h1_h","h1_h",50,-5,0);      double h_maxVal;
    auto h1_H = new TH1D("h1_H","h1_H",50,-5,0);      double H_maxVal;
    auto h1_c = new TH1D("h1_c","h1_c",50,-5,0);      double c_maxVal;

    int k = 0;
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        proton_tree->GetEntry(entry);    
        deuteron_tree->GetEntry(entry); 
        electron_tree->GetEntry(entry); 
        helium4_tree->GetEntry(entry);  
        helium3_tree->GetEntry(entry);  
        carbon_tree->GetEntry(entry);   

        for (size_t i = 22*k; i < 22*(k+1); i += 22)
        // for (size_t i = 0; i < 22; i += 22)
        // for (size_t i = 0; i < p_EnergyVec->size(); i += 22)
        {
            auto p_start = p_EnergyVec->begin() + i;  auto p_end = (i + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal = *std::max_element(p_start, p_end); if((*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 ) /*&& (*p_RMSVec)[1] < 45 && (*p_RMSVec)[2] < 45) /*if(p_First_Had_Layer==0)*/{ h1_p->Fill(log10(p_maxVal/p_E_total)); }
            auto d_start = d_EnergyVec->begin() + i;  auto d_end = (i + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end();  d_maxVal = *std::max_element(d_start, d_end); if((*d_RMSVec)[0]>15 && (*d_RMSVec)[1]>15 ) /*&& (*d_RMSVec)[1] < 45 && (*d_RMSVec)[2] < 45) /*if(d_First_Had_Layer==0)*/{ h1_d->Fill(log10(d_maxVal/d_E_total)); }
            auto e_start = e_EnergyVec->begin() + i;  auto e_end = (i + 22 < e_EnergyVec->size() ) ? e_start + 22 : e_EnergyVec->end();  e_maxVal = *std::max_element(e_start, e_end); if((*e_RMSVec)[0]>15 && (*e_RMSVec)[1]>15 ) /*&& (*e_RMSVec)[1] < 45 && (*e_RMSVec)[2] < 45) /*if(e_First_Had_Layer==0)*/{ h1_e->Fill(log10(e_maxVal/e_E_total)); }
            auto h_start = h_EnergyVec->begin() + i;  auto h_end = (i + 22 < h_EnergyVec->size() ) ? h_start + 22 : h_EnergyVec->end();  h_maxVal = *std::max_element(h_start, h_end); if((*h_RMSVec)[0]>15 && (*h_RMSVec)[1]>15 ) /*&& (*h_RMSVec)[1] < 45 && (*h_RMSVec)[2] < 45) /*if(h_First_Had_Layer==0)*/{ h1_h->Fill(log10(h_maxVal/h_E_total)); }
            auto H_start = H_EnergyVec->begin() + i;  auto H_end = (i + 22 < H_EnergyVec->size() ) ? H_start + 22 : H_EnergyVec->end();  H_maxVal = *std::max_element(H_start, H_end); if((*H_RMSVec)[0]>15 && (*H_RMSVec)[1]>15 ) /*&& (*H_RMSVec)[1] < 45 && (*H_RMSVec)[2] < 45) /*if(H_First_Had_Layer==0)*/{ h1_H->Fill(log10(H_maxVal/H_E_total)); }
            auto c_start = c_EnergyVec->begin() + i;  auto c_end = (i + 22 < c_EnergyVec->size() ) ? c_start + 22 : c_EnergyVec->end();  c_maxVal = *std::max_element(c_start, c_end); if((*c_RMSVec)[0]>15 && (*c_RMSVec)[1]>15 ) /*&& (*c_RMSVec)[1] < 45 && (*c_RMSVec)[2] < 45) /*if(c_First_Had_Layer==0)*/{ h1_c->Fill(log10(c_maxVal/c_E_total)); }
            
            // cout << " bar " << i << " Layer " << int(i/22) <<endl; 
            // std::copy(p_start, p_end, std::ostream_iterator<double>(std::cout, " "));
            // cout << " max = "<< p_maxVal <<  endl;
        }
    }

    h1_p->Sumw2(); /*h1_p->Scale(1.0/h1_p->Integral());*/ h1_p->SetLineColor(kRed);     h1_p->SetMarkerColor(kRed);     h1_p->SetLineWidth(2);
    h1_d->Sumw2(); /*h1_d->Scale(1.0/h1_d->Integral());*/ h1_d->SetLineColor(kBlue);    h1_d->SetMarkerColor(kBlue);    h1_d->SetLineWidth(2);
    h1_e->Sumw2(); /*h1_e->Scale(1.0/h1_e->Integral());*/ h1_e->SetLineColor(kOrange-3);h1_e->SetMarkerColor(kOrange-3);h1_e->SetLineWidth(2);
    h1_h->Sumw2(); /*h1_h->Scale(1.0/h1_h->Integral());*/ h1_h->SetLineColor(kGreen-3); h1_h->SetMarkerColor(kGreen-3); h1_h->SetLineWidth(2);
    h1_H->Sumw2(); /*h1_H->Scale(1.0/h1_H->Integral());*/ h1_H->SetLineColor(kGreen-3); h1_H->SetMarkerColor(kGreen-3); h1_H->SetLineWidth(2);
    h1_c->Sumw2(); /*h1_c->Scale(1.0/h1_c->Integral());*/ h1_c->SetLineColor(kMagenta); h1_c->SetMarkerColor(kGreen-3); h1_c->SetLineWidth(2);
    // h1_p->GetYaxis()->SetRangeUser(0,0.5);h1_p->SetTitle("1000 GeV incident; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L0/ Total Deposit);Normalized Count");
    h1_p->GetYaxis()->SetRangeUser(0,2000);h1_p->SetTitle(Form("1000 GeV incident; log_{10}(Rm) = log_{10}(Max Energy Deposit bar in L%d/ Total Deposit);Count",k));
    
    auto c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    c1->Clear();
    h1_p->Draw("hist");      h1_p->SetStats(kTRUE); 
    h1_d->Draw("histsame");h1_d->SetStats(kFALSE); 
    h1_e->Draw("same");h1_e->SetStats(kFALSE); 
    h1_h->Draw("same");h1_h->SetStats(kFALSE); 
    // h1_H->Draw("same");h1_H->SetStats(kFALSE); 
    h1_c->Draw("same");h1_c->SetStats(kFALSE); 

    auto legend1 = new TLegend(0.12, 0.68, 0.28, 0.88);
    legend1->AddEntry(h1_p, "Proton", "el");
    legend1->AddEntry(h1_d, "Deuteron", "el");
    legend1->AddEntry(h1_e, "Electron", "el");
    legend1->AddEntry(h1_h, "Helium4", "el");         
    // legend1->AddEntry(h1_H, "Helium3", "el");         
    legend1->AddEntry(h1_c, "Carbon", "el");         
    legend1->Draw();       
}