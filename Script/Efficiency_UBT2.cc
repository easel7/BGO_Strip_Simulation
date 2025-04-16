void CalculateEfficiency(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, double Eff_HET, double Eff_UBT, double Eff_MIT)
{
    Eff_HET = h2->Integral() / h1->Integral();
    Eff_UBT = h3->Integral() / h1->Integral();
    Eff_MIT = h4->Integral() / h1->Integral();
}

TGraphErrors* CreateGraphWithProperties(int n, double* x, double* y, double* xErr, double* yErr, int markerStyle, int markerColor, int lineColor, const char* title)
{
    auto graph = new TGraphErrors(n, x, y, xErr, yErr);
    graph->SetTitle(title);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    return graph;
}

void Efficiency_UBT2()
{
    double Proton_Eff_HET[15]={0};      double Deuteron_Eff_HET[15]={0};
    double Proton_Eff_UBT[15]={0};      double Deuteron_Eff_UBT[15]={0};
    double Proton_Eff_MIT[15]={0};      double Deuteron_Eff_MIT[15]={0};
    double Proton_Eff_RMS[15]={0};      double Deuteron_Eff_RMS[15]={0};

    int p_First_Had_Layer;  double p_Total_E;  double p_Energy;  std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_First_Had_Layer;  double d_Total_E;  double d_Energy;  std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr; double d_weight;
    double Energy[15]={0};
    double Energy_Err[15]={0};
    double Uncertainty[15]={0};
    double n_BGO = TMath::Na()*7.13/1245.8344; // cm-3
    TH1D *h1_p[15]; TH1D *h2_p[15]; TH1D *h3_p[15]; TH1D *h4_p[15]; TH1D *h5_p[15];
    TH1D *h1_d[15]; TH1D *h2_d[15]; TH1D *h3_d[15]; TH1D *h4_d[15]; TH1D *h5_d[15];

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);
    proton_tree->SetBranchAddress("weight"           ,&p_weight);
    proton_tree->SetBranchAddress("Energy"           ,&p_Energy);


    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    deuteron_tree->SetBranchAddress("Energy"           ,&d_Energy);

    int nbins = 15;
    double xmin = 1e1;
    double xmax = 1e4;

    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    std::vector<double> bin_edges(nbins + 1);
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * (logxmax - logxmin) / nbins);
    }

    auto H1_p = new TH1D("H1_p","H1_p",nbins, bin_edges.data());
    auto H2_p = new TH1D("H2_p","H2_p",nbins, bin_edges.data());
    auto H3_p = new TH1D("H3_p","H3_p",nbins, bin_edges.data());
    auto H4_p = new TH1D("H4_p","H4_p",nbins, bin_edges.data());


    auto H1_d = new TH1D("H1_d","H1_d",nbins, bin_edges.data());
    auto H2_d = new TH1D("H2_d","H2_d",nbins, bin_edges.data());
    auto H3_d = new TH1D("H3_d","H3_d",nbins, bin_edges.data());
    auto H4_d = new TH1D("H4_d","H4_d",nbins, bin_edges.data());



    for (int i = 0; i < 15; i++)
    {
        Energy[i]     = 1.1 + 0.2 * i;
        Energy_Err[i] = 0.1;
        h1_p[i] = new TH1D(Form("h1_p[%d]",i),Form("h1_p[%d]",i),100,0,1);              
        h2_p[i] = new TH1D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),100,0,1);              
        h3_p[i] = new TH1D(Form("h3_p[%d]",i),Form("h3_p[%d]",i),100,0,1);              
        h4_p[i] = new TH1D(Form("h4_p[%d]",i),Form("h4_p[%d]",i),100,0,1);              
        h5_p[i] = new TH1D(Form("h5_p[%d]",i),Form("h5_p[%d]",i),100,0,1);            
        h1_d[i] = new TH1D(Form("h1_d[%d]",i),Form("h1_d[%d]",i),100,0,1); 
        h2_d[i] = new TH1D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),100,0,1); 
        h3_d[i] = new TH1D(Form("h3_d[%d]",i),Form("h3_d[%d]",i),100,0,1); 
        h4_d[i] = new TH1D(Form("h4_d[%d]",i),Form("h4_d[%d]",i),100,0,1); 
        h5_d[i] = new TH1D(Form("h5_d[%d]",i),Form("h5_d[%d]",i),100,0,1);  
    }

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);   
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        h1_p[p_energy_index]->Fill(p_Total_E/p_Energy); 
        if ( (*p_EnergyVec)[0] > 0.0092 && (*p_EnergyVec)[1]  > 0.0092)    h3_p[p_energy_index]->Fill(p_Total_E/p_Energy); 
        if ( (*p_EnergyVec)[0] > 0.23   && (*p_EnergyVec)[1]  > 0.23   && (*p_EnergyVec)[2]  > 0.23     &&  (*p_EnergyVec)[3] > 0.046) h2_p[p_energy_index]->Fill(p_Total_E/p_Energy); 
        if (((*p_EnergyVec)[2] > 0.0092 && (*p_EnergyVec)[10] > 0.0092 && (*p_EnergyVec)[12] > 0.0092 ) || ((*p_EnergyVec)[3] > 0.0092 && (*p_EnergyVec)[11] > 0.0092 && (*p_EnergyVec)[13] > 0.0092 )) h4_p[p_energy_index]->Fill(p_Total_E/p_Energy);
        // if((*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 )                                       h5_p[p_energy_index]->Fill(p_Total_E/p_Energy);
        if((*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 && (*p_RMSVec)[2]<45 && (*p_RMSVec)[3]<45) h5_p[p_energy_index]->Fill(p_Total_E/p_Energy);
    }
           
    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
            deuteron_tree->GetEntry(entry); 
            int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
            if(d_energy_index < 0 || d_energy_index > 14) continue;
            h1_d[d_energy_index]->Fill(d_Total_E/d_Energy); 
            if ( (*d_EnergyVec)[0] > 0.0092 && (*d_EnergyVec)[1] > 0.0092) h3_d[d_energy_index]->Fill(d_Total_E/d_Energy); 
            if ( (*d_EnergyVec)[0] > 0.23   && (*d_EnergyVec)[1] > 0.23 && (*d_EnergyVec)[2] > 0.23 && (*d_EnergyVec)[3] > 0.046) h2_d[d_energy_index]->Fill(d_Total_E/d_Energy); 
            if (((*d_EnergyVec)[2] > 0.0092 && (*d_EnergyVec)[10] > 0.0092 && (*d_EnergyVec)[12] > 0.0092 ) || ((*d_EnergyVec)[3] > 0.0092 && (*d_EnergyVec)[11] > 0.0092 && (*d_EnergyVec)[13] > 0.0092 )) h4_d[d_energy_index]->Fill(d_Total_E/d_Energy);
            // if ((*d_RMSVec)[0]>15 && (*d_RMSVec)[1]>15 ) h5_d[d_energy_index]->Fill(d_Total_E/d_Energy);
            if((*d_RMSVec)[0]>15 && (*d_RMSVec)[1]>15 && (*d_RMSVec)[2]<45 && (*d_RMSVec)[3]<45) h5_d[d_energy_index]->Fill(d_Total_E/d_Energy);
    }


    for (int i = 0; i < 15; i++)
    {
        Proton_Eff_HET[i]   = h2_p[i]->Integral()/h3_p[i]->Integral();
        Proton_Eff_UBT[i]   = h3_p[i]->Integral()/h1_p[i]->Integral();
        Proton_Eff_MIT[i]   = h4_p[i]->Integral()/h1_p[i]->Integral();
        Proton_Eff_RMS[i]   = h5_p[i]->Integral()/h1_p[i]->Integral();

        Deuteron_Eff_HET[i] = h2_d[i]->Integral()/h3_d[i]->Integral();
        Deuteron_Eff_UBT[i] = h3_d[i]->Integral()/h1_d[i]->Integral();
        Deuteron_Eff_MIT[i] = h4_d[i]->Integral()/h1_d[i]->Integral();
        Deuteron_Eff_RMS[i] = h5_d[i]->Integral()/h1_d[i]->Integral();

        H1_p->SetBinContent(i+1,Proton_Eff_HET[i]);         H1_p->SetBinError(i+1,sqrt(h2_p[i]->Integral())/h3_p[i]->Integral());
        H2_p->SetBinContent(i+1,Proton_Eff_UBT[i]);         H2_p->SetBinError(i+1,sqrt(h3_p[i]->Integral())/h1_p[i]->Integral());
        H3_p->SetBinContent(i+1,Proton_Eff_MIT[i]);         H3_p->SetBinError(i+1,sqrt(h4_p[i]->Integral())/h1_p[i]->Integral());
        H4_p->SetBinContent(i+1,Proton_Eff_RMS[i]);         H4_p->SetBinError(i+1,sqrt(h5_p[i]->Integral())/h1_p[i]->Integral());

        H1_d->SetBinContent(i+1,Deuteron_Eff_HET[i]);       H1_d->SetBinError(i+1,sqrt(h2_d[i]->Integral())/h3_d[i]->Integral());
        H2_d->SetBinContent(i+1,Deuteron_Eff_UBT[i]);       H2_d->SetBinError(i+1,sqrt(h3_d[i]->Integral())/h1_d[i]->Integral());
        H3_d->SetBinContent(i+1,Deuteron_Eff_MIT[i]);       H3_d->SetBinError(i+1,sqrt(h4_d[i]->Integral())/h1_d[i]->Integral());
        H4_d->SetBinContent(i+1,Deuteron_Eff_RMS[i]);       H4_d->SetBinError(i+1,sqrt(h5_d[i]->Integral())/h1_d[i]->Integral());


        cout << "Energy = " << Energy[i] << " GeV !  ele Total eff : " << h1_p[i]->Integral() << " ele HET eff : " <<  h2_p[i]->Integral() <<" ele UBT eff : " << h3_p[i]->Integral() <<" ele MIP eff : " << h4_p[i]->Integral() <<endl;
    }

    H1_p->SetLineColor(kRed); H1_p->SetLineWidth(2);  H1_p->SetTitle("High-Energy Trigger; BGO Deposit Energy (GeV); Trigger Efficiency");
    H2_p->SetLineColor(kRed); H2_p->SetLineWidth(2);  H2_p->SetTitle("Unbiased Trigger;BGO Deposit Energy (GeV); Trigger Efficiency");
    H3_p->SetLineColor(kRed); H3_p->SetLineWidth(2);  H3_p->SetTitle("MIP Trigger;BGO Deposit Energy (GeV); Trigger Efficiency");
    H4_p->SetLineColor(kRed); H4_p->SetLineWidth(2);  H4_p->SetTitle("RMS_{0}>15 && RMS_{1}>15 && RMS_{2}<45 && RMS_{3} < 45;BGO Deposit Energy (GeV); Trigger Efficiency");

    H1_d->SetLineColor(kBlue); H1_d->SetLineWidth(2); 
    H2_d->SetLineColor(kBlue); H2_d->SetLineWidth(2); 
    H3_d->SetLineColor(kBlue); H3_d->SetLineWidth(2); 
    H4_d->SetLineColor(kBlue); H4_d->SetLineWidth(2); 

    auto c0 = new TCanvas("c0","c0",1200,1200);
    c0->Clear();
    c0->cd();
    gPad->SetLogx();gStyle->SetOptStat(0);
    
    H2_p->Draw("E1P");
    H2_p->GetYaxis()->SetRangeUser(0.998, 1.001);
    H2_p->GetYaxis()->SetNdivisions(505);
    H2_p->GetYaxis()->SetTitleOffset(1.45);
    H2_d->Draw("E1PSAME");

    auto legend0 = new TLegend(0.60, 0.68, 0.88, 0.88);

    legend0->AddEntry(H2_p, "UBT Proton", "ep");
    legend0->AddEntry(H2_d, "UBT Deuteron", "ep");
    legend0->Draw();

    auto c1 = new TCanvas("c1","c1",1200,1200);
    c1->Clear();
    c1->cd();
    gPad->SetLogx();gStyle->SetOptStat(0);

    H1_p->Draw("E1P");
    H1_d->Draw("E1PSAME");
    H1_p->GetYaxis()->SetRangeUser(0., 1.05);
    H1_p->GetYaxis()->SetNdivisions(505);

    auto legend1 = new TLegend(0.72, 0.78, 0.88, 0.88);

    legend1->AddEntry(H1_p, "HET Proton", "ep");
    legend1->AddEntry(H1_d, "HET Deuteron", "ep");
    legend1->Draw();

    auto c2 = new TCanvas("c2","c2",1200,1200);
    c2->Clear();
    c2->cd();
    gPad->SetLogx();gStyle->SetOptStat(0);
gStyle->SetOptStat(0);


    H3_p->Draw("E1P");
    H3_d->Draw("E1PSAME");

    H3_p->GetYaxis()->SetRangeUser(0., 1.05);
    H3_p->GetYaxis()->SetNdivisions(505);

    auto legend2 = new TLegend(0.72, 0.78, 0.88, 0.88);

    legend2->AddEntry(H3_p, "MIP Trigger Proton", "ep");
    legend2->AddEntry(H3_d, "MIP Trigger Deuteron", "ep");
    legend2->Draw();

    auto c3 = new TCanvas("c3","c3",1200,1200);
    c3->Clear();
    c3->cd();
    gPad->SetLogx();gStyle->SetOptStat(0);
gStyle->SetOptStat(0);


    H4_p->Draw("E1P");
    H4_d->Draw("E1PSAME");
    H4_p->GetYaxis()->SetRangeUser(0., 1.05);
    H4_p->GetYaxis()->SetNdivisions(505);

    auto legend3 = new TLegend(0.72, 0.78, 0.88, 0.88);

    legend3->AddEntry(H3_p, "Proton", "ep");
    legend3->AddEntry(H3_d, "Deuteron", "ep");
    legend3->Draw();

}