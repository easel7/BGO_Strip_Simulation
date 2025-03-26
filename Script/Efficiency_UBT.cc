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

void Efficiency_UBT()
{
    double Proton_Eff_HET[9]={0};
    double Proton_Eff_UBT[9]={0};
    double Proton_Eff_MIT[9]={0};

    double Deuteron_Eff_HET[9]={0};
    double Deuteron_Eff_UBT[9]={0};
    double Deuteron_Eff_MIT[9]={0};

    double Electron_Eff_HET[9]={0};
    double Electron_Eff_UBT[9]={0};
    double Electron_Eff_MIT[9]={0};

    double Helium4_Eff_HET[9]={0};
    double Helium4_Eff_UBT[9]={0};
    double Helium4_Eff_MIT[9]={0};

    double Helium3_Eff_HET[9]={0};
    double Helium3_Eff_UBT[9]={0};
    double Helium3_Eff_MIT[9]={0};

    double Carbon_Eff_HET[9]={0};
    double Carbon_Eff_UBT[9]={0};
    double Carbon_Eff_MIT[9]={0};

    std::vector<double>* p_EnergyVec = nullptr;    double p_Total_E; double p_Energy;
    std::vector<double>* d_EnergyVec = nullptr;    double d_Total_E; double d_Energy;
    std::vector<double>* e_EnergyVec = nullptr;    double e_Total_E; double e_Energy;
    std::vector<double>* h_EnergyVec = nullptr;    double h_Total_E; double h_Energy;
    std::vector<double>* H_EnergyVec = nullptr;    double H_Total_E; double H_Energy;
    std::vector<double>* c_EnergyVec = nullptr;    double c_Total_E; double c_Energy;
    double Energy[9]={0};
    double Energy_Err[9]={0};
    double Uncertainty[9]={0};
    double n_BGO = TMath::Na()*7.13/1245.8344; // cm-3
    TH1D *h1_p[9]; TH1D *h2_p[9]; TH1D *h3_p[9]; TH1D *h4_p[9];
    TH1D *h1_d[9]; TH1D *h2_d[9]; TH1D *h3_d[9]; TH1D *h4_d[9];
    TH1D *h1_e[9]; TH1D *h2_e[9]; TH1D *h3_e[9]; TH1D *h4_e[9];
    TH1D *h1_h[9]; TH1D *h2_h[9]; TH1D *h3_h[9]; TH1D *h4_h[9];
    TH1D *h1_H[9]; TH1D *h2_H[9]; TH1D *h3_H[9]; TH1D *h4_H[9];
    TH1D *h1_c[9]; TH1D *h2_c[9]; TH1D *h3_c[9]; TH1D *h4_c[9];

    for (int i = 0; i < 9; i++)
    {

        if(i<10)  {Energy[i] =  (i+1)*10; Energy_Err[i] = 5;}
        else   {Energy[i] =  i*100-800;Energy_Err[i] = 50;}
        h1_p[i] = new TH1D(Form("h1_p[%d]",i),Form("h1_p[%d]",i),100,0,1);
        h2_p[i] = new TH1D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),100,0,1);
        h3_p[i] = new TH1D(Form("h3_p[%d]",i),Form("h3_p[%d]",i),100,0,1);
        h4_p[i] = new TH1D(Form("h4_p[%d]",i),Form("h4_p[%d]",i),100,0,1);
        h1_d[i] = new TH1D(Form("h1_d[%d]",i),Form("h1_d[%d]",i),100,0,1);
        h2_d[i] = new TH1D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),100,0,1);
        h3_d[i] = new TH1D(Form("h3_d[%d]",i),Form("h3_d[%d]",i),100,0,1);
        h4_d[i] = new TH1D(Form("h4_d[%d]",i),Form("h4_d[%d]",i),100,0,1);
        h1_e[i] = new TH1D(Form("h1_e[%d]",i),Form("h1_e[%d]",i),100,0,1);
        h2_e[i] = new TH1D(Form("h2_e[%d]",i),Form("h2_e[%d]",i),100,0,1);
        h3_e[i] = new TH1D(Form("h3_e[%d]",i),Form("h3_e[%d]",i),100,0,1);
        h4_e[i] = new TH1D(Form("h4_e[%d]",i),Form("h4_e[%d]",i),100,0,1);
        h1_h[i] = new TH1D(Form("h1_h[%d]",i),Form("h1_h[%d]",i),100,0,1);
        h2_h[i] = new TH1D(Form("h2_h[%d]",i),Form("h2_h[%d]",i),100,0,1);
        h3_h[i] = new TH1D(Form("h3_h[%d]",i),Form("h3_h[%d]",i),100,0,1);
        h4_h[i] = new TH1D(Form("h4_h[%d]",i),Form("h4_h[%d]",i),100,0,1);
        h1_H[i] = new TH1D(Form("h1_H[%d]",i),Form("h1_H[%d]",i),100,0,1);
        h2_H[i] = new TH1D(Form("h2_H[%d]",i),Form("h2_H[%d]",i),100,0,1);
        h3_H[i] = new TH1D(Form("h3_H[%d]",i),Form("h3_H[%d]",i),100,0,1);
        h4_H[i] = new TH1D(Form("h4_H[%d]",i),Form("h4_H[%d]",i),100,0,1);
        h1_c[i] = new TH1D(Form("h1_c[%d]",i),Form("h1_c[%d]",i),100,0,1);
        h2_c[i] = new TH1D(Form("h2_c[%d]",i),Form("h2_c[%d]",i),100,0,1);
        h3_c[i] = new TH1D(Form("h3_c[%d]",i),Form("h3_c[%d]",i),100,0,1);
        h4_c[i] = new TH1D(Form("h4_c[%d]",i),Form("h4_c[%d]",i),100,0,1);

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
        proton_tree->SetBranchAddress("Total_E",&p_Total_E);
        proton_tree->SetBranchAddress("Energy",&p_Energy);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
        deuteron_tree->SetBranchAddress("Total_E",&d_Total_E);
        deuteron_tree->SetBranchAddress("Energy",&d_Energy);

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
        electron_tree->SetBranchAddress("Total_E",&e_Total_E);
        electron_tree->SetBranchAddress("Energy",&e_Energy);

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        helium4_tree->SetBranchAddress("LayerEnergyVector",&h_EnergyVec);
        helium4_tree->SetBranchAddress("Total_E",&h_Total_E);
        helium4_tree->SetBranchAddress("Energy",&h_Energy);

        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        helium3_tree->SetBranchAddress("LayerEnergyVector",&H_EnergyVec);
        helium3_tree->SetBranchAddress("Total_E",&H_Total_E);
        helium3_tree->SetBranchAddress("Energy",&H_Energy);

        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        carbon_tree->SetBranchAddress("LayerEnergyVector",&c_EnergyVec);
        carbon_tree->SetBranchAddress("Total_E",&c_Total_E);
        carbon_tree->SetBranchAddress("Energy",&c_Energy);

        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
        {
            proton_tree->GetEntry(entry);   
            h1_p[i]->Fill(p_Total_E/p_Energy); 
            if ( (*p_EnergyVec)[0] > 0.0092 && (*p_EnergyVec)[1]  > 0.0092) h3_p[i]->Fill(p_Total_E/p_Energy); 
            if ( (*p_EnergyVec)[0] > 0.23   && (*p_EnergyVec)[1]  > 0.23   && (*p_EnergyVec)[2]  > 0.23     &&  (*p_EnergyVec)[3] > 0.046) h2_p[i]->Fill(p_Total_E/p_Energy); 
            if (((*p_EnergyVec)[2] > 0.0092 && (*p_EnergyVec)[10] > 0.0092 && (*p_EnergyVec)[12] > 0.0092 ) || ((*p_EnergyVec)[3] > 0.0092 && (*p_EnergyVec)[11] > 0.0092 && (*p_EnergyVec)[13] > 0.0092 )) h4_p[i]->Fill(p_Total_E/p_Energy);
            
            deuteron_tree->GetEntry(entry); 
            h1_d[i]->Fill(d_Total_E/d_Energy); 
            if ( (*d_EnergyVec)[0] > 0.0092 && (*d_EnergyVec)[1] > 0.0092) h3_d[i]->Fill(d_Total_E/d_Energy); 
            if ( (*d_EnergyVec)[0] > 0.23   && (*d_EnergyVec)[1] > 0.23 && (*d_EnergyVec)[2] > 0.23 && (*d_EnergyVec)[3] > 0.046) h2_d[i]->Fill(d_Total_E/d_Energy); 
            if (((*d_EnergyVec)[2] > 0.0092 && (*d_EnergyVec)[10] > 0.0092 && (*d_EnergyVec)[12] > 0.0092 ) || ((*d_EnergyVec)[3] > 0.0092 && (*d_EnergyVec)[11] > 0.0092 && (*d_EnergyVec)[13] > 0.0092 )) h4_d[i]->Fill(d_Total_E/d_Energy);
            
            electron_tree->GetEntry(entry); 
            h1_e[i]->Fill(e_Total_E/e_Energy); 
            if ((*e_EnergyVec)[0] > 0.0092 && (*e_EnergyVec)[1] > 0.0092) h3_e[i]->Fill(e_Total_E/e_Energy); 
            if ((*e_EnergyVec)[0] > 0.23 && (*e_EnergyVec)[1] > 0.23 && (*e_EnergyVec)[2] > 0.23 && (*e_EnergyVec)[3] > 0.046) h2_e[i]->Fill(e_Total_E/e_Energy); 
            if (((*e_EnergyVec)[2] > 0.0092 && (*e_EnergyVec)[10] > 0.0092 && (*e_EnergyVec)[12] > 0.0092 ) || ((*e_EnergyVec)[3] > 0.0092 && (*e_EnergyVec)[11] > 0.0092 && (*e_EnergyVec)[13] > 0.0092 )) h4_e[i]->Fill(e_Total_E/e_Energy);
            
            helium4_tree->GetEntry(entry);  
            h1_h[i]->Fill(h_Total_E/h_Energy); 
            if ((*h_EnergyVec)[0] > 0.0092 && (*h_EnergyVec)[1] > 0.0092) h3_h[i]->Fill(h_Total_E/h_Energy); 
            if ((*h_EnergyVec)[0] > 0.23 && (*h_EnergyVec)[1] > 0.23 && (*h_EnergyVec)[2] > 0.23 && (*h_EnergyVec)[3] > 0.046) h2_h[i]->Fill(h_Total_E/h_Energy); 
            if (((*h_EnergyVec)[2] > 0.0092 && (*h_EnergyVec)[10] > 0.0092 && (*h_EnergyVec)[12] > 0.0092 ) || ((*h_EnergyVec)[3] > 0.0092 && (*h_EnergyVec)[11] > 0.0092 && (*h_EnergyVec)[13] > 0.0092 )) h4_h[i]->Fill(h_Total_E/h_Energy);
            helium3_tree->GetEntry(entry);  
            h1_H[i]->Fill(H_Total_E/H_Energy); 
            if ((*H_EnergyVec)[0] > 0.0092 && (*H_EnergyVec)[1] > 0.0092) h3_H[i]->Fill(H_Total_E/H_Energy); 
            if ((*H_EnergyVec)[0] > 0.23 && (*H_EnergyVec)[1] > 0.23 && (*H_EnergyVec)[2] > 0.23 && (*H_EnergyVec)[3] > 0.046) h2_H[i]->Fill(H_Total_E/H_Energy); 
            if (((*H_EnergyVec)[2] > 0.0092 && (*H_EnergyVec)[10] > 0.0092 && (*H_EnergyVec)[12] > 0.0092 ) || ((*H_EnergyVec)[3] > 0.0092 && (*H_EnergyVec)[11] > 0.0092 && (*H_EnergyVec)[13] > 0.0092 )) h4_H[i]->Fill(H_Total_E/H_Energy);
            carbon_tree->GetEntry(entry);   
            h1_c[i]->Fill(c_Total_E/c_Energy); 
            if ((*c_EnergyVec)[0] > 0.0092 && (*c_EnergyVec)[1] > 0.0092) h3_c[i]->Fill(c_Total_E/c_Energy); 
            if ((*c_EnergyVec)[0] > 0.23 && (*c_EnergyVec)[1] > 0.23 && (*c_EnergyVec)[2] > 0.23 && (*c_EnergyVec)[3] > 0.046) h2_c[i]->Fill(c_Total_E/c_Energy); 
            if (((*c_EnergyVec)[2] > 0.0092 && (*c_EnergyVec)[10] > 0.0092 && (*c_EnergyVec)[12] > 0.0092 ) || ((*c_EnergyVec)[3] > 0.0092 && (*c_EnergyVec)[11] > 0.0092 && (*c_EnergyVec)[13] > 0.0092 )) h4_c[i]->Fill(c_Total_E/c_Energy);
        }

        Proton_Eff_HET[i]   = h2_p[i]->Integral()/h1_p[i]->Integral();
        Proton_Eff_UBT[i]   = h3_p[i]->Integral()/h1_p[i]->Integral();
        Proton_Eff_MIT[i]   = h4_p[i]->Integral()/h1_p[i]->Integral();
        Deuteron_Eff_HET[i] = h2_d[i]->Integral()/h1_d[i]->Integral();
        Deuteron_Eff_UBT[i] = h3_d[i]->Integral()/h1_d[i]->Integral();
        Deuteron_Eff_MIT[i] = h4_d[i]->Integral()/h1_d[i]->Integral();
        Electron_Eff_HET[i] = h2_e[i]->Integral()/h1_e[i]->Integral();
        Electron_Eff_UBT[i] = h3_e[i]->Integral()/h1_e[i]->Integral();
        Electron_Eff_MIT[i] = h4_e[i]->Integral()/h1_e[i]->Integral();
        Helium4_Eff_HET[i]  = h2_h[i]->Integral()/h1_h[i]->Integral();
        Helium4_Eff_UBT[i]  = h3_h[i]->Integral()/h1_h[i]->Integral();
        Helium4_Eff_MIT[i]  = h4_h[i]->Integral()/h1_h[i]->Integral();
        Helium3_Eff_HET[i]  = h2_H[i]->Integral()/h1_H[i]->Integral();
        Helium3_Eff_UBT[i]  = h3_H[i]->Integral()/h1_H[i]->Integral();
        Helium3_Eff_MIT[i]  = h4_H[i]->Integral()/h1_H[i]->Integral();
        Carbon_Eff_HET[i]   = h2_c[i]->Integral()/h1_c[i]->Integral();
        Carbon_Eff_UBT[i]   = h3_c[i]->Integral()/h1_c[i]->Integral();
        Carbon_Eff_MIT[i]   = h4_c[i]->Integral()/h1_c[i]->Integral();

        cout << "Energy = " << int(Energy[i]) << " GeV !  ele Total eff : " << h1_e[i]->Integral() << " ele HET eff : " <<  h2_e[i]->Integral() <<" ele UBT eff : " << h3_e[i]->Integral() <<" ele MIP eff : " << h4_e[i]->Integral() <<endl;
    }

    auto gre_p_HET = CreateGraphWithProperties(9, Energy, Proton_Eff_HET, Energy_Err, Uncertainty, 20, kRed, kRed, "High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_p_UBT = CreateGraphWithProperties(9, Energy, Proton_Eff_UBT, Energy_Err, Uncertainty, 20, kRed, kRed, "Unbiased Trigger;Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_p_MIT = CreateGraphWithProperties(9, Energy, Proton_Eff_MIT, Energy_Err, Uncertainty, 20, kRed, kRed, "MIP Trigger;Kinetic Energy (GeV); Trigger Efficiency");

    auto gre_d_HET = CreateGraphWithProperties(9, Energy, Deuteron_Eff_HET, Energy_Err, Uncertainty, 21, kBlue, kBlue, "High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_d_UBT = CreateGraphWithProperties(9, Energy, Deuteron_Eff_UBT, Energy_Err, Uncertainty, 21, kBlue, kBlue, "Unbiased Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_d_MIT = CreateGraphWithProperties(9, Energy, Deuteron_Eff_MIT, Energy_Err, Uncertainty, 21, kBlue, kBlue, "MIP Trigger; Kinetic Energy (GeV); Trigger Efficiency");

    auto gre_e_HET = CreateGraphWithProperties(9, Energy, Electron_Eff_HET, Energy_Err, Uncertainty, 22, kOrange-3, kOrange-3, "High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_e_UBT = CreateGraphWithProperties(9, Energy, Electron_Eff_UBT, Energy_Err, Uncertainty, 22, kOrange-3, kOrange-3, "Unbiased Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_e_MIT = CreateGraphWithProperties(9, Energy, Electron_Eff_MIT, Energy_Err, Uncertainty, 22, kOrange-3, kOrange-3, "MIP Trigger; Kinetic Energy (GeV); Trigger Efficiency");

    auto gre_h_HET = CreateGraphWithProperties(9, Energy, Helium4_Eff_HET, Energy_Err, Uncertainty, 23, kGreen-3, kGreen-3, "High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_h_UBT = CreateGraphWithProperties(9, Energy, Helium4_Eff_UBT, Energy_Err, Uncertainty, 23, kGreen-3, kGreen-3, "Unbiased Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_h_MIT = CreateGraphWithProperties(9, Energy, Helium4_Eff_MIT, Energy_Err, Uncertainty, 23, kGreen-3, kGreen-3, "MIP Trigger; Kinetic Energy (GeV); Trigger Efficiency");

    auto gre_H_HET = CreateGraphWithProperties(9, Energy, Helium3_Eff_HET, Energy_Err, Uncertainty, 32, kGreen-3, kGreen-3, "High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_H_UBT = CreateGraphWithProperties(9, Energy, Helium3_Eff_UBT, Energy_Err, Uncertainty, 32, kGreen-3, kGreen-3, "Unbiased Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_H_MIT = CreateGraphWithProperties(9, Energy, Helium3_Eff_MIT, Energy_Err, Uncertainty, 32, kGreen-3, kGreen-3, "MIP Trigger; Kinetic Energy (GeV); Trigger Efficiency");

    auto gre_c_HET = CreateGraphWithProperties(9, Energy, Carbon_Eff_HET, Energy_Err, Uncertainty, 23, kMagenta, kMagenta, "High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_c_UBT = CreateGraphWithProperties(9, Energy, Carbon_Eff_UBT, Energy_Err, Uncertainty, 23, kMagenta, kMagenta, "Unbiased Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    auto gre_c_MIT = CreateGraphWithProperties(9, Energy, Carbon_Eff_MIT, Energy_Err, Uncertainty, 23, kMagenta, kMagenta, "MIP Trigger; Kinetic Energy (GeV); Trigger Efficiency");

    auto c0 = new TCanvas("c0","c0",1200,1200);
    c0->Clear();
    c0->cd();
    gPad->SetLogx();
    gre_p_UBT->Draw("AP");
    gre_p_UBT->GetYaxis()->SetRangeUser(0.998, 1.001);
    gre_p_UBT->GetYaxis()->SetNdivisions(505);
    gre_p_UBT->GetYaxis()->SetTitleOffset(1.45);
    gre_d_UBT->Draw("PSAME");
    gre_e_UBT->Draw("PSAME");
    gre_h_UBT->Draw("PSAME");
    gre_c_UBT->Draw("PSAME");


    auto legend0 = new TLegend(0.60, 0.68, 0.88, 0.88);
    legend0->SetNColumns(2);
    legend0->AddEntry(gre_p_UBT, "UBT Proton", "ep");
    legend0->AddEntry(gre_d_UBT, "UBT Deuteron", "ep");
    legend0->AddEntry(gre_e_UBT, "UBT Electron", "ep");
    legend0->AddEntry(gre_h_UBT, "UBT Helium4", "ep");
    legend0->AddEntry(gre_H_UBT, "UBT Helium3", "ep");
    legend0->AddEntry(gre_c_UBT, "UBT Carbon", "ep");
    legend0->Draw();

    auto c1 = new TCanvas("c1","c1",1200,1200);
    c1->Clear();
    c1->cd();
    gPad->SetLogx();
    gre_p_HET->Draw("AP");
    gre_d_HET->Draw("PSAME");
    gre_e_HET->Draw("PSAME");
    gre_h_HET->Draw("PSAME");
    gre_c_HET->Draw("PSAME");

    gre_p_HET->GetYaxis()->SetRangeUser(0., 1.05);
    gre_p_HET->GetYaxis()->SetNdivisions(505);

    auto legend1 = new TLegend(0.62, 0.58, 0.88, 0.78);
    legend1->SetNColumns(2);
    legend1->AddEntry(gre_p_HET, "HET Proton", "ep");
    legend1->AddEntry(gre_d_HET, "HET Deuteron", "ep");
    legend1->AddEntry(gre_e_HET, "HET Electron", "ep");
    legend1->AddEntry(gre_h_HET, "HET Helium4", "ep");
    legend1->AddEntry(gre_H_HET, "HET Helium3", "ep");
    legend1->AddEntry(gre_c_HET, "HET Carbon", "ep");
    legend1->Draw();

    auto c2 = new TCanvas("c2","c2",1200,1200);
    c2->Clear();
    c2->cd();
    gPad->SetLogx();
    gre_p_MIT->Draw("AP");
    gre_d_MIT->Draw("PSAME");
    gre_e_MIT->Draw("PSAME");
    gre_h_MIT->Draw("PSAME");
    gre_c_MIT->Draw("PSAME");

    gre_p_MIT->GetYaxis()->SetRangeUser(0., 1.05);
    gre_p_MIT->GetYaxis()->SetNdivisions(505);

    auto legend2 = new TLegend(0.62, 0.58, 0.88, 0.78);
    legend2->SetNColumns(2);
    legend2->AddEntry(gre_p_MIT, "MIP Trigger Proton", "ep");
    legend2->AddEntry(gre_d_MIT, "MIP Trigger Deuteron", "ep");
    legend2->AddEntry(gre_e_MIT, "MIP Trigger Electron", "ep");
    legend2->AddEntry(gre_h_MIT, "MIP Trigger Helium4", "ep");
    legend2->AddEntry(gre_H_MIT, "MIP Trigger Helium3", "ep");
    legend2->AddEntry(gre_c_MIT, "MIP Trigger Carbon", "ep");
    legend2->Draw();

}