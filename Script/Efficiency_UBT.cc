void Efficiency_UBT()
{
    double Proton_Eff_HET[19]={0};
    double Proton_Eff_UBT[19]={0};
    double Proton_Eff_MIT[19]={0};

    double Deuteron_Eff_HET[19]={0};
    double Deuteron_Eff_UBT[19]={0};
    double Deuteron_Eff_MIT[19]={0};

    double Electron_Eff_HET[19]={0};
    double Electron_Eff_UBT[19]={0};
    double Electron_Eff_MIT[19]={0};

    double Helium4_Eff_HET[19]={0};
    double Helium4_Eff_UBT[19]={0};
    double Helium4_Eff_MIT[19]={0};

    double Helium3_Eff_HET[19]={0};
    double Helium3_Eff_UBT[19]={0};
    double Helium3_Eff_MIT[19]={0};

    double Carbon_Eff_HET[19]={0};
    double Carbon_Eff_UBT[19]={0};
    double Carbon_Eff_MIT[19]={0};

    std::vector<double>* p_EnergyVec = nullptr;    double p_Total_E; double p_Energy;
    std::vector<double>* d_EnergyVec = nullptr;    double d_Total_E; double d_Energy;
    std::vector<double>* e_EnergyVec = nullptr;    double e_Total_E; double e_Energy;
    std::vector<double>* h_EnergyVec = nullptr;    double h_Total_E; double h_Energy;
    std::vector<double>* H_EnergyVec = nullptr;    double H_Total_E; double H_Energy;
    std::vector<double>* c_EnergyVec = nullptr;    double c_Total_E; double c_Energy;
    double Energy[19]={0};
    double Energy_Err[19]={0};
    double Uncertainty[19]={0};
    double n_BGO = TMath::Na()*7.13/1245.8344; // cm-3
    TH1D *h1_p[19]; TH1D *h2_p[19]; TH1D *h3_p[19]; TH1D *h4_p[19];
    TH1D *h1_d[19]; TH1D *h2_d[19]; TH1D *h3_d[19]; TH1D *h4_d[19];
    TH1D *h1_e[19]; TH1D *h2_e[19]; TH1D *h3_e[19]; TH1D *h4_e[19];
    TH1D *h1_h[19]; TH1D *h2_h[19]; TH1D *h3_h[19]; TH1D *h4_h[19];
    TH1D *h1_H[19]; TH1D *h2_H[19]; TH1D *h3_H[19]; TH1D *h4_H[19];
    TH1D *h1_c[19]; TH1D *h2_c[19]; TH1D *h3_c[19]; TH1D *h4_c[19];

    TCut UBT = "(L0_E>0.0092 && L1_E>0.0092)";
    TCut HET = "(L0_E>0.23 && L1_E >0.23 && L2_E>0.23 && L3_E>0.046)";
    TCut MIT = "(L2_E >0.0092 && L10_E >0.0092 && L12_E >0.0092)  || (L3_E >0.0092 && L11_E >0.0092 && L13_E >0.0092)";
    for (int i = 0; i < 19; i++)
    {
        h1_p[i] = new TH1D(Form("h1_p[%d]",i),Form("h1_p[%d]",i),100,0,1);        h1_d[i] = new TH1D(Form("h1_d[%d]",i),Form("h1_d[%d]",i),100,0,1);         h1_e[i] = new TH1D(Form("h1_e[%d]",i),Form("h1_e[%d]",i),100,0,1);         h1_h[i] = new TH1D(Form("h1_h[%d]",i),Form("h1_h[%d]",i),100,0,1);        h1_H[i] = new TH1D(Form("h1_H[%d]",i),Form("h1_H[%d]",i),100,0,1);        h1_c[i] = new TH1D(Form("h1_c[%d]",i),Form("h1_c[%d]",i),100,0,1);
        h2_p[i] = new TH1D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),100,0,1);        h2_d[i] = new TH1D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),100,0,1);         h2_e[i] = new TH1D(Form("h2_e[%d]",i),Form("h2_e[%d]",i),100,0,1);         h2_h[i] = new TH1D(Form("h2_h[%d]",i),Form("h2_h[%d]",i),100,0,1);        h2_H[i] = new TH1D(Form("h2_H[%d]",i),Form("h2_H[%d]",i),100,0,1);        h2_c[i] = new TH1D(Form("h2_c[%d]",i),Form("h2_c[%d]",i),100,0,1);
        h3_p[i] = new TH1D(Form("h3_p[%d]",i),Form("h3_p[%d]",i),100,0,1);        h3_d[i] = new TH1D(Form("h3_d[%d]",i),Form("h3_d[%d]",i),100,0,1);         h3_e[i] = new TH1D(Form("h3_e[%d]",i),Form("h3_e[%d]",i),100,0,1);         h3_h[i] = new TH1D(Form("h3_h[%d]",i),Form("h3_h[%d]",i),100,0,1);        h3_H[i] = new TH1D(Form("h3_H[%d]",i),Form("h3_H[%d]",i),100,0,1);        h3_c[i] = new TH1D(Form("h3_c[%d]",i),Form("h3_c[%d]",i),100,0,1);
        h4_p[i] = new TH1D(Form("h4_p[%d]",i),Form("h4_p[%d]",i),100,0,1);        h4_d[i] = new TH1D(Form("h4_d[%d]",i),Form("h4_d[%d]",i),100,0,1);         h4_e[i] = new TH1D(Form("h4_e[%d]",i),Form("h4_e[%d]",i),100,0,1);         h4_h[i] = new TH1D(Form("h4_h[%d]",i),Form("h4_h[%d]",i),100,0,1);        h4_H[i] = new TH1D(Form("h4_H[%d]",i),Form("h4_H[%d]",i),100,0,1);        h4_c[i] = new TH1D(Form("h4_c[%d]",i),Form("h4_c[%d]",i),100,0,1);
    }
    for (int i = 0; i < 19; i++)
    {
        if(i<10)  {Energy[i] =  (i+1)*10; Energy_Err[i] = 5;}
        else   {Energy[i] =  i*100-800;Energy_Err[i] = 50;}

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
        proton_tree->SetBranchAddress("Total_E",&p_Total_E);
        proton_tree->SetBranchAddress("Energy",&p_Energy);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
        deuteron_tree->SetBranchAddress("Total_E",&d_Total_E);

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
        electron_tree->SetBranchAddress("Total_E",&e_Total_E);

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        helium4_tree->SetBranchAddress("LayerEnergyVector",&h_EnergyVec);
        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        helium3_tree->SetBranchAddress("LayerEnergyVector",&H_EnergyVec);
        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        carbon_tree->SetBranchAddress("LayerEnergyVector",&c_EnergyVec);

        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
        {
            proton_tree->GetEntry(entry);   
            h1_p[i]->Fill(p_Total_E/p_Energy); 
            if ((*p_EnergyVec)[0] > 0.0092 && (*p_EnergyVec)[1] > 0.0092) h2_p[i]->Fill(p_Total_E/p_Energy); 
            if ((*p_EnergyVec)[0] > 0.23 && (*p_EnergyVec)[1] > 0.23 && (*p_EnergyVec)[2] > 0.23 && (*p_EnergyVec)[0] > 0.046) h3_p[i]->Fill(p_Total_E/p_Energy); 
            if (((*p_EnergyVec)[2] > 0.0092 && (*p_EnergyVec)[10] > 0.0092 && (*p_EnergyVec)[12] > 0.0092 ) || ((*p_EnergyVec)[3] > 0.0092 && (*p_EnergyVec)[11] > 0.0092 && (*p_EnergyVec)[13] > 0.0092 )) h4_p[i]->Fill(p_Total_E/p_Energy);
            deuteron_tree->GetEntry(entry); 
            h1_d[i]->Fill(d_Total_E/d_Energy); 
            if ((*d_EnergyVec)[0] > 0.0092 && (*d_EnergyVec)[1] > 0.0092) h2_d[i]->Fill(d_Total_E/d_Energy); 
            if ((*d_EnergyVec)[0] > 0.23 && (*d_EnergyVec)[1] > 0.23 && (*d_EnergyVec)[2] > 0.23 && (*d_EnergyVec)[0] > 0.046) h3_d[i]->Fill(d_Total_E/d_Energy); 
            if (((*d_EnergyVec)[2] > 0.0092 && (*d_EnergyVec)[10] > 0.0092 && (*d_EnergyVec)[12] > 0.0092 ) || ((*d_EnergyVec)[3] > 0.0092 && (*d_EnergyVec)[11] > 0.0092 && (*d_EnergyVec)[13] > 0.0092 )) h4_d[i]->Fill(d_Total_E/d_Energy);
            electron_tree->GetEntry(entry); 
            h1_e[i]->Fill(e_Total_E/e_Energy); 
            if ((*e_EnergyVec)[0] > 0.0092 && (*e_EnergyVec)[1] > 0.0092) h2_e[i]->Fill(e_Total_E/e_Energy); 
            if ((*e_EnergyVec)[0] > 0.23 && (*e_EnergyVec)[1] > 0.23 && (*e_EnergyVec)[2] > 0.23 && (*e_EnergyVec)[0] > 0.046) h3_e[i]->Fill(e_Total_E/e_Energy); 
            if (((*e_EnergyVec)[2] > 0.0092 && (*e_EnergyVec)[10] > 0.0092 && (*e_EnergyVec)[12] > 0.0092 ) || ((*e_EnergyVec)[3] > 0.0092 && (*e_EnergyVec)[11] > 0.0092 && (*e_EnergyVec)[13] > 0.0092 )) h4_e[i]->Fill(e_Total_E/e_Energy);
            helium4_tree->GetEntry(entry);  
            h1_h[i]->Fill(h_Total_E/h_Energy); 
            if ((*h_EnergyVec)[0] > 0.0092 && (*h_EnergyVec)[1] > 0.0092) h2_h[i]->Fill(h_Total_E/h_Energy); 
            if ((*h_EnergyVec)[0] > 0.23 && (*h_EnergyVec)[1] > 0.23 && (*h_EnergyVec)[2] > 0.23 && (*h_EnergyVec)[0] > 0.046) h3_h[i]->Fill(h_Total_E/h_Energy); 
            if (((*h_EnergyVec)[2] > 0.0092 && (*h_EnergyVec)[10] > 0.0092 && (*h_EnergyVec)[12] > 0.0092 ) || ((*h_EnergyVec)[3] > 0.0092 && (*h_EnergyVec)[11] > 0.0092 && (*h_EnergyVec)[13] > 0.0092 )) h4_h[i]->Fill(h_Total_E/h_Energy);
            helium3_tree->GetEntry(entry);  
            h1_H[i]->Fill(H_Total_E/H_Energy); 
            if ((*H_EnergyVec)[0] > 0.0092 && (*H_EnergyVec)[1] > 0.0092) h2_H[i]->Fill(H_Total_E/H_Energy); 
            if ((*H_EnergyVec)[0] > 0.23 && (*H_EnergyVec)[1] > 0.23 && (*H_EnergyVec)[2] > 0.23 && (*H_EnergyVec)[0] > 0.046) h3_H[i]->Fill(H_Total_E/H_Energy); 
            if (((*H_EnergyVec)[2] > 0.0092 && (*H_EnergyVec)[10] > 0.0092 && (*H_EnergyVec)[12] > 0.0092 ) || ((*H_EnergyVec)[3] > 0.0092 && (*H_EnergyVec)[11] > 0.0092 && (*H_EnergyVec)[13] > 0.0092 )) h4_H[i]->Fill(H_Total_E/H_Energy);
            carbon_tree->GetEntry(entry);   
            h1_c[i]->Fill(c_Total_E/c_Energy); 
            if ((*c_EnergyVec)[0] > 0.0092 && (*c_EnergyVec)[1] > 0.0092) h2_c[i]->Fill(c_Total_E/c_Energy); 
            if ((*c_EnergyVec)[0] > 0.23 && (*c_EnergyVec)[1] > 0.23 && (*c_EnergyVec)[2] > 0.23 && (*c_EnergyVec)[0] > 0.046) h3_c[i]->Fill(c_Total_E/c_Energy); 
            if (((*c_EnergyVec)[2] > 0.0092 && (*c_EnergyVec)[10] > 0.0092 && (*c_EnergyVec)[12] > 0.0092 ) || ((*c_EnergyVec)[3] > 0.0092 && (*c_EnergyVec)[11] > 0.0092 && (*c_EnergyVec)[13] > 0.0092 )) h4_c[i]->Fill(c_Total_E/c_Energy);
        }



    }

    Proton_Eff_HET[i] = h2_p[i]->Integral()/h1_p[i]->Integral();
    Proton_Eff_UBT[i] = h3_p[i]->Integral()/h1_p[i]->Integral();
    Proton_Eff_MIT[i] = h4_p[i]->Integral()/h1_p[i]->Integral();
    Deuteron_Eff_HET[i] = h2_d[i]->Integral()/h1_d[i]->Integral();
    Deuteron_Eff_UBT[i] = h3_d[i]->Integral()/h1_d[i]->Integral();
    Deuteron_Eff_MIT[i] = h4_d[i]->Integral()/h1_d[i]->Integral();
    Electron_Eff_HET[i] = h2_e[i]->Integral()/h1_e[i]->Integral();
    Electron_Eff_UBT[i] = h3_e[i]->Integral()/h1_e[i]->Integral();
    Electron_Eff_MIT[i] = h4_e[i]->Integral()/h1_e[i]->Integral();
    Helium4_Eff_HET[i] = h2_h[i]->Integral()/h1_h[i]->Integral();
    Helium4_Eff_UBT[i] = h3_h[i]->Integral()/h1_h[i]->Integral();
    Helium4_Eff_MIT[i] = h4_h[i]->Integral()/h1_h[i]->Integral();
    Helium3_Eff_HET[i] = h2_H[i]->Integral()/h1_H[i]->Integral();
    Helium3_Eff_UBT[i] = h3_H[i]->Integral()/h1_H[i]->Integral();
    Helium3_Eff_MIT[i] = h4_H[i]->Integral()/h1_H[i]->Integral();
    Carbon_Eff_HET[i] = h2_c[i]->Integral()/h1_c[i]->Integral();
    Carbon_Eff_UBT[i] = h3_c[i]->Integral()/h1_c[i]->Integral();
    Carbon_Eff_MIT[i] = h4_c[i]->Integral()/h1_c[i]->Integral();
    

    auto gre_p_HET = new TGraphErrors(19, Energy, Proton_Eff_HET, Energy_Err, Uncertainty);
    auto gre_p_UBT = new TGraphErrors(19, Energy, Proton_Eff_UBT, Energy_Err, Uncertainty);
    auto gre_p_MIT = new TGraphErrors(19, Energy, Proton_Eff_MIT, Energy_Err, Uncertainty);

    gre_p_HET->SetTitle("High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    gre_p_HET->SetMarkerStyle(20);
    gre_p_HET->SetMarkerColor(kRed);
    gre_p_HET->SetLineColor(kRed);
    gre_p_UBT->SetTitle("Unbiased Trigger;Kinetic Energy (GeV); Trigger Efficiency");
    gre_p_UBT->SetMarkerStyle(20);
    gre_p_UBT->SetMarkerColor(kRed);
    gre_p_UBT->SetLineColor(kRed);
    gre_p_MIT->SetTitle("MIP Trigger;Kinetic Energy (GeV); Trigger Efficiency");
    gre_p_MIT->SetMarkerStyle(20);
    gre_p_MIT->SetMarkerColor(kRed);
    gre_p_MIT->SetLineColor(kRed);

    auto gre_d_HET = new TGraphErrors(19, Energy, Deuteron_Eff_HET, Energy_Err, Uncertainty);
    auto gre_d_UBT = new TGraphErrors(19, Energy, Deuteron_Eff_UBT, Energy_Err, Uncertainty);
    auto gre_d_MIT = new TGraphErrors(19, Energy, Deuteron_Eff_MIT, Energy_Err, Uncertainty);

    gre_d_HET->SetMarkerStyle(21);
    gre_d_HET->SetMarkerColor(kBlue);
    gre_d_HET->SetLineColor(kBlue);
    gre_d_UBT->SetMarkerStyle(21);
    gre_d_UBT->SetMarkerColor(kBlue);
    gre_d_UBT->SetLineColor(kBlue);
    gre_d_MIT->SetMarkerStyle(21);
    gre_d_MIT->SetMarkerColor(kBlue);
    gre_d_MIT->SetLineColor(kBlue);

    auto gre_e_HET = new TGraphErrors(19, Energy, Electron_Eff_HET, Energy_Err, Uncertainty);
    auto gre_e_UBT = new TGraphErrors(19, Energy, Electron_Eff_UBT, Energy_Err, Uncertainty);
    auto gre_e_MIT = new TGraphErrors(19, Energy, Electron_Eff_MIT, Energy_Err, Uncertainty);


    gre_e_HET->SetMarkerStyle(22);
    gre_e_HET->SetMarkerColor(kOrange-3);
    gre_e_HET->SetLineColor(kOrange-3);
    gre_e_UBT->SetMarkerStyle(22);
    gre_e_UBT->SetMarkerColor(kOrange-3);
    gre_e_UBT->SetLineColor(kOrange-3);
    gre_e_MIT->SetMarkerStyle(22);
    gre_e_MIT->SetMarkerColor(kOrange-3);
    gre_e_MIT->SetLineColor(kOrange-3);


    auto gre_h_HET = new TGraphErrors(19, Energy, Helium4_Eff_HET, Energy_Err, Uncertainty);
    auto gre_h_UBT = new TGraphErrors(19, Energy, Helium4_Eff_UBT, Energy_Err, Uncertainty);
    auto gre_h_MIT = new TGraphErrors(19, Energy, Helium4_Eff_MIT, Energy_Err, Uncertainty);

    gre_h_HET->SetMarkerStyle(23);
    gre_h_HET->SetMarkerColor(kGreen-3);
    gre_h_HET->SetLineColor(kGreen-3);
    gre_h_UBT->SetMarkerStyle(23);
    gre_h_UBT->SetMarkerColor(kGreen-3);
    gre_h_UBT->SetLineColor(kGreen-3);
    gre_h_MIT->SetMarkerStyle(23);
    gre_h_MIT->SetMarkerColor(kGreen-3);
    gre_h_MIT->SetLineColor(kGreen-3);

    auto gre_H_HET = new TGraphErrors(19, Energy, Helium3_Eff_HET, Energy_Err, Uncertainty);
    auto gre_H_UBT = new TGraphErrors(19, Energy, Helium3_Eff_UBT, Energy_Err, Uncertainty);
    auto gre_H_MIT = new TGraphErrors(19, Energy, Helium3_Eff_MIT, Energy_Err, Uncertainty);

    gre_H_HET->SetMarkerStyle(32);
    gre_H_HET->SetMarkerColor(kGreen-3);
    gre_H_HET->SetLineColor(kGreen-3);
    gre_H_UBT->SetMarkerStyle(32);
    gre_H_UBT->SetMarkerColor(kGreen-3);
    gre_H_UBT->SetLineColor(kGreen-3);
    gre_H_MIT->SetMarkerStyle(32);
    gre_H_MIT->SetMarkerColor(kGreen-3);
    gre_H_MIT->SetLineColor(kGreen-3);

    auto gre_c_HET = new TGraphErrors(19, Energy, Carbon_Eff_HET, Energy_Err, Uncertainty);
    auto gre_c_UBT = new TGraphErrors(19, Energy, Carbon_Eff_UBT, Energy_Err, Uncertainty);
    auto gre_c_MIT = new TGraphErrors(19, Energy, Carbon_Eff_MIT, Energy_Err, Uncertainty);

    gre_c_HET->SetMarkerStyle(23);
    gre_c_HET->SetMarkerColor(kMagenta);
    gre_c_HET->SetLineColor(kMagenta);
    gre_c_UBT->SetMarkerStyle(32);
    gre_c_UBT->SetMarkerColor(kMagenta);
    gre_c_UBT->SetLineColor(kMagenta);
    gre_c_MIT->SetMarkerStyle(32);
    gre_c_MIT->SetMarkerColor(kMagenta);
    gre_c_MIT->SetLineColor(kMagenta);

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