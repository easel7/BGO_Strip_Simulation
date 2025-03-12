void Efficiency_UBT()
{
    double Proton_Eff_HET[19]={0};
    double Proton_Eff_UBT[19]={0};

    double Deuteron_Eff_HET[19]={0};
    double Deuteron_Eff_UBT[19]={0};

    double Electron_Eff_HET[19]={0};
    double Electron_Eff_UBT[19]={0};

    double Helium4_Eff_HET[19]={0};
    double Helium4_Eff_UBT[19]={0};

    double Helium3_Eff_HET[19]={0};
    double Helium3_Eff_UBT[19]={0};

    double Carbon_Eff_HET[19]={0};
    double Carbon_Eff_UBT[19]={0};

    double Energy[19]={0};
    double Energy_Err[19]={0};
    double Uncertainty[19]={0};
    double n_BGO = TMath::Na()*7.13/1245.8344; // cm-3
    TH1D *h1_p[19]; TH1D *h2_p[19]; TH1D *h3_p[19];
    TH1D *h1_d[19]; TH1D *h2_d[19]; TH1D *h3_d[19];
    TH1D *h1_e[19]; TH1D *h2_e[19]; TH1D *h3_e[19];
    TH1D *h1_h[19]; TH1D *h2_h[19]; TH1D *h3_h[19];
    TH1D *h1_H[19]; TH1D *h2_H[19]; TH1D *h3_H[19];
    TH1D *h1_c[19]; TH1D *h2_c[19]; TH1D *h3_c[19];

    TCut UBT = "(L0_E>0.0092 && L1_E>0.0092)";
    TCut HET = "(L0_E>0.23 && L1_E >0.23 &&   L2_E>0.23 && L3_E>0.046)";

    for (int i = 0; i < 19; i++)
    {
        if(i<10)  {Energy[i] =  (i+1)*10; Energy_Err[i] = 5;}
        else   {Energy[i] =  i*100-800;Energy_Err[i] = 50;}

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        h1_p[i] = new TH1D(Form("h1_p[%d]",i),Form("h1_p[%d]",i),100,0,1);
        h2_p[i] = new TH1D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),100,0,1);
        h3_p[i] = new TH1D(Form("h3_p[%d]",i),Form("h3_p[%d]",i),100,0,1);
        proton_tree->Draw(Form("Total_E/Energy>>h1_p[%d]",i),"","");
        proton_tree->Draw(Form("Total_E/Energy>>h2_p[%d]",i), HET,"");
        proton_tree->Draw(Form("Total_E/Energy>>h3_p[%d]",i), UBT,"");
        Proton_Eff_HET[i] = h2_p[i]->Integral()/h1_p[i]->Integral();
        Proton_Eff_UBT[i] = h3_p[i]->Integral()/h1_p[i]->Integral();

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        h1_d[i] = new TH1D(Form("h1_d[%d]",i),Form("h1_d[%d]",i),100,0,1);
        h2_d[i] = new TH1D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),100,0,1);
        h3_d[i] = new TH1D(Form("h3_d[%d]",i),Form("h3_d[%d]",i),100,0,1);
        deuteron_tree->Draw(Form("Total_E/Energy>>h1_d[%d]",i),"", "");
        deuteron_tree->Draw(Form("Total_E/Energy>>h2_d[%d]",i), HET,"");
        deuteron_tree->Draw(Form("Total_E/Energy>>h3_d[%d]",i), UBT,"");
        Deuteron_Eff_HET[i] = h2_d[i]->Integral()/h1_d[i]->Integral();
        Deuteron_Eff_UBT[i] = h3_d[i]->Integral()/h1_d[i]->Integral();

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        h1_e[i] = new TH1D(Form("h1_e[%d]",i),Form("h1_e[%d]",i),100,0,1);
        h2_e[i] = new TH1D(Form("h2_e[%d]",i),Form("h2_e[%d]",i),100,0,1);
        h3_e[i] = new TH1D(Form("h3_e[%d]",i),Form("h3_e[%d]",i),100,0,1);
        electron_tree->Draw(Form("Total_E/Energy>>h1_e[%d]",i),"","");
        electron_tree->Draw(Form("Total_E/Energy>>h2_e[%d]",i), HET,"");
        electron_tree->Draw(Form("Total_E/Energy>>h3_e[%d]",i), UBT,"");
        Electron_Eff_HET[i] = h2_e[i]->Integral()/h1_e[i]->Integral();
        Electron_Eff_UBT[i] = h3_e[i]->Integral()/h1_e[i]->Integral();

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        h1_h[i] = new TH1D(Form("h1_h[%d]",i),Form("h1_h[%d]",i),100,0,1);
        h2_h[i] = new TH1D(Form("h2_h[%d]",i),Form("h2_h[%d]",i),100,0,1);
        h3_h[i] = new TH1D(Form("h3_h[%d]",i),Form("h3_h[%d]",i),100,0,1);
        helium4_tree->Draw(Form("Total_E/Energy>>h1_h[%d]",i), "", "");
        helium4_tree->Draw(Form("Total_E/Energy>>h2_h[%d]",i),  HET,"");
        helium4_tree->Draw(Form("Total_E/Energy>>h3_h[%d]",i),  UBT,"");
        Helium4_Eff_HET[i] = h2_h[i]->Integral()/h1_h[i]->Integral();
        Helium4_Eff_UBT[i] = h3_h[i]->Integral()/h1_h[i]->Integral();
        cout << "Energy = " << int(Energy[i]) << " GeV !" << " He4 HET eff : " << h2_h[i]->Integral() <<" He4 UBT eff : " << h3_h[i]->Integral() <<endl;

        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        h1_H[i] = new TH1D(Form("h1_H[%d]",i),Form("h1_H[%d]",i),100,0,1);
        h2_H[i] = new TH1D(Form("h2_H[%d]",i),Form("h2_H[%d]",i),100,0,1);
        h3_H[i] = new TH1D(Form("h3_H[%d]",i),Form("h3_H[%d]",i),100,0,1);
        helium3_tree->Draw(Form("Total_E/Energy>>h1_H[%d]",i), "", "");
        helium3_tree->Draw(Form("Total_E/Energy>>h2_H[%d]",i),  HET,"");
        helium3_tree->Draw(Form("Total_E/Energy>>h3_H[%d]",i),  UBT,"");
        Helium3_Eff_HET[i] = h2_H[i]->Integral()/h1_H[i]->Integral();
        Helium3_Eff_UBT[i] = h3_H[i]->Integral()/h1_H[i]->Integral();
        // cout << "Energy = " << int(Energy[i]) << " GeV !" << " He3 HET eff : " << h2_H[i]->Integral() <<" He3 UBT eff : " << h3_H[i]->Integral() <<endl;

        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        h1_c[i] = new TH1D(Form("h1_c[%d]",i),Form("h1_c[%d]",i),100,0,1);
        h2_c[i] = new TH1D(Form("h2_c[%d]",i),Form("h2_c[%d]",i),100,0,1);
        h3_c[i] = new TH1D(Form("h3_c[%d]",i),Form("h3_c[%d]",i),100,0,1);
        carbon_tree->Draw(Form("Total_E/Energy>>h1_c[%d]",i), "", "");
        carbon_tree->Draw(Form("Total_E/Energy>>h2_c[%d]",i),  HET,"");
        carbon_tree->Draw(Form("Total_E/Energy>>h3_c[%d]",i),  UBT,"");
        Carbon_Eff_HET[i] = h2_c[i]->Integral()/h1_c[i]->Integral();
        Carbon_Eff_UBT[i] = h3_c[i]->Integral()/h1_c[i]->Integral();
    }
    

    auto gre_p_HET = new TGraphErrors(19, Energy, Proton_Eff_HET, Energy_Err, Uncertainty);
    auto gre_p_UBT = new TGraphErrors(19, Energy, Proton_Eff_UBT, Energy_Err, Uncertainty);
    
    gre_p_HET->SetTitle("High-Energy Trigger; Kinetic Energy (GeV); Trigger Efficiency");
    gre_p_HET->SetMarkerStyle(20);
    gre_p_HET->SetMarkerColor(kRed);
    gre_p_HET->SetLineColor(kRed);
    gre_p_UBT->SetTitle("Unbiased Trigger;Kinetic Energy (GeV); Trigger Efficiency");
    gre_p_UBT->SetMarkerStyle(20);
    gre_p_UBT->SetMarkerColor(kRed);
    gre_p_UBT->SetLineColor(kRed);

    auto gre_d_HET = new TGraphErrors(19, Energy, Deuteron_Eff_HET, Energy_Err, Uncertainty);
    auto gre_d_UBT = new TGraphErrors(19, Energy, Deuteron_Eff_UBT, Energy_Err, Uncertainty);

    gre_d_HET->SetMarkerStyle(21);
    gre_d_HET->SetMarkerColor(kBlue);
    gre_d_HET->SetLineColor(kBlue);
    gre_d_UBT->SetMarkerStyle(21);
    gre_d_UBT->SetMarkerColor(kBlue);
    gre_d_UBT->SetLineColor(kBlue);

    auto gre_e_HET = new TGraphErrors(19, Energy, Electron_Eff_HET, Energy_Err, Uncertainty);
    auto gre_e_UBT = new TGraphErrors(19, Energy, Electron_Eff_UBT, Energy_Err, Uncertainty);

    gre_e_HET->SetMarkerStyle(22);
    gre_e_HET->SetMarkerColor(kOrange-3);
    gre_e_HET->SetLineColor(kOrange-3);
    gre_e_UBT->SetMarkerStyle(22);
    gre_e_UBT->SetMarkerColor(kOrange-3);
    gre_e_UBT->SetLineColor(kOrange-3);


    auto gre_h_HET = new TGraphErrors(19, Energy, Helium4_Eff_HET, Energy_Err, Uncertainty);
    auto gre_h_UBT = new TGraphErrors(19, Energy, Helium4_Eff_UBT, Energy_Err, Uncertainty);
    gre_h_HET->SetMarkerStyle(23);
    gre_h_HET->SetMarkerColor(kGreen-3);
    gre_h_HET->SetLineColor(kGreen-3);
    gre_h_UBT->SetMarkerStyle(23);
    gre_h_UBT->SetMarkerColor(kGreen-3);
    gre_h_UBT->SetLineColor(kGreen-3);

    auto gre_H_HET = new TGraphErrors(19, Energy, Helium3_Eff_HET, Energy_Err, Uncertainty);
    auto gre_H_UBT = new TGraphErrors(19, Energy, Helium3_Eff_UBT, Energy_Err, Uncertainty);
    gre_H_HET->SetMarkerStyle(32);
    gre_H_HET->SetMarkerColor(kGreen-3);
    gre_H_HET->SetLineColor(kGreen-3);
    gre_H_UBT->SetMarkerStyle(32);
    gre_H_UBT->SetMarkerColor(kGreen-3);
    gre_H_UBT->SetLineColor(kGreen-3);

    auto gre_c_HET = new TGraphErrors(19, Energy, Carbon_Eff_HET, Energy_Err, Uncertainty);
    auto gre_c_UBT = new TGraphErrors(19, Energy, Carbon_Eff_UBT, Energy_Err, Uncertainty);
    gre_c_HET->SetMarkerStyle(23);
    gre_c_HET->SetMarkerColor(kMagenta);
    gre_c_HET->SetLineColor(kMagenta);
    gre_c_UBT->SetMarkerStyle(32);
    gre_c_UBT->SetMarkerColor(kMagenta);
    gre_c_UBT->SetLineColor(kMagenta);

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


    auto legend1 = new TLegend(0.60, 0.68, 0.88, 0.88);
    legend1->SetNColumns(2);
    legend1->AddEntry(gre_p_UBT, "UBT Proton", "ep");
    legend1->AddEntry(gre_d_UBT, "UBT Deuteron", "ep");
    legend1->AddEntry(gre_e_UBT, "UBT Electron", "ep");
    legend1->AddEntry(gre_h_UBT, "UBT Helium4", "ep");
    legend1->AddEntry(gre_H_UBT, "UBT Helium3", "ep");
    legend1->AddEntry(gre_c_UBT, "UBT Carbon", "ep");
    legend1->Draw();

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

    auto legend2 = new TLegend(0.62, 0.58, 0.88, 0.78);
    legend2->SetNColumns(2);
    legend2->AddEntry(gre_p_HET, "HET Proton", "ep");
    legend2->AddEntry(gre_d_HET, "HET Deuteron", "ep");
    legend2->AddEntry(gre_e_HET, "HET Electron", "ep");
    legend2->AddEntry(gre_h_HET, "HET Helium4", "ep");
    legend2->AddEntry(gre_H_HET, "HET Helium3", "ep");
    legend2->AddEntry(gre_c_HET, "HET Carbon", "ep");
    legend2->Draw();

}