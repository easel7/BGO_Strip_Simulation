void CrossSection_HI()
{
    double Proton_Section[19]={0};
    double Proton_Section_Err[19]={0};
    double Deuteron_Section[19]={0};
    double Deuteron_Section_Err[19]={0};
    double Electron_Section[19]={0};
    double Electron_Section_Err[19]={0};
    double Helium4_Section[19]={0};
    double Helium4_Section_Err[19]={0};
    double Helium3_Section[19]={0};
    double Helium3_Section_Err[19]={0};
    double Carbon_Section[19]={0};
    double Carbon_Section_Err[19]={0};
    double Energy[19]={0};
    double Energy_Err[19]={0};
    double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    TH1D *h1_p[19];    TF1  *fitFunc_p[19];
    TH1D *h1_d[19];    TF1  *fitFunc_d[19];
    TH1D *h1_e[19];    TF1  *fitFunc_e[19];
    TH1D *h1_h[19];    TF1  *fitFunc_h[19];
    TH1D *h1_H[19];    TF1  *fitFunc_H[19];
    TH1D *h1_c[19];    TF1  *fitFunc_c[19];

    TCut EM = "First_Type==0";
    TCut HD = "First_Had_Depth>=0";
    TCut HI = "First_Had_Depth>=0 && First_Had_Type==1";

    for (int i = 0; i < 19; i++)
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        cout << "Energy = " << int(Energy[i]) << " GeV !" << endl;

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        h1_p[i] = new TH1D(Form("h1_p[%d]",i),Form("h1_p[%d]",i),30,0,300);
        proton_tree->Draw(Form("First_Had_Depth>>h1_p[%d]",i),HI,"");
        fitFunc_p[i] = new TF1(Form("fitFunc_p[%d]",i), "[0]*exp(-x/[1])", 0, 200);
        fitFunc_p[i]->SetParameters(600, 15); 

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        h1_d[i] = new TH1D(Form("h1_d[%d]",i),Form("h1_d[%d]",i),35,0,350);
        deuteron_tree->Draw(Form("First_Had_Depth>>h1_d[%d]",i),HI, "");
        fitFunc_d[i] = new TF1(Form("fitFunc_d[%d]",i), "[0]*exp(-x/[1])", 0, 200);
        fitFunc_d[i]->SetParameters(100, 15); 

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        h1_e[i] = new TH1D(Form("h1_e[%d]",i),Form("h1_e[%d]",i),50,0,100);
        fitFunc_e[i] = new TF1(Form("fitFunc_e[%d]",i), "[0]*exp(-x/[1])", 0, 60);
        fitFunc_e[i]->SetParameters(100, 15); 
        electron_tree->Draw(Form("First_Had_Depth>>h1_e[%d]",i),HI, "");

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        h1_h[i] = new TH1D(Form("h1_h[%d]",i),Form("h1_h[%d]",i),35,0,350);
        helium4_tree->Draw(Form("First_Had_Depth>>h1_h[%d]",i),HI, "");
        fitFunc_h[i] = new TF1(Form("fitFunc_h[%d]",i), "[0]*exp(-x/[1])", 0, 200);
        fitFunc_h[i]->SetParameters(100, 15); 

        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        h1_H[i] = new TH1D(Form("h1_H[%d]",i),Form("h1_H[%d]",i),35,0,350);
        helium3_tree->Draw(Form("First_Had_Depth>>h1_H[%d]",i),HI, "");
        fitFunc_H[i] = new TF1(Form("fitFunc_H[%d]",i), "[0]*exp(-x/[1])", 0, 200);
        fitFunc_H[i]->SetParameters(100, 15); 

        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        h1_c[i] = new TH1D(Form("h1_c[%d]",i),Form("h1_c[%d]",i),35,0,350);
        carbon_tree->Draw(Form("First_Had_Depth>>h1_c[%d]",i),HI, "");
        fitFunc_c[i] = new TF1(Form("fitFunc_c[%d]",i), "[0]*exp(-x/[1])", 0, 200);
        fitFunc_c[i]->SetParameters(1000, 1); 

        auto c1 = new TCanvas("c1","c1",2000,1500);
        c1->Clear();
        c1->Divide(3,2);
        c1->cd(1);
        proton_file->cd();
        gPad->SetLogy();
        h1_p[i]->Draw("hist");
        h1_p[i]->SetTitle(Form("%d GeV Proton Depth Distribution; Depth (mm); Counts", int(Energy[i])));
        h1_p[i]->Fit(fitFunc_p[i], "RQ"); // 进行拟合
        fitFunc_p[i]->Draw("same");
        double proton_constant   = fitFunc_p[i]->GetParameter(0);
        double proton_lambda     = fitFunc_p[i]->GetParameter(1);
        double proton_lambda_err = fitFunc_p[i]->GetParError(1);
        Proton_Section[i] = 1 / (proton_lambda * n_BGO) * 1e25; // Convert to barn
        Proton_Section_Err[i] = Proton_Section[i] * proton_lambda_err / proton_lambda;
        cout << "Proton " << Energy[i] << " GeV: Sigma = " << Proton_Section[i] << " ± " << Proton_Section_Err[i] << " barn" << endl;

        // Process Deuteron Data
        c1->cd(2);
        deuteron_file->cd();
        gPad->SetLogy();
        h1_d[i]->Draw("hist");
        h1_d[i]->SetTitle(Form("%d GeV Deuteron Depth Distribution; Depth (mm); Counts", int(Energy[i])));
        h1_d[i]->Fit(fitFunc_d[i], "RQ");
        fitFunc_d[i]->Draw("same");
        double deuteron_constant   = fitFunc_d[i]->GetParameter(0);
        double deuteron_lambda     = fitFunc_d[i]->GetParameter(1);
        double deuteron_lambda_err = fitFunc_d[i]->GetParError(1);

        Deuteron_Section[i] = 1 / (deuteron_lambda * n_BGO) * 1e25;
        Deuteron_Section_Err[i] = Deuteron_Section[i] * deuteron_lambda_err / deuteron_lambda;
        cout << "Deuteron " << Energy[i] << " GeV: Sigma = " << Deuteron_Section[i] << " ± " << Deuteron_Section_Err[i] << " barn" << endl;

        c1->cd(3);
        electron_file->cd();
        gPad->SetLogy();
        h1_e[i]->Draw("hist");
        h1_e[i]->SetTitle(Form("%d GeV Electron Depth Distribution; Depth (mm); Counts", int(Energy[i])));
        h1_e[i]->Fit(fitFunc_e[i], "RQ");
        fitFunc_e[i]->Draw("same");
        double electron_constant   = fitFunc_e[i]->GetParameter(0);
        double electron_lambda     = fitFunc_e[i]->GetParameter(1);
        double electron_lambda_err = fitFunc_e[i]->GetParError(1);

        Electron_Section[i] = 1 / (electron_lambda * n_BGO) * 1e25;
        Electron_Section_Err[i] = Electron_Section[i] * electron_lambda_err / electron_lambda;
        cout << "Electron " << Energy[i] << " GeV: Sigma = " << Electron_Section[i] << " ± " << Electron_Section_Err[i] << " barn" << endl;
        
        c1->cd(4);
        helium4_file->cd();
        gPad->SetLogy();
        h1_h[i]->Draw("hist");
        h1_h[i]->SetTitle(Form("%d GeV Helium4 Depth Distribution; Depth (mm); Counts", int(Energy[i])));
        h1_h[i]->Fit(fitFunc_h[i], "RQ");
        fitFunc_h[i]->Draw("same");
        double helium4_constant   = fitFunc_h[i]->GetParameter(0);
        double helium4_lambda     = fitFunc_h[i]->GetParameter(1);
        double helium4_lambda_err = fitFunc_h[i]->GetParError(1);

        Helium4_Section[i] = 1 / (helium4_lambda * n_BGO) * 1e25;
        Helium4_Section_Err[i] = Helium4_Section[i] * helium4_lambda_err / helium4_lambda;
        cout << "Helium4 " << Energy[i] << " GeV: Sigma = " << Helium4_Section[i] << " ± " << Helium4_Section_Err[i] << " barn" << endl;

        c1->cd(5);
        helium3_file->cd();
        gPad->SetLogy();
        h1_H[i]->Draw("hist");
        h1_H[i]->SetTitle(Form("%d GeV Helium3 Depth Distribution; Depth (mm); Counts", int(Energy[i])));
        h1_H[i]->Fit(fitFunc_H[i], "RQ");
        fitFunc_H[i]->Draw("same");
        double helium3_constant   = fitFunc_H[i]->GetParameter(0);
        double helium3_lambda     = fitFunc_H[i]->GetParameter(1);
        double helium3_lambda_err = fitFunc_H[i]->GetParError(1);

        Helium3_Section[i] = 1 / (helium3_lambda * n_BGO) * 1e25;
        Helium3_Section_Err[i] = Helium3_Section[i] * helium3_lambda_err / helium3_lambda;
        cout << "Helium3 " << Energy[i] << " GeV: Sigma = " << Helium3_Section[i] << " ± " << Helium3_Section_Err[i] << " barn" << endl;

        c1->cd(6);
        carbon_file->cd();
        gPad->SetLogy();

        h1_c[i]->Draw("hist");
        h1_c[i]->SetTitle(Form("%d GeV Carbon Depth Distribution; Depth (mm); Counts", int(Energy[i])));
        h1_c[i]->Fit(fitFunc_c[i], "RQ");
        fitFunc_c[i]->Draw("same");
        double carbon_constant   = fitFunc_c[i]->GetParameter(0);
        double carbon_lambda     = fitFunc_c[i]->GetParameter(1);
        double carbon_lambda_err = fitFunc_c[i]->GetParError(1);

        Carbon_Section[i] = 1 / (carbon_lambda * n_BGO) * 1e25;
        Carbon_Section_Err[i] = Carbon_Section[i] * carbon_lambda_err / carbon_lambda;
        cout << "Carbon " << Energy[i] << " GeV: Sigma = " << Carbon_Section[i] << " ± " << Carbon_Section_Err[i] << " barn" << endl;


        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4c/Script/CrossS/CrossSection_HI_%dGeV.pdf",int(Energy[i])) );
    }
    
    auto c0 = new TCanvas("c0","c0",900,900);
    auto gre_p = new TGraphErrors(19, Energy, Proton_Section, Energy_Err, Proton_Section_Err);
    gre_p->SetTitle("Fitted Hadronic Interaction Cross Section vs Energy; Kinetic Energy (GeV); Cross Section per Atom (barn)");
    gre_p->SetMarkerStyle(20);
    gre_p->SetMarkerColor(kRed);
    gre_p->SetLineColor(kRed);

    auto gre_d = new TGraphErrors(19, Energy, Deuteron_Section, Energy_Err, Deuteron_Section_Err);
    gre_d->SetTitle("; Kinetic Energy (GeV); Cross Section per Atom (barn)");
    gre_d->SetMarkerStyle(21);
    gre_d->SetMarkerColor(kBlue);
    gre_d->SetLineColor(kBlue);

    auto gre_e = new TGraphErrors(19, Energy, Electron_Section, Energy_Err, Electron_Section_Err);
    gre_e->SetMarkerStyle(22);
    gre_e->SetMarkerColor(kOrange-3);
    gre_e->SetLineColor(kOrange-3);

    auto gre_h = new TGraphErrors(19, Energy, Helium4_Section, Energy_Err, Helium4_Section_Err);
    gre_h->SetMarkerStyle(23);
    gre_h->SetMarkerColor(kGreen-3);
    gre_h->SetLineColor(kGreen-3);

    auto gre_H = new TGraphErrors(19, Energy, Helium3_Section, Energy_Err, Helium3_Section_Err);
    gre_H->SetMarkerStyle(32);
    gre_H->SetMarkerColor(kGreen-3);
    gre_H->SetLineColor(kGreen-3);

    auto gre_c = new TGraphErrors(19, Energy, Carbon_Section, Energy_Err, Carbon_Section_Err);
    gre_c->SetMarkerStyle(33);
    gre_c->SetMarkerColor(kMagenta);
    gre_c->SetLineColor(kMagenta);

    c0->Clear();
    c0->cd();
    gPad->SetGrid(1, 1);      // Vertical ,Horizontal grid
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    gre_p->Draw("AP");

    gre_p->GetYaxis()->SetRangeUser(0, 1.7);
    gre_p->GetYaxis()->SetNdivisions(505);
    gre_p->GetYaxis()->SetLabelSize(0.05);
    gre_p->GetYaxis()->SetTitleSize(0.05);
    gre_p->GetYaxis()->SetTitleOffset(0.9);

    gre_d->Draw("PSAME");
    // gre_e->Draw("PSAME");
    gre_h->Draw("PSAME");
    gre_H->Draw("PSAME");
    gre_c->Draw("PSAME");

    auto legend1 = new TLegend(0.60, 0.12, 0.88, 0.32);
    legend1->SetNColumns(2);
    legend1->AddEntry(gre_p, "Proton", "ep");
    legend1->AddEntry(gre_d, "Deuteron", "ep");
    // legend1->AddEntry(gre_e, "Electron", "ep");
    legend1->AddEntry(gre_h, "Helium4", "ep");
    legend1->AddEntry(gre_H, "Helium3", "ep");
    legend1->AddEntry(gre_c, "Carbon", "ep");
    legend1->Draw();
    c0->SaveAs("/Users/xiongzheng/software/B4/B4c/Script/CrossS/CrossSection_HI.pdf");


}