void CrossSection_HE_Extract()
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

    double proton_constant[19]={0};     
    double proton_lambda[19]={0};      
    double proton_lambda_err[19]={0};   
    double deuteron_constant[19]={0};  
    double deuteron_lambda[19]={0};    
    double deuteron_lambda_err[19]={0};
    double electron_constant[19]={0};
    double electron_lambda[19]={0};
    double electron_lambda_err[19]={0};
    double helium4_constant[19]={0};
    double helium4_lambda[19]={0};
    double helium4_lambda_err[19]={0};
    double helium3_constant[19]={0};
    double helium3_lambda[19]={0};
    double helium3_lambda_err[19]={0};
    double carbon_constant[19]={0};
    double carbon_lambda[19]={0};
    double carbon_lambda_err[19]={0};

    double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    TH1D *h1_p[19]; TH1D *hC_p[19]; TF1  *fitFunc_p[19];
    TH1D *h1_d[19]; TH1D *hC_d[19]; TF1  *fitFunc_d[19];
    TH1D *h1_e[19]; TH1D *hC_e[19]; TF1  *fitFunc_e[19];
    TH1D *h1_h[19]; TH1D *hC_h[19]; TF1  *fitFunc_h[19];
    TH1D *h1_H[19]; TH1D *hC_H[19]; TF1  *fitFunc_H[19];
    TH1D *h1_c[19]; TH1D *hC_c[19]; TF1  *fitFunc_c[19];

    TCut EM = "First_Type==0";
    TCut HD = "First_Had_Depth>=0";
    TCut HI = "First_Had_Depth>=0 && First_Had_Type==1";
    TCut HE = "First_Had_Depth>=0 && First_Had_Type==2";

    for (int i =0; i < 19; i++)
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        cout << "Energy = " << int(Energy[i]) << " GeV !" << endl;

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        h1_p[i] = new TH1D(Form("h1_p[%d]",i),Form("h1_p[%d]",i),35,0,35);
        hC_p[i] = new TH1D(Form("hC_p[%d]",i),Form("hC_p[%d]",i),35,0,35);
        proton_tree->Draw(Form("First_Had_Depth>>h1_p[%d]",i),HE,"");
        fitFunc_p[i] = new TF1(Form("fitFunc_p[%d]",i), "[0]*exp(-x/[1])",0, 10);
        fitFunc_p[i]->SetParameters(1e4, 15); 

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        h1_d[i] = new TH1D(Form("h1_d[%d]",i),Form("h1_d[%d]",i),35,0,35);
        hC_d[i] = new TH1D(Form("hC_d[%d]",i),Form("hC_d[%d]",i),35,0,35);
        deuteron_tree->Draw(Form("First_Had_Depth>>h1_d[%d]",i),HE, "");
        fitFunc_d[i] = new TF1(Form("fitFunc_d[%d]",i), "[0]*exp(-x/[1])",0, 10);
        fitFunc_d[i]->SetParameters(1e4, 15); 

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        h1_e[i] = new TH1D(Form("h1_e[%d]",i),Form("h1_e[%d]",i),35,0,35);
        hC_e[i] = new TH1D(Form("hC_e[%d]",i),Form("hC_e[%d]",i),35,0,35);

        fitFunc_e[i] = new TF1(Form("fitFunc_e[%d]",i), "[0]*exp(-x/[1])",0, 10);
        fitFunc_e[i]->SetParameters(1e4, 15); 
        electron_tree->Draw(Form("First_Had_Depth>>h1_e[%d]",i),HE, "");

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        h1_h[i] = new TH1D(Form("h1_h[%d]",i),Form("h1_h[%d]",i),35,0,35);
        hC_h[i] = new TH1D(Form("hC_h[%d]",i),Form("hC_h[%d]",i),35,0,35);
        helium4_tree->Draw(Form("First_Had_Depth>>h1_h[%d]",i),HE, "");
        fitFunc_h[i] = new TF1(Form("fitFunc_h[%d]",i), "[0]*exp(-x/[1])",0, 10);
        fitFunc_h[i]->SetParameters(1e4, 15); 

        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        h1_H[i] = new TH1D(Form("h1_H[%d]",i),Form("h1_H[%d]",i),35,0,35);
        hC_H[i] = new TH1D(Form("hC_H[%d]",i),Form("hC_H[%d]",i),35,0,35);
        helium3_tree->Draw(Form("First_Had_Depth>>h1_H[%d]",i),HE, "");
        fitFunc_H[i] = new TF1(Form("fitFunc_H[%d]",i), "[0]*exp(-x/[1])",0, 10);
        fitFunc_H[i]->SetParameters(1e4, 15); 

        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        h1_c[i] = new TH1D(Form("h1_c[%d]",i),Form("h1_c[%d]",i),35,0,35);
        hC_c[i] = new TH1D(Form("hC_c[%d]",i),Form("hC_c[%d]",i),35,0,35);
        carbon_tree->Draw(Form("First_Had_Depth>>h1_c[%d]",i),HE, "");
        fitFunc_c[i] = new TF1(Form("fitFunc_c[%d]",i), "[0]*exp(-x/[1])",0, 10);
        fitFunc_c[i]->SetParameters(1e4, 1); 

        for (int jj=1 ; jj<=35 ;jj++)
        {
            hC_p[i]->SetBinContent(jj,1e4-h1_p[i]->Integral(0,jj));
            hC_d[i]->SetBinContent(jj,1e4-h1_d[i]->Integral(0,jj));
            hC_e[i]->SetBinContent(jj,1e4-h1_e[i]->Integral(0,jj));
            hC_h[i]->SetBinContent(jj,1e4-h1_h[i]->Integral(0,jj));
            hC_H[i]->SetBinContent(jj,1e4-h1_H[i]->Integral(0,jj));
            hC_c[i]->SetBinContent(jj,1e4-h1_c[i]->Integral(0,jj));
        }

        auto c1 = new TCanvas("c1","c1",2000,1500);
        c1->Clear();
        c1->Divide(3,2);
        c1->cd(1);
        proton_file->cd();
        //gPad->SetLogy();
        hC_p[i]->Draw("hist");
        hC_p[i]->SetTitle(Form("%d GeV Proton Depth Distribution; Depth (mm);N_{survive}", int(Energy[i])));
        hC_p[i]->Fit(fitFunc_p[i], "RQ"); // 进行拟合
        gStyle->SetOptFit(1111);

        fitFunc_p[i]->Draw("same");
        proton_constant[i]   = fitFunc_p[i]->GetParameter(0);
        proton_lambda[i]     = fitFunc_p[i]->GetParameter(1);
        proton_lambda_err[i] = fitFunc_p[i]->GetParError(1);
        Proton_Section[i] = 1 / (proton_lambda[i] * n_BGO) * 1e25; // Convert to barn
        Proton_Section_Err[i] = Proton_Section[i] * proton_lambda_err[i] / proton_lambda[i];
        cout << "Proton " << Energy[i] << " GeV: Sigma = " << Proton_Section[i] << " ± " << Proton_Section_Err[i] << " barn" << endl;

        // Process Deuteron Data
        c1->cd(2);
        deuteron_file->cd();
        //gPad->SetLogy();
        hC_d[i]->Draw("hist");
        hC_d[i]->SetTitle(Form("%d GeV Deuteron Depth Distribution; Depth (mm);N_{survive}", int(Energy[i])));
        hC_d[i]->Fit(fitFunc_d[i], "RQ");
        gStyle->SetOptFit(1111);

        fitFunc_d[i]->Draw("same");
        deuteron_constant[i]   = fitFunc_d[i]->GetParameter(0);
        deuteron_lambda[i]     = fitFunc_d[i]->GetParameter(1);
        deuteron_lambda_err[i] = fitFunc_d[i]->GetParError(1);
        Deuteron_Section[i] = 1 / (deuteron_lambda[i] * n_BGO) * 1e25;
        Deuteron_Section_Err[i] = Deuteron_Section[i] * deuteron_lambda_err[i] / deuteron_lambda[i];
        cout << "Deuteron " << Energy[i] << " GeV: Sigma = " << Deuteron_Section[i] << " ± " << Deuteron_Section_Err[i] << " barn" << endl;

        c1->cd(3);
        electron_file->cd();
        //gPad->SetLogy();
        hC_e[i]->Draw("hist");
        hC_e[i]->SetTitle(Form("%d GeV Electron Depth Distribution; Depth (mm);N_{survive}", int(Energy[i])));
        hC_e[i]->Fit(fitFunc_e[i], "RQ");
        gStyle->SetOptFit(1111);

        fitFunc_e[i]->Draw("same");
        electron_constant[i]   = fitFunc_e[i]->GetParameter(0);
        electron_lambda[i]     = fitFunc_e[i]->GetParameter(1);
        electron_lambda_err[i] = fitFunc_e[i]->GetParError(1);
        Electron_Section[i] = 1 / (electron_lambda[i] * n_BGO) * 1e25;
        Electron_Section_Err[i] = Electron_Section[i] * electron_lambda_err[i] / electron_lambda[i];
        cout << "Electron " << Energy[i] << " GeV: Sigma = " << Electron_Section[i] << " ± " << Electron_Section_Err[i] << " barn" << endl;
        
        c1->cd(4);
        helium4_file->cd();
        //gPad->SetLogy();
        hC_h[i]->Draw("hist");
        hC_h[i]->SetTitle(Form("%d GeV Helium4 Depth Distribution; Depth (mm);N_{survive}", int(Energy[i])));
        hC_h[i]->Fit(fitFunc_h[i], "RQ");
        gStyle->SetOptFit(1111);

        fitFunc_h[i]->Draw("same");
        helium4_constant[i]   = fitFunc_h[i]->GetParameter(0);
        helium4_lambda[i]     = fitFunc_h[i]->GetParameter(1);
        helium4_lambda_err[i] = fitFunc_h[i]->GetParError(1);

        Helium4_Section[i] = 1 / (helium4_lambda[i] * n_BGO) * 1e25;
        Helium4_Section_Err[i] = Helium4_Section[i] * helium4_lambda_err[i] / helium4_lambda[i];
        cout << "Helium4 " << Energy[i] << " GeV: Sigma = " << Helium4_Section[i] << " ± " << Helium4_Section_Err[i] << " barn" << endl;

        c1->cd(5);
        helium3_file->cd();
        //gPad->SetLogy();
        hC_H[i]->Draw("hist");
        hC_H[i]->SetTitle(Form("%d GeV Helium3 Depth Distribution; Depth (mm);N_{survive}", int(Energy[i])));
        hC_H[i]->Fit(fitFunc_H[i], "RQ");
        gStyle->SetOptFit(1111);
        fitFunc_H[i]->Draw("same");
        helium3_constant[i]   = fitFunc_H[i]->GetParameter(0);
        helium3_lambda[i]     = fitFunc_H[i]->GetParameter(1);
        helium3_lambda_err[i] = fitFunc_H[i]->GetParError(1);        
        Helium3_Section[i] = 1 / (helium3_lambda[i] * n_BGO) * 1e25;
        Helium3_Section_Err[i] = Helium3_Section[i] * helium3_lambda_err[i] / helium3_lambda[i];
        cout << "Helium3 " << Energy[i] << " GeV: Sigma = " << Helium3_Section[i] << " ± " << Helium3_Section_Err[i] << " barn" << endl;

        c1->cd(6);
        carbon_file->cd();
        //gPad->SetLogy();
        hC_c[i]->Draw("hist");
        hC_c[i]->SetTitle(Form("%d GeV Carbon Depth Distribution; Depth (mm);N_{survive}", int(Energy[i])));
        hC_c[i]->Fit(fitFunc_c[i], "RQ");
        gStyle->SetOptFit(1111);

        fitFunc_c[i]->Draw("same");
        carbon_constant[i]   = fitFunc_c[i]->GetParameter(0);
        carbon_lambda[i]     = fitFunc_c[i]->GetParameter(1);
        carbon_lambda_err[i] = fitFunc_c[i]->GetParError(1);

        Carbon_Section[i] = 1 / (carbon_lambda[i] * n_BGO) * 1e25;
        Carbon_Section_Err[i] = Carbon_Section[i] * carbon_lambda_err[i] / carbon_lambda[i];
        cout << "Carbon " << Energy[i] << " GeV: Sigma = " << Carbon_Section[i] << " ± " << Carbon_Section_Err[i] << " barn" << endl;
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4c/Script/Theory/CrossSection_HE_%dGeV.pdf",int(Energy[i])) );
    }

    auto file1 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/proton_BGO.root");
    auto hist1 = (TH1D*)file1->Get("h6");  //
    auto hist6 = (TH1D*)file1->Get("h3");  

    auto file2 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/deuteron_BGO.root");
    auto hist2 = (TH1D*)file2->Get("h6");  
    auto hist7 = (TH1D*)file2->Get("h3");  

    auto file3 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/alpha_BGO.root");
    auto hist3 = (TH1D*)file3->Get("h6");  
    auto hist8 = (TH1D*)file3->Get("h3");  


    auto file4 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/He3_BGO.root");
    auto hist4 = (TH1D*)file4->Get("h6");  
    auto hist9 = (TH1D*)file4->Get("h3");  


    auto file5 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/Carbon_BGO.root");
    auto hist5 = (TH1D*)file5->Get("h6");  
    auto hist0 = (TH1D*)file5->Get("h3");  

    int Nbins = hist1->GetNbinsX();
    double KN_Energy[60] = {0};
    double Proton[60] = {0};
    double Deuteron[60] = {0};
    double Helium3[60] = {0};
    double Helium4[60] = {0};
    double Carbon[60] = {0};

    double KN_Proton[60] = {0};
    double KN_Deuteron[60] = {0};
    double KN_Helium3[60] = {0};
    double KN_Helium4[60] = {0};
    double KN_Carbon[60] = {0};
    for (int ii = 0 ; ii < Nbins ; ii++)
    {
        KN_Energy[ii]   = pow(10,hist1->GetBinCenter(ii+1)-3);
        file1->cd();Proton[ii]   = hist1->GetBinContent(ii+1); KN_Proton[ii]   = hist6->GetBinContent(ii+1)* 10.;
        file2->cd();Deuteron[ii] = hist2->GetBinContent(ii+1); KN_Deuteron[ii] = hist7->GetBinContent(ii+1)* 10.;
        file3->cd();Helium3[ii]  = hist3->GetBinContent(ii+1); KN_Helium3[ii]  = hist8->GetBinContent(ii+1)* 10.;
        file4->cd();Helium4[ii]  = hist4->GetBinContent(ii+1); KN_Helium4[ii]  = hist9->GetBinContent(ii+1)* 10.;
        file5->cd();Carbon[ii]   = hist5->GetBinContent(ii+1); KN_Carbon[ii]   = hist0->GetBinContent(ii+1)* 10.;
    }
    auto gre1 = new TGraph(Nbins,KN_Energy, Proton);
    auto gre2 = new TGraph(Nbins,KN_Energy, Deuteron);
    auto gre3 = new TGraph(Nbins,KN_Energy, Helium3);
    auto gre4 = new TGraph(Nbins,KN_Energy, Helium4);
    auto gre5 = new TGraph(Nbins,KN_Energy, Carbon);

    auto gre6 = new TGraph(Nbins,KN_Energy, KN_Proton);    
    auto gre7 = new TGraph(Nbins,KN_Energy, KN_Deuteron);   
    auto gre8 = new TGraph(Nbins,KN_Energy, KN_Helium3);  
    auto gre9 = new TGraph(Nbins,KN_Energy, KN_Helium4);  
    auto gre0 = new TGraph(Nbins,KN_Energy, KN_Carbon);      
    
    gre1->SetLineColor(kRed);
    gre2->SetLineColor(kBlue);
    gre3->SetLineColor(kGreen+3);
    gre4->SetLineColor(kGreen+3);
    gre4->SetLineStyle(2);
    gre5->SetLineColor(kMagenta);

    gre6->SetLineColor(kRed);
    gre7->SetLineColor(kBlue);
    gre8->SetLineColor(kGreen+3);
    gre9->SetLineColor(kGreen+3);
    gre9->SetLineStyle(2);
    gre0->SetLineColor(kMagenta);

    auto gr_p = new TGraphErrors(19, Energy, proton_lambda, Energy_Err, proton_lambda_err);        gr_p->SetMarkerStyle(20);gr_p->SetMarkerColor(kRed);gr_p->SetLineColor(kRed);
    auto gr_d = new TGraphErrors(19, Energy, deuteron_lambda, Energy_Err, deuteron_lambda_err);    gr_d->SetMarkerStyle(21);gr_d->SetMarkerColor(kBlue);gr_d->SetLineColor(kBlue);
    auto gr_h = new TGraphErrors(19, Energy, helium4_lambda, Energy_Err, helium4_lambda_err);      gr_h->SetMarkerStyle(23);gr_h->SetMarkerColor(kGreen-3);gr_h->SetLineColor(kGreen-3);
    auto gr_H = new TGraphErrors(19, Energy, helium3_lambda, Energy_Err, helium3_lambda_err);      gr_H->SetMarkerStyle(32);gr_H->SetMarkerColor(kGreen-3);gr_H->SetLineColor(kGreen-3);
    auto gr_c = new TGraphErrors(19, Energy, carbon_lambda, Energy_Err, carbon_lambda_err);        gr_c->SetMarkerStyle(33);gr_c->SetMarkerColor(kMagenta);gr_c->SetLineColor(kMagenta);
 
    gr_p->SetTitle("Fitted Hadronic Collision Length vs Energy; Kinetic Energy (GeV); Length (mm)");


    auto gre_p = new TGraphErrors(19, Energy, Proton_Section, Energy_Err, Proton_Section_Err);
    gre_p->SetTitle("Fitted Hadronic Collision Cross Section vs Energy; Kinetic Energy (GeV); Cross Section per Atom (barn)");
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

    auto c0 = new TCanvas("c0","c0",900,900);
    c0->Clear();
    c0->cd();
    gPad->SetGrid(1, 1);      // Vertical ,Horizontal grid
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    gre_p->Draw("AP");
    gre_p->GetXaxis()->SetLimits(1,1e5);
    gre_p->GetYaxis()->SetRangeUser(0, 2);
    gre_p->GetYaxis()->SetNdivisions(505);
    gre_p->GetYaxis()->SetLabelSize(0.05);
    gre_p->GetYaxis()->SetTitleSize(0.05);
    gre_p->GetYaxis()->SetTitleOffset(0.9);

    gre_d->Draw("PSAME");
    gre_h->Draw("PSAME");
    gre_H->Draw("PSAME");
    gre_c->Draw("PSAME");

    gre1->Draw("LSAME");
    gre2->Draw("LSAME");
    gre3->Draw("LSAME");
    gre4->Draw("LSAME");
    gre5->Draw("LSAME");

    auto legend1 = new TLegend(0.60, 0.12, 0.88, 0.32);
    legend1->SetNColumns(2);
    legend1->AddEntry(gre_p, "Proton", "ep");
    legend1->AddEntry(gre_d, "Deuteron", "ep");
    legend1->AddEntry(gre_h, "Helium4", "ep");
    legend1->AddEntry(gre_H, "Helium3", "ep");
    legend1->AddEntry(gre_c, "Carbon", "ep");
    legend1->Draw();

    auto legend2 = new TLegend(0.12, 0.12, 0.88, 0.32);
    legend2->SetNColumns(2);
    legend2->AddEntry(gre1, "GEANT4 Proton", "l");
    legend2->AddEntry(gre2, "GEANT4 Deuteron", "l");
    legend2->AddEntry(gre3, "GEANT4 Helium4", "l");
    legend2->AddEntry(gre4, "GEANT4 Helium3", "l");
    legend2->AddEntry(gre5, "GEANT4 Carbon", "l");
    legend2->Draw();
    c0->SaveAs("/Users/xiongzheng/software/B4/B4c/Script/Theory/CrossSection_HE.pdf");

    auto c2 = new TCanvas("c2","c2",900,900);
    c2->cd();
    gPad->SetGrid(1, 1);      // Vertical ,Horizontal grid
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    gr_p->Draw("AP");
    gr_p->GetXaxis()->SetLimits(1,1e5);
    gr_p->GetYaxis()->SetRangeUser(0, 250);
    gr_p->GetYaxis()->SetNdivisions(505);
    gr_p->GetYaxis()->SetLabelSize(0.05);
    gr_p->GetYaxis()->SetTitleSize(0.05);
    gr_p->GetYaxis()->SetTitleOffset(0.9);
    gr_d->Draw("PSAME");
    gr_h->Draw("PSAME");
    gr_H->Draw("PSAME");
    gr_c->Draw("PSAME");

    gre6->Draw("LSAME");
    gre7->Draw("LSAME");
    gre8->Draw("LSAME");
    gre9->Draw("LSAME");
    gre0->Draw("LSAME");

    auto legend3 = new TLegend(0.60, 0.12, 0.88, 0.32);
    legend3->SetNColumns(2);
    legend3->AddEntry(gr_p, "Proton", "ep");
    legend3->AddEntry(gr_d, "Deuteron", "ep");
    legend3->AddEntry(gr_h, "Helium4", "ep");
    legend3->AddEntry(gr_H, "Helium3", "ep");
    legend3->AddEntry(gr_c, "Carbon", "ep");
    legend3->Draw();

    auto legend4 = new TLegend(0.12, 0.12, 0.52, 0.32);
    legend4->SetNColumns(2);
    legend4->AddEntry(gre6, "GEANT4 Proton", "l");
    legend4->AddEntry(gre7, "GEANT4 Deuteron", "l");
    legend4->AddEntry(gre8, "GEANT4 Helium4", "l");
    legend4->AddEntry(gre9, "GEANT4 Helium3", "l");
    legend4->AddEntry(gre0, "GEANT4 Carbon", "l");
    legend4->Draw();

    c0->SaveAs("/Users/xiongzheng/software/B4/B4c/Script/Theory/FreePath_HE.pdf");

}