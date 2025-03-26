void Extract()
{
    auto file1 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/proton_BGO.root");
    auto hist1 = (TH1D*)file1->Get("h4");  
    auto hist6 = (TH1D*)file1->Get("h7");  

    auto file2 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/deuteron_BGO.root");
    auto hist2 = (TH1D*)file2->Get("h4");  
    auto hist7 = (TH1D*)file2->Get("h7");  

    auto file3 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/alpha_BGO.root");
    auto hist3 = (TH1D*)file3->Get("h4");  
    auto hist8 = (TH1D*)file3->Get("h7");  


    auto file4 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/He3_BGO.root");
    auto hist4 = (TH1D*)file4->Get("h4");  
    auto hist9 = (TH1D*)file4->Get("h7");  

    auto file5 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/Carbon_BGO.root");
    auto hist5 = (TH1D*)file5->Get("h4");  
    auto hist0 = (TH1D*)file5->Get("h7");  

    int Nbins = hist1->GetNbinsX();
    double Energy[60] = {0};
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
        Energy[ii]   = pow(10,hist1->GetBinCenter(ii+1)-3);
        file1->cd();Proton[ii]   = hist1->GetBinContent(ii+1); KN_Proton[ii]   = hist6->GetBinContent(ii+1);
        file2->cd();Deuteron[ii] = hist2->GetBinContent(ii+1); KN_Deuteron[ii] = hist7->GetBinContent(ii+1);
        file3->cd();Helium3[ii]  = hist3->GetBinContent(ii+1); KN_Helium3[ii]  = hist8->GetBinContent(ii+1);
        file4->cd();Helium4[ii]  = hist4->GetBinContent(ii+1); KN_Helium4[ii]  = hist9->GetBinContent(ii+1);
        file5->cd();Carbon[ii]   = hist5->GetBinContent(ii+1); KN_Carbon[ii]   = hist0->GetBinContent(ii+1);
    }
    auto gre1 = new TGraph(Nbins,Energy, Proton);
    auto gre2 = new TGraph(Nbins,Energy, Deuteron);
    auto gre3 = new TGraph(Nbins,Energy, Helium3  );
    auto gre4 = new TGraph(Nbins,Energy, Helium4  );
    auto gre5 = new TGraph(Nbins,Energy, Carbon  );

    auto gre6 = new TGraph(Nbins,Energy, KN_Proton);
    auto gre7 = new TGraph(Nbins,Energy, KN_Deuteron);
    auto gre8 = new TGraph(Nbins,Energy, KN_Helium3  );
    auto gre9 = new TGraph(Nbins,Energy, KN_Helium4  );
    auto gre0 = new TGraph(Nbins,Energy, KN_Carbon  );

    auto c1 = new TCanvas("c1","c1",900,600);
    c1->Clear();
    c1->cd();
    gPad->SetLogx();

    gre1->SetTitle("Hadronic Interaction Cross Section; Kinetic Energy (GeV); Cross Section (barn)");
    gre1->SetLineColor(kRed);
    gre2->SetLineColor(kBlue);
    gre3->SetLineColor(kGreen+3);
    gre4->SetLineColor(kGreen+3);
    gre4->SetLineStyle(2);
    gre5->SetLineColor(kMagenta);
    gre1->GetYaxis()->SetRangeUser(0.5,1.7);

    gre1->Draw("AL");
    gre2->Draw("LSAME");
    gre3->Draw("LSAME");
    gre4->Draw("LSAME");
    gre5->Draw("LSAME");

    auto legend1 = new TLegend(0.60, 0.52, 0.88, 0.72);
    legend1->SetNColumns(2);
    legend1->AddEntry(gre1, "Proton", "l");
    legend1->AddEntry(gre2, "Deuteron", "l");
    legend1->AddEntry(gre3, "Helium4", "l");
    legend1->AddEntry(gre4, "Helium3", "l");
    legend1->AddEntry(gre5, "Carbon", "l");
    legend1->Draw();

    auto c2 = new TCanvas("c2","c2",900,600);
    c2->Clear();
    c2->cd();
    gPad->SetLogx();

    gre6->SetTitle("Hadronic Interaction Cross Section; Kinetic Energy per nucleon (GeV/n); Cross Section (barn)");
    gre6->SetLineColor(kRed);
    gre7->SetLineColor(kBlue);
    gre8->SetLineColor(kGreen+3);
    gre9->SetLineColor(kGreen+3);
    gre9->SetLineStyle(2);
    gre0->SetLineColor(kMagenta);
    gre6->GetYaxis()->SetRangeUser(0.5,1.7);

    

    gre6->Draw("AL");
    gre7->Draw("LSAME");
    gre8->Draw("LSAME");
    gre9->Draw("LSAME");
    gre0->Draw("LSAME");

    auto legend2 = new TLegend(0.60, 0.52, 0.88, 0.72);
    legend2->SetNColumns(2);
    legend2->AddEntry(gre6, "Proton", "l");
    legend2->AddEntry(gre7, "Deuteron", "l");
    legend2->AddEntry(gre8, "Helium4", "l");
    legend2->AddEntry(gre9, "Helium3", "l");
    legend2->AddEntry(gre0, "Carbon", "l");
    legend2->Draw();

}