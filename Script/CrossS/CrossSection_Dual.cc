void CrossSection_Dual()
{

    TCut EM  = "First_Type==0";
    TCut HD  = "First_Had_Depth>=0";
    TCut HI  = "First_Had_Depth>=0 && First_Had_Type==1";

    auto file1 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
    auto tree1 = (TTree*)file1->Get("B4");
    auto h1 = new TH1D("h1","h1",35,0,350);
    auto hC1 = new TH1D("hC1","hC1",35,0,350);

    tree1->Draw("First_Had_Depth>>h1",HI,"");
    hC1->SetTitle("1000 GeV Inelastic Hadronic Depth Distribution;Depth(mm);N_{survive}");

    auto file2 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_1000GeV.root");
    auto tree2 = (TTree*)file2->Get("B4");
    auto h2 = new TH1D("h2","h2",35,0,350);
    auto hC2 = new TH1D("hC2","hC2",35,0,350);
    tree2->Draw("First_Had_Depth>>h2",HI,"");
    for(int ii = 1 ; ii <= 35 ; ii++)
    {
        hC1->SetBinContent(ii, (1e4 -h1->Integral(0,ii) )*(0.5) );
        hC2->SetBinContent(ii, ( 1e4 -h2->Integral(0,ii) )*(0.5) );
    }

    TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-x/[1])", 0, 40); fitFunc1->SetParameters(100, 10); fitFunc1->SetLineColor(kRed);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 0, 40); fitFunc2->SetParameters(100, 10); fitFunc2->SetLineColor(kBlue);
    TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*exp(-x/[1])", 0, 40); fitFunc3->SetParameters(100, 10); fitFunc3->SetLineColor(kBlack);


    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top

    auto c1 = new TCanvas("c1","c1",1200,1200);

    c1->Clear();
    gPad->SetLogy();
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    // Set styles and colors for clarity
    hC1->GetYaxis()->SetRangeUser(1e1,2e4);
    hC1->SetLineColor(kRed);   // Proton in red
    hC1->SetLineWidth(2);
    hC2->SetLineColor(kBlue);  // Deuteron in blue
    hC2->SetLineWidth(2);

    hC1->Draw();
    hC2->Draw("same");
    hC1->Fit(fitFunc1, "R"); // 进行拟合
    double constant1   = fitFunc1->GetParameter(0);
    double lambda1     = fitFunc1->GetParameter(1);
    double lambda1_err = fitFunc1->GetParError(1);

    hC2->Fit(fitFunc2, "R"); // 进行拟合
    double constant2   = fitFunc2->GetParameter(0);
    double lambda2     = fitFunc2->GetParameter(1);
    double lambda2_err = fitFunc2->GetParError(1);

    auto hC3 = (TH1D*)hC2->Clone();
    hC3->Add(hC1);
    hC3->SetLineColor(kBlack);  // Deuteron in blue
    hC3->SetLineWidth(2);
    hC3->Draw("histsame");
    hC3->Fit(fitFunc3, "R"); // 进行拟合
    double constant3   = fitFunc3->GetParameter(0);
    double lambda3     = fitFunc3->GetParameter(1);
    double lambda3_err = fitFunc3->GetParError(1);

    double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    double hd_section = 1 / lambda1 / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hd_section_err = hd_section * lambda1_err/lambda1; // barn
    double hi_section = 1 / lambda2 / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hi_section_err = hi_section * lambda2_err/lambda2; // barn
    double hm_section = 1 / lambda3 / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hm_section_err = hm_section * lambda3_err/lambda3; // barn

    latex.DrawLatex(5,pow(10,3.3),"Fitting Function: N_{leave} =N_{total} *exp(-x/#lambda)");
    latex.DrawLatex(5,pow(10,3.0),Form("Proton Fitting #lambda: %.2f mm",lambda1));
    latex.DrawLatex(5,pow(10,2.7),Form("Deuteron Fitting #lambda: %.2f mm",lambda2));
    // latex.DrawLatex(5,pow(10,2.4),Form("Mixed Fitting #lambda: %.2f mm",lambda3));

    cout << "Proton Constant: " << constant1  << ",HD Lambda: " << lambda1 << " ± " << lambda1_err << " mm " << endl;
    cout << "Deuteron Constant: " << constant2  << ",HI Lambda: " << lambda2 << " ± " << lambda2_err << " mm " << endl;

    cout << "Proton Nuclear Interaction Section : " << hd_section << " ± " << hd_section_err << " barn" << endl;
    cout << "Deuteron Nuclear Interaction Section : " << hi_section << " ± " << hi_section_err << " barn" << endl;

    c1->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/CrossS/CrossSection_Dual_1000GeV.pdf");
}