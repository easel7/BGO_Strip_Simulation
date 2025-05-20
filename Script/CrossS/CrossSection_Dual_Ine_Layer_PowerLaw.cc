void CrossSection_Dual_Ine_Layer_PowerLaw()
{

    TCut EN  = "Total_E > 64 && Total_E < 1e2";
    TCut HI  = "First_Ine_Depth>=0 && Total_E > 64 && Total_E < 1e2";

    auto file1 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_PowerLaw.root");
    auto tree1 = (TTree*)file1->Get("B4");
    auto h1_0 = new TH1D("h1_0","h1_0",14,0,14);
    auto h1 = new TH1D("h1","h1",14,0,14);
    auto hC1 = new TH1D("hC1","hC1",14,0,14);
    auto h1_energy = new TH1D("h1_energy","h1_energy",20,0,4);

    tree1->Draw("First_Had_Layer>>h1_0",EN,"");
    tree1->Draw("First_Had_Layer>>h1",HI,"");
    tree1->Draw("log10(Energy)>>h1_energy",HI,"");

    auto file2 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_PowerLaw.root");
    auto tree2 = (TTree*)file2->Get("B4");
    auto h2_0 = new TH1D("h2_0","h2_0",14,0,14);
    auto h2 = new TH1D("h2","h2",14,0,14);
    auto hC2 = new TH1D("hC2","hC2",14,0,14);
    auto h2_energy = new TH1D("h2_energy","h2_energy",20,0,4);

    tree2->Draw("First_Had_Layer>>h2_0",EN,"");
    tree2->Draw("First_Had_Layer>>h2",HI,"");
    tree2->Draw("log10(Energy)>>h2_energy",HI,"");

    for(int ii = 1 ; ii <= 14 ; ii++)
    {
        hC1->SetBinContent(ii, ( 1e4 -h1->Integral(0,ii) ) * 0.5 );
        hC2->SetBinContent(ii, ( 1e4 -h2->Integral(0,ii) ) * 0.5 );
        // cout << h2->GetBinContent(ii) << endl;
    }

    TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-x/[1])", 0, 3); fitFunc1->SetParameters(100, 10); fitFunc1->SetLineColor(kRed);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 0, 3); fitFunc2->SetParameters(100, 10); fitFunc2->SetLineColor(kBlue);
    TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*exp(-x/[1])", 0, 3); fitFunc3->SetParameters(100, 10); fitFunc3->SetLineColor(kBlack);

    // Fit for Interaction 
    TF1 *fitFunc4 = new TF1("fitFunc4", "[0]/8.75*exp(-x/8.75)+[1]/6.73*exp(-x/6.73)", 0, 10); fitFunc4->SetParameters(1e4, 1e4); fitFunc4->SetLineColor(kBlack);

    // Fit for Survive
    TF1 *fitFunc5 = new TF1("fitFunc5", "[0]*(1-[1])*exp(-x/8.75)+ [0]*[1]*exp(-x/6.73)", 0, 10); fitFunc5->SetParameters(1e4, 0.5); fitFunc5->SetLineColor(kBlack);

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top

    auto c0 = new TCanvas("c0","c0",1200,1200);

    h1_0->SetLineColor(kRed);   // Proton in red
    h2_0->SetLineColor(kBlue);  // Deuteron in blue
    h1_0->Scale(0.5);
    h2_0->Scale(0.5);
    auto h3_0 = (TH1D*)h2_0->Clone();
    h3_0->Add(h1_0);
    h3_0->SetLineColor(kBlack);  // Deuteron in blue
    h3_0->SetLineWidth(2);
    h3_0->SetTitle("N_{interaction};BGO Layer;Counts");
    h3_0->Draw("hist");

    h1_0->Draw("histsame");
    h2_0->Draw("histsame");
    h3_0->Fit(fitFunc4, "R"); // 进行拟合

    fitFunc4->Draw("same");
    double proton_num   = fitFunc4->GetParameter(0);
    double proton_num_err = fitFunc4->GetParError(0);
    double deuteron_num = fitFunc4->GetParameter(1);
    double deuteron_num_err = fitFunc4->GetParError(1);
    latex.DrawLatex(0,3500,"Fitting Function: N_{int} = N_{p}/#lambda_{p} *exp(-x/#lambda_{p}) + N_{d}/#lambda_{d} *exp(-x/#lambda_{d})");
    latex.DrawLatex(5,3000,Form("Fitting Proton Number: %.1f",proton_num));
    latex.DrawLatex(5,2500,Form("Fitting Deuteron Number: %.1f",deuteron_num));
    latex.DrawLatex(5,2000,Form("N_{deuteron}/(N_{proton}+N_{deuteron}): %.2f #pm %.2f",deuteron_num / (proton_num + deuteron_num), deuteron_num_err / sqrt(pow(proton_num_err,2)+pow(deuteron_num_err,2))));

    auto c1 = new TCanvas("c1","c1",1200,1200);
    c1->Clear();
    gPad->SetLogy();
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    // Set styles and colors for clarity
    hC1->SetTitle("64 - 100 GeV Deposit Energy, Inelastic Hadronic Layer Distribution;Layer;N_{survive}");
    hC1->GetYaxis()->SetRangeUser(1e1,5e5);
    hC1->SetLineColor(kRed);   // Proton in red
    hC1->SetLineWidth(2);
    hC2->SetLineColor(kBlue);  // Deuteron in blue
    hC2->SetLineWidth(2);

    hC1->Draw("hist");
    hC2->Draw("histsame");
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
    double hd_section = 1 / (lambda1*25) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hd_section_err = hd_section * lambda1_err/lambda1; // barn
    double hi_section = 1 / (lambda2*25) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hi_section_err = hi_section * lambda2_err/lambda2; // barn
    double hm_section = 1 / (lambda3*25) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hm_section_err = hm_section * lambda3_err/lambda3; // barn

    latex.DrawLatex(0,pow(10,3.3),"Fitting Function: N_{leave} =N_{total} *exp(-x/#lambda)");
    latex.DrawLatex(0,pow(10,3.0),Form("Proton Fitting #lambda: %.2f mm",lambda1*25));
    latex.DrawLatex(0,pow(10,2.7),Form("Deuteron Fitting #lambda: %.2f mm",lambda2*25));
    // latex.DrawLatex(5,pow(10,2.4),Form("Mixed Fitting #lambda: %.2f mm",lambda3));

    cout << "Proton Constant: " << constant1  << ",HD Lambda: " << lambda1 << " ± " << lambda1_err << " mm " << endl;
    cout << "Deuteron Constant: " << constant2  << ",HI Lambda: " << lambda2 << " ± " << lambda2_err << " mm " << endl;
    cout << "Proton Nuclear Interaction Section : " << hd_section << " ± " << hd_section_err << " barn" << endl;
    cout << "Deuteron Nuclear Interaction Section : " << hi_section << " ± " << hi_section_err << " barn" << endl;

    // c1->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/CrossS/CrossSection_Dual_1000GeV.pdf");

    auto c2 = new TCanvas("c2","c2",1200,1200);
    c2->Clear();
    gPad->SetLogy();
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    h1_energy->SetTitle("64 - 100 GeV Deposit Energy,log10(Incident Energy);N_{count}");
    h1_energy->GetYaxis()->SetRangeUser(1e1,5e5);
    h1_energy->SetLineColor(kRed);   // Proton in red
    h1_energy->SetLineWidth(2);
    h2_energy->SetLineColor(kBlue);  // Deuteron in blue
    h2_energy->SetLineWidth(2);
    h1_energy->Draw("hist");
    h2_energy->Draw("histsame");
}