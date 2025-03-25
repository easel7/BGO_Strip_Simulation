void CrossSection_Electron()
{
    auto file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Gamma_100GeV.root");
    auto tree = (TTree*)file->Get("B4");
    auto c1 = new TCanvas("c1","c1",1800,600);
    auto *h0 = new TH1D("h0","h0",30,0,30);
    auto *h1 = new TH1D("h1","h1",30,0,30);
    auto *h2 = new TH1D("h2","h2",30,0,30);

    TCut EM  = "First_Type==0";
    TCut HD  = "First_Had_Depth>=0";
    TCut HI  = "First_Had_Depth>=0 && First_Had_Type==1";

    tree->Draw("First_Depth>>h0",EM,"");
    tree->Draw("First_Had_Depth>>h1",HD,"");
    tree->Draw("First_Had_Depth>>h2",HI,"");

    h0->SetTitle("100 GeV Gamma - Electromagnetic Depth Distribution;Depth(mm);Counts");
    h1->SetTitle("100 GeV Gamma - Hadronic Depth Distribution;Depth(mm);Counts");
    h2->SetTitle("100 GeV Gamma - Inelastic Hadronic Depth Distribution;Depth(mm);Counts");

    TH1D *hC0 = (TH1D*)h0->Clone(); 
    TH1D *hC1 = (TH1D*)h1->Clone(); 
    TH1D *hC2 = (TH1D*)h2->Clone(); 

    hC0->SetTitle("100 GeV Gamma - Electromagnetic Depth Distribution;Depth(mm);N_{Survive}");
    hC1->SetTitle("100 GeV Gamma - Hadronic Depth Distribution;Depth(mm);N_{Survive}");
    hC2->SetTitle("100 GeV Gamma - Inelastic Hadronic Depth Distribution;Depth(mm);N_{Survive}");

    cout << "total h0 = " << h0->Integral() << endl;
    for(int ii=1 ;ii<=30 ; ii++)
    {
        hC0->SetBinContent(ii,1e4 - h0->Integral(0,ii));
        hC1->SetBinContent(ii,1e4 - h1->Integral(0,ii));
        hC2->SetBinContent(ii,1e4 - h2->Integral(0,ii));

        cout << "h0 = " << h0->GetBinContent(ii) << " , hC0 = " << hC0->GetBinContent(ii) << endl;
    }

    TF1 *fitFunc0 = new TF1("fitFunc0", "[0]*exp(-x/[1])", 0, 10);
    TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-x/[1])", 0, 30);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 0, 30);

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top

    c1->Clear();
    c1->Divide(3,1);
    c1->cd(1);
    gPad->SetLogy();
    gStyle->SetOptFit(1111);

    h0->Draw();
    fitFunc0->SetParameters(100, 10); // 初始参数：振幅=100, λ=10 cm
    h0->Fit(fitFunc0, "R"); // 进行拟合
    double constant0   = fitFunc0->GetParameter(0);
    double lambda0     = fitFunc0->GetParameter(1);
    double lambda0_err = fitFunc0->GetParError(1);
    latex.DrawLatex(5,100,"Fitting Function: P=[0]*exp(-x/[1])");
    
    c1->cd(2);
    gPad->SetLogy();
    gStyle->SetOptFit(1111);

    hC1->Draw();
    fitFunc1->SetParameters(100, 10); // 初始参数：振幅=100, λ=10 cm
    hC1->Fit(fitFunc1, "R"); // 进行拟合
    double constant1   = fitFunc1->GetParameter(0);
    double lambda1     = fitFunc1->GetParameter(1);
    double lambda1_err = fitFunc1->GetParError(1);
    latex.DrawLatex(5,1000,"Fitting Function: P=[0]*exp(-x/[1])");

    c1->cd(3);
    gPad->SetLogy();
    gStyle->SetOptFit(1111);

    hC2->Draw();
    fitFunc2->SetParameters(100, 10); // 初始参数：振幅=100, λ=10 cm
    hC2->Fit(fitFunc2, "R"); // 进行拟合
    double constant2   = fitFunc2->GetParameter(0);
    double lambda2     = fitFunc2->GetParameter(1);
    double lambda2_err = fitFunc2->GetParError(1);
    latex.DrawLatex(5,1000,"Fitting Function: P=[0]*exp(-x/[1])");

    double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    double em_section = 1 / lambda0 / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double em_section_err = em_section * lambda0_err/lambda0; // barn
    double hd_section = 1 / lambda1 / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hd_section_err = hd_section * lambda1_err/lambda1; // barn
    double hi_section = 1 / lambda2 / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    double hi_section_err = hi_section * lambda2_err/lambda2; // barn

    cout << "EM Interaction Constant: " << constant0  << ",EM Lambda: " << lambda0 << " ± " << lambda0_err << " mm " << endl;
    cout << "Nuclear Collsion Constant: " << constant1  << ",HD Lambda: " << lambda1 << " ± " << lambda1_err << " mm " << endl;
    cout << "Nuclear Interaction Constant: " << constant2  << ",HI Lambda: " << lambda2 << " ± " << lambda2_err << " mm " << endl;

    cout << "EM Section : " << em_section << " ± " << em_section_err << " barn" << endl;
    cout << "Nuclear Collsion Section : " << hd_section << " ± " << hd_section_err << " barn" << endl;
    cout << "Nuclear Interaction Section : " << hi_section << " ± " << hi_section_err << " barn" << endl;

    c1->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/CrossS/CrossSection_Gamma.pdf");
}