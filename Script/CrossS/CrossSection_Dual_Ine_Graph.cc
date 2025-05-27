void CrossSection_Dual_Ine_Graph()
{
    TCut HI  = "First_Ine_Depth>=0";
    auto file1 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
    auto tree1 = (TTree*)file1->Get("B4");
    auto h1_0 = new TH1D("h1_0","h1_0",16,-25,375);
    auto h1   = new TH1D("h1"  ,"h1"  ,16,-25,375);
    auto hC1  = new TH1D("hC1" ,"hC1" ,16,-25,375);
    tree1->Draw("First_Ine_Depth>>h1_0","","");
    tree1->Draw("First_Ine_Depth>>h1","","");

    auto file2 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_1000GeV.root");
    auto tree2 = (TTree*)file2->Get("B4");
    auto h2_0 = new TH1D("h2_0","h2_0",16,-25,375);
    auto h2   = new TH1D("h2"  ,"h2"  ,16,-25,375);
    auto hC2  = new TH1D("hC2" ,"hC2" ,16,-25,375);
    tree2->Draw("First_Ine_Depth>>h2_0","","");
    tree2->Draw("First_Ine_Depth>>h2","","");

    auto p_int = new TGraphErrors();
    auto d_int = new TGraphErrors();
    auto t_int = new TGraphErrors();

    auto p_sur = new TGraphErrors();
    auto d_sur = new TGraphErrors();
    auto t_sur = new TGraphErrors();

    auto h3 = (TH1D*)h2->Clone();    h3->Add(h1);
    h1->Sumw2();
    h2->Sumw2();
    h3->Sumw2();

    cout << h1_0->Integral() << " , " <<  h2_0->Integral()  << endl;
    for(int ii = 0 ; ii <= h1->GetNbinsX() ; ii++)
    {
        hC1->SetBinContent(ii+1, ( h1_0->Integral() - h1->Integral(0,ii) ) );  // N survive
        hC2->SetBinContent(ii+1, ( h2_0->Integral() - h2->Integral(0,ii) ) );  // N survive
        cout << ii << " , " << hC1->GetBinContent(ii+1) << " , " << h1->Integral(0,ii) << " , " << h1->GetBinContent(ii+1) <<  endl;
        p_int->SetPoint(ii,25.5*(ii+1), h1->GetBinContent(ii+1) );         p_int->SetPointError(ii,0, h1->GetBinError(ii+1) );
        d_int->SetPoint(ii,25.5*(ii+1), h2->GetBinContent(ii+1) );         d_int->SetPointError(ii,0, h2->GetBinError(ii+1) );
        t_int->SetPoint(ii,25.5*(ii+1), h3->GetBinContent(ii+1) );         t_int->SetPointError(ii,0, h3->GetBinError(ii+1) );
    }

    auto hC3 = (TH1D*)hC2->Clone();    hC3->Add(hC1);
    hC1->Sumw2();
    hC2->Sumw2();
    hC3->Sumw2();

    for(int ii = 0 ; ii <= h1->GetNbinsX() ; ii++)
    {
        p_sur->SetPoint(ii,25.5*ii    , hC1->GetBinContent(ii+1));       p_sur->SetPointError(ii,0    , hC1->GetBinError(ii+1));
        d_sur->SetPoint(ii,25.5*ii    , hC2->GetBinContent(ii+1));       d_sur->SetPointError(ii,0    , hC2->GetBinError(ii+1));
        t_sur->SetPoint(ii,25.5*ii    , hC3->GetBinContent(ii+1));       t_sur->SetPointError(ii,0    , hC3->GetBinError(ii+1));
        cout << hC1->GetBinContent(ii+1) << " , " << h1->Integral(0,ii) << " , " << h1->GetBinContent(ii+1) <<  endl;
    }

    // Fit for Survive
    TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*(1-[1])*exp(-x/208)+ [0]*[1]*exp(-x/173)", 50,300); fitFunc3->SetParameters(1e4, 0.5); fitFunc3->SetLineColor(kBlack);
    TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-x/[1])", 50,300); fitFunc1->SetParameters(100, 10); fitFunc1->SetLineColor(kRed);  fitFunc1->FixParameter(0,1e4);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 50,300); fitFunc2->SetParameters(100, 10); fitFunc2->SetLineColor(kBlue); fitFunc2->FixParameter(0,1e4);

    // Fit for Interaction 
    TF1 *fitFunc4 = new TF1("fitFunc4", "[0]*[2] *(1-[1])/208*exp(-x/208)+[0]*25 *[1]/173*exp(-x/173)", 0,300); fitFunc4->SetParameters(1e4, 0.5); fitFunc4->SetLineColor(kBlack);
    TF1 *fitFunc5 = new TF1("fitFunc5", "[0]*[2] /[1]*exp(-x/[1])", 0,300); fitFunc5->SetParameters(1e3,8); fitFunc5->SetLineColor(kRed);  fitFunc5->FixParameter(0,1e4);
    TF1 *fitFunc6 = new TF1("fitFunc6", "[0]*[2] /[1]*exp(-x/[1])", 0,300); fitFunc6->SetParameters(1e2,8); fitFunc6->SetLineColor(kBlue); fitFunc6->FixParameter(0,1e4);
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(72);
    latex.SetTextAlign(13);  //align at top

    // auto c0 = new TCanvas("c0","c0",1200,900);
    // p_int->SetLineColor(kRed);       p_int->SetMarkerColor(kRed);   
    // d_int->SetLineColor(kBlue);      d_int->SetMarkerColor(kBlue);  
    // t_int->SetLineColor(kBlack);     t_int->SetMarkerColor(kBlack); 
    // t_int->SetLineWidth(2);
    // t_int->SetTitle("N_{interaction};Depth(mm);Counts");
    // t_sur->GetXaxis()->SetLimits(-10,360);
    // t_int->GetYaxis()->SetRangeUser(0,1.2*h3->GetMaximum());
    // t_int->Draw("AP");
    // p_int->Draw("Psame");
    // d_int->Draw("Psame");

    // t_int->Fit(fitFunc4,"R"); // 进行拟合
    // p_int->Fit(fitFunc5,"R");
    // d_int->Fit(fitFunc6,"R");
    // fitFunc4->Draw("same");
    // fitFunc5->Draw("same");
    // fitFunc6->Draw("same");
    // double total_num1   = fitFunc4->GetParameter(0);
    // double total_num1_err = fitFunc4->GetParError(0);
    // double deuteron_ratio1 = fitFunc4->GetParameter(1);
    // double deuteron_ratio1_err = fitFunc4->GetParError(1);
    // latex.DrawLatex(0,500,"Fitting Function: N_{int} = N_{tot}*(1-^{}r^{}_{d})/#lambda_{p} *exp(-x/#lambda_{p}) + N_{tot}*r_{d}/#lambda_{d} *exp(-x/#lambda_{d})");
    // latex.DrawLatex(0,400,Form("Fitting Total Number: %.2f#pm %.2f",total_num1 , total_num1_err ));
    // latex.DrawLatex(50,300,Form("Fitting Deuteron Ratio: %.2f#pm %.2f",deuteron_ratio1 , deuteron_ratio1_err ));
    // latex.DrawLatex(0,208,Form("Proton Fitting Num: %.2f ",fitFunc5->GetParameter(0)));
    // latex.DrawLatex(0,100,Form("Deuteron Fitting Num: %.2f ",fitFunc6->GetParameter(0)));

    auto c1 = new TCanvas("c1","c1",1200,900);
    c1->Clear();
    // gPad->SetLogy();
    // gStyle->SetOptFit(0);
    // gStyle->SetOptStat(0);
    t_sur->SetTitle("1000 GeV Inelastic Hadronic Depth Distribution;Depth(mm);N_{survive}");
    t_sur->GetYaxis()->SetRangeUser(1e1,3.6e4);
    t_sur->GetXaxis()->SetLimits(-10,360);

    p_sur->SetLineColor(kRed);       p_sur->SetMarkerColor(kRed);   
    d_sur->SetLineColor(kBlue);      d_sur->SetMarkerColor(kBlue);  
    t_sur->SetLineColor(kBlack);     t_sur->SetMarkerColor(kBlack); 
    t_sur->Draw("AP");
    p_sur->Draw("Psame");
    d_sur->Draw("Psame");
    p_sur->Fit(fitFunc1, "R"); // 进行拟合
    d_sur->Fit(fitFunc2, "R"); // 进行拟合
    t_sur->Fit(fitFunc3, "R"); // 进行拟合
    fitFunc1->Draw("same");
    fitFunc2->Draw("same");
    fitFunc3->Draw("same");

    // double constant1   = fitFunc1->GetParameter(0);
    // double lambda1     = fitFunc1->GetParameter(1);
    // double lambda1_err = fitFunc1->GetParError(1);
    // double constant2   = fitFunc2->GetParameter(0);
    // double lambda2     = fitFunc2->GetParameter(1);
    // double lambda2_err = fitFunc2->GetParError(1);
    // double total_num   = fitFunc3->GetParameter(0);
    // double total_num_err = fitFunc3->GetParError(0);
    // double deuteron_ratio = fitFunc3->GetParameter(1);
    // double deuteron_ratio_err = fitFunc3->GetParError(1);

    // double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
    // double hd_section = 1 / (lambda1) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    // double hd_section_err = hd_section * lambda1_err/lambda1; // barn
    // double hi_section = 1 / (lambda2) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
    // double hi_section_err = hi_section * lambda2_err/lambda2; // barn
    // latex.DrawLatex(0,pow(10,4.3),"Fitting Function: N_{sur} = N_{tot}*(1-^{}r^{}_{d})*exp(-x/#lambda_{p}) + N_{tot}*r_{d}/#lambda_{d} *exp(-x/#lambda_{d})");
    // latex.DrawLatex(0,pow(10,4.0),Form("Fitting Total Number: %.2f#pm %.2f",total_num , total_num_err ));
    // latex.DrawLatex(0,pow(10,3.7),Form("Fitting Deuteron Ratio: %.2f#pm %.2f",deuteron_ratio , deuteron_ratio_err ));
    // latex.DrawLatex(0,pow(10,3.4),Form("Proton Fitting Num: %.2f ",constant1));
    // latex.DrawLatex(0,pow(10,3.1),Form("Deuteron Fitting Num: %.2f ",constant2));

    // cout << "Proton Ratio: " << total_num  << " ± " << total_num_err << " " << endl;
    // cout << "Deuteron Ratio: " << deuteron_ratio  << " ± " << deuteron_ratio_err << " " << endl;
    // cout << "Proton Constant: " << constant1  << ",HD Lambda: " << lambda1 << " ± " << lambda1_err << " mm " << endl;
    // cout << "Deuteron Constant: " << constant2  << ",HI Lambda: " << lambda2 << " ± " << lambda2_err << " mm " << endl;
    // cout << "Proton Nuclear Interaction Section : " << hd_section << " ± " << hd_section_err << " barn" << endl;
    // cout << "Deuteron Nuclear Interaction Section : " << hi_section << " ± " << hi_section_err << " barn" << endl;

    // c1->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/CrossS/CrossSection_Dual_1000GeV.pdf");
}