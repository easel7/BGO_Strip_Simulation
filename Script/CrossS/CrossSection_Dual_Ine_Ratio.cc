void CrossSection_Dual_Ine_Ratio()
{
    double Ratio[18] = {0};

    TCut HI  = "First_Ine_Depth>=0 && Nhits >= 10";
    auto file1 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
    auto tree1 = (TTree*)file1->Get("B4");
    auto file2 = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_1000GeV.root");
    auto tree2 = (TTree*)file2->Get("B4");

    // Fit for deuteron ratio and error;
    // Fit from N_int
    // Fit from N_sur
    auto gre_int = new TGraphErrors();  
    auto gre_sur = new TGraphErrors();

    for (int i =0; i < 18; i++)
    {
        if (i < 10)      Ratio[i] = (i + 1) * 0.01;     
        else             Ratio[i] = (i - 9 + 1) * 0.1;  

        cout << Ratio[i] << " , " << 1-Ratio[i] <<  endl;
        file1->cd();
        if (gDirectory->FindObject("h1_0")) delete gDirectory->FindObject("h1_0");
        if (gDirectory->FindObject("h1"  )) delete gDirectory->FindObject("h1"  );
        if (gDirectory->FindObject("hC1" )) delete gDirectory->FindObject("hC1" );
        auto h1_0 = new TH1D("h1_0","h1_0",16,-25,375);
        auto h1   = new TH1D("h1","h1",14,0,357);
        auto hC1  = new TH1D("hC1","hC1",14,0,357);
        tree1->Draw("First_Ine_Depth>>h1_0","","");
        tree1->Draw("First_Ine_Depth>>h1",HI,"");
        
        file2->cd();
        if (gDirectory->FindObject("h2_0")) delete gDirectory->FindObject("h2_0");
        if (gDirectory->FindObject("h2"  )) delete gDirectory->FindObject("h2"  );
        if (gDirectory->FindObject("hC2" )) delete gDirectory->FindObject("hC2" );
        auto h2_0 = new TH1D("h2_0","h2_0",16,-25,375);
        auto h2   = new TH1D("h2","h2",14,0,357);
        auto hC2  = new TH1D("hC2","hC2",14,0,357);
        tree2->Draw("First_Ine_Depth>>h2_0","","");
        tree2->Draw("First_Ine_Depth>>h2",HI,"");
        h1->Scale(1-Ratio[i]);  h1_0->Scale(1-Ratio[i]);
        h2->Scale(Ratio[i]);    h2_0->Scale(Ratio[i]);
        for(int k = 1 ; k <= h1->GetNbinsX() ; k++)
        {
            hC1->SetBinContent(k, ( h1_0->Integral() - h1->Integral(1,k) ) );
            hC2->SetBinContent(k, ( h2_0->Integral() - h2->Integral(1,k) ) );
            // cout << h1->Integral(1,k) << " , " <<  h1->GetBinContent(k) << endl;
        }
        auto h3 = (TH1D*)h2->Clone();    h3->Add(h1);
        auto hC3 = (TH1D*)hC2->Clone();  hC3->Add(hC1);

        // Fit for Survive
        TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*(1-[1])*exp(-x/208)+ [0]*[1]*exp(-x/173)", 0, 300); fitFunc3->SetParameters(1e4, Ratio[i]); fitFunc3->SetLineColor(kBlack); // fitFunc3->SetParLimits(1, 1e-6 ,1);
        TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-x/[1])", 0, 300); fitFunc1->SetParameters(100, 10); fitFunc1->SetLineColor(kRed);
        TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-x/[1])", 0, 300); fitFunc2->SetParameters(100, 10); fitFunc2->SetLineColor(kBlue);

        // Fit for Interaction 
        TF1 *fitFunc4 = new TF1("fitFunc4", "[0]*[2] *(1-[1])/208*exp(-x/208)+[0]*25 *[1]/173*exp(-x/173)", 0, 300); fitFunc4->SetParameters(1e4, Ratio[i]); fitFunc4->SetLineColor(kBlack); //fitFunc4->SetParLimits(1, 1e-6 ,1);
        TF1 *fitFunc5 = new TF1("fitFunc5", "[0]*[2] /[1]*exp(-x/[1])", 0, 300); fitFunc5->SetParameters(1e3,8); fitFunc5->SetLineColor(kRed);
        TF1 *fitFunc6 = new TF1("fitFunc6", "[0]*[2] /[1]*exp(-x/[1])", 0, 300); fitFunc6->SetParameters(1e2,8); fitFunc6->SetLineColor(kBlue);
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetTextFont(72);
        latex.SetTextAlign(13);  //align at top

        auto c0 = new TCanvas("c0","c0",1800,900);
        c0->Divide(2,1);
        c0->cd(1);
        h1->SetTitle("N_{interaction};Depth(mm);Counts");
        h1->SetLineColor(kRed);   // Proton in red
        h2->SetLineColor(kBlue);  // Deuteron in blue


        h3->SetLineColor(kBlack);  // Deuteron in blue
        h3->SetLineWidth(2);
        h1->GetYaxis()->SetRangeUser(0,1.2*h3->GetMaximum());
        h1->Draw("hist");
        h2->Draw("histsame");
        h3->Draw("histsame");
        h3->Fit(fitFunc4, "RSQ"); // 进行拟合
        h1->Fit(fitFunc5,"RSQ");
        h2->Fit(fitFunc6,"RSQ");
        fitFunc4->Draw("same");
        fitFunc5->Draw("same");
        fitFunc6->Draw("same");
        double total_num1   = fitFunc4->GetParameter(0);
        double total_num1_err = fitFunc4->GetParError(0);
        double deuteron_ratio1 = fitFunc4->GetParameter(1);
        double deuteron_ratio1_err = fitFunc4->GetParError(1);

        gre_int->SetPoint(i,Ratio[i],deuteron_ratio1);
        gre_int->SetPointError(i,0,deuteron_ratio1_err);

        c0->cd(2);
        // gPad->SetLogy();
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

        // Set styles and colors for clarity
        hC1->GetYaxis()->SetRangeUser(1e1,2e4);
        hC1->SetLineColor(kRed);   // Proton in red
        hC1->SetLineWidth(2);
        hC2->SetLineColor(kBlue);  // Deuteron in blue
        hC2->SetLineWidth(2);


        hC3->SetTitle("1000 GeV Inelastic Hadronic Depth Distribution;Depth(mm);N_{survive}");
        hC3->Draw("hist");
        hC1->Draw("hist");
        hC2->Draw("histsame");
        hC1->Fit(fitFunc1, "RSQ"); // 进行拟合
        hC2->Fit(fitFunc2, "RSQ"); // 进行拟合
        double constant1   = fitFunc1->GetParameter(0);
        double lambda1     = fitFunc1->GetParameter(1);
        double lambda1_err = fitFunc1->GetParError(1);
        double constant2   = fitFunc2->GetParameter(0);
        double lambda2     = fitFunc2->GetParameter(1);
        double lambda2_err = fitFunc2->GetParError(1);


        hC3->SetLineColor(kBlack);  // Deuteron in blue
        hC3->SetLineWidth(2);
        hC3->Draw("histsame");
        hC3->Fit(fitFunc3, "RSQ"); // 进行拟合
        fitFunc3->Draw("same");
        double total_num   = fitFunc3->GetParameter(0);
        double total_num_err = fitFunc3->GetParError(0);
        double deuteron_ratio = fitFunc3->GetParameter(1);
        double deuteron_ratio_err = fitFunc3->GetParError(1);

        double n_BGO = TMath::Na()*7.13/ (1245.8344/19.); // cm-3
        double hd_section = 1 / (lambda1) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
        double hd_section_err = hd_section * lambda1_err/lambda1; // barn
        double hi_section = 1 / (lambda2) / n_BGO * 1e25; // barn, mm = 1e-1 cm, 1e24 barn = 1 cm^2
        double hi_section_err = hi_section * lambda2_err/lambda2; // barn

        gre_sur->SetPoint(i,Ratio[i],deuteron_ratio);
        gre_sur->SetPointError(i,0,deuteron_ratio_err);

        c0->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/CrossS/Ratio_%.2f_1000GeV.pdf",Ratio[i]));
    }

    auto c2 = new TCanvas("c2","c2",1900,900);
    gre_int->SetLineColor(kRed);
    gre_sur->SetLineColor(kBlue);
    gre_int->SetLineWidth(2);
    gre_sur->SetLineWidth(2);
    gre_int->SetMarkerColor(kRed);
    gre_sur->SetMarkerColor(kBlue);
    gre_int->SetMarkerStyle(20);
    gre_sur->SetMarkerStyle(21);

    auto line_ref = new TLine(5e-3,5e-3,2,2);
    line_ref->SetLineColor(kBlack);
    line_ref->SetLineWidth(2);
    line_ref->SetLineStyle(2);

    auto line_vet = new TLine(0.5,5e-3,0.5,2);
    line_vet->SetLineColor(kBlack);
    line_vet->SetLineWidth(2);
    line_vet->SetLineStyle(2);

    c2->cd();
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gre_int->GetXaxis()->SetLimits(5e-3,2);
    gre_int->GetYaxis()->SetRangeUser(5e-3,2);
    gre_int->SetTitle(";True Deuteron Ratio;Fitted Deuteron Ratio");
    gre_int->Draw("AP");
    gre_sur->Draw("PSAME");
    line_ref->Draw("same");
    line_vet->Draw("same");

    auto lg2 = new TLegend(0.68,0.12,0.88,0.22);
    lg2->AddEntry(gre_int,"Fited Ratio From N_{int}","p");
    lg2->AddEntry(gre_sur,"Fited Ratio From N_{sur}","p");
    lg2->AddEntry(line_ref,"y = x","l");
    lg2->AddEntry(line_vet,"x = 0.5","l");
    lg2->Draw();
    c2->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/CrossS/FittedRatio.pdf");
}