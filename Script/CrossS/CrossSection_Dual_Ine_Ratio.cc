void CrossSection_Dual_Ine_Ratio()
{
    int Energy_Name[28]={0};
    int Energy_Name_Err[28]={0};
    double Ratio[28] = {0};

    double Ratio_Edge[29];
    double Energy_Edge[29];

    for (int j = 0; j < 29; j++)
    {
        if (j <= 10)                 { Energy_Edge[j] = 5 + 10 * j;             Ratio_Edge[j] = 0.0005 + 0.001 * j;     }   // Center: 10~100 → edges: 5~105
        else if (j > 10 && j <= 19)  { Energy_Edge[j] = 150 + 100 * (j - 10);   Ratio_Edge[j] = 0.015 + 0.01 * (j - 10);      }   // Center: 10~100 → edges: 5~105  
        else                         { Energy_Edge[j] = 1500 + 1000 * (j - 19); Ratio_Edge[j] = 0.15 + 0.1 * (j - 19); } // Center: 2000~10000 → edges: 1500~10500     
    }

    auto hist_int = new TH2D("hist_int","hist_int",28,Energy_Edge,28,Ratio_Edge);
    auto hist_sur = new TH2D("hist_sur","hist_sur",28,Energy_Edge,28,Ratio_Edge);
    auto lamb_int = new TH2D("lamb_int","lamb_int",28,Energy_Edge,28,Ratio_Edge);
    auto lamb_sur = new TH2D("lamb_sur","lamb_sur",28,Energy_Edge,28,Ratio_Edge);

    hist_int->SetTitle("N_{int};Incident Energy (GeV);True r_{d}; Fitted r_{d}");
    hist_sur->SetTitle("N_{sur};Incident Energy (GeV);True r_{d}; Fitted r_{d}");

    lamb_int->SetTitle("N_{int};Incident Energy (GeV);True r_{d}; Fitted #lambda_{d}");
    lamb_sur->SetTitle("N_{sur};Incident Energy (GeV);True r_{d}; Fitted #lambda_{d}");

    for (int k =18; k < 19; k++) // Energy
    // for (int k =0; k < 28; k++)
    {
        if (k < 10)      {Energy_Name[k] = (k + 1) * 10;        }       // 10 ~ 100
        else if (k < 19) {Energy_Name[k] = (k - 9 + 1) * 100;   }        // 200 ~ 1000
        else             {Energy_Name[k] = (k - 18 + 1) * 1000; }        // 2000 ~ 10000

        TCut HI  = "First_Ine_Depth>=0 && Nhits >= 10";
        auto file_p  = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",Energy_Name[k]));
        auto tree_p  = (TTree*)file_p->Get("B4");
        auto file_d  = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",Energy_Name[k]));
        auto tree_d  = (TTree*)file_d->Get("B4");

        // Fit for deuteron ratio and error;
        // Fit from N_int
        // Fit from N_sur
        auto gre_int = new TGraphErrors();  
        auto gre_sur = new TGraphErrors();

        auto gre_lambda_d_int = new TGraphErrors();  
        auto gre_lambda_d_sur = new TGraphErrors();

        for (int i =5; i < 6; i++) // Ratio
        // for (int i =0; i < 27; i++)
        {
            if      (i < 10)          Ratio[i] = (i + 1) * 0.001;     
            else if (i < 19)          Ratio[i] = (i - 9 + 1) * 0.01;
            else                      Ratio[i] = (i - 18 + 1) * 0.1;  
            cout << Ratio[i] << " , " << 1-Ratio[i] <<  endl;
            if (gDirectory->FindObject("c0")) delete gDirectory->FindObject("c0");

            file_p->cd();
            if (gDirectory->FindObject("h_p_tot")) delete gDirectory->FindObject("h_p_tot");
            if (gDirectory->FindObject("h_p_int")) delete gDirectory->FindObject("h_p_int");
            if (gDirectory->FindObject("h_p_sur")) delete gDirectory->FindObject("h_p_sur");
            auto h_p_tot = new TH1D("h_p_tot","h_p_tot",15,-25.5,357);  h_p_tot->Sumw2();
            auto h_p_int = new TH1D("h_p_int","h_p_int",14,0,357);      h_p_int->Sumw2();
            auto h_p_sur = new TH1D("h_p_sur","h_p_sur",14,0,357);      
            tree_p->Draw("First_Ine_Depth>>h_p_tot","","");
            tree_p->Draw("First_Ine_Depth>>h_p_int",HI,"");
            
            file_d->cd();
            if (gDirectory->FindObject("h_d_tot")) delete gDirectory->FindObject("h_d_tot");
            if (gDirectory->FindObject("h_d_int")) delete gDirectory->FindObject("h_d_int");
            if (gDirectory->FindObject("h_d_sur")) delete gDirectory->FindObject("h_d_sur");
            auto h_d_tot = new TH1D("h_d_tot","h_d_tot",15,-25.5,357);  h_d_tot->Sumw2();
            auto h_d_int = new TH1D("h_d_int","h_d_int",14,0,357);      h_d_int->Sumw2();
            auto h_d_sur = new TH1D("h_d_sur","h_d_sur",14,0,357);     
            tree_d->Draw("First_Ine_Depth>>h_d_tot","","");
            tree_d->Draw("First_Ine_Depth>>h_d_int",HI,"");

            h_p_int->Scale( 1-Ratio[i]);  h_p_tot->Scale(1-Ratio[i]);
            h_d_int->Scale( Ratio[i]  );    h_d_tot->Scale(Ratio[i]);
            
            for(int ii = 1 ; ii <= h_p_int->GetNbinsX() ; ii++)
            {
                h_p_sur->SetBinContent(ii, ( h_p_tot->Integral() - h_p_int->Integral(1,ii) ) );
                h_d_sur->SetBinContent(ii, ( h_d_tot->Integral() - h_d_int->Integral(1,ii) ) );
            }
            auto h_2_int = (TH1D*)h_p_int->Clone(); h_2_int->Add(h_d_int);
            auto h_2_sur = (TH1D*)h_p_sur->Clone(); h_2_sur->Add(h_d_sur);
            h_2_int->Sumw2();    h_p_int->Sumw2();
            h_2_sur->Sumw2();    h_d_sur->Sumw2();

            // Fit for Survive
            TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x+12.75)/[1])", 0,300); fitFunc1->SetParameters(1e4, 200); fitFunc1->SetLineColor(kRed);  fitFunc1->FixParameter(0,h_p_tot->Integral());
            TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x+12.75)/[1])", 0,300); fitFunc2->SetParameters(1e4, 200); fitFunc2->SetLineColor(kBlue); fitFunc2->FixParameter(0,h_d_tot->Integral());
            TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*(1-[1])*exp(-(x+12.75)/[2])+ [0]*[1]*exp(-(x+12.75)/ [3] )", 0,300); fitFunc3->SetParameters(2e4,Ratio[i],200,170); fitFunc3->SetLineColor(kBlack); 
            fitFunc3->FixParameter(0,(h_p_tot->Integral()+h_d_tot->Integral()));
            fitFunc3->FixParameter(2,fitFunc1->GetParameter(1));
            // fitFunc3->FixParameter(3,fitFunc2->GetParameter(1));
            fitFunc3->SetParLimits(1,1e-3,1);
            // fitFunc3->SetParLimits(2,fitFunc1->GetParameter(1)/2,fitFunc1->GetParameter(1)*2);
            fitFunc3->SetParLimits(3,fitFunc1->GetParameter(1)/2,fitFunc1->GetParameter(1));
            // fitFunc3->SetParLimits(3,-100,0);
            // Fit for Interaction 
            TF1 *fitFunc5 = new TF1("fitFunc5", "[0]*25.5/[1]*exp(-x/[1])", 50,300); fitFunc5->SetParameters(1e5,200); fitFunc5->SetLineColor(kRed);  fitFunc5->FixParameter(0,(h_p_tot->Integral()));
            TF1 *fitFunc6 = new TF1("fitFunc6", "[0]*25.5/[1]*exp(-x/[1])", 50,300); fitFunc6->SetParameters(1e5,200); fitFunc6->SetLineColor(kBlue); fitFunc6->FixParameter(0,(h_d_tot->Integral()));
            TF1 *fitFunc4 = new TF1("fitFunc4", "[0]*25.5*(1-[1]) /[2] *exp(-x/[2])+ [0]*25.5*[1]/[3] *exp(-x/ [3] )", 50,300); fitFunc4->SetParameters(1e5, Ratio[i],200,170); fitFunc4->SetLineColor(kBlack); 
            fitFunc4->FixParameter(0,(h_p_tot->Integral()+h_d_tot->Integral()));
            fitFunc4->FixParameter(2,fitFunc5->GetParameter(1));
            // fitFunc4->FixParameter(3,fitFunc6->GetParameter(1));
            fitFunc4->SetParLimits(1,1e-3,1);
            // fitFunc4->SetParLimits(2,fitFunc5->GetParameter(1)/2,fitFunc5->GetParameter(1)*2);
            fitFunc4->SetParLimits(3,fitFunc5->GetParameter(1)/2,fitFunc5->GetParameter(1));
            // fitFunc4->SetParLimits(3,-100,0);
            TLatex latex;
            latex.SetTextSize(0.04);
            latex.SetTextFont(72);
            latex.SetTextAlign(13);  //align at top

            auto c0 = new TCanvas("c0","c0",2400,1000);
            c0->Divide(2,1);
            c0->cd(1);
            gStyle->SetOptStat(0);
            h_p_sur->SetLineColor(kRed);     h_p_sur->SetLineWidth(2);
            h_d_sur->SetLineColor(kBlue);    h_d_sur->SetLineWidth(2);
            h_2_sur->SetLineColor(kBlack);   h_2_sur->SetLineWidth(2);
            h_2_sur->SetTitle("N_{survive};Depth(mm);Counts");
            h_2_sur->GetYaxis()->SetRangeUser(0,1.2*h_2_sur->GetMaximum());
            h_2_sur->Draw("hist");
            h_p_sur->Draw("histsame");
            h_d_sur->Draw("histsame");
            h_p_sur->Fit(fitFunc1, "SR"); 
            h_d_sur->Fit(fitFunc2, "SR"); 
            h_2_sur->Fit(fitFunc3, "SR"); 
            fitFunc1->Draw("same");
            fitFunc2->Draw("same");
            fitFunc3->Draw("same");

            double proton_lambda       = fitFunc3->GetParameter(2);
            double proton_lambda_err   = fitFunc3->GetParError(2);
            double deuteron_lambda     = fitFunc3->GetParameter(3);
            double deuteron_lambda_err = fitFunc3->GetParError(3);
            double deuteron_ratio      = fitFunc3->GetParameter(1);
            double deuteron_ratio_err  = fitFunc3->GetParError(1);
            double lambda1     = fitFunc1->GetParameter(1);
            double lambda1_err = fitFunc1->GetParError(1);
            double lambda2     = fitFunc2->GetParameter(1);
            double lambda2_err = fitFunc2->GetParError(1);
            latex.DrawLatex(0,h_2_sur->GetMaximum()-1000,"Function: N_{sur} = N_{tot} #upoint [ (1-^{}r^{}_{d}) #upoint e^{-x/#lambda_{p}} + r_{d} #upoint e^{-x/#lambda_{d}} ] ");
            latex.DrawLatex(0,h_2_sur->GetMaximum()-2000,Form("Fitted r_{d}: %.2f#pm %.2f",deuteron_ratio , deuteron_ratio_err ));
            latex.DrawLatex(0,h_2_sur->GetMaximum()-3000,Form("Fitted #lambda_{d}: %.2f#pm %.2f mm",deuteron_lambda , deuteron_lambda_err ));
            latex.DrawLatex(0,h_2_sur->GetMaximum()-4000,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda1 , lambda1_err));
            latex.DrawLatex(0,h_2_sur->GetMaximum()-5000,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda2 , lambda2_err));

            c0->cd(2);
            h_p_int->SetLineColor(kRed);     h_p_int->SetLineWidth(2);
            h_d_int->SetLineColor(kBlue);    h_d_int->SetLineWidth(2);
            h_2_int->SetLineColor(kBlack);   h_2_int->SetLineWidth(2);
            h_2_int->SetTitle("N_{interaction};Depth(mm);Counts");
            h_2_int->GetYaxis()->SetRangeUser(0,1.2*h_2_int->GetMaximum());
            h_2_int->Draw("hist");
            h_p_int->Draw("histsame");
            h_d_int->Draw("histsame");
            h_p_int->Fit(fitFunc5,"R");
            h_d_int->Fit(fitFunc6,"R");
            h_2_int->Fit(fitFunc4,"R"); 
            fitFunc4->Draw("same");
            fitFunc5->Draw("same");
            fitFunc6->Draw("same");

            double proton_lambda1       = fitFunc4->GetParameter(2);
            double proton_lambda1_err   = fitFunc4->GetParError(2);
            double deuteron_lambda1     = fitFunc4->GetParameter(3);
            double deuteron_lambda1_err = fitFunc4->GetParError(3);
            double deuteron_ratio1      = fitFunc4->GetParameter(1);
            double deuteron_ratio1_err  = fitFunc4->GetParError(1);
            double lambda3     = fitFunc5->GetParameter(1);
            double lambda3_err = fitFunc5->GetParError(1);
            double lambda4     = fitFunc6->GetParameter(1);
            double lambda4_err = fitFunc6->GetParError(1);
            latex.DrawLatex(0,h_2_int->GetMaximum()-50,"Function: N_{int} = N_{tot} #upoint [  #frac{ 1-^{}r^{}_{d} }{#lambda_{p} } #upoint e^{-x/#lambda_{p} } + #frac{r_{d}}{#lambda_{d}} #upoint e^{-x/#lambda_{d}} ]");
            latex.DrawLatex(0,h_2_int->GetMaximum()-200,Form("Fitted r_{d}: %.2f#pm %.2f",deuteron_ratio1 , deuteron_ratio1_err ));
            latex.DrawLatex(0,h_2_int->GetMaximum()-300,Form("Fitted #lambda_{d}: %.2f#pm %.2f mm",deuteron_lambda1 , deuteron_lambda1_err ));
            latex.DrawLatex(0,h_2_int->GetMaximum()-400,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda3 , lambda3_err));
            latex.DrawLatex(0,h_2_int->GetMaximum()-500,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda4 , lambda4_err));


            gre_lambda_d_int->SetPoint(i,0.99*Ratio[i],deuteron_lambda);
            gre_lambda_d_int->SetPointError(i,0,deuteron_lambda_err);
            gre_lambda_d_sur->SetPoint(i,1.01*Ratio[i],deuteron_lambda1);
            gre_lambda_d_sur->SetPointError(i,0,deuteron_lambda1_err);

            gre_sur->SetPoint(i,0.99*Ratio[i],deuteron_ratio);
            gre_sur->SetPointError(i,0,deuteron_ratio_err);
            gre_int->SetPoint(i,1.01*Ratio[i],deuteron_ratio1);
            gre_int->SetPointError(i,0,deuteron_ratio1_err);

            hist_sur->SetBinContent(k+1,i+1, deuteron_ratio);
            hist_int->SetBinContent(k+1,i+1, deuteron_ratio1);
            hist_sur->SetBinError(k+1,i+1, deuteron_ratio_err);
            hist_int->SetBinError(k+1,i+1, deuteron_ratio1_err);
            lamb_sur->SetBinContent(k+1,i+1, deuteron_lambda);
            lamb_int->SetBinContent(k+1,i+1, deuteron_lambda1);
            lamb_sur->SetBinError(k+1,i+1, deuteron_lambda_err);
            lamb_int->SetBinError(k+1,i+1, deuteron_lambda1_err);

            c0->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/CrossS/RatioEnergySearch/Ratio_%.2f_%dGeV.pdf",Ratio[i],Energy_Name[k]));
        } // i Ratio

        auto c2 = new TCanvas("c2","c2",1900,900);
        gre_int->SetLineColor(kRed);
        gre_sur->SetLineColor(kBlue);
        gre_int->SetLineWidth(2);
        gre_sur->SetLineWidth(2);
        gre_int->SetMarkerColor(kRed);
        gre_sur->SetMarkerColor(kBlue);
        gre_int->SetMarkerStyle(20);
        gre_sur->SetMarkerStyle(21);

        auto line_ref = new TLine(5e-4,5e-4,2,2);
        line_ref->SetLineColor(kBlack);
        line_ref->SetLineWidth(2);
        line_ref->SetLineStyle(2);

        auto line_vet = new TLine(0.5,5e-4,0.5,2);
        line_vet->SetLineColor(kOrange-3);
        line_vet->SetLineWidth(2);
        line_vet->SetLineStyle(2);

        c2->cd();
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gre_int->GetXaxis()->SetLimits(5e-4,2);
        gre_int->GetYaxis()->SetRangeUser(5e-4,2);
        gre_int->SetTitle(";True r_{d};Fitted r_{d}");
        gre_int->Draw("AP");
        gre_sur->Draw("PSAME");
        line_ref->Draw("same");
        line_vet->Draw("same");

        auto lg2 = new TLegend(0.68,0.12,0.88,0.32);
        lg2->AddEntry(gre_int,"Fitted From N_{int}","ep");
        lg2->AddEntry(gre_sur,"Fitted From N_{sur}","ep");
        lg2->AddEntry(line_ref,"y = x","l");
        lg2->AddEntry(line_vet,"x = 0.5","l");
        lg2->Draw();
        c2->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/CrossS/RatioEnergySearch/FittedRatio_%dGeV.pdf",Energy_Name[k]));

        auto c3 = new TCanvas("c3","c3",1900,900);
        gre_lambda_d_int->SetLineColor(kRed);
        gre_lambda_d_sur->SetLineColor(kBlue);
        gre_lambda_d_int->SetLineWidth(2);
        gre_lambda_d_sur->SetLineWidth(2);
        gre_lambda_d_int->SetMarkerColor(kRed);
        gre_lambda_d_sur->SetMarkerColor(kBlue);
        gre_lambda_d_int->SetMarkerStyle(20);
        gre_lambda_d_sur->SetMarkerStyle(21);

        c3->cd();
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gre_lambda_d_int->GetXaxis()->SetLimits(5e-4,2e0);
        gre_lambda_d_int->GetYaxis()->SetRangeUser(0,300);
        gre_lambda_d_int->SetTitle(";True r_{d};Fitted #lambda_{d} (mm)");
        gre_lambda_d_int->Draw("AP");
        gre_lambda_d_sur->Draw("PSAME");
        auto lg3 = new TLegend(0.68,0.12,0.88,0.32);
        lg3->AddEntry(gre_int,"Fitted From N_{int}","ep");
        lg3->AddEntry(gre_sur,"Fitted From N_{sur}","ep");
        lg3->Draw();
    }// k Energy


    auto file2 = TFile::Open("/Users/xiongzheng/software/Hadr00/build/deuteron_BGO.root");
    auto hist2 = (TH1D*)file2->Get("h4");
    auto hist7 = (TH1D*)file2->Get("h1");
    int Nbins = hist2->GetNbinsX();
    double KN_Energy[60] = {0};
    double Deuteron[60] = {0};
    double KN_Deuteron[60] = {0};
    for (int ii = 0 ; ii < Nbins ; ii++)
    {
        KN_Energy[ii]   = pow(10,hist2->GetBinCenter(ii+1)-3);
        file2->cd();Deuteron[ii] = hist2->GetBinContent(ii+1); KN_Deuteron[ii] = hist7->GetBinContent(ii+1)* 10.;
    }
    auto gre7 = new TGraph(Nbins,KN_Energy, KN_Deuteron);   gre7->SetLineColor(kBlue);

    auto c4 = new TCanvas("c4","c4",2500,2500);
    TMultiGraph* mg1 = new TMultiGraph();
    TMultiGraph* mg2 = new TMultiGraph();
    TMultiGraph* mg3 = new TMultiGraph();
    TMultiGraph* mg4 = new TMultiGraph();

    for (int j = 0; j < 18; j++) {
        TH1D* proj = hist_sur->ProjectionX(Form("proj_%d", j), j+1, j+1); proj->SetLineColor(j+11);  proj->SetLineWidth(2);
        TH1D* prok = hist_int->ProjectionX(Form("prok_%d", j), j+1, j+1); prok->SetLineColor(j+11);  prok->SetLineWidth(2);
        TH1D* prol = lamb_sur->ProjectionX(Form("prol_%d", j), j+1, j+1); prol->SetLineColor(j+11);  prol->SetLineWidth(2);
        TH1D* prom = lamb_int->ProjectionX(Form("prom_%d", j), j+1, j+1); prom->SetLineColor(j+11);  prom->SetLineWidth(2);

        TGraph* g1 = new TGraph(proj);  g1->SetLineColor(j+11); g1->SetTitle(Form("r_{d}=%.2f", hist_sur->GetYaxis()->GetBinCenter(j+1))); mg1->Add(g1); // graphs.push_back(g1);
        TGraph* g2 = new TGraph(prok);  g2->SetLineColor(j+11); g2->SetTitle(Form("r_{d}=%.2f", hist_int->GetYaxis()->GetBinCenter(j+1))); mg2->Add(g2); // graphs.push_back(g2);
        TGraph* g3 = new TGraph(prol);  g3->SetLineColor(j+11); g3->SetTitle(Form("r_{d}=%.2f", lamb_sur->GetYaxis()->GetBinCenter(j+1))); mg3->Add(g3); // graphs.push_back(g3);
        TGraph* g4 = new TGraph(prom);  g4->SetLineColor(j+11); g4->SetTitle(Form("r_{d}=%.2f", lamb_int->GetYaxis()->GetBinCenter(j+1))); mg4->Add(g4); // graphs.push_back(g4);
    }

    c4->Divide(2,2);
    c4->cd(1);
    mg1->GetYaxis()->SetRangeUser(5e-4,2);
    mg1->GetXaxis()->SetLimits(9e0,1.1e4);
    mg1->GetYaxis()->SetTitleSize(0.05);
    mg1->GetXaxis()->SetTitleSize(0.05);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    mg1->GetHistogram()->SetTitle("Fitted r_{d} vs Energy at fixed true r_{d};Energy (GeV);Fitted r_{d} based on N_{sur}");
    mg1->Draw("AL");

    c4->cd(2);
    mg2->GetYaxis()->SetRangeUser(5e-4,2);
    mg2->GetXaxis()->SetLimits(9e0,1.1e4);
    mg2->GetYaxis()->SetTitleSize(0.05);
    mg2->GetXaxis()->SetTitleSize(0.05);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    mg2->GetHistogram()->SetTitle("Fitted r_{d} vs Energy at fixed true r_{d};Energy (GeV);Fitted r_{d} based on N_{int}");
    mg2->Draw("AL");

    c4->cd(3);
    mg3->GetYaxis()->SetRangeUser(20,300);
    mg3->GetXaxis()->SetLimits(9e0,1.1e4);
    mg3->GetYaxis()->SetTitleSize(0.05);
    mg3->GetXaxis()->SetTitleSize(0.05);
    gPad->SetLogx(1);
    mg3->GetHistogram()->SetTitle("Fitted #lambda_{d} vs Energy at fixed true r_{d};Energy (GeV);Fitted #lambda_{d} based on N_{sur}");
    mg3->Draw("AL");
    gre7->Draw("LSAME");

    auto lg4 = new TLegend(0.68,0.12,0.88,0.28);
    lg4->AddEntry(mg3,"Fitted #lambda_{d}","l");
    lg4->AddEntry(gre7,"GEANT4 #lambda_{d}","l");
    lg4->Draw();

    c4->cd(4);
    mg4->GetYaxis()->SetRangeUser(20,300);
    mg4->GetXaxis()->SetLimits(9e0,1.1e4);
    mg4->GetYaxis()->SetTitleSize(0.05);
    mg4->GetXaxis()->SetTitleSize(0.05);
    gPad->SetLogx(1);
    mg4->GetHistogram()->SetTitle("Fitted #lambda_{d} vs Energy at fixed true r_{d};Energy (GeV);Fitted #lambda_{d} based on N_{int}");
    mg4->Draw("AL");
    gre7->Draw("LSAME");
    lg4->Draw();

}