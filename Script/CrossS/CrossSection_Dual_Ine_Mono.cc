void CrossSection_Dual_Ine_Mono()
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

    for (int k =18; k < 19; k++)
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
        auto gre_p_int = new TGraphErrors();  
        auto gre_p_sur = new TGraphErrors();
        auto gre_d_int = new TGraphErrors();  
        auto gre_d_sur = new TGraphErrors();

        // for (int i =27; i < 28; i++)
        for (int i =0; i < 27; i++)
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

            h_p_int->Scale( 1-Ratio[i]);    h_p_tot->Scale(1-Ratio[i]);
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

            for(int ii = 0; ii<= h_2_int-> GetNbinsX() ; ii++)
            {
                cout << "Check Error " << h_2_int->GetBinError(ii) << " Proton " << h_p_int->GetBinError(ii) << " Deuteron " << h_d_int->GetBinError(ii) << endl;
            }

            // Fit for Survive
            TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x+12.75)/[1])", 60,300); fitFunc1->SetParameters(1e4,200); fitFunc1->SetLineColor(kRed);  fitFunc1->FixParameter(0,h_p_tot->Integral());
            TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x+12.75)/[1])", 60,300); fitFunc2->SetParameters(1e4,200); fitFunc2->SetLineColor(kBlue); fitFunc2->FixParameter(0,h_d_tot->Integral());
            TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*exp(-(x+12.75)/[1])", 60,300); fitFunc3->SetParameters(1e4,200); fitFunc3->SetLineColor(kBlack); fitFunc3->FixParameter(0,(h_p_tot->Integral()+h_d_tot->Integral()));

            // Fit for Interaction 
            TF1 *fitFunc5 = new TF1("fitFunc5", "[0]/[1]*exp(-x/[1])", 60,300); fitFunc5->SetParameters(1e5,200); fitFunc5->SetLineColor(kRed);  fitFunc5->FixParameter(0,(h_p_tot->Integral()*h_p_int->GetBinWidth(1)));
            TF1 *fitFunc6 = new TF1("fitFunc6", "[0]/[1]*exp(-x/[1])", 60,300); fitFunc6->SetParameters(1e5,200); fitFunc6->SetLineColor(kBlue); fitFunc6->FixParameter(0,(h_d_tot->Integral()*h_d_int->GetBinWidth(1)));
            TF1 *fitFunc4 = new TF1("fitFunc4", "[0]/[1]*exp(-x/[1])", 60,300); fitFunc4->SetParameters(1e5,200); fitFunc4->SetLineColor(kBlack); fitFunc4->FixParameter(0,(h_p_tot->Integral()+h_d_tot->Integral())*h_p_int->GetBinWidth(1));

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
            h_p_sur->Fit(fitFunc1, "R"); 
            h_d_sur->Fit(fitFunc2, "R"); 
            h_2_sur->Fit(fitFunc3, "R"); 
            fitFunc1->Draw("same");
            fitFunc2->Draw("same");
            fitFunc3->Draw("same");

            double mixture_length      = fitFunc3->GetParameter(1);
            double mixture_length_err  = fitFunc3->GetParError(1);
            double lambda1     = fitFunc1->GetParameter(1);
            double lambda1_err = fitFunc1->GetParError(1);
            double lambda2     = fitFunc2->GetParameter(1);
            double lambda2_err = fitFunc2->GetParError(1);
            latex.DrawLatex(0,h_2_sur->GetMaximum()-1000,"Function: N_{sur} = N_{tot} #upoint e^{-x/#lambda_{x}}  ");
            latex.DrawLatex(0,h_2_sur->GetMaximum()-2000,Form("Fitted #lambda_{d}: %.2f#pm %.2f mm",mixture_length , mixture_length_err ));
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
            
            // 
            // fitFunc4->FixParameter(1,fitFunc5->GetParameter(1));
            h_2_int->Fit(fitFunc4,"R"); 
            fitFunc4->Draw("same");
            fitFunc5->Draw("same");
            fitFunc6->Draw("same");

            double mixture_length1      = fitFunc4->GetParameter(1);
            double mixture_length1_err  = fitFunc4->GetParError(1);
            double lambda3     = fitFunc5->GetParameter(1);
            double lambda3_err = fitFunc5->GetParError(1);
            double lambda4     = fitFunc6->GetParameter(1);
            double lambda4_err = fitFunc6->GetParError(1);
            latex.DrawLatex(0,h_2_int->GetMaximum()-50,"Function: N_{int} = #frac{ N_{tot} #upoint #Delta x }{#lambda_{x} }  #upoint e^{-x/#lambda_{x} } ");
            latex.DrawLatex(0,h_2_int->GetMaximum()-200,Form("Fitted r_{d}: %.2f#pm %.2f mm",mixture_length1 , mixture_length1_err ));
            latex.DrawLatex(0,h_2_int->GetMaximum()-400,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda3 , lambda3_err));
            latex.DrawLatex(0,h_2_int->GetMaximum()-500,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda4 , lambda4_err));

            gre_sur->SetPoint(i,Ratio[i],mixture_length);
            gre_sur->SetPointError(i,0,mixture_length_err);
            gre_int->SetPoint(i,Ratio[i],mixture_length1);
            gre_int->SetPointError(i,0,mixture_length1_err);

            gre_p_sur->SetPoint(i,Ratio[i],lambda1);
            gre_p_sur->SetPointError(i,0,lambda1_err);
            gre_p_int->SetPoint(i,Ratio[i],lambda3);
            gre_p_int->SetPointError(i,0,lambda3_err);

            gre_d_sur->SetPoint(i,Ratio[i],lambda2);
            gre_d_sur->SetPointError(i,0,lambda2_err);
            gre_d_int->SetPoint(i,Ratio[i],lambda4);
            gre_d_int->SetPointError(i,0,lambda4_err);

            hist_sur->SetBinContent(k+1,i+1, mixture_length);
            hist_int->SetBinContent(k+1,i+1, mixture_length1);
            hist_sur->SetBinError(k+1,i+1, mixture_length_err);
            hist_int->SetBinError(k+1,i+1, mixture_length1_err);

            // c0->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/CrossS/RatioEnergySearch/Ratio_%.2f_%dGeV.pdf",Ratio[i],Energy_Name[k]));
        } // i Ratio

        auto c2 = new TCanvas("c2","c2",2100,900);
        c2->Divide(2,1);
        c2->cd(1);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gre_int->SetLineColor(kBlack);
        gre_int->SetLineWidth(2);
        gre_int->SetMarkerColor(kBlack);
        gre_int->SetMarkerStyle(20);

        gre_p_int->SetLineColor(kRed);
        gre_p_int->SetLineWidth(2);
        gre_p_int->SetMarkerColor(kRed);
        gre_p_int->SetMarkerStyle(20);

        gre_d_int->SetLineColor(kBlue);
        gre_d_int->SetLineWidth(2);
        gre_d_int->SetMarkerColor(kBlue);
        gre_d_int->SetMarkerStyle(20);

        gre_int->GetXaxis()->SetLimits(5e-4,2);
        gre_int->GetYaxis()->SetRangeUser(50,250);
        gre_int->SetTitle("Fitted From N_{int};True r_{d};#lambda (mm)");
        gre_int->Draw("AP");
        gre_p_int->Draw("PSAME");
        gre_d_int->Draw("PSAME");

        auto lg2_1 = new TLegend(0.68,0.12,0.88,0.32);
        lg2_1->AddEntry(gre_int,"Mixture","ep");
        lg2_1->AddEntry(gre_d_int,"Deuteron","ep");
        lg2_1->AddEntry(gre_p_int,"Proton","ep");
        lg2_1->Draw();

        c2->cd(2);
        gPad->SetGrid(1,1);
        gPad->SetLogx();

        gre_sur->SetLineColor(kBlack);
        gre_sur->SetLineWidth(2);
        gre_sur->SetMarkerColor(kBlack);
        gre_sur->SetMarkerStyle(21);

        gre_p_sur->SetLineColor(kRed);
        gre_p_sur->SetLineWidth(2);
        gre_p_sur->SetMarkerColor(kRed);
        gre_p_sur->SetMarkerStyle(20);

        gre_d_sur->SetLineColor(kBlue);
        gre_d_sur->SetLineWidth(2);
        gre_d_sur->SetMarkerColor(kBlue);
        gre_d_sur->SetMarkerStyle(20);

        gre_sur->GetXaxis()->SetLimits(5e-4,2);
        gre_sur->GetYaxis()->SetRangeUser(50,250);
        gre_sur->SetTitle("Fitted From N_{sur};True r_{d};#lambda (mm)");
        gre_sur->Draw("AP");
        gre_p_sur->Draw("PSAME");
        gre_d_sur->Draw("PSAME");

        auto lg2_2 = new TLegend(0.68,0.12,0.88,0.32);
        lg2_2->AddEntry(gre_sur,"Mixture","ep");
        lg2_2->AddEntry(gre_d_sur,"Deuteron","ep");
        lg2_2->AddEntry(gre_p_sur,"Proton","ep");
        lg2_2->Draw();
        // c2->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/CrossS/RatioEnergySearch/FittedRatio_%dGeV.pdf",Energy_Name[k]));

    }// k Energy
}