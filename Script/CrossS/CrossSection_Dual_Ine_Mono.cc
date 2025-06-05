void CrossSection_Dual_Ine_Mono()
{
    int Energy_Name[28]={0};
    int Energy_Name_Err[28]={0};
    double Ratio[28] = {0};
    double CL95_Sur[28]  = {0};    double CL95_Sur_Err[28]  = {0};    
    double CL90_Sur[28]  = {0};    double CL90_Sur_Err[28]  = {0};    
    double CL95_Int[28]  = {0};    double CL95_Int_Err[28]  = {0};    
    double CL90_Int[28]  = {0};    double CL90_Int_Err[28]  = {0};    
    double Energy[28]= {0};    double Energy_Err[28]= {0};

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

    // for (int k =18; k < 19; k++)
    for (int k =0; k < 28; k++)
    {
        if (k < 10)      {Energy_Name[k] = (k + 1) * 10;         Energy[k] = (k + 1) * 10.;       }       // 10 ~ 100
        else if (k < 19) {Energy_Name[k] = (k - 9 + 1) * 100;    Energy[k] = (k - 9 + 1) * 100.;  }        // 200 ~ 1000
        else             {Energy_Name[k] = (k - 18 + 1) * 1000;  Energy[k] = (k - 18 + 1) * 1000.;}        // 2000 ~ 10000

        cout << " Processing Energy " << Energy[k] << " GeV !!!! " << endl; 
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

        auto chi2_int = new TGraph();
        auto chi2_sur = new TGraph();

        CL95_Sur[k]=1e-3; int Sur_95_Tag = -1;
        CL90_Sur[k]=1e-3; int Sur_90_Tag = -1;
        CL95_Int[k]=1e-3; int Int_95_Tag = -1;
        CL90_Int[k]=1e-3; int Int_90_Tag = -1;
        auto range = ROOT::Fit::DataRange();
        range.SetRange(60, 300);
        auto opt = ROOT::Fit::DataOptions();

        // for (int i =16; i < 17; i++)
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
            if (gDirectory->FindObject("C_p_int")) delete gDirectory->FindObject("C_p_int");
            if (gDirectory->FindObject("C_p_sur")) delete gDirectory->FindObject("C_p_sur");
            auto h_p_tot = new TH1D("h_p_tot","h_p_tot",15,-25.5,357);  h_p_tot->Sumw2();
            auto h_p_int = new TH1D("h_p_int","h_p_int",14,0,357);      h_p_int->Sumw2();
            auto h_p_sur = new TH1D("h_p_sur","h_p_sur",14,0,357);      
            auto C_p_int = new TH1D("C_p_int","C_p_int",14,0,357);      C_p_int->Sumw2();
            auto C_p_sur = new TH1D("C_p_sur","C_p_sur",14,0,357);      
            tree_p->Draw("First_Ine_Depth>>h_p_tot","","");
            tree_p->Draw("First_Ine_Depth>>h_p_int",HI,"");
            tree_p->Draw("First_Ine_Depth>>C_p_int",HI,"");

            file_d->cd();
            if (gDirectory->FindObject("h_d_tot")) delete gDirectory->FindObject("h_d_tot");
            if (gDirectory->FindObject("h_d_int")) delete gDirectory->FindObject("h_d_int");
            if (gDirectory->FindObject("h_d_sur")) delete gDirectory->FindObject("h_d_sur");
            auto h_d_tot = new TH1D("h_d_tot","h_d_tot",15,-25.5,357);  h_d_tot->Sumw2();
            auto h_d_int = new TH1D("h_d_int","h_d_int",14,0,357);      h_d_int->Sumw2();
            auto h_d_sur = new TH1D("h_d_sur","h_d_sur",14,0,357);     
            tree_d->Draw("First_Ine_Depth>>h_d_tot","","");
            tree_d->Draw("First_Ine_Depth>>h_d_int",HI,"");

            cout << h_p_tot->Integral() << endl;
            for(int ii = 1 ; ii <= C_p_int->GetNbinsX() ; ii++)
            {
                C_p_sur->SetBinContent(ii , ( h_p_tot->Integral() - C_p_int->Integral(1,ii) ));
                // cout << ii << " , " << C_p_sur->GetBinContent(ii) << endl;
            }
            C_p_sur->Sumw2();

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

            // for(int ii = 0; ii<= h_2_int->GetNbinsX() ; ii++)
            // {
            //     cout << "Check Error " << h_2_int->GetBinError(ii) << " Proton " << h_p_int->GetBinError(ii) << " Deuteron " << h_d_int->GetBinError(ii) << endl;
            // }

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
            gStyle->SetOptStat(0);
            c0->Divide(2,1);
            c0->cd(1);
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


            cout << (h_p_tot->Integral()+h_d_tot->Integral()) << endl;
            auto data_sur = ROOT::Fit::BinData(opt, range);
            ROOT::Fit::FillData(data_sur, C_p_sur); 
            ROOT::Fit::FillData(data_sur, h_2_sur); 
            auto *f1_sur = new TF1("f1_sur",Form("%.2f*exp(-(x+12.75)/[0])", (h_p_tot->Integral()+h_d_tot->Integral()) ),60,300);
            f1_sur->SetParameter(0,200);
            auto wf_sur = ROOT::Math::WrappedTF1(*f1_sur);
            auto fitter_sur = ROOT::Fit::Fitter() ;
            fitter_sur.SetFunction(wf_sur);
            fitter_sur.Fit(data_sur);
            auto result_sur = fitter_sur.Result();
            result_sur.Print(std::cout);
            cout << result_sur.Chi2() << endl;
            cout << fitFunc1->GetChisquare() << endl;
            cout << fitFunc3->GetChisquare() << endl;
            double Chi2_combine_sur   = result_sur.Chi2();
            double Chi2_proton_sur    = fitFunc1->GetChisquare();
            double Chi2_deuteron_sur  = fitFunc2->GetChisquare();
            double Chi2_mixture_sur   = fitFunc3->GetChisquare();

            auto c0_1 = new TCanvas("c0_1","c0_1",2400,1000);
            c0_1->cd();
            h_2_sur->Draw("hist");
            C_p_sur->Draw("histsame");
            f1_sur->Draw("same");
    
            auto lg01 = new TLegend(0.72,0.72,0.88,0.88);
            lg01->AddEntry(C_p_sur,"Proton Sample","l");
            lg01->AddEntry(h_2_sur,"Mixture","l");
            lg01->AddEntry(f1_sur ,"Function share same #lambda","l");
            lg01->Draw();

            if( (Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur) ) >=2.71 && Sur_90_Tag < 0) 
            { 
                CL90_Sur[k] = Ratio[i];
                Sur_90_Tag = 1;
            }
            if( (Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur) ) >=3.84 && Sur_95_Tag < 0) 
            { 
                CL95_Sur[k] = Ratio[i];
                Sur_95_Tag = 1;
            }

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

            auto data_int = ROOT::Fit::BinData(opt, range);
            ROOT::Fit::FillData(data_int, C_p_int); // NDF = 13
            ROOT::Fit::FillData(data_int, h_2_int); //
            auto *f1_int = new TF1("f1_int",Form("%.2f/[0]*exp(-x/[0])", (h_p_tot->Integral()+h_d_tot->Integral())*h_p_int->GetBinWidth(1) ),60,300);
            f1_int->SetParameter(0,200);
            auto wf_int = ROOT::Math::WrappedTF1(*f1_int);
            auto fitter_int = ROOT::Fit::Fitter() ;
            fitter_int.SetFunction(wf_int);
            fitter_int.Fit(data_int);
            auto result_int = fitter_int.Result();
            result_int.Print(std::cout);
            cout << result_int.Chi2() << endl;
            cout << fitFunc4->GetChisquare() << endl;
            cout << fitFunc5->GetChisquare() << endl;
            double Chi2_combine_int   = result_int.Chi2();
            double Chi2_proton_int    = fitFunc5->GetChisquare();
            double Chi2_deuteron_int  = fitFunc6->GetChisquare();
            double Chi2_mixture_int   = fitFunc4->GetChisquare();

            if( (Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int) ) >=2.71 && Int_90_Tag < 0) 
            { 
                CL90_Int[k] = Ratio[i];
                Int_90_Tag = 1;
            }
            if( (Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int) ) >=3.84 && Int_95_Tag < 0) 
            { 
                CL95_Int[k] = Ratio[i];
                Int_95_Tag = 1;
            }

            chi2_int->SetPoint(i,Ratio[i],Chi2_combine_int - (Chi2_mixture_int + Chi2_proton_int));
            chi2_sur->SetPoint(i,Ratio[i],Chi2_combine_sur - (Chi2_mixture_sur + Chi2_proton_sur));

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

        auto c2 = new TCanvas("c2","c2",2100,2100);
        c2->Divide(2,2);
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

        c2->cd(3);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        chi2_int->SetLineColor(kRed);
        chi2_int->SetLineWidth(2);
        chi2_int->GetXaxis()->SetLimits(5e-4,2);
        chi2_int->GetYaxis()->SetRangeUser(0,20);
        chi2_int->SetTitle("Fitted From N_{int};True r_{d};#Delta#chi^{2}");
        chi2_int->Draw("AL");

        c2->cd(4);
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        chi2_sur->SetLineColor(kRed);
        chi2_sur->SetLineWidth(2);
        chi2_sur->GetXaxis()->SetLimits(5e-4,2);
        chi2_sur->GetYaxis()->SetRangeUser(0,20);
        chi2_sur->SetTitle("Fitted From N_{sur};True r_{d};#Delta#chi^{2}");
        chi2_sur->Draw("AL");
        // c2->SaveAs(Form("/Users/xiongzheng/software/B4/B4e/Script/CrossS/RatioEnergySearch/FittedRatio_%dGeV.pdf",Energy_Name[k]));

        cout << "Survival   90 CL " << CL90_Sur[k] << " , 95 CL = " << CL95_Sur[k] << endl;
        cout << "Interaction 90 CL " << CL90_Int[k] << " , 95 CL = " << CL95_Int[k] << endl;


    }// k Energy
    auto gre_90_sur = new TGraph(28,Energy,CL90_Sur); //,Energy_Name_Err,CL90_Sur_Err
    auto gre_95_sur = new TGraph(28,Energy,CL95_Sur); //,Energy_Name_Err,CL95_Sur_Err

    auto gre_90_int = new TGraph(28,Energy,CL90_Int); //,Energy_Name_Err,CL90_Sur_Err
    auto gre_95_int = new TGraph(28,Energy,CL95_Int); //,Energy_Name_Err,CL95_Sur_Err

    gre_90_sur->SetLineColor(kBlue);
    gre_95_sur->SetLineColor(kRed);
    gre_90_sur->SetLineWidth(2);
    gre_95_sur->SetLineWidth(2);

    gre_90_int->SetLineColor(kBlue);
    gre_95_int->SetLineColor(kRed);
    gre_90_int->SetLineWidth(2);
    gre_95_int->SetLineWidth(2);

    auto c4 = new TCanvas("c4","c4",2500,1200);
    c4->Divide(2,1);
    c4->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    gPad->SetLogx();
    gre_90_sur->GetYaxis()->SetRangeUser(9e-4,2);
    gre_90_sur->GetXaxis()->SetLimits(9e0,2e4);
    gre_90_sur->SetTitle("Sensitivity Curve from N_{sur};Energy(GeV); True r_{d}");
    gre_90_sur->Draw("ALP");
    gre_95_sur->Draw("LPSAME");
    auto lg4 = new TLegend(0.72,0.12,0.88,0.22);
    lg4->AddEntry(gre_90_sur,"90% C.L.","l");
    lg4->AddEntry(gre_95_sur,"95% C.L.","l");
    lg4->Draw();

    c4->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    gPad->SetLogx();
    gre_90_int->GetYaxis()->SetRangeUser(9e-4,2);
    gre_90_int->GetXaxis()->SetLimits(9e0,2e4);
    gre_90_int->SetTitle("Sensitivity Curve from N_{int};Energy(GeV); True r_{d}");
    gre_90_int->Draw("ALP");
    gre_95_int->Draw("LPSAME");
    lg4->Draw();

}