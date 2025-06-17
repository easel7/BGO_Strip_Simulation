void CrossSection_Dual_Ine_PowerLaw()
{
    auto file_p  = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_PowerLaw.root");
    auto tree_p  = (TTree*)file_p->Get("B4");
    auto h_p_tot = new TH1D("h_p_tot","h_p_tot",15,-25.5,357);  h_p_tot->Sumw2();
    auto h_p_int = new TH1D("h_p_int","h_p_int",14,0,357);      h_p_int->Sumw2();
    auto h_p_sur = new TH1D("h_p_sur","h_p_sur",14,0,357);      
    tree_p->Draw("First_Ine_Depth>>h_p_tot","(Total_E>40 && Total_E<49)","");
    tree_p->Draw("First_Ine_Depth>>h_p_int","(Total_E>40 && Total_E<49)","");

    auto file_d  = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_PowerLaw.root");
    auto tree_d  = (TTree*)file_d->Get("B4");
    auto h_d_tot = new TH1D("h_d_tot","h_d_tot",15,-25.5,357);  h_d_tot->Sumw2();
    auto h_d_int = new TH1D("h_d_int","h_d_int",14,0,357);      h_d_int->Sumw2();
    auto h_d_sur = new TH1D("h_d_sur","h_d_sur",14,0,357);     
    tree_d->Draw("First_Ine_Depth>>h_d_tot"," (Total_E>40 && Total_E<49 )","");
    tree_d->Draw("First_Ine_Depth>>h_d_int"," (Total_E>40 && Total_E<49 )","");

    // h_p_int->Scale(0.5);    h_p_tot->Scale(0.5);
    // h_d_int->Scale(0.5);    h_d_tot->Scale(0.5);

    cout << h_p_tot->Integral() << " , " <<  h_d_tot->Integral()  << endl;
    for(int ii = 1 ; ii <= h_p_int->GetNbinsX() ; ii++)
    {
        h_p_sur->SetBinContent(ii, ( h_p_tot->Integral() - h_p_int->Integral(1,ii) ) );
        h_d_sur->SetBinContent(ii, ( h_d_tot->Integral() - h_d_int->Integral(1,ii) ) );
        cout << ii << " , Sum N_{int} " << h_p_int->Integral(1,ii) << " ,  N_{int} " <<  h_p_int->GetBinContent(ii) << " , N_{sur} " << h_p_sur->GetBinContent(ii)  << endl;
    }
    auto h_2_int = (TH1D*)h_p_int->Clone(); h_2_int->Add(h_d_int);
    auto h_2_sur = (TH1D*)h_p_sur->Clone(); h_2_sur->Add(h_d_sur);
    h_2_int->Sumw2();    h_p_int->Sumw2();
    h_2_sur->Sumw2();    h_d_sur->Sumw2();

    // Fit for Survive
    TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*exp(-(x+12.75)/[1])", 0,250); fitFunc1->SetParameters(1e4, 200); fitFunc1->SetLineColor(kRed);  // fitFunc1->FixParameter(0,h_p_tot->Integral());
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*exp(-(x+12.75)/[1])", 0,250); fitFunc2->SetParameters(1e4, 200); fitFunc2->SetLineColor(kBlue); // fitFunc2->FixParameter(0,h_d_tot->Integral());

    // Fit for Interaction 
    TF1 *fitFunc5 = new TF1("fitFunc5", "[0]*25.5 / [1]*exp(-x/[1])", 0,250); fitFunc5->SetParameters(1e5,200); fitFunc5->SetLineColor(kRed);   fitFunc5->FixParameter(0,h_p_tot->Integral());
    TF1 *fitFunc6 = new TF1("fitFunc6", "[0]*25.5 / [1]*exp(-x/[1])", 0,250); fitFunc6->SetParameters(1e5,200); fitFunc6->SetLineColor(kBlue);  fitFunc6->FixParameter(0,h_d_tot->Integral());
    
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

    TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*exp(-(x+12.75)/[1])", 0,250); 
    fitFunc3->SetParameters(1e4,200); 
    fitFunc3->SetLineColor(kBlack); 
    // fitFunc3->FixParameter(0,(h_p_tot->Integral()+h_d_tot->Integral()));

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
    latex.DrawLatex(0,h_2_sur->GetMaximum()-10000,"Function: N_{sur} = N_{tot} #upoint e^{-x/#lambda_{x}}  ");
    latex.DrawLatex(0,h_2_sur->GetMaximum()-20000,Form("Fitted #lambda_{d}: %.2f#pm %.2f mm",mixture_length , mixture_length_err ));
    latex.DrawLatex(0,h_2_sur->GetMaximum()-40000,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda1 , lambda1_err));
    latex.DrawLatex(0,h_2_sur->GetMaximum()-50000,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda2 , lambda2_err));


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
    TF1 *fitFunc4 = new TF1("fitFunc4", "[0]*25.5/[1]*exp(-x/[1])", 0,250); 
    fitFunc4->SetParameters(1e5,200); 
    fitFunc4->SetLineColor(kBlack); 
    fitFunc4->FixParameter(0,(h_p_tot->Integral()+h_d_tot->Integral()));
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
    latex.DrawLatex(0,h_2_int->GetMaximum()-1000,"Function: N_{int} = #frac{ N_{tot} #upoint #Delta x }{#lambda_{x} }  #upoint e^{-x/#lambda_{x} } ");
    latex.DrawLatex(0,h_2_int->GetMaximum()-2000,Form("Fitted r_{d}: %.2f#pm %.2f mm",mixture_length1 , mixture_length1_err ));
    latex.DrawLatex(0,h_2_int->GetMaximum()-4000,Form("#color[2]{Fitted #lambda_{p} alone: %.2f#pm %.2f mm}",lambda3 , lambda3_err));
    latex.DrawLatex(0,h_2_int->GetMaximum()-5000,Form("#color[4]{Fitted #lambda_{d} alone: %.2f#pm %.2f mm}",lambda4 , lambda4_err));

    c0->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/CrossS/CrossSection_Dual_PowerLaw.pdf");
}