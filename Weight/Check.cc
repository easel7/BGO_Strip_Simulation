void Check()
{

    TCut UBT = "(L0_E>0.0092 && L1_E>0.0092)";
    TCut HET = "(L0_E>0.23 && L1_E >0.23 &&   L2_E>0.23 && L3_E>0.046)";

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Root/Carbon_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");

    auto h1_p = new TH1D("h1_p", "h1_p",40,0,4);  
    auto h2_p = new TH1D("h2_p", "h2_p",40,0,4);  

    proton_tree->Draw("log10(Energy)>>h1_p","","");
    proton_tree->Draw("log10(Total_E)>>h2_p","","");

    h1_p->SetTitle("Before Weighted;log_{10}(Energy/GeV);Counts");
    h2_p->SetTitle("Before Weighted;log_{10}(E_dep/GeV);Counts");

    auto proton2_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Carbon_PowerLaw.root");
    auto proton2_tree = (TTree*)proton2_file->Get("B4");

    auto h1_w = new TH1D("h1_w", "h1_w",40,0,4);  
    auto h2_w = new TH1D("h2_w", "h2_w",40,0,4);  

    proton2_tree->Draw("log10(Energy)>>h1_w","weight","");
    proton2_tree->Draw("log10(Total_E)>>h2_w","weight","");

    h1_w->SetTitle("After Weighted;log_{10}(Energy/GeV);Counts");
    h2_w->SetTitle("After Weighted;log_{10}(E_dep/GeV);Counts");

    auto c1 = new TCanvas("c1","c1",900,900);
    c1->Clear();
    c1->cd();
    gStyle->SetOptStat(0);
    h1_p->Draw();

    auto c2 = new TCanvas("c2","c2",900,900);
    c2->Clear();
    c2->cd();
    gStyle->SetOptStat(0);
    h2_p->Draw();




    auto c4 = new TCanvas("c4","c4",900,900);
    c4->Clear();
    c4->cd();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    h2_w->Draw();

    auto c3 = new TCanvas("c3","c3",900,900);
    c3->Clear();
    c3->cd();
    gPad->SetLogy();
    gStyle->SetOptStat(0);



    h1_w->Draw();
    auto tf1 = new TF1("tf1","[0]*TMath::Power(pow(10,x), [1])",1,3.5);
    h1_w->Fit(tf1,"R");

    // Get slope parameter from fit
    double slope     = tf1->GetParameter(1);
    double slope_err = tf1->GetParError(1);

    // Add slope text to the canvas
    TLatex *tex = new TLatex();
    tex->SetNDC();  // Use normalized coordinates (0-1)
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.2, 0.3, Form("Fitted index = %.2f #pm %.2f", slope,slope_err));
}

