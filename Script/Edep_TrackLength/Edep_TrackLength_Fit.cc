void Edep_TrackLength_Fit()
{
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Root/Electron_100GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    auto c1 = new TCanvas("c1","c1",1200,1200);

    auto L0 = new TH2D("L0","L0",50,0,25,50,0,25);
    auto L1 = new TProfile("L1","L1",50,0,25,0,25,"s");
    
    TCut UBT = "(L0_E>0.0092 && L1_E>0.0092)";
    TCut HET = "(L0_E>0.23 && L1_E >0.23 &&   L2_E>0.23 && L3_E>0.046)";

    proton_tree->Draw("L6_E:L6_L  >>L0",UBT,"colz");
    proton_tree->Draw("L6_E:L6_L   >>L1",UBT,"profs");


    c1->Clear();

    gStyle->SetOptStat(0);
    L0->SetTitle("100 GeV ELectron Edep vs TrackLength in L6;TrackLength(m);Energy(GeV)");
    L0->Draw("colz");
    L1->SetLineColor(kRed);
    L1->SetLineWidth(2);
    L1->Draw("profsame");
    auto fitFunc = new TF1("fitFunc", "pol1", 5, 15);
    L1->Fit(fitFunc, "", "", 5, 15);

    // Get slope parameter from fit
    double slope = fitFunc->GetParameter(1);
    double slope_err = fitFunc->GetParError(1);


    // Add slope text to the canvas
    TLatex *tex = new TLatex();
    tex->SetNDC();  // Use normalized coordinates (0-1)
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.3, 0.7, Form("Fitted Slope = %.2f #pm %.2f", slope,slope_err));


}