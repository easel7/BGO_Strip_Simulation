void Edep_TrackLength()
{
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Root/Helium4_100GeV.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    auto c1 = new TCanvas("c1","c1",1200,1200);

    auto L0 = new TH2D("L0","L0",50,0,25,50,0,25);
    auto L1 = new TH2D("L1","L1",50,0,25,50,0,25);
    auto L2 = new TH2D("L2","L2",50,0,25,50,0,25);
    auto L3 = new TH2D("L3","L3",50,0,25,50,0,25);

    TCut UBT = "(L0_E>0.0092 && L1_E>0.0092)";
    // TCut HET = "(L0_E>0.23 && L1_E >0.23 &&   L2_E>0.23 && L3_E>0.046)";
    TCut HET = "";


    proton_tree->Draw("L0_E:L0_L  >>L0",HET,"colz");
    proton_tree->Draw("L4_E:L4_L  >>L1",HET,"colz");
    proton_tree->Draw("L8_E:L8_L  >>L2",HET,"colz");
    proton_tree->Draw("L12_E:L12_L>>L3",HET,"colz");

    c1->Clear();
    // gStyle->SetOptStat(0);
    c1->Divide(2,2);
    
    c1->cd(1);
    L0->SetTitle("Edep vs TrackLength in L0;TrackLength(m);Energy(GeV)");
    L0->Draw("colz");

    c1->cd(2);
    L1->SetTitle("Edep vs TrackLength in L4;TrackLength(m);Energy(GeV)");
    L1->Draw("colz");

    c1->cd(3);
    L2->SetTitle("Edep vs TrackLength in L8;TrackLength(m);Energy(GeV)");
    L2->Draw("colz");

    c1->cd(4);
    L3->SetTitle("Edep vs TrackLength in L12;TrackLength(m);Energy(GeV)");
    L3->Draw("colz");

}