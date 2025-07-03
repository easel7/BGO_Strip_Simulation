void EnergyMatrix()
{
    std::vector<double>* p_EnergyVec = nullptr;    double p_weight; double p_Total_E; double p_Energy; std::vector<double>* p_RMSVec = nullptr;
    std::vector<double>* d_EnergyVec = nullptr;    double d_weight; double d_Total_E; double d_Energy; std::vector<double>* d_RMSVec = nullptr;
    std::vector<double>* e_EnergyVec = nullptr;    double e_weight; double e_Total_E; double e_Energy; std::vector<double>* e_RMSVec = nullptr;
    std::vector<double>* h_EnergyVec = nullptr;    double h_weight; double h_Total_E; double h_Energy; std::vector<double>* h_RMSVec = nullptr;
    std::vector<double>* H_EnergyVec = nullptr;    double H_weight; double H_Total_E; double H_Energy; std::vector<double>* H_RMSVec = nullptr;
    std::vector<double>* c_EnergyVec = nullptr;    double c_weight; double c_Total_E; double c_Energy; std::vector<double>* c_RMSVec = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Total_E",&p_Total_E);
    proton_tree->SetBranchAddress("Energy",&p_Energy);
    proton_tree->SetBranchAddress("weight",&p_weight);

    proton_tree->SetBranchAddress("RMS",&p_RMSVec);
    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Total_E",&d_Total_E);
    deuteron_tree->SetBranchAddress("Energy",&d_Energy);
    deuteron_tree->SetBranchAddress("RMS",&d_RMSVec);
    deuteron_tree->SetBranchAddress("weight",&d_weight);

    auto electron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Electron_PowerLaw.root");
    auto electron_tree = (TTree*)electron_file->Get("B4");
    electron_tree->SetBranchAddress("LayerEnergyVector",&e_EnergyVec);
    electron_tree->SetBranchAddress("Total_E",&e_Total_E);
    electron_tree->SetBranchAddress("Energy",&e_Energy);
    electron_tree->SetBranchAddress("RMS",&e_RMSVec);
    electron_tree->SetBranchAddress("weight",&e_weight);

    auto helium4_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Helium4_PowerLaw.root");
    auto helium4_tree = (TTree*)helium4_file->Get("B4");
    helium4_tree->SetBranchAddress("LayerEnergyVector",&h_EnergyVec);
    helium4_tree->SetBranchAddress("Total_E",&h_Total_E);
    helium4_tree->SetBranchAddress("Energy",&h_Energy);
    helium4_tree->SetBranchAddress("RMS",&h_RMSVec);
    helium4_tree->SetBranchAddress("weight",&h_weight);

    auto helium3_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Helium3_PowerLaw.root");
    auto helium3_tree = (TTree*)helium3_file->Get("B4");
    helium3_tree->SetBranchAddress("LayerEnergyVector",&H_EnergyVec);
    helium3_tree->SetBranchAddress("Total_E",&H_Total_E);
    helium3_tree->SetBranchAddress("Energy",&H_Energy);
    helium3_tree->SetBranchAddress("RMS",&H_RMSVec);
    helium3_tree->SetBranchAddress("weight",&H_weight);

    auto carbon_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Carbon_PowerLaw.root");
    auto carbon_tree = (TTree*)carbon_file->Get("B4");
    carbon_tree->SetBranchAddress("LayerEnergyVector",&c_EnergyVec);
    carbon_tree->SetBranchAddress("Total_E",&c_Total_E);
    carbon_tree->SetBranchAddress("Energy",&c_Energy);
    carbon_tree->SetBranchAddress("RMS",&c_RMSVec);
    carbon_tree->SetBranchAddress("weight",&c_weight);

    
    // TCut UBT = "weight*(L0_E>0.0092 && L1_E>0.0092)";
    // TCut HET = "weight*(L0_E>0.23 && L1_E >0.23 && L2_E>0.23 && L3_E>0.046)";
    // TCut HET = "weight";

    auto *h2_p = new TH2D("h2_p","h2_p",40,0,4.0,40,0,4.0); auto *h1_p = new TH2D("h1_p","h1_p",40,0,4.0,40,0,4.0);
    auto *h2_d = new TH2D("h2_d","h2_d",40,0,4.0,40,0,4.0); auto *h1_d = new TH2D("h1_d","h1_d",40,0,4.0,40,0,4.0);
    auto *h2_e = new TH2D("h2_e","h2_e",40,0,4.0,40,0,4.0); auto *h1_e = new TH2D("h1_e","h1_e",40,0,4.0,40,0,4.0);
    auto *h2_h = new TH2D("h2_h","h2_h",40,0,4.0,40,0,4.0); auto *h1_h = new TH2D("h1_h","h1_h",40,0,4.0,40,0,4.0);
    auto *h2_H = new TH2D("h2_H","h2_H",40,0,4.0,40,0,4.0); auto *h1_H = new TH2D("h1_H","h1_H",40,0,4.0,40,0,4.0);
    auto *h2_c = new TH2D("h2_c","h2_c",40,0,4.0,40,0,4.0); auto *h1_c = new TH2D("h1_c","h1_c",40,0,4.0,40,0,4.0);

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
            proton_tree->GetEntry(entry);   h2_p->Fill(log10(p_Energy),log10(p_Total_E), p_weight);
            deuteron_tree->GetEntry(entry); h2_d->Fill(log10(d_Energy),log10(d_Total_E), d_weight);
            electron_tree->GetEntry(entry); h2_e->Fill(log10(e_Energy),log10(e_Total_E), e_weight);
            helium4_tree->GetEntry(entry);  h2_h->Fill(log10(h_Energy),log10(h_Total_E), h_weight);
            helium3_tree->GetEntry(entry);  h2_H->Fill(log10(H_Energy),log10(H_Total_E), H_weight);
            carbon_tree->GetEntry(entry);   h2_c->Fill(log10(c_Energy),log10(c_Total_E), c_weight);           
    }

    auto tex = new TLatex(3.8,4.1,"log_{10} Entries");
    tex->SetTextSize(0.03);
    for (int ii=0;ii<40;ii++)
    {
        for(int jj=0;jj<40;jj++)
        {   
            h1_p->SetBinContent(ii+1,jj+1,log10(h2_p->GetBinContent(ii+1,jj+1)));
            h1_d->SetBinContent(ii+1,jj+1,log10(h2_d->GetBinContent(ii+1,jj+1)));
            h1_e->SetBinContent(ii+1,jj+1,log10(h2_e->GetBinContent(ii+1,jj+1)));
            h1_h->SetBinContent(ii+1,jj+1,log10(h2_h->GetBinContent(ii+1,jj+1)));
            h1_H->SetBinContent(ii+1,jj+1,log10(h2_H->GetBinContent(ii+1,jj+1)));
            h1_c->SetBinContent(ii+1,jj+1,log10(h2_c->GetBinContent(ii+1,jj+1)));

        }
    }

    auto c1 = new TCanvas("c1","c1",2400,1600);
    
    c1->Divide(3,2);
    c1->cd(1);
    // gPad->SetLogz();  // 设置 Z 轴为对数刻度
    gStyle->SetOptStat(0);
    h1_p->SetTitle("Proton Energy Response Matrix;log_{10}Energy Incident(GeV);log_{10}Energy Deposit(GeV)");
    h1_p->Draw("colz");
    tex->Draw();

    c1->cd(2);
    // gPad->SetLogz();  // 设置 Z 轴为对数刻度
    gStyle->SetOptStat(0);
    h1_d->SetTitle("Deuteron Energy Response Matrix;log_{10}Energy Incident(GeV);log_{10}Energy Deposit(GeV)");
    h1_d->Draw("colz");
    tex->Draw();

    c1->cd(3);
    // gPad->SetLogz();  // 设置 Z 轴为对数刻度
    gStyle->SetOptStat(0);
    h1_e->SetTitle("Electron Energy Response Matrix;log_{10}Energy Incident(GeV);log_{10}Energy Deposit(GeV)");
    h1_e->Draw("colz");
    tex->Draw();

    c1->cd(4);
    // gPad->SetLogz();  // 设置 Z 轴为对数刻度
    gStyle->SetOptStat(0);
    h1_h->SetTitle("Helium4 Energy Response Matrix;log_{10}Energy Incident(GeV);log_{10}Energy Deposit(GeV)");
    h1_h->Draw("colz");
    tex->Draw();

    c1->cd(5);
    // gPad->SetLogz();  // 设置 Z 轴为对数刻度
    gStyle->SetOptStat(0);
    h1_H->SetTitle("Helium3 Energy Response Matrix;log_{10}Energy Incident(GeV);log_{10}Energy Deposit(GeV)");
    h1_H->Draw("colz");
    tex->Draw();


    c1->cd(6);
    // gPad->SetLogz();  // 设置 Z 轴为对数刻度
    gStyle->SetOptStat(0);
    h1_c->SetTitle("Carbon Energy Response Matrix;log_{10}Energy Incident(GeV);log_{10}Energy Deposit(GeV)");
    h1_c->Draw("colz");
    tex->Draw();

    c1->SaveAs("/Users/xiongzheng/software/B4/B4e/Script/Matrix/EnergyMatrix.pdf");

}