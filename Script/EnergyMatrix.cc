void EnergyMatrix()
{
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    auto electron_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Electron_PowerLaw.root");
    auto electron_tree = (TTree*)electron_file->Get("B4");
    auto helium4_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Helium4_PowerLaw.root");
    auto helium4_tree = (TTree*)helium4_file->Get("B4");
    auto helium3_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Helium3_PowerLaw.root");
    auto helium3_tree = (TTree*)helium3_file->Get("B4");
    auto carbon_file = TFile::Open("/Users/xiongzheng/software/B4/B4c/Weight/Carbon_PowerLaw.root");
    auto carbon_tree = (TTree*)carbon_file->Get("B4");
    
    // TCut UBT = "weight*(L0_E>0.0092 && L1_E>0.0092)";
    TCut HET = "weight*(L0_E>0.23 && L1_E >0.23 && L2_E>0.23 && L3_E>0.046)";
    // TCut HET = "weight";

    auto *h2_p = new TH2D("h2_p","h2_p",36,0,3.6,36,0,3.6); auto *h1_p = new TH2D("h1_p","h1_p",36,0,3.6,36,0,3.6);
    auto *h2_d = new TH2D("h2_d","h2_d",36,0,3.6,36,0,3.6); auto *h1_d = new TH2D("h1_d","h1_d",36,0,3.6,36,0,3.6);
    auto *h2_e = new TH2D("h2_e","h2_e",36,0,3.6,36,0,3.6); auto *h1_e = new TH2D("h1_e","h1_e",36,0,3.6,36,0,3.6);
    auto *h2_h = new TH2D("h2_h","h2_h",36,0,3.6,36,0,3.6); auto *h1_h = new TH2D("h1_h","h1_h",36,0,3.6,36,0,3.6);
    auto *h2_H = new TH2D("h2_H","h2_H",36,0,3.6,36,0,3.6); auto *h1_H = new TH2D("h1_H","h1_H",36,0,3.6,36,0,3.6);
    auto *h2_c = new TH2D("h2_c","h2_c",36,0,3.6,36,0,3.6); auto *h1_c = new TH2D("h1_c","h1_c",36,0,3.6,36,0,3.6);

    proton_tree->Draw("log10(Total_E):log10(Energy)>>h2_p",HET,"");   
    deuteron_tree->Draw("log10(Total_E):log10(Energy)>>h2_d",HET,""); 
    electron_tree->Draw("log10(Total_E):log10(Energy)>>h2_e",HET,""); 
    helium4_tree->Draw("log10(Total_E):log10(Energy)>>h2_h",HET,"");  
    helium3_tree->Draw("log10(Total_E):log10(Energy)>>h2_H",HET,"");  
    carbon_tree ->Draw("log10(Total_E):log10(Energy)>>h2_c",HET,"");  

    auto tex = new TLatex(3.4,3.7,"log_{10} Entries");
    tex->SetTextSize(0.03);
    for (int ii=0;ii<36;ii++)
    {
        for(int jj=0;jj<36;jj++)
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

    c1->SaveAs("/Users/xiongzheng/software/B4/B4c/Script/EnergyMatrix.pdf");

}