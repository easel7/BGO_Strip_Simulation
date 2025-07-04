void EnergyMatrix3()
{
    int p_FH_Lay; int p_FH_Type; double p_Total_E;      int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_FH_Lay; int d_FH_Type; double d_Total_E;      int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
    int p_FI_Lay;    double p_FI_Dep;    int p_particle; double p_Energy;
    int d_FI_Lay;    double d_FI_Dep;    int d_particle; double d_Energy;

    int ii;
    double Energy[30]             = {0.};
    double P_Bias[30]             = {0.};    
    double P_Reso[30]             = {0.};    
    double D_Bias[30]             = {0.};      
    double D_Reso[30]             = {0.};    


    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("Particle"         ,&p_particle);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type"   ,&p_FH_Type);
    proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
    proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);
    proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
    proton_tree->SetBranchAddress("weight"          ,&p_weight);
    proton_tree->SetBranchAddress("Energy",&p_Energy);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);


    proton_tree->SetBranchAddress("RMS",&p_RMSVec);
    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("Particle"         ,&d_particle);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_FH_Lay);
    deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_FH_Type);
    deuteron_tree->SetBranchAddress("First_Ine_Depth", &d_FI_Dep);
    deuteron_tree->SetBranchAddress("First_Ine_Layer", &d_FI_Lay);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);
    deuteron_tree->SetBranchAddress("Energy"         , &d_Energy);
    deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);
    
    auto *h2_p = new TH2D("h2_p","h2_p",30,1,4.0,30,1,4.0); auto *h1_p = new TH2D("h1_p","h1_p",30,1,4.0,30,1,4.0);
    auto *h2_d = new TH2D("h2_d","h2_d",30,1,4.0,30,1,4.0); auto *h1_d = new TH2D("h1_d","h1_d",30,1,4.0,30,1,4.0);

    auto *hR_p = new TH2D("hR_p","hR_p",30,1,4.0,30,-1.5,1.5);  auto *hr_p = new TH2D("hr_p","hr_p",30,1,4.0,30,-1.5,1.5); 
    auto *hR_d = new TH2D("hR_d","hR_d",30,1,4.0,30,-1.5,1.5);  auto *hr_d = new TH2D("hr_d","hr_d",30,1,4.0,30,-1.5,1.5);

    auto *prof_p = new TProfile("prof_p","prof_p",30,1,4,1,4,"s");
    auto *prof_d = new TProfile("prof_d","prof_d",30,1,4,1,4,"s");


    TH1D *h0_p[30];
    TH1D *h0_d[30];
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry); 
        if (entry%100000==0) cout << " Proton : " << entry << endl;
        int p_energy_index = int(floor((log10(p_Total_E) ) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 20) continue;
        if(p_FI_Dep < 0) continue;
        if (p_Nhits < 10 ) continue;
        h2_p->Fill(  log10(p_Total_E), log10(p_Energy), p_weight); 
        // hR_p->Fill( log10(p_Energy) , (p_Total_E-p_Energy)/p_Energy ,p_weight);
        prof_p->Fill( log10(p_Total_E), log10(p_Energy), p_weight);
    }
    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry); 
        if (entry%100000==0) cout << " Deuteron : " << entry << endl;
        int d_energy_index = int(floor((log10(d_Total_E) ) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 20) continue;
        if(d_FI_Dep < 0) continue;
        if (d_Nhits < 10 ) continue;
        h2_d->Fill( log10(d_Total_E), log10(d_Energy), d_weight); 
        // hR_d->Fill( log10(d_Energy) , (d_Total_E-d_Energy)/d_Energy ,d_weight);
        prof_d->Fill(log10(d_Total_E), log10(d_Energy), d_weight);
    }
    auto tex = new TLatex(3.8,4.1,"log_{10} Entries");
    tex->SetTextSize(0.03);

    auto c2 = new TCanvas("c2","c2",2500,1200); c2->Divide(8,5);
    auto c3 = new TCanvas("c3","c3",2500,1200); c3->Divide(8,5);


    TF1 *fit_d = new TF1("fit_d", "pol1", 1, 3.5);
    prof_d->Fit(fit_d, "R");
    fit_d->SetLineColor(kBlack);
    fit_d->SetLineWidth(2);
    TF1 *fit_p = new TF1("fit_p", "pol1", 1, 3.5);
    prof_p->Fit(fit_p, "R");
    fit_p->SetLineColor(kBlack);
    fit_p->SetLineWidth(2);
    std::cout << "Proton fit (pol1): " << fit_p->GetParameter(0) << ", "  << fit_p->GetParameter(1) <<  std::endl;
    std::cout << "Deuteron fit (pol1): " << fit_d->GetParameter(0) << ", " << fit_d->GetParameter(1) <<  std::endl;
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry); 
        if (entry%100000==0) cout << " Proton : " << entry << endl;
        int p_energy_index = int(floor((log10(p_Total_E) ) / 0.2));
        double energy_rec_p =  pow(10, log10(p_Total_E) * fit_p->GetParameter(1) + fit_p->GetParameter(0)) ;
        if(p_energy_index < 0 || p_energy_index > 20) continue;
        if(p_FI_Dep < 0) continue;
        if (p_Nhits < 10 ) continue;
        hR_p->Fill( log10(p_Energy) , (energy_rec_p-p_Energy)/p_Energy ,p_weight);
    }
    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry); 
        if (entry%100000==0) cout << " Deuteron : " << entry << endl;
        int d_energy_index = int(floor((log10(d_Total_E) ) / 0.2));
        double energy_rec_d =  pow(10, log10(d_Total_E) * fit_d->GetParameter(1) + fit_d->GetParameter(0)) ;
        if(d_energy_index < 0 || d_energy_index > 20) continue;
        if(d_FI_Dep < 0) continue;
        if (d_Nhits < 10 ) continue;
        hR_d->Fill( log10(d_Energy) , (energy_rec_d-d_Energy)/d_Energy ,d_weight);
    }

    for (int ii=0;ii<30;ii++)
    {
        for(int jj=0;jj<30;jj++)
        {   
            h1_p->SetBinContent(ii+1,jj+1,log10(h2_p->GetBinContent(ii+1,jj+1)));
            h1_d->SetBinContent(ii+1,jj+1,log10(h2_d->GetBinContent(ii+1,jj+1)));
            hr_p->SetBinContent(ii+1,jj+1,log10(hR_p->GetBinContent(ii+1,jj+1)));
            hr_d->SetBinContent(ii+1,jj+1,log10(hR_d->GetBinContent(ii+1,jj+1)));
        }

        Energy[ii]     = 1.05 + 0.1 * ii;
        h0_p[ii] = hR_p->ProjectionY(Form("h0_p[%d]",ii),ii,ii+1,"");
        c2->cd(ii+1);
        h0_p[ii]->Scale(1./h0_p[ii]->Integral());
        h0_p[ii]->SetTitle(Form("Proton [%.1f TeV, %.1f TeV];#frac{E_{rec}-E_{true}}{E_{true}};Fraction",pow(10,0.1 * ii+0),pow(10,0.1*ii+0.1)));
        h0_p[ii]->GetYaxis()->CenterTitle();
        h0_p[ii]->GetYaxis()->SetTitleOffset(1.1);
        h0_p[ii]->GetYaxis()->SetRangeUser(0,h0_p[ii]->GetMaximum()*1.2);
        h0_p[ii]->GetYaxis()->SetLabelSize(0.03);
        h0_p[ii]->GetXaxis()->CenterTitle();
        h0_p[ii]->GetXaxis()->SetLabelSize(0.03);
        h0_p[ii]->SetLineColor(kRed);
        h0_p[ii]->SetLineWidth(2);
        h0_p[ii]->Draw("hist"); 
        P_Bias[ii]      = h0_p[ii]->GetMean();
        P_Reso[ii]      = h0_p[ii]->GetStdDev();

        h0_d[ii] = hR_d->ProjectionY(Form("h0_d[%d]",ii),ii,ii+1,"");
        c3->cd(ii+1);
        h0_d[ii]->Scale(1./h0_d[ii]->Integral());
        h0_d[ii]->SetTitle(Form("Dueteron [%.1f TeV, %.1f TeV];#frac{E_{rec}-E_{true}}{E_{true}};Fraction",pow(10,0.1 * ii+0),pow(10,0.1*ii+0.1)));
        h0_d[ii]->GetYaxis()->CenterTitle();
        h0_d[ii]->GetYaxis()->SetTitleOffset(1.1);
        h0_d[ii]->GetYaxis()->SetRangeUser(0,h0_d[ii]->GetMaximum()*1.2);
        h0_d[ii]->GetYaxis()->SetLabelSize(0.03);
        h0_d[ii]->GetXaxis()->CenterTitle();
        h0_d[ii]->GetXaxis()->SetLabelSize(0.03);
        h0_d[ii]->SetLineColor(kBlue);
        h0_d[ii]->SetLineWidth(2);
        h0_d[ii]->Draw("hist"); 
        D_Bias[ii]      = h0_d[ii]->GetMean();
        D_Reso[ii]      = h0_d[ii]->GetStdDev();
    }
    auto P_bias_gre = new TGraph(30,Energy,P_Bias);
    auto P_reso_gre = new TGraph(30,Energy,P_Reso);
    auto D_bias_gre = new TGraph(30,Energy,D_Bias);
    auto D_reso_gre = new TGraph(30,Energy,D_Reso);

    auto c1 = new TCanvas("c1","c1",1200,1200);
    c1->Divide(3,2);
    c1->cd(1);
    gStyle->SetOptStat(0);
    h1_p->SetTitle("Proton Energy Response Matrix;log_{10}Energy Deposit(GeV);log_{10}Energy Incident(GeV)");
    h1_p->Draw("colz");
    prof_p->Draw("same");
    tex->Draw();

    // Fit the proton profile with a quadratic polynomial (pol2)
   
    fit_p->Draw("same");


    c1->cd(2);
    gStyle->SetOptStat(0);
    hr_p->SetTitle(";log_{10}Energy Incident(GeV);#frac{E_{rec}-E_{true}}{E_{true}}");
    hr_p->Draw("colz");
    tex->Draw();

    c1->cd(3);
    P_bias_gre->SetMarkerColor(kRed);
    D_bias_gre->SetMarkerColor(kBlue);
    P_bias_gre->SetMarkerStyle(20);
    D_bias_gre->SetMarkerStyle(21);
    P_bias_gre->SetTitle(";log_{10}Energy Incident(GeV);Bias");
    P_bias_gre->Draw("AP");
    D_bias_gre->Draw("PSAME");

    c1->cd(4);
    gStyle->SetOptStat(0);
    h1_d->SetTitle("Deuteron Energy Response Matrix;log_{10}Energy Deposit(GeV);log_{10}Energy Incident(GeV)");
    h1_d->Draw("colz");
    prof_d->Draw("same");
    tex->Draw();

    // Fit the deuteron profile with a quadratic polynomial (pol2)

    fit_d->Draw("same");


    c1->cd(5);
    gStyle->SetOptStat(0);
    hr_d->SetTitle(";log_{10}Energy Incident(GeV);#frac{E_{rec}-E_{true}}{E_{true}}");
    hr_d->Draw("colz");
    tex->Draw();

    c1->cd(6);
    P_reso_gre->SetMarkerColor(kRed);
    D_reso_gre->SetMarkerColor(kBlue);
    P_reso_gre->SetMarkerStyle(20);
    D_reso_gre->SetMarkerStyle(21);
    P_reso_gre->SetTitle(";log_{10}Energy Incident(GeV);Resolution");
    P_reso_gre->Draw("AP");
    D_reso_gre->Draw("PSAME");

    auto write_file = new TFile("/Users/xiongzheng/software/B4/B4e/Weight/RecE.root", "RECREATE");
    write_file->cd();
    c1->Write();
    fit_p->Write();
    fit_d->Write();
    write_file->Close();
}