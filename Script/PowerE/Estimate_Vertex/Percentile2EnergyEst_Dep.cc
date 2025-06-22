#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 \hat{x_ine} =  - 4.789 * sigma + x_mid

void Percentile2EnergyEst_Dep()
{
    double Ratio[18] = {0};
    double Ratio_Edge[19];
    double Energy_Dep[8] = {8,19,43,100,320,1e3,3.2e3,1e4};

    auto mean_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/DepthEst_Fitting.root");
    auto gre_p = (TGraphErrors*)mean_file->Get("gre_p");
    auto gre_d = (TGraphErrors*)mean_file->Get("gre_d");

    for (int j = 0; j < 19; j++)
    {
        if (j <= 10)       Ratio_Edge[j] = 0.005 + 0.01 * j;        // Center: 10~100 → edges: 5~105
        else               Ratio_Edge[j] = 0.15 + 0.1 * (j - 10); // Center: 2000~10000 → edges: 1500~10500    
        // cout << Ratio_Edge[j] << endl;

    }
    int nbins = 60;
    double xmin = 1e1;
    double xmax = 1e4;
    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    std::vector<double> bin_edges(nbins + 1);
    for (int j = 0; j <= nbins; j++) {
        bin_edges[j] = pow(10, logxmin + j * (logxmax - logxmin) / nbins);
        // cout << bin_edges[j] << endl;
    }

    int p_FH_Lay; int p_FH_Type; double p_energy_res;  double p_Energy;  double p_Total_E;  int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_FH_Lay; int d_FH_Type; double d_energy_res;  double d_Energy;  double d_Total_E;  int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
    int p_FI_Lay;    double p_FI_Dep;    int p_particle;
    int d_FI_Lay;    double d_FI_Dep;    int d_particle;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("Particle"         ,&p_particle);
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
    proton_tree->SetBranchAddress("BarEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_FH_Lay);
    proton_tree->SetBranchAddress("First_Had_Type"   ,&p_FH_Type);
    proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
    proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);
    proton_tree->SetBranchAddress("energy_res"     ,&p_energy_res);
    proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
    proton_tree->SetBranchAddress("Energy"         , &p_Energy);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);
    proton_tree->SetBranchAddress("weight"          ,&p_weight);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("Particle"         ,&d_particle);
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
    deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_FH_Lay);
    deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_FH_Type);
    deuteron_tree->SetBranchAddress("First_Ine_Depth", &d_FI_Dep);
    deuteron_tree->SetBranchAddress("First_Ine_Layer", &d_FI_Lay);
    deuteron_tree->SetBranchAddress("energy_res"     ,&d_energy_res);
    deuteron_tree->SetBranchAddress("Total_E"         ,&d_Total_E);
    deuteron_tree->SetBranchAddress("Energy"         , &d_Energy);
    deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);


    // Deposit
    auto h_DepE_p = new TH2D("h_DepE_p","h_DepE_p",nbins,bin_edges.data(),7,Energy_Dep);
    auto h_DepE_d = new TH2D("h_DepE_d","h_DepE_d",nbins,bin_edges.data(),7,Energy_Dep);
    auto hist_p_2 = new TH1D("hist_p_2","hist_p_2",7,Energy_Dep);
    auto hist_d_2 = new TH1D("hist_d_2","hist_d_2",7,Energy_Dep);
    auto Energy_Mid_p = new TH1D("Energy_Mid_p","Energy_Mid_p",7,Energy_Dep);
    auto Energy_Mid_d = new TH1D("Energy_Mid_d","Energy_Mid_d",7,Energy_Dep);

    TH1D *E_br_p[7]; TF1 *langaus_p[7];
    TH1D *E_br_d[7]; TF1 *langaus_d[7];

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    // for (Long64_t entry = 0; entry < 10000; entry++)
    {        
        proton_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Proton : " << entry << endl;
        int p_energy_index = int(floor((log10(p_energy_res) ) / 0.2));
        if (p_energy_index < 0 || p_energy_index > 20) continue;
        if (p_FI_Dep < 0) {p_FI_Lay = 14;}
        if (p_Nhits < 10 ) continue;
        h_DepE_p->Fill(p_Energy,p_Total_E,p_weight);
    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries() ; ++entry)
    // for (Long64_t entry = 0; entry < 10000 ; ++entry)
    {
        deuteron_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Deuteron : " << entry << endl;
        int d_energy_index = int(floor((log10(d_energy_res) ) / 0.2));
        if (d_energy_index < 0 || d_energy_index > 20) continue;
        if (d_FI_Dep < 0) {d_FI_Lay = 14;}
        if (d_Nhits < 10 ) continue;
        h_DepE_d->Fill(d_Energy,d_Total_E,d_weight);
    }

    auto c0 = new TCanvas("c0", "c0", 2500, 2500);
    c0->Divide(2, 2); // 上面两格放 colz 图，下面两格放投影图
    gStyle->SetOptStat(0);
    // 画 Deposit vs True Energy（colz）
    c0->cd(1);
    gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz();
    h_DepE_p->SetTitle("Proton;True Energy (GeV);Deposit Energy (GeV)");
    h_DepE_p->Draw("colz");
    
    c0->cd(2);
    gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz();
    h_DepE_d->SetTitle("Deuteron;True Energy (GeV);Deposit Energy (GeV)");
    h_DepE_d->Draw("colz");
    
    // Legend
    auto lg0 = new TLegend(0.12, 0.7, 0.88, 0.88);
    lg0->SetNColumns(2);
    
    // 投影图
    for (int j = 0; j < h_DepE_p->GetNbinsY(); j++)
    // for (int j = 0; j < 2; j++) 
    {
        E_br_p[j] = (TH1D*)h_DepE_p->ProjectionX(Form("E_br_p_%d", j), j+1, j+1);
        E_br_d[j] = (TH1D*)h_DepE_d->ProjectionX(Form("E_br_d_%d", j), j+1, j+1);
        E_br_p[j]->SetLineColor(j+1); E_br_p[j]->SetLineWidth(2); E_br_p[j]->Scale(1/E_br_p[j]->Integral());
        E_br_d[j]->SetLineColor(j+1); E_br_d[j]->SetLineWidth(2); E_br_d[j]->Scale(1/E_br_d[j]->Integral());
    
        lg0->AddEntry(E_br_p[j], Form("Dep %.1f-%.1f GeV", h_DepE_p->GetYaxis()->GetBinLowEdge(j+1), h_DepE_p->GetYaxis()->GetBinUpEdge(j+1)), "l");
    
        if (j == 0) 
        {
            c0->cd(3); E_br_p[0]->SetTitle("Proton: True Energy Distribution;True Energy (GeV);Count");   gPad->SetLogx(); E_br_p[0]->GetYaxis()->SetRangeUser(0,0.5);  E_br_p[0]->Draw("hist");
            c0->cd(4); E_br_d[0]->SetTitle("Deuteron: True Energy Distribution;True Energy (GeV);Count"); gPad->SetLogx(); E_br_d[0]->GetYaxis()->SetRangeUser(0,0.5);  E_br_d[0]->Draw("hist");
        } 
        else 
        {
            c0->cd(3); E_br_p[j]->Draw("histsame");
            c0->cd(4); E_br_d[j]->Draw("histsame");
        }


        c0->cd(3); 
        langaus_p[j] = new TF1(Form("langaus_p[%d]",j) ,langaufun, Energy_Dep[j]  , 7*Energy_Dep[j+1], 4); 
        langaus_p[j]->SetLineColor(j+1); 
        langaus_p[j]->SetParameters(E_br_p[j]->GetStdDev(),E_br_p[j]->GetMean(),E_br_p[j]->Integral(),0.01);
        langaus_p[j]->SetParLimits(0,0,E_br_p[j]->GetStdDev()*3);
        langaus_p[j]->SetParLimits(1,0,E_br_p[j]->GetMean()*3);
        langaus_p[j]->SetParLimits(2,0,1e4);
        langaus_p[j]->SetParLimits(3,0,1e1);
        if(j==0) E_br_p[j]->Draw("hist");
        else E_br_p[j]->Draw("histsame");
        E_br_p[j]->Fit(langaus_p[j],"R");        
        langaus_p[j]->Draw("lsame");
        Energy_Mid_p->SetBinContent(j+1,langaus_p[j]->GetParameter(1)); //  Get True Energy Peak
        Energy_Mid_p->SetBinError(j+1,  langaus_p[j]->GetParameter(0)); //  Get True Energy Width
        mean_file->cd();
        hist_p_2->SetBinContent(j+1, gre_p->Eval(langaus_p[j]->GetParameter(1)));

        c0->cd(4); 
        langaus_d[j] = new TF1(Form("langaus_d[%d]",j) , langaufun, Energy_Dep[j]  , 7*Energy_Dep[j+1],4);
        langaus_d[j]->SetLineColor(j+1); 
        langaus_d[j]->SetParameters(E_br_d[j]->GetStdDev(),E_br_d[j]->GetMean(),E_br_d[j]->Integral(),0.01);
        langaus_d[j]->SetParLimits(0,0,E_br_d[j]->GetStdDev()*3);
        langaus_d[j]->SetParLimits(1,0,E_br_d[j]->GetMean()*3);
        langaus_d[j]->SetParLimits(2,0,1e4);
        langaus_d[j]->SetParLimits(3,0,1e1);
        if(j==0) E_br_d[j]->Draw("hist");
        else E_br_d[j]->Draw("histsame");
        E_br_d[j]->Fit(langaus_d[j],"R");        
        langaus_d[j]->Draw("lsame");
        Energy_Mid_d->SetBinContent(j+1, langaus_d[j]->GetParameter(1)); //  Get True Energy Peak
        Energy_Mid_d->SetBinError(j+1,   langaus_d[j]->GetParameter(0)); //  Get True Energy Width
        mean_file->cd();
        hist_d_2->SetBinContent(j+1, gre_d->Eval(langaus_d[j]->GetParameter(1)));
    }
    
    c0->cd(3); lg0->Draw();
    c0->cd(4); lg0->Draw();
    hist_p_2->SetTitle("Proton;Deposit Energy (GeV); BEST FIT MEAN");
    hist_d_2->SetTitle("Deuteron;Deposit Energy (GeV); BEST FIT MEAN");
    Energy_Mid_p->SetTitle("Proton;Deposit Energy (GeV); True Energy(GeV)");
    Energy_Mid_d->SetTitle("Deuteron;Deposit Energy (GeV); True Energy(GeV)");

    auto write_file = new TFile("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Estimate_Vertex/EnergyEst_Dep.root", "RECREATE");
    write_file->cd();
    c0->Write();

    hist_p_2->Write();
    hist_d_2->Write();
    Energy_Mid_p->Write();
    Energy_Mid_d->Write();
    write_file->Close();
    mean_file->Close();
}