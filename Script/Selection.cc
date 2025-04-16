void Selection()
{
    int p_First_Had_Layer;  double p_Total_E;  double p_Energy; double p_weight; int p_First_Had_Type;
    int d_First_Had_Layer;  double d_Total_E;  double d_Energy; double d_weight; int d_First_Had_Type;
    std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr; 
    std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr; 

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("First_Had_Type"   ,&d_First_Had_Type);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);
    proton_tree->SetBranchAddress("weight"           ,&p_weight);
    proton_tree->SetBranchAddress("Energy"           ,&p_Energy);


    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_First_Had_Type);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    deuteron_tree->SetBranchAddress("Energy"           ,&d_Energy);

    double Energy[15] = {0};
    double Energy_LL[15] = {0};
    double Energy_UL[15] = {0};

    TH1I *h1_p[15]; TH1I *h2_p[15]; TH1I *h3_p[15];
    TH1I *h1_d[15]; TH1I *h2_d[15]; TH1I *h3_d[15];

    for (int i = 0; i < 15; i++)
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;
        h1_p[i] = new TH1I(Form("h1_p[%d]",i),Form("h1_p[%d]",i),7,0,7);
        h2_p[i] = new TH1I(Form("h2_p[%d]",i),Form("h2_p[%d]",i),7,0,7); 
        h3_p[i] = new TH1I(Form("h3_p[%d]",i),Form("h3_p[%d]",i),7,0,7);
        h1_d[i] = new TH1I(Form("h1_d[%d]",i),Form("h1_d[%d]",i),7,0,7); 
        h2_d[i] = new TH1I(Form("h2_d[%d]",i),Form("h2_d[%d]",i),7,0,7); 
        h3_d[i] = new TH1I(Form("h3_d[%d]",i),Form("h3_d[%d]",i),7,0,7);
    }
    
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);   
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        h1_p[p_energy_index]->Fill(p_First_Had_Layer); 
        if ((*p_EnergyVec)[0] > 0.23 && (*p_EnergyVec)[1] > 0.23 && (*p_EnergyVec)[2]  > 0.23 && (*p_EnergyVec)[3] > 0.046) h2_p[p_energy_index]->Fill(p_First_Had_Layer); 
        if ((*p_EnergyVec)[0] > 0.23 && (*p_EnergyVec)[1] > 0.23 && (*p_RMSVec)[2]<40 && (*p_RMSVec)[3]<40)                               h3_p[p_energy_index]->Fill(p_First_Had_Layer);
    }


    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry); 
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 14) continue;
        h1_d[d_energy_index]->Fill(d_First_Had_Layer); 
        if ((*d_EnergyVec)[0] > 0.23 && (*d_EnergyVec)[1] > 0.23 && (*d_EnergyVec)[2] > 0.23 && (*d_EnergyVec)[3] > 0.046) h2_d[d_energy_index]->Fill(d_First_Had_Layer); 
        if ((*d_EnergyVec)[0] > 0.23 && (*d_EnergyVec)[1] > 0.23 && (*d_RMSVec)[2]<40 && (*d_RMSVec)[3]<40)                                h3_d[d_energy_index]->Fill(d_First_Had_Layer);
    }

    auto c1 = new TCanvas("c1","c1",2500,1500);
    c1->Clear();
    c1->Divide(5,3);
    gStyle->SetOptStat(0);

    auto c2 = new TCanvas("c2","c2",2500,1500);
    c2->Clear();
    c2->Divide(5,3);
    gStyle->SetOptStat(0);
    for (int i = 0; i < 14; i++)
    {
        c1->cd(i+1);
        h1_p[i]->SetMinimum(0);
        h1_p[i]->SetTitle(Form("Proton Deposit Energy[%.2fGeV, %.2fGeV]; Hadronic Interaction Layer;Events",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        h1_p[i]->SetLineColor(kBlack); h1_p[i]->SetLineWidth(2); h1_p[i]->Draw("histE1");
        h2_p[i]->SetLineColor(kRed  ); h2_p[i]->SetLineWidth(2); h2_p[i]->Draw("histsame");
        h3_p[i]->SetLineColor(kBlue ); h3_p[i]->SetLineWidth(2); h3_p[i]->Draw("histsame");


        c2->cd(i+1);
        h1_d[i]->SetMinimum(0);
        h1_d[i]->SetTitle(Form("Deuteron Deposit Energy[%.2fGeV, %.2fGeV]; Hadronic Interaction Layer;Events",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        h1_d[i]->SetLineColor(kBlack); h1_d[i]->SetLineWidth(2); h1_d[i]->Draw("histE1");
        h2_d[i]->SetLineColor(kRed  ); h2_d[i]->SetLineWidth(2); h2_d[i]->Draw("histsame");
        h3_d[i]->SetLineColor(kBlue ); h3_d[i]->SetLineWidth(2); h3_d[i]->Draw("histsame");
    }

    c1->cd(15);
    auto legend1 = new TLegend(0.12, 0.62, 0.88, 0.88);
    legend1->AddEntry(h1_p[0], "No Trigger", "l");
    legend1->AddEntry(h2_p[0], "HET", "l");
    legend1->AddEntry(h3_p[0], "BGO_{0}>10 MIP && BGO_{1}>10 MIP && RMS_{2}<45 && RMS_{3} < 45", "l");
    legend1->Draw();       
    c2->cd(15);
    legend1->Draw();       

}