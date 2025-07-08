#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 (x_ine - x_mid)

void Percentile2Energy_Rec()
{
    int p_FH_Lay; int p_FH_Type; double p_Total_E;      int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; 
    int d_FH_Lay; int d_FH_Type; double d_Total_E;      int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; 
    int p_FI_Lay;    double p_FI_Dep;    double p_energy;  double p_energy_res;     int p_particle;  double p_weight;
    int d_FI_Lay;    double d_FI_Dep;    double d_energy;  double d_energy_res;     int d_particle;  double d_weight;
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_200GeV.root");
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
    proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
    proton_tree->SetBranchAddress("Energy"         , &p_energy);
    proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_200GeV.root");
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
    deuteron_tree->SetBranchAddress("Total_E"         ,&d_Total_E);
    deuteron_tree->SetBranchAddress("Energy"         , &d_energy);
    deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);
    double Energy[15]={0};
    double Energy_LL[15]={0};      
    double Energy_UL[15]={0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    //
    auto h2_p = new TH2D("h2_p","Proton;SSP-exit distance/BGO L;Energy Deposit Ratio",28,0,1,25,0,1);
    auto h2_d = new TH2D("h2_d","Deuteron;SSP-exit distance/BGO L;Energy Deposit Ratio",28,0,1,25,0,1); 

    auto h1_p = new TH2D("h1_p","Proton;End Fraction;Energy Deposit Ratio",28,0,1,25,0,1);
    auto h1_d = new TH2D("h1_d","Deuteron;End Fraction;Energy Deposit Ratio",28,0,1,25,0,1); 

    for(int i =0 ; i<15 ; i++)  // Deposit Energy Bin
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;
    }

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); entry++)
    {        
        proton_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Proton : " << entry << endl;
        // int p_energy_index = int(floor((log10(p_Total_E) ) / 0.2));
        // if(p_energy_index < 0 || p_energy_index > 20) continue;
        if(p_FI_Dep < 0) continue;
        // if (p_Nhits < 10 ) continue;
        // if(p_Total_E / p_energy < 0.2 ) continue;
        h2_p->Fill((355-p_FI_Dep)/355, p_Total_E / p_energy);
        h1_p->Fill((*p_Efrac)[13], p_Total_E / p_energy);

    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry);
        if (entry%10000==0) cout << " Deuteron : " << entry << endl;
        // int d_energy_index = int(floor((log10(d_Total_E) ) / 0.2));
        // if(d_energy_index < 0 || d_energy_index > 20) continue;
        if(d_FI_Dep < 0) continue;
        // if (d_Nhits < 10 ) continue;
        // if(d_Total_E / d_energy < 0.2 ) continue;
        h2_d->Fill((355-d_FI_Dep)/355, d_Total_E / d_energy);
        h1_d->Fill((*d_Efrac)[13], d_Total_E / d_energy);
    }

    
    auto c3 = new TCanvas("c3","c3",2000,2000);
    c3->cd();
    gStyle->SetOptFit(1111);
    c3->Divide(2,2);

    c3->cd(1);
    gPad->SetLogz();
    h2_p->Draw("colz");
    cout << h2_p->GetCorrelationFactor() << endl;
    c3->cd(2);
    gPad->SetLogz();
    h2_d->Draw("colz");
    cout << h2_d->GetCorrelationFactor() << endl;

    c3->cd(3);
    gPad->SetLogz();
    h1_p->Draw("colz");
    cout << h1_p->GetCorrelationFactor() << endl;

    c3->cd(4);
    gPad->SetLogz();
    h1_d->Draw("colz");
    cout << h1_d->GetCorrelationFactor() << endl;
    
}