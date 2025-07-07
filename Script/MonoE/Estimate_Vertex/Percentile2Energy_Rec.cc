#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti_hist.hh"
// 比较拟合差距 (x_ine - x_mid)

void Percentile2Energy_Rec()
{
    int Energy_Name[28]={0};
    int Energy_Name_Err[28]={0};
    auto gre_p = new TGraphErrors();
    auto gre_d = new TGraphErrors();
    if (gDirectory->FindObject("hist_p")) delete gDirectory->FindObject("hist_p");
    if (gDirectory->FindObject("hist_d")) delete gDirectory->FindObject("hist_d");
    auto hist_p = new TH1D("hist_p","hist_p",28,0,28);
    auto hist_d = new TH1D("hist_d","hist_d",28,0,28);

    // for (int k =0; k < 1; k++)
    // for (int k =0; k < 28; k++)
    for (int k =18; k < 19; k++)
    {
        if (k < 10)      {Energy_Name[k] = (k + 1) * 10;        }       // 10 ~ 100
        else if (k < 19) {Energy_Name[k] = (k - 9 + 1) * 100;   }        // 200 ~ 1000
        else             {Energy_Name[k] = (k - 18 + 1) * 1000; }        // 2000 ~ 10000
        int p_FH_Lay; int p_FH_Type; double p_Total_E;      int p_Nhits;std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_L_EnergyVec = nullptr;   std::vector<double>* p_EnergyVec = nullptr;   std::vector<double>* p_Efrac = nullptr; double p_weight;
        int d_FH_Lay; int d_FH_Type; double d_Total_E;      int d_Nhits;std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_L_EnergyVec = nullptr;   std::vector<double>* d_EnergyVec = nullptr;   std::vector<double>* d_Efrac = nullptr; double d_weight;
        int p_FI_Lay;    double p_FI_Dep;    double p_energy;
        int d_FI_Lay;    double d_FI_Dep;    double d_energy;
        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",Energy_Name[k]));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("Energy"         ,&p_energy);
        proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
        proton_tree->SetBranchAddress("BarEnergyVector",&p_EnergyVec);
        proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
        proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_FH_Lay);
        proton_tree->SetBranchAddress("First_Had_Type"   ,&p_FH_Type);
        proton_tree->SetBranchAddress("First_Ine_Depth", &p_FI_Dep);
        proton_tree->SetBranchAddress("First_Ine_Layer", &p_FI_Lay);
        proton_tree->SetBranchAddress("Total_E"         ,&p_Total_E);
        // proton_tree->SetBranchAddress("weight"          ,&p_weight);
        proton_tree->SetBranchAddress("Nhits"          , &p_Nhits);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",Energy_Name[k]));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("Energy"         ,&d_energy);
        deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
        deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
        deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
        deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_FH_Lay);
        deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_FH_Type);
        deuteron_tree->SetBranchAddress("First_Ine_Depth", &d_FI_Dep);
        deuteron_tree->SetBranchAddress("First_Ine_Layer", &d_FI_Lay);
        deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
        // deuteron_tree->SetBranchAddress("weight"           ,&d_weight);
        deuteron_tree->SetBranchAddress("Nhits"          , &d_Nhits);

        double Energy[15]={0};
        double Energy_LL[15]={0};      
        double Energy_UL[15]={0};

        double Layer[14]={0};
        double Layer_Err[14]={0};

        //
        auto h2_p = new TH2D("h2_p",";SSP-exit distance/BGO L;Energy Deposit Ratio",28,0,1,25,0,1);
        auto h2_d = new TH2D("h2_d",";SSP-exit distance/BGO L;Energy Deposit Ratio",28,0,1,25,0,1); 

        auto h1_p = new TH2D("h1_p",";End Fraction;Energy Deposit Ratio",28,0,1,25,0,1);
        auto h1_d = new TH2D("h1_d",";End Fraction;Energy Deposit Ratio",28,0,1,25,0,1); 

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
            int p_energy_index = int(floor((log10(p_Total_E) ) / 0.2));
            if(p_energy_index < 0 || p_energy_index > 20) continue;
            if(p_FI_Dep < 0) continue;
            if (p_Nhits < 10 ) continue;
            h2_p->Fill((355-p_FI_Dep)/355, p_Total_E / p_energy);
            h1_p->Fill((*p_Efrac)[13], p_Total_E / p_energy);

        }

        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
        {
            deuteron_tree->GetEntry(entry);
            if (entry%10000==0) cout << " Deuteron : " << entry << endl;
            int d_energy_index = int(floor((log10(d_Total_E) ) / 0.2));
            if(d_energy_index < 0 || d_energy_index > 20) continue;
            if(d_FI_Dep < 0) continue;
            if (d_Nhits < 10 ) continue;
            h2_d->Fill((355-d_FI_Dep)/355, d_Total_E / d_energy);
            h1_d->Fill((*d_Efrac)[13], d_Total_E / d_energy);
        }

      
        auto c3 = new TCanvas("c3","c3",2000,2000);
        c3->cd();
        gStyle->SetOptFit(1111);
        c3->Divide(2,2);

        c3->cd(1);
        h2_p->Draw("colz");
        c3->cd(2);
        h2_d->Draw("colz");
        c3->cd(3);
        h1_p->Draw("colz");
        c3->cd(4);
        h1_d->Draw("colz");
    }
}