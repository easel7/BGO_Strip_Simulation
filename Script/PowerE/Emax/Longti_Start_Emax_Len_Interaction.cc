#include "/Users/xiongzheng/software/B4/B4e/Script/Ulti.hh"

void Longti_Start_Emax_Len_Interaction()
{
    int p_First_Had_Layer; int p_First_Had_Type; double p_Total_E;  std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr; double p_weight;
    int d_First_Had_Layer; int d_First_Had_Type; double d_Total_E;  std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr; double d_weight;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("First_Had_Type"  ,&p_First_Had_Type);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);
    proton_tree->SetBranchAddress("weight"           ,&p_weight);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_First_Had_Type);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    deuteron_tree->SetBranchAddress("weight"           ,&d_weight);
    
    double Energy[15]={0};
    double Energy_LL[15]={0};      
    double Energy_UL[15]={0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[15][14];
    TH1D *h1_d[15][14];

    TH1D *h1_p_inter[15];
    TH1D *h1_d_inter[15];


    for(int i =0 ; i<15 ; i++)  // Deposit Energy Bin
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;

        h1_p_inter[i] =new TH1D(Form("h1_p_inter[%d]",i),Form("h1_p_inter[%d]",i), 14,0,14);  
        h1_d_inter[i] =new TH1D(Form("h1_d_inter[%d]",i),Form("h1_d_inter[%d]",i), 14,0,14);  

        for( int j= 0; j<14 ;j++)
        {
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),14,0,14);  
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),14,0,14);    
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
        }
    }

    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    {
        proton_tree->GetEntry(entry);   
        double sum_p = 0;
        double bar_Change_info[14] = {0};
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;    
        if(p_First_Had_Type < 1 ) continue;
        bar_Change_info[0] = log10((*p_EnergyVec)[0] / 0.023);
        for (int layer = 1; layer < 14; ++layer) 
        {
            if ((*p_EnergyVec)[layer - 1] == 0 || (*p_EnergyVec)[layer] == 0) bar_Change_info[layer] = -5; 
            else  bar_Change_info[layer] = log10((*p_EnergyVec)[layer] / (*p_EnergyVec)[layer - 1]);
            // cout << bar_Change_info[layer] << endl;
        }
        FindMaxPositiveSegment(bar_Change_info,14,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(bar_Change_info,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);


        auto p_start = p_EnergyVec->begin();  auto p_end = p_EnergyVec->end();  
        auto p_maxIt = std::max_element(p_start, p_end);
        auto p_minIt = std::min_element(p_start, p_end);
        int p_maxIndex = std::distance(p_start, p_maxIt);
        int p_minIndex = std::distance(p_start, p_minIt);
        double p_maxVal = *p_maxIt;
        double p_minVal = *p_minIt;
        // cout << seg_peak_value << endl;
        h1_p[p_energy_index][p_First_Had_Layer]->Fill(p_maxIndex-seg_start_idx);
        h1_p_inter[p_energy_index]->Fill(p_maxIndex-seg_start_idx) ;
    }

    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    {
        deuteron_tree->GetEntry(entry);
        double sum_d = 0;
        double bar_Change_info[14] = {0};
        double rate_max_min      = 0;
        double seg_sum           = 0;   // 总增长和
        int    seg_len           = 0;   // 连续正增长长度
        int    seg_start_idx     = 0;   // 连续正增长起点索引
        double seg_peak_value    = 0;   // 正段增长最大值
        int    seg_peak_idx      = 0;   // 正段增长最大值的索引
        double seg_sum_to_peak   = 0;   // 从起点到增长最大值的增长和
        int    seg_len_to_peak   = 0;   // 从起点到增长最大值的索引
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 14) continue;
        if(d_First_Had_Type < 1 ) continue;
        bar_Change_info[0] = log10((*d_EnergyVec)[0] / 0.023);
        for(int layer = 1 ; layer<14 ; layer++)
        {   
            if ((*d_EnergyVec)[layer - 1] == 0 || (*d_EnergyVec)[layer] == 0) bar_Change_info[layer] = -5; 
            else  bar_Change_info[layer] = log10((*d_EnergyVec)[layer] / (*d_EnergyVec)[layer - 1]);
        }
        FindMaxPositiveSegment(bar_Change_info,14,seg_sum,seg_len,seg_start_idx);
        FindMaxValueInPositiveSegment(bar_Change_info,seg_start_idx,seg_len,seg_peak_value,seg_peak_idx);

        auto d_start = d_EnergyVec->begin();  auto d_end = d_EnergyVec->end();  
        auto d_maxIt = std::max_element(d_start, d_end);
        auto d_minIt = std::min_element(d_start, d_end);
        int d_maxIndex = std::distance(d_start, d_maxIt);
        int d_minIndex = std::distance(d_start, d_minIt);
        double d_maxVal = *d_maxIt;
        double d_minVal = *d_minIt;

        h1_d[d_energy_index][d_First_Had_Layer]->Fill(d_maxIndex-seg_start_idx);
        h1_d_inter[d_energy_index]->Fill(d_maxIndex-seg_start_idx) ;
    }

    for (int i = 9; i < 10; i++) // Deposit Energy Bin
    {
        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);

        for (int j = 0; j < 14; j++) // layer
        {
            h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
            h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
            
            c1->cd(j + 1);

            h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral());
            h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral());
            h1_p[i][j]->GetYaxis()->SetRangeUser(0,h1_p[i][j]->GetMaximum()*1.2);

            h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV]Interaction happened in L%d;Layer Distance (Peak - Start);Normalized Count", pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j ));
            h1_p[i][j]->Draw("hist");
            h1_d[i][j]->Draw("histsame");

        }
        c1->cd(15);
        TLatex *tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));tex->SetNDC();tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[i][0], "Proton", "l");
        legend1->AddEntry(h1_d[i][0], "Deuteron", "l");     
        legend1->Draw();       

        auto c2 = new TCanvas("c2","c2",1000,1000);
        c2->cd();
        h1_p_inter[i]->SetLineColor(kRed);   h1_p_inter[i]->SetMarkerColor(kRed);  h1_p_inter[i]->SetLineWidth(2);   h1_p_inter[i]->Sumw2();
        h1_d_inter[i]->SetLineColor(kBlue);  h1_d_inter[i]->SetMarkerColor(kBlue); h1_d_inter[i]->SetLineWidth(2);   h1_d_inter[i]->Sumw2();
        h1_p_inter[i]->Scale(1.0/h1_p_inter[i]->Integral()); 
        h1_d_inter[i]->Scale(1.0/h1_d_inter[i]->Integral()); 
        h1_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
        h1_p_inter[i]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] Stack Multi Layer;Layer Distance (Peak - Start);Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        h1_p_inter[i]->Draw("hist");
        h1_d_inter[i]->Draw("histsame");
    }
}
