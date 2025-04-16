void Y2()
{
    int p_First_Had_Layer;  double p_Total_E;  std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr;   std::vector<double>* p_L_EnergyVec = nullptr;
    int d_First_Had_Layer;  double d_Total_E;  std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr;   std::vector<double>* d_L_EnergyVec = nullptr;

    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_PowerLaw.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
    proton_tree->SetBranchAddress("BarEnergyVector"  ,&p_EnergyVec);
    proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
    proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
    proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("Total_E"          ,&p_Total_E);

    auto deuteron_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_PowerLaw.root");
    auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
    deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
    deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
    deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
    deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
    deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
    deuteron_tree->SetBranchAddress("Total_E"          ,&d_Total_E);
    
    double Energy[15] = {0};
    double Energy_LL[15] = {0};
    double Energy_UL[15] = {0};

    double Layer[14]={0};
    double Layer_Err[14]={0};

    // Depsit and Layer
    TH1D *h1_p[15][14];  TH1D *hC_p[15][14];       TH2D *h2_p[15];      TGraph *gr_proton[15];          int point_counter_p = 0;
    TH1D *h1_d[15][14];  TH1D *hC_d[15][14];       TH2D *h2_d[15];      TGraph *gr_deuteron[15];        int point_counter_d = 0;

    for(int i =0 ; i<15 ; i++)
    {
        Energy[i]    = 1.1 + 0.2 * i;
        Energy_LL[i] = 1.0 + 0.2 * i;
        Energy_UL[i] = 1.2 + 0.2 * i;
        gr_proton[i] = new TGraph();
        gr_deuteron[i] = new TGraph();
        h2_p[i] = new TH2D(Form("h2_p[%d]",i),Form("h2_p[%d]",i),14,0,14,100,-2,0);
        h2_d[i] = new TH2D(Form("h2_d[%d]",i),Form("h2_d[%d]",i),14,0,14,100,-2,0);
        for( int j= 0; j<14 ;j++)
        {
            h1_p[i][j] = new TH1D(Form("h1_p[%d][%d]",i,j), Form("h1_p[%d][%d]",i,j),100,0,1); 
            h1_d[i][j] = new TH1D(Form("h1_d[%d][%d]",i,j), Form("h1_d[%d][%d]",i,j),100,0,1); 
            hC_p[i][j] = new TH1D(Form("hC_p[%d][%d]",i,j), Form("hC_p[%d][%d]",i,j),100,0,1);  
            hC_d[i][j] = new TH1D(Form("hC_d[%d][%d]",i,j), Form("hC_d[%d][%d]",i,j),100,0,1);  
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
            
        }
    }
    
    // proton_file->cd();
    for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        proton_tree->GetEntry(entry);   
        double p_maxVal=0;
        int p_energy_index = int(floor((log10(p_Total_E) - 1) / 0.2));
        if(p_energy_index < 0 || p_energy_index > 14) continue;
        if( (*p_RMSVec)[0]>15 && (*p_RMSVec)[1]>15 && (*p_RMSVec)[2]<45 && (*p_RMSVec)[3]<45 )  // 
        {
            auto p_start = p_Efrac->begin();  auto p_end = p_Efrac->end();  
            p_maxVal = *std::max_element(p_start, p_end); 
            for (int k = 0; k < p_Efrac->size(); k ++)
            {
                h1_p[p_energy_index][k]->Fill((*p_Efrac)[k]/p_maxVal);
            }
        }
    }

    deuteron_file->cd();
    for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
    // for (Long64_t entry = 0; entry < 1; ++entry)
    {
        deuteron_tree->GetEntry(entry); 
        double d_maxVal=0;
        int d_energy_index = int(floor((log10(d_Total_E) - 1) / 0.2));
        if(d_energy_index < 0 || d_energy_index > 14) continue;
        // cout << d_Total_E << "," <<  d_energy_index << endl;
        // cout << d_EnergyVec->size() << endl;
        if((*d_RMSVec)[0]>15 && (*d_RMSVec)[1]>15 && (*d_RMSVec)[2]<45 && (*d_RMSVec)[3]<45) // 
        {
            auto d_start = d_Efrac->begin();  auto d_end = d_Efrac->end();  
            d_maxVal = *std::max_element(d_start, d_end); 
            for (int k = 0; k < d_Efrac->size(); k ++)
            {
                h1_d[d_energy_index][k]->Fill((*d_Efrac)[k]/d_maxVal);
            }
        }
        // cout << "Next " << endl;
    }

    for(int i = 0 ; i<15 ; i++) // Energy
    {
        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
    
        auto c2 = new TCanvas("c2","c2",2500,1500);
        c2->Clear();
        c2->Divide(5,3);
        gStyle->SetOptStat(0);
        double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};     // double SUM_Proton_Ratio[14]={0};     double SUM_Deuteron_Ratio[14]={0};    
        double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0};  // double SUM_Proton_Ratio_LL[14]={0};  double SUM_Deuteron_Ratio_LL[14]={0}; 
        double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0};  // double SUM_Proton_Ratio_UL[14]={0};  double SUM_Deuteron_Ratio_UL[14]={0}; 
        for( int j= 0; j<14 ;j++) // layer
        {
            h1_p[i][j]->SetLineColor(kRed);   h1_p[i][j]->SetMarkerColor(kRed);  // h1_p[i][j]->SetLineWidth(2);   h1_p[i][j]->Sumw2();
            h1_d[i][j]->SetLineColor(kBlue);  h1_d[i][j]->SetMarkerColor(kBlue); // h1_d[i][j]->SetLineWidth(2);   h1_d[i][j]->Sumw2();
            hC_p[i][j]->SetLineColor(kRed);   hC_p[i][j]->SetMarkerColor(kRed);  // hC_p[i][j]->SetLineWidth(2);   hC_p[i][j]->Sumw2(); 
            hC_d[i][j]->SetLineColor(kBlue);  hC_d[i][j]->SetMarkerColor(kBlue); // hC_d[i][j]->SetLineWidth(2);   hC_d[i][j]->Sumw2(); 
            
            c1->cd(j+1);
            h1_p[i][j]->Scale(1.0/h1_p[i][j]->Integral()); 
            h1_d[i][j]->Scale(1.0/h1_d[i][j]->Integral()); 

            h1_p[i][j]->GetYaxis()->SetRangeUser(0,0.14);
            h1_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] in L%d;log10(Jm) = log10(Max Energy Deposit bar/ Total Deposit);Normalized Count",pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j));
            h1_p[i][j]->Draw();
            h1_p[i][j]->Draw("hist");
            h1_d[i][j]->Draw("histsame");

            double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles

            double p_values[3];      h1_p[i][j]->GetQuantiles(3, p_values, quantiles);      Proton_Ratio[j]       = p_values[1];             Proton_Ratio_LL[j] = p_values[1] - p_values[0];               Proton_Ratio_UL[j]   =  p_values[2] - p_values[1];
            double d_values[3];      h1_d[i][j]->GetQuantiles(3, d_values, quantiles);      Deuteron_Ratio[j]     = d_values[1];           Deuteron_Ratio_LL[j] = d_values[1] - d_values[0];               Deuteron_Ratio_UL[j] =  d_values[2] - d_values[1];
            // double SUM_p_values[3];  hC_p[i][j]->GetQuantiles(3, SUM_p_values, quantiles);  SUM_Proton_Ratio[j]   = SUM_p_values[1];   SUM_Proton_Ratio_LL[j]   = SUM_p_values[1] - SUM_p_values[0];   SUM_Proton_Ratio_UL[j]   = SUM_p_values[2] - SUM_p_values[1];
            // double SUM_d_values[3];  hC_d[i][j]->GetQuantiles(3, SUM_d_values, quantiles);  SUM_Deuteron_Ratio[j] = SUM_d_values[1];   SUM_Deuteron_Ratio_LL[j] = SUM_d_values[1] - SUM_d_values[0];   SUM_Deuteron_Ratio_UL[j] = SUM_d_values[2] - SUM_d_values[1];
        
            // c2->cd(j + 1);
            // hC_p[i][j]->Scale(1.0/hC_p[i][j]->Integral());
            // hC_d[i][j]->Scale(1.0/hC_d[i][j]->Integral());
            // hC_p[i][j]->SetTitle(Form("Deposit Energy[%.2fGeV, %.2fGeV] in L%d;#sum_{0}^{%d} (dE/dx) / Edep",pow(10,Energy_LL[i]),pow(10,Energy_UL[i]),j,j));
            // hC_p[i][j]->Draw("hist");
            // hC_d[i][j]->Draw("histsame");
        }
        auto tex = new TLatex(0.1,0.9,Form("Deposit Energy[%.2fGeV, %.2fGeV]",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        c1->cd(15);
        tex->SetNDC();
        tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[i][0], "Proton", "l");
        legend1->AddEntry(h1_d[i][0], "Deuteron", "l");
        legend1->Draw();       
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Y_PowerE/PDF/Jm_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        // c2->cd(15);
        // tex->Draw();
        // legend1->Draw();       
        // c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Y_PowerE/CDF/Jm_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        // auto c3 = new TCanvas("c3","c3",1000,1000);
        // auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_Ratio  ,Layer_Err,Layer_Err,Proton_Ratio_LL  ,Proton_Ratio_UL);
        // auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_Ratio,Layer_Err,Layer_Err,Deuteron_Ratio_LL,Deuteron_Ratio_UL);
        // gre_p->SetMarkerStyle(20);gre_p->SetMarkerColor(kRed);     gre_p->SetLineColor(kRed);     gre_p->SetLineWidth(2);
        // gre_d->SetMarkerStyle(21);gre_d->SetMarkerColor(kBlue);    gre_d->SetLineColor(kBlue);    gre_d->SetLineWidth(2);  
        // gre_d->SetTitle(Form("Deposit Energy [%.2fGeV, %.2fGeV]; BGO Layer; Jm = (dE/dx) / Edep",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));

        // gre_d->GetYaxis()->SetRangeUser(0,1);
        // gre_d->GetXaxis()->SetLimits(0,14);
        // gre_d->Draw("AP");
        // gre_p->Draw("PSAME");


        // auto SUM_axis_top = new TGaxis(0, 1, 14, 1, 0, 32, 510, "-L"); 
        // SUM_axis_top->SetTitle("Radiation Length"); 
        // SUM_axis_top->SetLabelSize(0.02); 
        // SUM_axis_top->SetTitleSize(0.02); 
        // SUM_axis_top->SetTitleOffset(1.2);
        // SUM_axis_top->Draw();

        // auto legend3 = new TLegend(0.70, 0.80, 0.88, 0.88);
        // legend3->SetNColumns(2);
        // legend3->AddEntry(gre_p, "Proton", "ep");
        // legend3->AddEntry(gre_d, "Deuteron", "ep");
        // legend3->Draw();
        // c3->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Y_PowerE/PDF/Jm_BGOLayer_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        // auto c4 = new TCanvas("c4","c4",1500,1000);
        // // gPad->SetLogy(1);
        // auto SUM_gre_p = new TGraphAsymmErrors(14,Layer,SUM_Proton_Ratio  ,Layer_Err,Layer_Err,SUM_Proton_Ratio_LL  ,SUM_Proton_Ratio_UL);
        // auto SUM_gre_d = new TGraphAsymmErrors(14,Layer,SUM_Deuteron_Ratio,Layer_Err,Layer_Err,SUM_Deuteron_Ratio_LL,SUM_Deuteron_Ratio_UL);
        // SUM_gre_p->SetMarkerStyle(20);    SUM_gre_p->SetMarkerColor(kRed);     SUM_gre_p->SetLineColor(kRed);     SUM_gre_p->SetLineWidth(2);
        // SUM_gre_d->SetMarkerStyle(21);    SUM_gre_d->SetMarkerColor(kBlue);    SUM_gre_d->SetLineColor(kBlue);    SUM_gre_d->SetLineWidth(2);  
        // SUM_gre_d->SetTitle(Form("Deposit Energy [%.2fGeV, %.2fGeV]; BGO Layer; #sum (dE/dx) / Edep",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
    
        // SUM_gre_d->GetYaxis()->SetRangeUser(1e-1,1);
        // SUM_gre_d->GetXaxis()->SetLimits(0,14);
        // SUM_gre_d->Draw("AP");
        // SUM_gre_p->Draw("PSAME");

        // SUM_axis_top->Draw();
        // auto legend4 = new TLegend(0.75, 0.12, 0.88, 0.32);
        // legend4->SetNColumns(2);
        // legend4->AddEntry(SUM_gre_p, "Proton", "ep");
        // legend4->AddEntry(SUM_gre_d, "Deuteron", "ep");
        // legend4->Draw();
        // c4->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Y_PowerE/CDF/Jm_BGOLayer_SUM_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        // auto c5 = new TCanvas("c5", "c5", 1500, 1000);
        // // gPad->SetLogy(1);
        // gr_proton[i]->SetTitle(Form("Deposit Energy [%.2fGeV, %.2fGeV]; BGO Layer; #sum (dE/dx) / Edep",pow(10,Energy_LL[i]),pow(10,Energy_UL[i])));
        // gr_proton[i]->GetYaxis()->SetRangeUser(1e-1,1);

        // gr_proton[i]->SetMarkerStyle(20);  // proton: circle
        // gr_proton[i]->SetMarkerColorAlpha(kRed, 0.1);  // 0.0 = fully transparent, 1.0 = fully opaque
        // gr_proton[i]->SetMarkerSize(0.8);

        // gr_deuteron[i]->SetMarkerStyle(21);  // deuteron: square
        // gr_deuteron[i]->SetMarkerColorAlpha(kBlue, 0.1);
        // gr_deuteron[i]->SetMarkerSize(0.8);

        // gr_proton[i]->Draw("AP");
        // gr_deuteron[i]->Draw("P SAME");
        // SUM_axis_top->Draw();
        // auto legend5 = new TLegend(0.75, 0.12, 0.88, 0.32);
        // legend5->AddEntry(gr_proton[i], "Proton", "p");
        // legend5->AddEntry(gr_deuteron[i], "Deuteron", "p");
        // legend5->Draw();
        // c5->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Y_PowerE/CDF/Jm_BGOLayer_SCAT_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));

        // auto c6 = new TCanvas("c6", "c6", 1500, 1000);
        // h2_p[i]->SetTitle(Form("Deposit Energy [%.2fGeV, %.2fGeV]; BGO Layer; #sum (dE/dx) / Edep",Energy_LL[i],Energy_UL[i]));
        // h2_p[i]->GetYaxis()->SetRangeUser(-1.2,0.2);
        // h2_p[i]->SetBarWidth(0.4);
        // h2_p[i]->SetBarOffset(-0.25);
        // h2_p[i]->SetFillColorAlpha(kRed,0.5);
        // h2_p[i]->SetLineColorAlpha(kRed,0.5);
      
        // h2_d[i]->SetBarWidth(0.4);
        // h2_d[i]->SetBarOffset(0.25);
        // h2_d[i]->SetFillColorAlpha(kBlue,0.5);
        // h2_d[i]->SetLineColorAlpha(kBlue,0.5);
        // h2_p[i]->Draw("CANDLEX(112111)");
        // h2_d[i]->Draw("CANDLEX(112111) same");
        // auto axis_top = new TGaxis(0, 0, 14, 0, 0, 32, 510, "-L"); 
        // axis_top->SetTitle("Radiation Length"); // 设置顶部X轴的标题
        // axis_top->SetLabelSize(0.03);  // 设置刻度数字大小
        // axis_top->SetTitleSize(0.03);  // 设置标题大小
        // axis_top->SetTitleOffset(1.2); // 设置标题偏移
        // axis_top->Draw();
        // legend5->Draw();
        // c6->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/PowerE/Y_PowerE/CDF/Jm_BGOLayer_CANDLE_DP_%.2f_%.2f.pdf",Energy_LL[i],Energy_UL[i]));
    }
}
