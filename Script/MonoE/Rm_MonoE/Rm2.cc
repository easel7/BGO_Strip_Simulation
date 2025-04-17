void Rm2()
{
    double Energy[19]={0};
    for (int i = 0; i < 19; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else if(i>=9 && i < 19)   {Energy[i] =  i*100-800;}
        else   {Energy[i] = i*1000 - 18000;}

        int p_First_Had_Layer; double p_E_total;   std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr;     std::vector<double>* p_L_EnergyVec = nullptr;
        int d_First_Had_Layer; double d_E_total;   std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr;     std::vector<double>* d_L_EnergyVec = nullptr;

        TH1D *h1_p[14];  TH1D *hC_p[14];  TF1  *fitFunc_p[14];         double Layer[14]={0};
        TH1D *h1_d[14];  TH1D *hC_d[14];  TF1  *fitFunc_d[14];          double Layer_Err[14]={0};

        TH2D *h2_p = new TH2D("h2_p","h2_p",14,0,14,50,-5,0);
        TH2D *h2_d = new TH2D("h2_d","h2_d",14,0,14,50,-5,0);
    
        double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};            double SUM_Proton_Ratio[14]={0};     double SUM_Deuteron_Ratio[14]={0};    
        double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0};         double SUM_Proton_Ratio_LL[14]={0};  double SUM_Deuteron_Ratio_LL[14]={0}; 
        double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0};         double SUM_Proton_Ratio_UL[14]={0};  double SUM_Deuteron_Ratio_UL[14]={0}; 


        auto *gr_proton = new TGraph();
        auto *gr_deuteron = new TGraph();
        int point_counter_p = 0;
        int point_counter_d = 0;

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/Root/Proton_1000GeV.root");
        // auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/build/Test_3.root");
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
        proton_tree->SetBranchAddress("BarEnergyVector",  &p_EnergyVec);
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_L_EnergyVec);
        proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
        proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
        proton_tree->SetBranchAddress("Total_E"          ,&p_E_total);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
        deuteron_tree->SetBranchAddress("BarEnergyVector"  ,&d_EnergyVec);
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_L_EnergyVec);
        deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
        deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
        deuteron_tree->SetBranchAddress("Total_E"          ,&d_E_total);

        auto c1 = new TCanvas("c1","c1",2500,1500);
        // auto c1 = new TCanvas("c1","c1",1200,900);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);

        auto c2 = new TCanvas("c2","c2",2500,1500);
        c2->Clear();
        c2->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j < 14; j++)
        {
            h1_p[j] = new TH1D(Form("h1_p[%d]",j), Form("h1_p[%d]",j),40,0,0.4);  
            h1_d[j] = new TH1D(Form("h1_d[%d]",j), Form("h1_d[%d]",j),40,0,0.4);  
            hC_p[j] = new TH1D(Form("hC_p[%d]",j), Form("hC_p[%d]",j),100,0,1);  
            hC_d[j] = new TH1D(Form("hC_d[%d]",j), Form("hC_d[%d]",j),100,0,1);  
        }

        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
        // for (Long64_t entry = 0; entry < 1; ++entry)
        {
            proton_tree->GetEntry(entry);   
            double sum_p = 0;
            double p_maxVal[14]={0}; 
            if((*p_L_EnergyVec)[0]>0.23 && (*p_L_EnergyVec)[1]>0.23 && (*p_RMSVec)[2]<40 && (*p_RMSVec)[3]<40  ) // && (*p_RMSVec)[0]<40 && (*p_RMSVec)[1]<40
            {
                // for (size_t k = 0; k < p_EnergyVec->size(); k ++) {sum_d += (*p_EnergyVec)[k];}
                for (size_t k = 0; k < p_EnergyVec->size(); k += 22)
                {
                    int index = int(k / 22);  // Get the Layer
                    auto p_start = p_EnergyVec->begin() + k;  auto p_end = (k + 22 < p_EnergyVec->size() ) ? p_start + 22 : p_EnergyVec->end();  p_maxVal[index] = *std::max_element(p_start, p_end); 
                    h1_p[index]->Fill(p_maxVal[index]/p_E_total);
                    sum_p += p_maxVal[index];
                    if (p_E_total > 0 && sum_p > 0)
                    {
                        hC_p[index]->Fill(sum_p / p_E_total);
                        gr_proton->SetPoint(point_counter_p++, index+0.9, sum_p / p_E_total);
                        h2_p->Fill(index,log10(sum_p / p_E_total));
                    }
                }
                // cout << " bar " << k << " Layer " << index << " p_max" << p_maxVal[index] << " Sum p_max " << sum_p << endl; 
                // std::copy(p_start, p_end, std::ostream_iterator<double>(std::cout, ", "));
                // cout << " max = "<< p_maxVal[index] <<  endl;
            }
            // cout <<  "sum p = " << sum_p <<  endl;
            // cout << "Total " << p_E_total << endl;
            // cout << p_EnergyVec->size() << endl;
        }

        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
        {
            deuteron_tree->GetEntry(entry); 
            double sum_d = 0;
            double d_maxVal[14]={0};
            if((*d_L_EnergyVec)[0]>0.23 && (*d_L_EnergyVec)[1]>0.23 && (*d_RMSVec)[2]<40 && (*d_RMSVec)[3]<40  ) //  && (*d_RMSVec)[0]<40 && (*d_RMSVec)[1]<40
            {
                for (size_t k = 0; k < d_EnergyVec->size(); k += 22)
                {
                    int index = int(k / 22);  // Get the Layer
                    auto d_start = d_EnergyVec->begin() + k;  auto d_end = (k + 22 < d_EnergyVec->size() ) ? d_start + 22 : d_EnergyVec->end();  d_maxVal[index] = *std::max_element(d_start, d_end); 
                    h1_d[index]->Fill(d_maxVal[index]/d_E_total);
                    sum_d += d_maxVal[index];
                    if (d_E_total > 0 && sum_d > 0)
                    {
                        hC_d[index]->Fill(sum_d / d_E_total);
                        gr_deuteron->SetPoint(point_counter_d++, index+1.1, sum_d / d_E_total);
                        h2_d->Fill(index,log10(sum_d / d_E_total));
                    }
                }
            }
            // cout <<  "sum d = " << sum_d <<  endl;
        }

        for (int j = 0; j < 14; j++)
        {
            h1_p[j]->SetLineColor(kRed);     h1_p[j]->SetMarkerColor(kRed);     h1_p[j]->SetLineWidth(2);   h1_p[j]->Sumw2();
            h1_d[j]->SetLineColor(kBlue);    h1_d[j]->SetMarkerColor(kBlue);    h1_d[j]->SetLineWidth(2);   h1_d[j]->Sumw2();
            hC_p[j]->SetLineColor(kRed);     hC_p[j]->SetMarkerColor(kRed);     hC_p[j]->SetLineWidth(2);   hC_p[j]->Sumw2(); 
            hC_d[j]->SetLineColor(kBlue);    hC_d[j]->SetMarkerColor(kBlue);    hC_d[j]->SetLineWidth(2);   hC_d[j]->Sumw2(); 

            c2->cd(j + 1);
            hC_p[j]->Scale(1.0/hC_p[j]->Integral());
            hC_d[j]->Scale(1.0/hC_d[j]->Integral());
            hC_p[j]->SetTitle(Form("%.f GeV incident in L%d;#sum Rm = #sum_{L0}^{L%d} Max Energy Deposit bar / Total Deposit",Energy[i],j,j));
            hC_p[j]->Draw("hist");
            hC_d[j]->Draw("histsame");
            
            c1->cd(j + 1);
            h1_p[j]->Scale(1.0/h1_p[j]->Integral()); 
            h1_d[j]->Scale(1.0/h1_d[j]->Integral()); 
            h1_p[j]->GetYaxis()->SetRangeUser(0,0.6);
            h1_p[j]->SetTitle(Form("%.f GeV incident in L%d;Rm = Max Energy Deposit bar/ Total Deposit;Normalized Count",Energy[i],j));
            h1_p[j]->Draw();
            h1_p[j]->Draw("hist");
            h1_d[j]->Draw("histsame");

            double quantiles[3] = {0.16, 0.50, 0.84};  // Percentiles
            double p_values[3];  h1_p[j]->GetQuantiles(3, p_values, quantiles);  Proton_Ratio[j]   = p_values[1];     Proton_Ratio_LL[j] = p_values[1] - p_values[0];   Proton_Ratio_UL[j]   = p_values[2] - p_values[1];
            double d_values[3];  h1_d[j]->GetQuantiles(3, d_values, quantiles);  Deuteron_Ratio[j] = d_values[1];   Deuteron_Ratio_LL[j] = d_values[1] - d_values[0];   Deuteron_Ratio_UL[j] = d_values[2] - d_values[1];
            double SUM_p_values[3];  hC_p[j]->GetQuantiles(3, SUM_p_values, quantiles);  SUM_Proton_Ratio[j]   = SUM_p_values[1];   SUM_Proton_Ratio_LL[j]   = SUM_p_values[1] - SUM_p_values[0];   SUM_Proton_Ratio_UL[j]   = SUM_p_values[2] - SUM_p_values[1];
            double SUM_d_values[3];  hC_d[j]->GetQuantiles(3, SUM_d_values, quantiles);  SUM_Deuteron_Ratio[j] = SUM_d_values[1];   SUM_Deuteron_Ratio_LL[j] = SUM_d_values[1] - SUM_d_values[0];   SUM_Deuteron_Ratio_UL[j] = SUM_d_values[2] - SUM_d_values[1];
            
            cout << SUM_Proton_Ratio[j]  << endl;
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
        }        
        c1->cd(15);
        TLatex *tex = new TLatex(0.1,0.9,Form("IncidentParicle %dGeV",int(Energy[i])));tex->SetNDC();tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(hC_p[0], "Proton", "l");
        legend1->AddEntry(hC_d[0], "Deuteron", "l");
        legend1->Draw();       
        c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Rm_MonoE/PDF/Rm_DP_%dGeV.pdf",int(Energy[i])));

        c2->cd(15);
        tex->Draw();
        legend1->Draw();       

        c2->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Rm_MonoE/CDF/Rm_DP_%dGeV.pdf",int(Energy[i])));

        auto c3 = new TCanvas("c3","c3",1000,1000);
        c3->cd();
        auto gre_p = new TGraphAsymmErrors(14,Layer,Proton_Ratio  ,Layer_Err,Layer_Err,Proton_Ratio_LL  ,Proton_Ratio_UL);
        auto gre_d = new TGraphAsymmErrors(14,Layer,Deuteron_Ratio,Layer_Err,Layer_Err,Deuteron_Ratio_LL,Deuteron_Ratio_UL);

        gre_p->SetMarkerStyle(20);gre_p->SetMarkerColor(kRed);     gre_p->SetLineColor(kRed);     gre_p->SetLineWidth(2);
        gre_d->SetMarkerStyle(21);gre_d->SetMarkerColor(kBlue);    gre_d->SetLineColor(kBlue);    gre_d->SetLineWidth(2);  
        gre_d->SetTitle(Form("Incident Energy %d GeV ; BGO Layer; Rm = (dE/dx) / Edep",int(Energy[i])));

        gre_d->GetYaxis()->SetRangeUser(0,0.2);
        gre_d->GetXaxis()->SetLimits(0,14);
        gre_d->Draw("AP");
        gre_p->Draw("PSAME");

        auto axis_top = new TGaxis(0, 0.2, 14, 0.2, 0, 32, 510, "-L"); 
        axis_top->SetTitle("Radiation Length"); 
        axis_top->SetLabelSize(0.03);  
        axis_top->SetTitleSize(0.03);  
        axis_top->SetTitleOffset(1.2); 
        axis_top->Draw();
        auto legend3 = new TLegend(0.70, 0.80, 0.88, 0.88);
        legend3->SetNColumns(2);
        legend3->AddEntry(gre_p, "Proton", "ep");
        legend3->AddEntry(gre_d, "Deuteron", "ep");
        legend3->Draw();
        c3->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Rm_MonoE/PDF/Rm_BGOLayer_DP_%dGeV.pdf",int(Energy[i])));


        auto c4 = new TCanvas("c4","c4",1500,1000);
        c4->cd();
        gPad->SetLogy(1);

        auto SUM_gre_p = new TGraphAsymmErrors(14,Layer,SUM_Proton_Ratio  ,Layer_Err,Layer_Err,SUM_Proton_Ratio_LL  ,SUM_Proton_Ratio_UL);
        auto SUM_gre_d = new TGraphAsymmErrors(14,Layer,SUM_Deuteron_Ratio,Layer_Err,Layer_Err,SUM_Deuteron_Ratio_LL,SUM_Deuteron_Ratio_UL);
        SUM_gre_p->SetMarkerStyle(20);    SUM_gre_p->SetMarkerColor(kRed);     SUM_gre_p->SetLineColor(kRed);     SUM_gre_p->SetLineWidth(2);
        SUM_gre_d->SetMarkerStyle(21);    SUM_gre_d->SetMarkerColor(kBlue);    SUM_gre_d->SetLineColor(kBlue);    SUM_gre_d->SetLineWidth(2);  
        SUM_gre_d->SetTitle(Form("Incident Energy %d GeV ; BGO Layer; #sum Rm = #sum Max Energy Deposit bar/ Total Deposit ",int(Energy[i])));
    
        SUM_gre_d->GetYaxis()->SetRangeUser(1e-5,1);
        SUM_gre_d->GetXaxis()->SetLimits(0,14);
        SUM_gre_d->Draw("AP");
        SUM_gre_p->Draw("PSAME");

        auto SUM_axis_top = new TGaxis(0, 1, 14, 1, 0, 32, 510, "-L"); 
        SUM_axis_top->SetTitle("Radiation Length"); // 设置顶部X轴的标题
        SUM_axis_top->SetLabelSize(0.02);  // 设置刻度数字大小
        SUM_axis_top->SetTitleSize(0.02);  // 设置标题大小
        SUM_axis_top->SetTitleOffset(1.2); // 设置标题偏移
        SUM_axis_top->Draw();
        auto legend4 = new TLegend(0.75, 0.12, 0.88, 0.32);
        legend4->SetNColumns(2);
        legend4->AddEntry(SUM_gre_p, "Proton", "ep");
        legend4->AddEntry(SUM_gre_d, "Deuteron", "ep");
        legend4->Draw();
        c4->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Rm_MonoE/CDF/Rm_BGOLayer_SUM_DP_%dGeV.pdf",int(Energy[i])));


        auto c5 = new TCanvas("c5", "c5", 1500, 1000);
        c5->cd();
        gPad->SetLogy(1);
        gr_proton->SetTitle(Form("Incident Energy %d GeV ; BGO Layer; #sum Rm = #sum Max Energy Deposit bar/ Total Deposit ",int(Energy[i])));
        gr_proton->GetYaxis()->SetRangeUser(1e-5,1);

        gr_proton->SetMarkerStyle(20);  // proton: circle
        gr_proton->SetMarkerColorAlpha(kRed, 0.1);  // 0.0 = fully transparent, 1.0 = fully opaque
        gr_proton->SetMarkerSize(0.8);

        gr_deuteron->SetMarkerStyle(21);  // deuteron: square
        gr_deuteron->SetMarkerColorAlpha(kBlue, 0.1);
        gr_deuteron->SetMarkerSize(0.8);

        gr_proton->Draw("AP");
        gr_deuteron->Draw("P SAME");
        SUM_axis_top->Draw();
        auto legend5 = new TLegend(0.75, 0.12, 0.88, 0.32);
        legend5->AddEntry(gr_proton, "Proton", "p");
        legend5->AddEntry(gr_deuteron, "Deuteron", "p");
        legend5->Draw();
        c5->SaveAs( Form("/Users/xiongzheng/software/B4/B4e/Script/MonoE/Rm_MonoE/CDF/Rm_BGOLayer_SCAT_DP_%dGeV.pdf",int(Energy[i])));

        auto c6 = new TCanvas("c6", "c6", 1500, 1000);
        c6->cd();
        h2_p->SetTitle(Form("Incident Energy %d GeV ; BGO Layer; #sum Rm = #sum Max Energy Deposit bar/ Total Deposit ",int(Energy[i])));
        h2_p->GetYaxis()->SetRangeUser(-5,0.2);
        h2_p->SetBarWidth(0.4);
        h2_p->SetBarOffset(-0.25);
        h2_p->SetFillColorAlpha(kRed,0.5);
        h2_p->SetLineColorAlpha(kRed,0.5);
      
        h2_d->SetBarWidth(0.4);
        h2_d->SetBarOffset(0.25);
        h2_d->SetFillColorAlpha(kBlue,0.5);
        h2_d->SetLineColorAlpha(kBlue,0.5);
        h2_p->Draw("CANDLEX(112111)");
        h2_d->Draw("CANDLEX(112111) same");
        axis_top->Draw();
        legend5->Draw();

    }
}