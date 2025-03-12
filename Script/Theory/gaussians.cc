void optimize_cross_section_limits(TH1D* h_blue, TH1D* h_red, TH1D* h_green, double &best_s, double &best_c, double &best_sta, double &best_low, double &best_high) {
    int nbins = h_blue->GetNbinsX();

    // 获取搜索起点和终点
    int bin_start = h_blue->FindFirstBinAbove(0);  // 蓝色直方图最左边有值的 bin
    int bin_end = h_blue->FindLastBinAbove(0);     // 蓝色直方图最右边有值的 bin

    best_s = 0;
    best_c = 0;
    best_sta = 0;
    best_low = h_blue->GetXaxis()->GetBinCenter(bin_start);
    best_high = h_blue->GetXaxis()->GetBinCenter(bin_end);

    // 遍历新的搜索范围
    for (int bin_low = bin_start; bin_low <= bin_end; bin_low++) {
        for (int bin_high = bin_low + 1; bin_high <= bin_end; bin_high++) {
            double low_edge = h_blue->GetXaxis()->GetBinCenter(bin_low);
            double high_edge = h_blue->GetXaxis()->GetBinCenter(bin_high);

            // 计算积分
            double sum_blue = h_blue->Integral(bin_low, bin_high);
            double sum_red = h_red->Integral(bin_low, bin_high);
            double sum_green = h_green->Integral(bin_low, bin_high);

            double total = sum_blue + sum_red + sum_green;
            if (total < 50) continue; // 统计量必须大于 1000

            double contamination = (sum_red + sum_green) / total;
            if (contamination >= 0.05) continue; // 污染率必须小于 0.05

            if (sum_blue > 0) {
                double significance = sum_blue / sqrt(total);

                // 记录最高显著性的范围
                if (significance > best_s) {
                    best_s = significance;
                    best_c = contamination;
                    best_sta = sum_blue;
                    best_low = low_edge;
                    best_high = high_edge;
                }
            }
        }
    }

    if (best_s > 0) {
        cout << "最佳范围: [" << best_low << ", " << best_high << "]" << endl;
        cout << "最大显著性: " << best_s << endl;
        cout << "污染比例: " << best_c << endl;
        cout << "统计量: " << best_sta << endl;
    } else {
        cout << "未找到满足条件的区间" << endl;
    }
}

void gaussians() {
    // Create a random number generator
    TRandom3 randGen(0); 

    double total_proton[7]={4.12E+08,1.85E+08,8.30E+07,3.67E+07,1.66E+07,7.48E+06,3.52E+06};
    double total_deuteron[7]={1.57E+05,1.01E+05,4.57E+04,1.92E+04,7.99E+03,3.40E+03,1.40E+03};
    double total_helium[7]={1.93E+06,1.15E+06,6.17E+05,3.05E+05,1.20E+05,5.30E+04,2.31E+04};

    double section_proton[7]={680,659,659,645,645,659,659};
    double section_deuteron[7]={886.16,886.84,886.57,887.36,889.31,892.04,898.63};
    double section_helium[7]={1105,1102,1102,1137,1137,1166,1166};

    double section_proton_reso[7]={3.12E-02,4.49E-02,4.49E-02,5.53E-02,5.53E-02,8.40E-02,8.40E-02};
    double section_helium_reso[7]={3.10E-02,3.76E-02,3.76E-02,4.48E-02,4.48E-02,5.83E-02,5.83E-02};

    double Energy[7]={0};
    // Define histograms
    TH1D* h1[7];  
    TH1D* h2[7];  
    TH1D* h3[7];  
    TLine* UL[7];
    TLine* LL[7];

    double best_s[7] =  {0};
    double best_c[7] =  {0};
    double best_sta[7] =  {0};
    double best_low[7] =  {0};
    double best_high[7] =  {0};

    for (int ii = 0 ; ii <7 ; ii++)
    {
        h1[ii] = new TH1D(Form("h1[%d]",ii), Form("h1[%d]",ii), 600, 0, 1800);
        h2[ii] = new TH1D(Form("h2[%d]",ii), Form("h2[%d]",ii), 600, 0, 1800);
        h3[ii] = new TH1D(Form("h3[%d]",ii), Form("h3[%d]",ii), 600, 0, 1800);
        
        for (int jj = 0; jj < total_proton[ii]; jj++) {
            h1[ii]->Fill(randGen.Gaus(section_proton[ii], section_proton[ii] * section_proton_reso[ii]  ));  
        }
    
        // Second Gaussian: Mean 900, Resolution 10% (Sigma = 90)
        for (int jj = 0; jj < total_deuteron[ii]; jj++) {
            h2[ii]->Fill(randGen.Gaus(section_deuteron[ii], section_deuteron[ii] * section_proton_reso[ii] ));  
        }
    
        // Second Gaussian: Mean 900, Resolution 10% (Sigma = 90)
        for (int jj = 0; jj < total_helium[ii]; jj++) {
            h3[ii]->Fill(randGen.Gaus(section_helium[ii], section_helium[ii] * section_helium_reso[ii] ));  
        }
        // Create a canvas
        auto c1 = new TCanvas("c1", "Gaussian Distributions", 800, 600);
        gPad->SetLogy();
        gStyle->SetOptStat(0);
        // Set histogram styles
        h1[ii]->SetTitle(Form("Gaussian Distributions log_{10}(Energy/GeV)[%.1f, %.1f];Cross Section (mbarn);Counts",  1.6 + 0.2 * ii , 1.8 + 0.2 * ii ));
        h1[ii]->GetYaxis()->SetRangeUser(1,h1[ii]->GetMaximum() * 2 );
        h1[ii]->SetLineColor(kRed);
        h1[ii]->SetLineWidth(2);
        h2[ii]->SetLineColor(kBlue);
        h2[ii]->SetLineWidth(2);
        h3[ii]->SetLineColor(kGreen);
        h3[ii]->SetLineWidth(2);

        // Draw histograms
        h1[ii]->Draw("HIST");
        h2[ii]->Draw("HIST SAME");
        h3[ii]->Draw("HIST SAME");

        // Add legend
        auto legend = new TLegend(0.5, 0.75, 0.88, 0.88);
        legend->AddEntry(h1[ii], Form("Proton Mean %.0f, Sigma %.2f%%",section_proton[ii], section_proton_reso[ii]*100), "l");
        legend->AddEntry(h2[ii], Form("Deuteron Mean %.0f, Sigma %.2f%%",section_deuteron[ii], section_proton_reso[ii]*100), "l");
        legend->AddEntry(h3[ii], Form("Helium4 Mean %.0f, Sigma %.2f%%",section_helium[ii], section_helium_reso[ii]*100), "l");
        legend->Draw();
        // Update canvas

        cout << " log_{10}(E/GeV) = [ " <<  1.6 + 0.2 * ii << " , " << 1.8 + 0.2 * ii  << " ] " << endl;
        optimize_cross_section_limits(h2[ii],h1[ii],h3[ii], best_s[ii], best_c[ii], best_sta[ii], best_low[ii], best_high[ii]);

        LL[ii] = new TLine(best_low[ii],1,best_low[ii],h2[ii]->GetMaximum() * 2);   LL[ii]->SetLineWidth(2); LL[ii]->Draw();
        UL[ii] = new TLine(best_high[ii],1,best_high[ii],h2[ii]->GetMaximum() * 2); UL[ii]->SetLineWidth(2); UL[ii]->Draw();

        c1->SaveAs(Form("/Users/xiongzheng/software/B4/B4c/Script/Theory/gaussians[%.2f,%.2f].pdf",1.6 + 0.2 * ii,1.8 + 0.2 * ii));


    }

    
    
}