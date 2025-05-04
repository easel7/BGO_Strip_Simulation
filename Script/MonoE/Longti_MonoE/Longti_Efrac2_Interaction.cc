void Longti_Efrac2_Interaction()
{
    double Energy[19]={0};
    for (int i = 18; i < 19; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
    int p_First_Had_Layer; int p_First_Had_Type; double p_Total_E;     std::vector<double>* p_RMSVec = nullptr;    std::vector<double>* p_EnergyVec = nullptr;    std::vector<double>* p_Efrac = nullptr;
    int d_First_Had_Layer; int d_First_Had_Type; double d_Total_E;     std::vector<double>* d_RMSVec = nullptr;    std::vector<double>* d_EnergyVec = nullptr;    std::vector<double>* d_Efrac = nullptr;

        TH1D *h1_p[14];           TH1D *hC_p[14];
        TH1D *h1_d[14];           TH1D *hC_d[14];

        TH1D *h1_p_inter[15];     TH1D *hC_p_inter[15];
        TH1D *h1_d_inter[15];     TH1D *hC_d_inter[15];

        TH2D *h2_p = new TH2D("h2_p","h2_p",14,0,14,50,-5,0);
        TH2D *h2_d = new TH2D("h2_d","h2_d",14,0,14,50,-5,0);

        double Proton_Ratio[14]={0};     double Deuteron_Ratio[14]={0};    
        double Proton_Ratio_LL[14]={0};  double Deuteron_Ratio_LL[14]={0}; 
        double Proton_Ratio_UL[14]={0};  double Deuteron_Ratio_UL[14]={0}; 

        double SUM_Proton_Ratio[14]={0};     double SUM_Deuteron_Ratio[14]={0};    
        double SUM_Proton_Ratio_LL[14]={0};  double SUM_Deuteron_Ratio_LL[14]={0}; 
        double SUM_Proton_Ratio_UL[14]={0};  double SUM_Deuteron_Ratio_UL[14]={0}; 

        double Layer[14]={0};
        double Layer_Err[14]={0};

        auto *gr_proton = new TGraph();
        auto *gr_deuteron = new TGraph();
        int point_counter_p = 0;
        int point_counter_d = 0;

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        proton_tree->SetBranchAddress("RMS"              ,&p_RMSVec);
        proton_tree->SetBranchAddress("LayerEnergyVector",&p_EnergyVec);
        proton_tree->SetBranchAddress("Efrac"            ,&p_Efrac);
        proton_tree->SetBranchAddress("First_Had_Layer"  ,&p_First_Had_Layer);
    proton_tree->SetBranchAddress("First_Had_Type"  ,&p_First_Had_Type);

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4e/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        deuteron_tree->SetBranchAddress("RMS"              ,&d_RMSVec);
        deuteron_tree->SetBranchAddress("LayerEnergyVector",&d_EnergyVec);
        deuteron_tree->SetBranchAddress("Efrac"            ,&d_Efrac);
        deuteron_tree->SetBranchAddress("First_Had_Layer"  ,&d_First_Had_Layer);
        deuteron_tree->SetBranchAddress("First_Had_Type"   ,&d_First_Had_Type);
        auto c1 = new TCanvas("c1","c1",2500,1500);
        c1->Clear();
        c1->Divide(5,3);
        gStyle->SetOptStat(0);
        auto c2 = new TCanvas("c2","c2",2500,1500);
        c2->Clear();
        c2->Divide(5,3);
        gStyle->SetOptStat(0);
        for (int j = 0; j < 14; j++)
        {

            h1_p[j] = new TH1D(Form("h1_p[%d]",j), Form("h1_p[%d]",j),50,-5,0);  
            h1_d[j] = new TH1D(Form("h1_d[%d]",j), Form("h1_d[%d]",j),50,-5,0);  
            hC_p[j] = new TH1D(Form("hC_p[%d]",j), Form("hC_p[%d]",j),50,-5,0);  
            hC_d[j] = new TH1D(Form("hC_d[%d]",j), Form("hC_d[%d]",j),50,-5,0);  

            h1_p_inter[j] =new TH1D(Form("h1_p_inter[%d]",j),Form("h1_p_inter[%d]",j), 50,-5,0);  
            h1_d_inter[j] =new TH1D(Form("h1_d_inter[%d]",j),Form("h1_d_inter[%d]",j), 50,-5,0);
            hC_p_inter[j] =new TH1D(Form("hC_p_inter[%d]",j),Form("hC_p_inter[%d]",j), 50,-5,0);  
            hC_d_inter[j] =new TH1D(Form("hC_d_inter[%d]",j),Form("hC_d_inter[%d]",j), 50,-5,0); 
        }


        for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
        {
            double sum_p = 0;
            proton_tree->GetEntry(entry);   
            if(p_First_Had_Type!=1) continue;
            h1_p[p_First_Had_Layer]->Fill(log10((*p_Efrac)[p_First_Had_Layer]));    
            h1_p_inter[p_First_Had_Layer]->Fill(log10((*p_Efrac)[p_First_Had_Layer]));    
            for(int i=0 ;i<=p_First_Had_Layer; i++)
            {
                sum_p += (*p_Efrac)[i];
            }
            hC_p[p_First_Had_Layer]->Fill(log10(sum_p));
            hC_p_inter[p_First_Had_Layer]->Fill(log10 (sum_p) );
        }
        for (Long64_t entry = 0; entry < deuteron_tree->GetEntries(); ++entry)
        {
            deuteron_tree->GetEntry(entry);
            double sum_d = 0;
            if(d_First_Had_Type!=1) continue;
            h1_d[d_First_Had_Layer]->Fill(log10((*d_Efrac)[d_First_Had_Layer]));
            h1_d_inter[d_First_Had_Layer]->Fill(log10 ((*d_Efrac)[d_First_Had_Layer]) );
            for(int i=0 ;i<=d_First_Had_Layer; i++)
            {
                sum_d += (*d_Efrac)[i];
            }
            hC_d[d_First_Had_Layer]->Fill(log10(sum_d));
            hC_d_inter[d_First_Had_Layer]->Fill(log10 (sum_d) );
        }



        for (int j = 0; j < 14; j++)
        {
            h1_p[j]->Sumw2(); h1_p[j]->Scale(1.0/h1_p[j]->Integral());h1_p[j]->SetLineColor(kRed);     h1_p[j]->SetMarkerColor(kRed);     h1_p[j]->SetLineWidth(2);
            h1_d[j]->Sumw2(); h1_d[j]->Scale(1.0/h1_d[j]->Integral());h1_d[j]->SetLineColor(kBlue);    h1_d[j]->SetMarkerColor(kBlue);    h1_d[j]->SetLineWidth(2);
            hC_p[j]->Sumw2(); hC_p[j]->Scale(1.0/hC_p[j]->Integral());hC_p[j]->SetLineColor(kRed);     hC_p[j]->SetMarkerColor(kRed);     hC_p[j]->SetLineWidth(2);
            hC_d[j]->Sumw2(); hC_d[j]->Scale(1.0/hC_d[j]->Integral());hC_d[j]->SetLineColor(kBlue);    hC_d[j]->SetMarkerColor(kBlue);    hC_d[j]->SetLineWidth(2);
            c1->cd(j + 1);
          
            h1_p[j]->GetYaxis()->SetRangeUser(0,0.25);h1_p[j]->SetTitle(Form("Distrubution of EdepRatio in L%d;log_{10}(Energy Deposit Ratio);Normalized Count",j));
            h1_p[j]->Draw("hist");
            h1_d[j]->Draw("histsame");

            c2->cd(j + 1);
            hC_p[j]->GetYaxis()->SetRangeUser(0,0.25);hC_p[j]->SetTitle(Form("Distrubution of #sum_{0}^{L%d} EdepRatio;log_{10}(#sum(Energy Deposit Ratio));Normalized Count",j));
            hC_p[j]->Draw("hist");
            hC_d[j]->Draw("histsame");
        
            Layer[j] = 0.5 + j;
            Layer_Err[j] = 0.5;
        }
        c1->cd(15);
        TLatex *tex = new TLatex(0.1,0.9,Form("IncidentParicle %dGeV",int(Energy[i])));tex->SetNDC();tex->Draw(); 
        auto legend1 = new TLegend(0.12, 0.12, 0.88, 0.88);
        legend1->AddEntry(h1_p[0], "Proton", "l");
        legend1->AddEntry(h1_d[0], "Deuteron", "l");  
        legend1->Draw();       

        c2->cd(15);
        tex->Draw();
        legend1->Draw();       

        auto c3 = new TCanvas("c3","c3",1000,1000);
        c3->cd();
        h1_p_inter[i]->SetLineColor(kRed);   h1_p_inter[i]->SetMarkerColor(kRed);  h1_p_inter[i]->SetLineWidth(2);   h1_p_inter[i]->Sumw2();
        h1_d_inter[i]->SetLineColor(kBlue);  h1_d_inter[i]->SetMarkerColor(kBlue); h1_d_inter[i]->SetLineWidth(2);   h1_d_inter[i]->Sumw2();
        h1_p_inter[i]->Scale(1.0/h1_p_inter[i]->Integral()); 
        h1_d_inter[i]->Scale(1.0/h1_d_inter[i]->Integral()); 
        h1_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
        h1_p_inter[i]->SetTitle(Form("IncidentParicle %dGeV Stack Multi Layer;log_{10}(Energy Deposit in Layer/Total Deposit);Normalized Count",int(Energy[i])));
        h1_p_inter[i]->Draw("hist");
        h1_d_inter[i]->Draw("histsame");

        auto c4 = new TCanvas("c4","c4",1000,1000);
        c4->cd();
        hC_p_inter[i]->SetLineColor(kRed);   hC_p_inter[i]->SetMarkerColor(kRed);  hC_p_inter[i]->SetLineWidth(2);   hC_p_inter[i]->Sumw2();
        hC_d_inter[i]->SetLineColor(kBlue);  hC_d_inter[i]->SetMarkerColor(kBlue); hC_d_inter[i]->SetLineWidth(2);   hC_d_inter[i]->Sumw2();
        hC_p_inter[i]->Scale(1.0/hC_p_inter[i]->Integral()); 
        hC_d_inter[i]->Scale(1.0/hC_d_inter[i]->Integral()); 
        hC_p_inter[i]->GetYaxis()->SetRangeUser(0,h1_p_inter[i]->GetMaximum()*1.2);
        hC_p_inter[i]->SetTitle(Form("IncidentParicle %dGeV Stack Multi Layer;#sum_{0}^{Interaction Layer} log_{10}(Energy Deposit in Layer/Total Deposit);Normalized Count",int(Energy[i])));
        hC_p_inter[i]->Draw("hist");
        hC_d_inter[i]->Draw("histsame");

        
    }
    
}