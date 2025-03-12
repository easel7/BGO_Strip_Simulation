void Vertex()
{
    double Energy[19]={0};
    TH1D* h1_p[19]; TH1D* hC_p[19];
    TH1D* h1_d[19]; TH1D* hC_d[19];
    TH1D* h1_e[19]; TH1D* hC_e[19];
    TH1D* h1_h[19]; TH1D* hC_h[19];
    TH1D* h1_H[19]; TH1D* hC_H[19];
    TH1D* h1_c[19]; TH1D* hC_c[19];

    TCut EM = "First_Type==0";
    TCut HD = "First_Had_Depth>=0";
    TCut HI = "First_Had_Depth>=0 && First_Had_Type==1";

    double N_proton[19][10];
    double N_deuteron[19][10];
    double N_electron[19][10];
    double N_helium4[19][10];
    double N_helium3[19][10];
    double N_carbon[19][10];


    for (int i = 0; i < 1; i++) // Energy
    {
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}

        auto proton_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Proton_%dGeV.root",int(Energy[i])));
        auto proton_tree = (TTree*)proton_file->Get("B4");
        h1_p[i] = new TH1D(Form("h1_p[%d]",i), Form("h1_p[%d]",i),15,-1,14);  
        proton_tree->Draw(Form("First_Had_Layer>>h1_p[%d]",i),"","");
        hC_p[i] = new TH1D(Form("h1_p[%d]",i), Form("h1_p[%d]",i),14,0,14);  

        auto deuteron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Deuteron_%dGeV.root",int(Energy[i])));
        auto deuteron_tree = (TTree*)deuteron_file->Get("B4");
        h1_d[i] = new TH1D(Form("h1_d[%d]",i), Form("h1_d[%d]",i),15,-1,14);  
        deuteron_tree->Draw(Form("First_Had_Layer>>h1_d[%d]",i),"", "");
        hC_d[i] = new TH1D(Form("h1_d[%d]",i), Form("h1_d[%d]",i),14,0,14);  

        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        h1_e[i] = new TH1D(Form("h1_e[%d]",i), Form("h1_e[%d]",i),15,-1,14);  
        electron_tree->Draw(Form("First_Had_Layer>>h1_e[%d]",i),"", "");
        hC_e[i] = new TH1D(Form("h1_e[%d]",i), Form("h1_e[%d]",i),14,0,14);  

        auto helium4_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium4_%dGeV.root",int(Energy[i])));
        auto helium4_tree = (TTree*)helium4_file->Get("B4");
        h1_h[i] = new TH1D(Form("h1_h[%d]",i), Form("h1_h[%d]",i),15,-1,14);  
        helium4_tree->Draw(Form("First_Had_Layer>>h1_h[%d]",i),"", "");
        hC_h[i] = new TH1D(Form("h1_h[%d]",i), Form("h1_h[%d]",i),14,0,14);  

        auto helium3_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Helium3_%dGeV.root",int(Energy[i])));
        auto helium3_tree = (TTree*)helium3_file->Get("B4");
        h1_H[i] = new TH1D(Form("h1_H[%d]",i), Form("h1_H[%d]",i),15,-1,14);  
        helium3_tree->Draw(Form("First_Had_Layer>>h1_H[%d]",i),"", "");
        hC_H[i] = new TH1D(Form("h1_H[%d]",i), Form("h1_H[%d]",i),14,0,14);  

        auto carbon_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Carbon_%dGeV.root",int(Energy[i])));
        auto carbon_tree = (TTree*)carbon_file->Get("B4");
        h1_c[i] = new TH1D(Form("h1_c[%d]",i), Form("h1_c[%d]",i),15,-1,14);  
        carbon_tree->Draw(Form("First_Had_Layer>>h1_c[%d]",i),"", "");
        hC_c[i] = new TH1D(Form("h1_c[%d]",i), Form("h1_c[%d]",i),14,0,14);  



        for (int jj=1 ; jj<=14 ;jj++)
        {
                hC_p[i]->SetBinContent(jj,1e4-h1_p[i]->Integral(2,jj+1));  
                hC_d[i]->SetBinContent(jj,1e4-h1_d[i]->Integral(2,jj+1));  
                hC_e[i]->SetBinContent(jj,1e4-h1_e[i]->Integral(2,jj+1));  
                hC_h[i]->SetBinContent(jj,1e4-h1_h[i]->Integral(2,jj+1));  
                hC_H[i]->SetBinContent(jj,1e4-h1_H[i]->Integral(2,jj+1));  
                hC_c[i]->SetBinContent(jj,1e4-h1_c[i]->Integral(2,jj+1));  
                // cout << hC_p[i] ->GetBinContent(jj) <<  endl;
                // cout << h1_p[i]->FindBin(0)<<  endl; // 找到 x=0 对应的 bin 索引
        }
        for (int jj=0 ; jj<10 ;jj++)
        {                
            N_proton[i][jj] = h1_p[i]->GetBinContent(jj+2);
            N_deuteron[i][jj] = h1_d[i]->GetBinContent(jj+2);
            if(jj==9)
            {
                N_proton[i][jj] = h1_p[i]->GetBinContent(15)+h1_p[i]->GetBinContent(14)+h1_p[i]->GetBinContent(13)+h1_p[i]->GetBinContent(12)+h1_p[i]->GetBinContent(11);
                N_deuteron[i][jj] = h1_d[i]->GetBinContent(15)+h1_d[i]->GetBinContent(14)+h1_d[i]->GetBinContent(13)+h1_d[i]->GetBinContent(12)+h1_d[i]->GetBinContent(11);
            }
            cout << N_proton[i][jj] << ", " <<  N_deuteron[i][jj] << endl;
        }

        h1_p[i]->SetLineWidth(2);h1_p[i]->SetLineColor(kRed);
        h1_d[i]->SetLineWidth(2);h1_d[i]->SetLineColor(kBlue);
        h1_e[i]->SetLineWidth(2);h1_e[i]->SetLineColor(kGreen+3);
        h1_h[i]->SetLineWidth(2);h1_h[i]->SetLineColor(kGreen+3);
        h1_H[i]->SetLineWidth(2);h1_H[i]->SetLineColor(kOrange-3);h1_H[i]->SetLineStyle(2);
        h1_c[i]->SetLineWidth(2);h1_c[i]->SetLineColor(kMagenta);

        hC_p[i]->SetLineWidth(2);hC_p[i]->SetLineColor(kRed);
        hC_d[i]->SetLineWidth(2);hC_d[i]->SetLineColor(kBlue);
        hC_e[i]->SetLineWidth(2);hC_e[i]->SetLineColor(kGreen+3);
        hC_h[i]->SetLineWidth(2);hC_h[i]->SetLineColor(kGreen+3);
        hC_H[i]->SetLineWidth(2);hC_H[i]->SetLineColor(kOrange-3);hC_H[i]->SetLineStyle(2);
        hC_c[i]->SetLineWidth(2);hC_c[i]->SetLineColor(kMagenta);

        auto c1 = new TCanvas("c1","c1",1200,900);
        c1->Clear();
        gStyle->SetOptStat(0);
        h1_p[i]->SetTitle("Inelastic Interaction Layer;Layer;N_{counts}");
        h1_p[i]->GetYaxis()->SetRangeUser(0,2e3);
        h1_p[i]->Draw("");
        h1_d[i]->Draw("same");

        auto c2 = new TCanvas("c2","c2",1200,900);
        c2->Clear();
        c2->cd();
        gStyle->SetOptStat(0);
        hC_p[i]->SetTitle("Inelastic Interaction Layer;Layer;N_{survive}");
        hC_p[i]->GetYaxis()->SetRangeUser(0,1e4);
        hC_p[i]->Draw("");
        hC_d[i]->Draw("same");


        // h1_d[i]->Draw("same");
        // c1->SaveAs( Form("/Users/xiongzheng/software/B4/B4c/Script/Vertex/HD_Layer_%dGeV.pdf",int(Energy[i])));
        
    

    }
}