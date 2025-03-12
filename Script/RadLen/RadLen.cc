void RadLen()
{
    double Energy[19]={0};
    auto c1 = new TCanvas("c1","c1",1800,900);
    c1->Divide(6,3);
    auto c2 = new TCanvas("c2","c2",900,900);
    c2->cd();
    for (int i =0; i < 19; i++)
    {
        
        if(i<9)  {Energy[i] =  (i+1)*10;}
        else   {Energy[i] =  i*100-800;}
        cout << "Energy = " << int(Energy[i]) << " GeV !" << endl;
        auto h1_e = new TH1D("h1_e","h1_e",14,0,357); // 14 bins, 0 - 357 mm, 25.5mm layer  
        auto electron_file = TFile::Open(Form("/Users/xiongzheng/software/B4/B4c/Root/Electron_%dGeV.root",int(Energy[i])));
        auto electron_tree = (TTree*)electron_file->Get("B4");
        int nEntries = electron_tree->GetEntries();

        const int numLayers = 14;  // Assuming 14 layers
        double L_E[numLayers];     // Array to store branch values
        double Layer_energy[numLayers];// Total energy sum
        double total_energy;// Total energy sum
        

        // Set branch addresses
        for (int j = 0; j < numLayers; j++)
        {
            electron_tree->SetBranchAddress(Form("L%d_E", j), &L_E[j]);
        }

        for (int entry = 0; entry < nEntries; entry++)
        {
            electron_tree->GetEntry(entry);
            for (int j = 0; j < numLayers; j++)
            {
                Layer_energy[j] += L_E[j];
                total_energy += L_E[j];
            }
        }
        double inte_energy = 0;
        for (int j = 0; j < numLayers; j++)
        {
            cout << "Layer " << j << " , Layer_energy[j] = " << Layer_energy[j] << endl; 
            inte_energy += Layer_energy[j]/1e4;
            h1_e->SetBinContent(j+1,  Energy[i] - inte_energy );   
        }
        cout << "Total_energy = " << total_energy / 1e4 << "GeV / events"<< endl;    
        
        c1->cd(i+1);
        gPad->SetLogy();
        gStyle->SetOptFit(1111);
        h1_e->SetTitle("Radiation Length; Energy ;");
        h1_e->Draw();
        auto fitFunc_e = new TF1("fitFunc_e","[0]*exp(-x/[1])", 0, 40);
        fitFunc_e->SetParameters(10000, 2); // 初始参数：振幅=100, λ=10 cm
        h1_e->Fit(fitFunc_e, "R"); // 进行拟合
    }
}