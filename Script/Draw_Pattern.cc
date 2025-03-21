void Draw_Pattern()
{
    std::vector<double>* energyVec = nullptr;
    
    auto proton_file = TFile::Open("/Users/xiongzheng/software/B4/B4e/build/Test_BK.root");
    auto proton_tree = (TTree*)proton_file->Get("B4");
    // Set the branch address to the vector pointer
    proton_tree->SetBranchAddress("BarEnergyVector", &energyVec);
    cout  << proton_tree->GetEntries() << endl;
    auto c1    = new TCanvas("c1","c1",1800,600);
    auto hXZ   = new TH2D("hXZ","BGO X-Z Plane",24,-12,12,14,0,14);
    auto hYZ   = new TH2D("hYZ","BGO Y-Z Plane",24,-12,12,14,0,14);
    auto h3box = new TH3D("h3box","3D View 100 GeV Proton Cascade",24,-12,12,24,-12,12,14,0,14);

    double minVal = -3;  // 设置最小值
    double maxVal = 2;  // 设置最大值
    // for (Long64_t entry = 0; entry < proton_tree->GetEntries(); ++entry)
    // {
        proton_tree->GetEntry(1);
        for(size_t i = 0; i < energyVec->size(); ++i) 
        {
            int layer = i / 24; 
            int bar = i % 24;
            double energy = 0;
            if ((*energyVec)[i] < 1e-3)  {energy = -5;              }
            else                         {energy = log10((*energyVec)[i]);}
            // std::cout << "ECEnergyVector[" << i << "] = " << (*energyVec)[i] << " , layer " <<  layer  << " , bar " <<  bar << " Calculated Energy " << energy << " , "<< layer % 2 << std::endl;
            if(layer % 2 == 0)
            {
                hXZ->SetBinContent(bar+1, 14-layer, energy);
                hYZ->SetBinContent(bar+1, 14-layer, -5);
                for(int j = 0 ; j<24 ; ++j)
                {
                    h3box->SetBinContent(bar+1 , j+1  , 14-layer , (*energyVec)[i]);
                }
                
            }
            else
            {
                hYZ->SetBinContent(bar+1, 14-layer, energy);
                hXZ->SetBinContent(bar+1, 14-layer, -5);
                for(int j = 0 ; j<24 ; ++j)
                {
                    h3box->SetBinContent(j+1 , bar+1, 14-layer ,(*energyVec)[i]);
                }
            }

        }


        hXZ->SetMinimum(minVal); 
        hXZ->SetMaximum(maxVal); 
        hYZ->SetMinimum(minVal); 
        hYZ->SetMaximum(maxVal); 
        h3box->SetMinimum(1e-3); 
        h3box->SetMaximum(12); 
        c1->Clear();
        TLatex latex;
        latex.SetTextSize(0.03);
        latex.SetTextColor(kBlack);
        c1->Divide(3,1);
        c1->cd(1);
        gStyle->SetOptStat(0);
        hXZ->Draw("colz");
        gPad->Update();
        auto *palette1 = (TPaletteAxis*) hXZ->GetListOfFunctions()->FindObject("palette");
        if(palette1)
        {
            palette1->SetX1NDC(0.90);
            palette1->SetX2NDC(0.93);
            gPad->Modified();
            gPad->Update();
        }
        latex.DrawLatexNDC(0.82, 0.92, "log_{10}(E/GeV)");

        c1->cd(2);
        gStyle->SetOptStat(0);
        hYZ->Draw("colz");
        gPad->Update();
        auto  *palette2 = (TPaletteAxis*) hYZ->GetListOfFunctions()->FindObject("palette");
        if(palette2)
        {
            palette2->SetX1NDC(0.90);
            palette2->SetX2NDC(0.93);
            gPad->Modified();
            gPad->Update();
        }
        latex.DrawLatexNDC(0.82, 0.92, "log_{10}(E/GeV)");

        c1->cd(3);
        h3box->GetZaxis()->SetNdivisions(10);
        h3box->Draw("BOX2Z");
        auto  *palette3 = (TPaletteAxis*) h3box->GetListOfFunctions()->FindObject("palette");
        latex.DrawLatexNDC(0.8, 0.92, "(Edep/GeV)");
        // }

}
