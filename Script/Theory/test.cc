void test(){
    TH1* h = new TH1D("h", "h", 100, -5., 5.);
    gRandom->SetSeed();
    h->FillRandom("gaus", 1u << 16);
    // get the cumulative of h
    TH1* hc = h->GetCumulative(kFALSE);
    // check that c has the "right" contents
    Double_t* integral = h->GetIntegral();
    for (Int_t i = 1; i <= hc->GetNbinsX(); ++i) {
       assert(std::abs(integral[i] * h->GetEntries() - hc->GetBinContent(i)) < 1e-7);
    }
    // draw histogram together with its cumulative distribution
    TCanvas* c = new TCanvas;
    c->Divide(1,2);
    c->cd(1);
    h->Draw();
    c->cd(2);
    hc->Draw();
    c->Update();

 }