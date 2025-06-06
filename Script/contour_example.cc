void contour_example() {
    // Dummy data
    const int n = 10;
    double x[n] = {1,2,3,4,5,6,7,8,9,10};
    double y[n] = {2.1,4.0,6.1,8.3,10.1,12.1,14.0,16.2,18.2,20.1};

    // Define chi2 function
    auto chi2 = [&](int& npar, double* grad, double& fval, double* par, int flag) {
        double a = par[0]; // slope
        double b = par[1]; // intercept
        fval = 0;
        for (int i = 0; i < n; ++i) {
            double yi = a * x[i] + b;
            fval += pow(y[i] - yi, 2); // chi2
        }
        return fval;
    };

    TMinuit* minuit = new TMinuit(2);
    minuit->SetFCN(chi2);
    minuit->DefineParameter(0, "a", 1.0, 0.1, 0, 0);
    minuit->DefineParameter(1, "b", 0.0, 0.1, 0, 0);
    minuit->Migrad();

    // Create contour plot
    const int npoints = 100;
    TGraph* gr = (TGraph*)minuit->Contour(npoints, 0, 1); // param 0 vs param 1
    TCanvas* c1 = new TCanvas("c1", "Contour", 800, 600);
    gr->SetTitle("Contour of #chi^{2};a; b");
    gr->Draw("AL");
}