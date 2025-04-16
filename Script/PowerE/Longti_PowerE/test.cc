void test(){
    int nbins = 80;
    double xmin = 1e-4;
    double xmax = 1;

    double logxmin = log10(xmin);
    double logxmax = log10(xmax);

    // 构造 log 等距 bin 边界
    std::vector<double> bin_edges(nbins + 1);
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * (logxmax - logxmin) / nbins);
    }

    // 创建直方图，传入 bin 边界数组
    TH1D* h = new TH1D("h", "Log-binned histogram;X axis;Counts", nbins, bin_edges.data());

    // // 设置画布 logx
    // gPad->SetLogx();

    // 随便填点数据测试（你也可以用你的实际数据）
    for (int i = 0; i < 10000; ++i) {
        double x = pow(10, gRandom->Uniform(logxmin, logxmax)); // log-uniform 随机数
        h->Fill(x);
    }

    h->Draw();
}