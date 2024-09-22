void draw()
{
    int run = 199;
    //int run = 303;
    //int run = 113;
    //int run = 167;
    //int run = 168;
    //int run = 7728;
    //int run = 169; // 0, 1, 10, 11, 12
    //int run = 170; // 7, 8, 9
    //int run = 171; // 1, 2, 3

    auto top = new LKDrawingGroup();
    top -> AddFile(Form("data/stark_%04d.hist_c3.root",run));
    //top -> AddFile("data/compare_168_with_199.hist_c3.root");
    top -> Draw("v");
}
