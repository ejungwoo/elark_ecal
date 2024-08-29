void ha()
{
    auto file = new TFile("data/stark_0113.hist.root");
    file -> ls();
    auto set = (LKDrawingSet*) file -> Get("DrawingSet");

    //auto hist = set -> FindHist("hist_d31_j_s3_rpos_esum");
    //auto hist = set -> FindHist("hist_j_det_j_esum");
    //auto hist = set -> FindHist("hist_d0_j_s0_g1_left_right");
    auto hist = set -> FindHist("hist_d1_j_s0_left_right");
    cout << hist << endl;
    //hist -> Draw();

    //auto drawing = set -> FindDrawing("drawing_hist_d0_j_s0_g1_left_right");
    //auto drawing = set -> FindDrawing("drawing_hist_d28_j_s5_left_right");
    //cout << drawing << endl;

    (new LKDataViewer(set)) -> Draw();
    //auto drawing = set -> FindDrawing("drawing_hist_d31_j_s3_rpos_esum");
    //auto drawing = set -> FindDrawing("drawing_hist_j_det_j_esum");
    //drawing -> GetMainHist() -> Draw();
    //cout << drawing<< endl;
    //drawing -> Get
    //set -> Print();
    //auto a = set -> FindObject("hist_d31_j_s6_left_right");
    //cout << a << endl;
}
