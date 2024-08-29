#include "si_calibration.h"
#include "energy_calibration.h"

TH1D* histEnergy        [fNumDetectors][2][fNumStrips];
TH1D* histEnergySum     [fNumDetectors][2][fNumStrips];
TH2D* histLeftRight     [fNumDetectors][2][fNumStrips];
TH2D* histEnergyPosition[fNumDetectors][2][fNumStrips];
TH2D* histEnergyDetector[2][2];
TH2D* histEnergyDetectorAll[2];

double g0Array[40][2][8];
double g1Array[40][2][8];
double g2Array[40][2][8];
double b0Array[40][2][8];
double b1Array[40][2][8];
double b2Array[40][2][8];
double resolutionArray[40][2][8];

void check(int det=0, bool drawFigures=true);
void get_parameters();
void write_parameters();

void summary(int run1=199, int run2=303)
{
    TString calName = "c2_";

    for (auto det=0; det<40; ++det)
        for (auto side=0; side<2; ++side)
            for (auto strip=0; strip<8; ++strip) {
                g0Array[det][side][strip] = 0;
                g1Array[det][side][strip] = 0;
                g2Array[det][side][strip] = 0;
                b0Array[det][side][strip] = 0;
                b1Array[det][side][strip] = 0;
                b2Array[det][side][strip] = 0;
            }

    /////////////////////////////////////////////////////////////////////
    // 199
    /////////////////////////////////////////////////////////////////////
    if (run1>0)
    {
        MakeRun(run1);
        get_parameters();
        auto file1 = new TFile(Form("data/stark_0%d.hist_c2.root",fRun),"read");
        auto set1 = (LKDrawingSet*) file1 -> Get("DrawingSet");
        for (auto dss : fStripArrayS)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            histEnergy[dss.det][dss.side][dss.strip] = (TH1D*) set1 -> FindHist((MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
        }
        for (auto dss : fStripArrayR)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            histEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) set1 -> FindHist((MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
            histLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) set1 -> FindHist((MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
            histEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) set1 -> FindHist((MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
        }
        histEnergyDetector[0][0] = (TH2D*) set1 -> FindHist((MakeHistName("det_j",calName+"esum",-1,0,-1,-1));
        histEnergyDetector[0][1] = (TH2D*) set1 -> FindHist((MakeHistName("det_o",calName+"esum",-1,1,-1,-1));
    }

    /////////////////////////////////////////////////////////////////////
    // 303
    /////////////////////////////////////////////////////////////////////
    if (run2>0)
    {
        MakeRun(run2);
        get_parameters();
        auto file2 = new TFile(Form("data/stark_0%d.hist_c2.root",fRun),"read");
        auto set2 = (LKDrawingSet*) file2 -> Get("DrawingSet");
        for (auto dss : fStripArrayS)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            histEnergy[dss.det][dss.side][dss.strip] = (TH1D*) set2 -> FindHist((MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
        }
        for (auto dss : fStripArrayR)
        {
            if (ContinueRegardingToDataType(dss.det)) continue;
            histEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) set2 -> FindHist((MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
            histLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) set2 -> FindHist((MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
            histEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) set2 -> FindHist((MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
        }
        histEnergyDetector[1][0] = (TH2D*) set2 -> FindHist((MakeHistName("det_j",calName+"esum",-1,0,-1,-1));
        histEnergyDetector[1][1] = (TH2D*) set2 -> FindHist((MakeHistName("det_o",calName+"esum",-1,1,-1,-1));
    }

    /////////////////////////////////////////////////////////////////////
    // merge hit pattern
    /////////////////////////////////////////////////////////////////////
    histEnergyDetectorAll[0] = (TH2D*) histEnergyDetector[0][0] -> Clone("h0");
    histEnergyDetectorAll[1] = (TH2D*) histEnergyDetector[0][1] -> Clone("h1");
    histEnergyDetectorAll[0] -> Add(histEnergyDetector[1][0]);
    histEnergyDetectorAll[1] -> Add(histEnergyDetector[1][1]);
    histEnergyDetectorAll[0] -> SetStats(0);
    histEnergyDetectorAll[1] -> SetStats(0);
    histEnergyDetectorAll[0] -> GetXaxis() -> SetTitleOffset(2.7);
    histEnergyDetectorAll[1] -> GetXaxis() -> SetTitleOffset(2.7);
    for (auto det=0; det<fNumDetectors; ++det)
    {
        auto detector = fStark -> GetSiDetector(det);
        auto name = detector -> GetDetTypeName();
        auto ring = detector -> GetLayer();
        TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E"; 
        histEnergyDetectorAll[0] -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        histEnergyDetectorAll[1] -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
    }
    auto cvs = MakeCanvas("cvs_energy_detector",2);
    auto pad1 = cvs -> cd(1); pad1 -> SetMargin(0.1,0.12,0.18,0.1); histEnergyDetectorAll[0] -> Draw("colz");
    auto pad2 = cvs -> cd(2); pad2 -> SetMargin(0.1,0.12,0.18,0.1); histEnergyDetectorAll[1] -> Draw("colz");
    cvs -> SaveAs(Form("figures/%s.png",cvs->GetName()));

    /////////////////////////////////////////////////////////////////////
    // Write parameters
    /////////////////////////////////////////////////////////////////////
    for (auto i=0; i<40; ++i) check(i,false);
    write_parameters();
    //for (auto i=0; i<10; ++i) check(i,true);
    //for (auto i=32; i<38; ++i) check(i,true);
}

void check(int det, bool drawFigures)
{
    /////////////////////////////////////////////////////////////////////
    // create group
    /////////////////////////////////////////////////////////////////////
    strip_group dssArrayR;
    strip_group dssArrayS;
    auto detector = fStark -> GetSiDetector(det);
    auto numJStrips = detector -> GetNumJunctionStrips();
    auto numOStrips = detector -> GetNumOhmicStrips();
    auto detType = detector -> GetDetType();
    for (auto side=0; side<2; ++side)
    {
        if (IsPositionSensitiveStrip(det,side))
        {
            auto numStrips = (side==0?numJStrips:numOStrips);
            for (auto strip=0; strip<numStrips; ++strip) {
                dssArrayR.array.push_back(strip_info(detType,det,side,strip));
            }
        }
        else {
            auto numStrips = (side==0?numJStrips:numOStrips);
            for (auto strip=0; strip<numStrips; ++strip) {
                dssArrayS.array.push_back(strip_info(detType,det,side,strip));
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // create canvas
    /////////////////////////////////////////////////////////////////////
    auto spectrum = new TSpectrum(fNumGates*2);

    auto fitAndDraw = [drawFigures,spectrum](strip_info dss, TH1D* hist)
    {
        auto lg = new TLegend(0.35,0.75,0.80,0.88);
        lg -> SetBorderSize(0);
        lg -> SetFillStyle(0);
        TF1 *fits[2];
        fits[0] = nullptr;
        fits[1] = nullptr;
        double x1=-999, x2=-999;
        double m1=-999, m2=-999;
        double s1=-999, s2=-999;
        if (hist->GetEntries()>fEntriesCut)
        {
            auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
            double* xPeaks = spectrum -> GetPositionX();
            if (xPeaks[1]<xPeaks[0]) {
                auto xx = xPeaks[1];
                xPeaks[1] = xPeaks[0];
                xPeaks[0] = xx;
            }
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
            {
                if (iPeak+1>numPeaks) {
                    e_warning << hist -> GetName() << " (" << dss.det << ", " << dss.side << ", " << dss.strip << ") : " << " #peaks = " << numPeaks << " #entries = " << hist -> GetEntries() << endl;
                    continue;
                }
                double xPeak = xPeaks[iPeak];
                fits[iPeak] = new TF1(Form("fitGaus_%d_%d_%d_%d",dss.det,dss.side,dss.strip,iPeak),"gaus(0)",0,6000);
                auto fit = fits[iPeak];
                fit -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
                fit -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
                hist -> Fit(fit,"Q0NR");
                auto amp = fit -> GetParameter(0);
                auto mean = fit -> GetParameter(1);
                auto sigma = fit -> GetParameter(2);
                fit -> SetRange(mean-2.5*sigma,mean+2.5*sigma);
                if (dss.det==34||dss.det==35||dss.det==36)
                    fit -> SetRange(mean-.8*sigma,mean+2.5*sigma);
                hist -> Fit(fit,"Q0NR");
                amp = fit -> GetParameter(0);
                mean = fit -> GetParameter(1);
                sigma = fit -> GetParameter(2);
                if (iPeak==0) { x1 = mean - 20*sigma; m1 = mean; s1 = sigma; }
                if (iPeak==1) { x2 = mean + 20*sigma; m2 = mean; s2 = sigma; }
                //fit -> SetRange(mean-2.5*sigma,mean+2.5*sigma);
                lg -> AddEntry((TObject*)nullptr,Form("mean_{%d} = %.3f",iPeak,mean),"");
            }
            resolutionArray[dss.det][dss.side][dss.strip] = 100*s1/m1;
        }
        if (drawFigures)
        {
            if (x1>-998&&x2>-998)
                hist -> GetXaxis() -> SetRangeUser(x1,x2);
            hist -> Draw();
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak) {
                if (fits[iPeak]!=nullptr)
                    fits[iPeak] -> Draw("samel");
            }
            lg -> Draw("same");
        }
    };

    /////////////////////////////////////////////////////////////////////
    // create canvas
    /////////////////////////////////////////////////////////////////////
    int icvs = 1;
    TCanvas *cvs = nullptr;
    if (drawFigures) {
        cvs = dssArrayS.MakeGroupCanvas(Form("summary_det%d",det),20);
        for (auto dss : dssArrayR.array) { cvs -> cd(icvs++); histEnergyPosition[dss.det][dss.side][dss.strip] -> Draw("colz"); }
    }

    for (auto dss : dssArrayR.array)
    {
        if (drawFigures) cvs -> cd(icvs++);
        auto hist = histEnergySum[dss.det][dss.side][dss.strip];
        fitAndDraw(dss,hist);
    }

    for (auto dss : dssArrayS.array)
    {
        if (drawFigures) cvs -> cd(icvs++);
        auto hist = histEnergy[dss.det][dss.side][dss.strip];
        fitAndDraw(dss,hist);
    }

    if (drawFigures) cvs -> SaveAs(Form("figures/%s.png",cvs->GetName()));
}

void get_parameters()
{
    int det, side, strip;
    double entries, g0, g1, g2, b0, b1, b2;

    auto fpar0 = new TFile(fC0ParFileName,"read");
    auto tree0 = (TTree*) fpar0 -> Get("parameters");
    tree0 -> SetBranchAddress("det"    ,&det    );
    tree0 -> SetBranchAddress("side"   ,&side   );
    tree0 -> SetBranchAddress("strip"  ,&strip  );
    tree0 -> SetBranchAddress("g0"     ,&g0     );
    auto n0 = tree0 -> GetEntries();
    for (auto i=0; i<n0; ++i) {
        tree0 -> GetEntry(i);
        g0Array[det][side][strip] = g0;
    }
    fpar0 -> Close();

    auto fpar1 = new TFile(fC1ParFileName,"read");
    auto tree1 = (TTree*) fpar1 -> Get("parameters");
    tree1 -> SetBranchAddress("det"    ,&det    );
    tree1 -> SetBranchAddress("side"   ,&side   );
    tree1 -> SetBranchAddress("strip"  ,&strip  );
    tree1 -> SetBranchAddress("g1"     ,&g1     );
    tree1 -> SetBranchAddress("g2"     ,&g2     );
    auto n1 = tree1 -> GetEntries();
    for (auto i=0; i<n1; ++i) {
        tree1 -> GetEntry(i);
        g1Array[det][side][strip] = g1;
        g2Array[det][side][strip] = g2;
    }
    fpar1 -> Close();

    auto fpar2 = new TFile(fC2ParFileName,"read");
    auto tree2 = (TTree*) fpar2 -> Get("parameters");
    tree2 -> SetBranchAddress("det"    ,&det    );
    tree2 -> SetBranchAddress("side"   ,&side   );
    tree2 -> SetBranchAddress("strip"  ,&strip  );
    tree2 -> SetBranchAddress("b0"     ,&b0     );
    tree2 -> SetBranchAddress("b1"     ,&b1     );
    tree2 -> SetBranchAddress("b2"     ,&b2     );
    auto n2 = tree2 -> GetEntries();
    for (auto i=0; i<n2; ++i) {
        tree2 -> GetEntry(i);
        b0Array[det][side][strip] = b0;
        b1Array[det][side][strip] = b1;
        b2Array[det][side][strip] = b2;
    }
    fpar2 -> Close();
}

void write_parameters()
{
    int det, side, strip;
    double entries, g0, g1, g2, b0, b1, b2, resolution;
    auto fpar = new TFile(fAllParFileName,"recreate");

    (new TParameter<double>("calibration_energy",f241AmAlphaEnergy1)) -> Write();;

    auto tree = new TTree("parameters","all (c0,c1,c2)");
    tree -> Branch("det"    ,&det    );
    tree -> Branch("side"   ,&side   );
    tree -> Branch("strip"  ,&strip  );
    tree -> Branch("g1"     ,&g1     );
    tree -> Branch("g2"     ,&g2     );
    tree -> Branch("g0"     ,&g0     );
    tree -> Branch("b0"     ,&b0     );
    tree -> Branch("b1"     ,&b1     );
    tree -> Branch("b2"     ,&b2     );
    tree -> Branch("resolution",&resolution);

    for (auto vv : {fAllGroupArrayR, fAllGroupArrayS})
    {
        for (auto dssGroup : vv)
        {
            for (auto dss : dssGroup.array)
            {
                det = dss.det;
                side = dss.side;
                strip = dss.strip;
                //////////////////////////////
                double user_gain = 1.;
                g0 = user_gain * g0Array[det][side][strip];
                g1 = user_gain * g1Array[det][side][strip];
                g2 = user_gain * g2Array[det][side][strip];
                b0 = user_gain * b0Array[det][side][strip];
                b1 = user_gain * b1Array[det][side][strip];
                b2 = user_gain * b2Array[det][side][strip];
                //////////////////////////////
                resolution = resolutionArray[det][side][strip];
                tree -> Fill();
            }
        }
    }

    fpar -> cd();
    tree -> Write();
    tree -> Print("toponly");
    cout << fpar -> GetName() << endl;
    fpar -> Close();
}
