#include "si_calibration.h"
#include "energy_calibration.h"

void run_energy_calibration(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    /////////////////////////////////////////////////////////////////////
    // 1) Get histograms
    /////////////////////////////////////////////////////////////////////
    auto fileHist = new TFile(fHistFileName,"read");
    auto setInput = (LKDrawingSet*) fileHist -> Get("DrawingSet");
    for (auto dss : fStripArrayS) {
        //fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) fileHist -> Get(MakeHistName("energy","",dss.det,dss.side,dss.strip));
        fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) setInput -> FindHist(MakeHistName("energy","",dss.det,dss.side,dss.strip));
    }

    /////////////////////////////////////////////////////////////////////
    // 2) Calculate energy calibration parameter
    /////////////////////////////////////////////////////////////////////
    int det, side, strip;
    double entries, a1, m1, s1, a2, m2, s2, g0;
    auto fout = new TFile(fC0ParFileName,"recreate");
    auto tree = new TTree("parameters","energy calibration parameter g0: adc * g0 = energy");
    tree -> Branch("det"    ,&det    );
    tree -> Branch("side"   ,&side   );
    tree -> Branch("strip"  ,&strip  );
    tree -> Branch("entries",&entries);
    tree -> Branch("a1"     ,&a1     );
    tree -> Branch("m1"     ,&m1     );
    tree -> Branch("s1"     ,&s1     );
    tree -> Branch("a2"     ,&a2     );
    tree -> Branch("m2"     ,&m2     );
    tree -> Branch("s2"     ,&s2     );
    tree -> Branch("g0"     ,&g0     );
    double gausParameters[40][2][8][3][3] = {0}; // det, strip, gate, range
    auto spectrum = new TSpectrum(fMaxPeaks);
    TF1 *fitGaus = new TF1("fitGaus","gaus(0)",0,6000);
    for (auto dss : fStripArrayS)
    {
        det = dss.det;
        side = dss.side;
        strip = dss.strip;
        auto hist = fHistEnergy[det][side][strip];
        entries = hist -> GetEntries();
        a1 = 0;
        m1 = 0;
        s1 = 0;
        a2 = 0;
        m2 = 0;
        s2 = 0;
        g0 = 0;
        if (entries<fEntriesCut)
        {
            e_warning << hist->GetName() << " entries = " << entries << endl;
            tree -> Fill();
            continue;
        }
        auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
        double* xPeaks = spectrum -> GetPositionX();
        if (numPeaks<fMaxPeaks) {
            e_warning << hist->GetName() << " #peaks =" << numPeaks << endl;
            tree -> Fill();
            continue;
        }
        if (xPeaks[1]<xPeaks[0]) {
            auto xx = xPeaks[1];
            xPeaks[1] = xPeaks[0];
            xPeaks[0] = xx;
        }
        for (auto iPeak=0; iPeak<fMaxPeaks; ++iPeak)
        {
            double xPeak = xPeaks[iPeak];
            fitGaus -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
            fitGaus -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
            hist -> Fit(fitGaus,"Q0NR");
            auto amp =   fitGaus -> GetParameter(0);
            auto mean =  fitGaus -> GetParameter(1);
            auto sigma = fitGaus -> GetParameter(2);
            fitGaus -> SetRange(mean-1.0*sigma, mean+2.5*sigma);
            hist -> Fit(fitGaus,"Q0NR");
            amp =   fitGaus -> GetParameter(0);
            mean =  fitGaus -> GetParameter(1);
            sigma = fitGaus -> GetParameter(2);
            gausParameters[det][side][strip][iPeak][0] = amp;
            gausParameters[det][side][strip][iPeak][1] = mean;
            gausParameters[det][side][strip][iPeak][2] = sigma;
            if (iPeak==0) {
                a1 = amp;
                m1 = mean;
                s1 = sigma;
            }
            if (iPeak==1) {
                a2 = amp;
                m2 = mean;
                s2 = sigma;
                g0 = f241AmAlphaEnergy1/mean;
            }
        }
        tree -> Fill();
    }
    fout -> cd();
    tree -> Write();
    fout -> Close();
    cout << fC0ParFileName << endl;

    /////////////////////////////////////////////////////////////////////
    // 3) Draw examples
    /////////////////////////////////////////////////////////////////////
    auto set = new LKDrawingSet();
    auto group = set -> CreateGroup("energy");
    //for (auto dssGroup : fExampleGroupArrayS)
    for (auto dssGroup : fAllGroupArrayS)
    {
        auto sub = group -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            auto drawing = sub -> CreateDrawing("");
            drawing -> Add(fHistEnergy[dss.det][dss.side][dss.strip]);
            //sub -> AddHist(fHistEnergy[dss.det][dss.side][dss.strip]);
            //auto lg = new TLegend(0.4,0.55,0.9,0.88);
            auto lg = new TLegend(0.35,0.60,0.85,0.88);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            lg -> SetNColumns(fMaxPeaks);
            for (auto iPeak=0; iPeak<fMaxPeaks; ++iPeak)
            {
                auto amp = gausParameters[dss.det][dss.side][dss.strip][iPeak][0];
                auto mean = gausParameters[dss.det][dss.side][dss.strip][iPeak][1];
                auto sigma = gausParameters[dss.det][dss.side][dss.strip][iPeak][2];
                if (amp==0)
                    continue;
                TF1 *fitGaus = new TF1(Form("fit%d%d%d%d",dss.det,dss.side,dss.strip,iPeak),"gaus(0)",0,6000);
                fitGaus -> SetParameter(0,amp);
                fitGaus -> SetParameter(1,mean);
                fitGaus -> SetParameter(2,sigma);
                fitGaus -> SetRange(mean-1.0*sigma, mean+2.5*sigma);
                drawing -> Add(fitGaus,"","samel");
                //fitGaus -> Draw("samel");
                lg -> AddEntry((TObject*)nullptr,Form("a_{%d} = %.1f",iPeak,amp),"");
                lg -> AddEntry((TObject*)nullptr,Form("m_{%d} = %.1f",iPeak,mean),"");
                lg -> AddEntry((TObject*)nullptr,Form("s_{%d} = %.1f",iPeak,sigma),"");
                lg -> AddEntry((TObject*)nullptr,Form("x_{%d} = %.1f",iPeak,fitGaus->GetChisquare()),"");
            }
            //lg -> Draw();
            drawing -> Add(lg);
            //icvs++;
        }
    }
    set -> Draw("viewer");
}
