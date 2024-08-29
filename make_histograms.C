#include "si_calibration.h"
#include "energy_calibration.h"

void make_histograms(int calibration=0, int run=-1, bool drawExample=true)
{
    MakeRun(run);
    int ne = fNBinA;
    double e1 = fBinA1;
    double e2 = fBinA2;

    /////////////////////////////////////////////////////////////////////
    // 1) Get calibration parameter if calibration is 1 (slope correction
    /////////////////////////////////////////////////////////////////////
    TString calName = "";
    if (calibration==1) {
        calName = "c1_";
        GetC1Parameters();
        GetC0Parameters();
        ne = fNBinE;
        e1 = fBinE1;
        e2 = fBinE2;
    }

    /////////////////////////////////////////////////////////////////////
    // 2) Create histograms
    /////////////////////////////////////////////////////////////////////
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = MakeHist1(calName+"energy","",dss.det,dss.side,dss.strip,-1);
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = MakeHist1(calName+"esum",             "",dss.det,dss.side,dss.strip,-1,ne,e1,e2);
        fHistLeftRight     [dss.det][dss.side][dss.strip] = MakeHist2(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1,ne,e1,e2,ne,e1,e2);
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1,fNBinX,fBinX1,fBinX2,ne,e1,1.5*e2);
        for (auto gate=0; gate<fNumGates; ++gate) {
            fHistLeftRightGate     [dss.det][dss.side][dss.strip][gate] = MakeHist2(calName+"left",calName+"right",dss.det,dss.side,dss.strip,gate,ne,e1,e2,ne,e1,e2);
            fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,gate,fNBinX,fBinX1,fBinX2,ne,e1,1.5*e2);
        }
    }
    fHistEnergyDetector[0] = MakeHist2("det_j",calName+"esum",-1,0,-1,-1,40*10,0,40*10,ne,e1,e2);
    fHistEnergyDetector[1] = MakeHist2("det_o",calName+"esum",-1,1,-1,-1,40*10,0,40*10,ne,e1,e2);

    /////////////////////////////////////////////////////////////////////
    // 3) Fill histogram from stark event tree
    /////////////////////////////////////////////////////////////////////
    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);
    //auto numEvents = tree -> GetEntries();
    auto numEvents = 100000;
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0) cout << "Filling raw histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            auto strip = channel -> GetStrip();
            if (ContinueRegardingToDataType(det)) continue;
            if (!IsPositionSensitiveStrip(channel))
            {
                double g0 = (calibration==1) ? fC0Parameters[det][side][strip] : 1.;
                auto energy = channel -> GetEnergy() * g0;
                fHistEnergy[det][side][strip] -> Fill(energy);
                fHistEnergyDetector[side] -> Fill(det*10+strip,energy);
            }
            else
            {
                double g1 = (calibration==1) ? fC1Parameters[det][0][strip][0] : 1.;
                double g2 = (calibration==1) ? fC1Parameters[det][0][strip][1] : 1.;
                auto energy1 = channel -> GetEnergy() * g1;
                auto energy2 = channel -> GetEnergy2() * g2;
                if (energy2>0) {
                    auto sum = energy1 + energy2;
                    auto pos = (energy1 - energy2) / sum;
                    fHistEnergyDetector[side] -> Fill(det*10+strip,sum);
                    fHistEnergySum     [det][side][strip] -> Fill(sum);
                    fHistLeftRight     [det][side][strip] -> Fill(energy1, energy2);
                    fHistEnergyPosition[det][side][strip] -> Fill(pos,sum);
                    if (calibration==1)
                    {
                        for (auto gate=0; gate<fNumGates; ++gate) {
                            double eAlpha = (gate==0?f148GdAlphaEnergy:f241AmAlphaEnergy1);
                            double range1 = eAlpha - 6*eAlpha*fExpectedResolution;
                            double range2 = eAlpha + 6*eAlpha*fExpectedResolution;
                            if (sum>range1 && sum<range2) {
                                fHistLeftRightGate     [det][side][strip][gate] -> Fill(energy1, energy2);
                                fHistEnergyPositionGate[det][side][strip][gate] -> Fill(pos,sum);
                            }
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 5) Calculate energy gate (0:Gd, 1:Am) for alpha source
    /////////////////////////////////////////////////////////////////////
    double gatingRange[40][8][fNumGates][2] = {0}; // det, strip, gate, range
    if (calibration==0)
    {
        auto spectrum = new TSpectrum(fNumGates*2);
        TF1 *fitGaus = new TF1("fitGaus","gaus(0)",0,6000);
        for (auto dss : fStripArrayR)
        {
            auto hist = fHistEnergySum[dss.det][dss.side][dss.strip];
            auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
            double* xPeaks = spectrum -> GetPositionX();
            if (numPeaks<fMaxPeaks) {
                e_warning << hist->GetName() << " #peaks =" << numPeaks << endl;
                //continue;
            }
            if (xPeaks[1]<xPeaks[0]) {
                auto xx = xPeaks[1];
                xPeaks[1] = xPeaks[0];
                xPeaks[0] = xx;
            }
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
            {
                if (iPeak+1>numPeaks) {
                    e_warning << hist -> GetName() << " (" << dss.det << ", " << dss.side << ", " << dss.strip << ") : " << " #peaks = " << numPeaks << " #entries = " << hist -> GetEntries() << endl;
                    gatingRange[dss.det][dss.strip][iPeak][0] = 0;
                    gatingRange[dss.det][dss.strip][iPeak][1] = 0;
                    continue;
                }
                double xPeak = xPeaks[iPeak];
                fitGaus -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
                fitGaus -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
                hist -> Fit(fitGaus,"Q0NR");
                gatingRange[dss.det][dss.strip][iPeak][0] = fitGaus -> GetParameter(1) - 5 * fitGaus -> GetParameter(2);
                gatingRange[dss.det][dss.strip][iPeak][1] = fitGaus -> GetParameter(1) + 7 * fitGaus -> GetParameter(2);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 6) Fill histogram with energy gate
    /////////////////////////////////////////////////////////////////////
    if (calibration==0)
    {
        for (auto iEvent=0; iEvent<numEvents; ++iEvent)
        {
            tree -> GetEntry(iEvent);
            auto numChannels = array -> GetEntries();
            if (iEvent%20000==0) cout << "Filling gate histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
            if (numChannels!=3) continue;
            for (auto iChannel=0; iChannel<numChannels; ++iChannel)
            {
                auto channel = (LKSiChannel*) array -> At(iChannel);
                auto det = channel -> GetDetID();
                auto side = channel -> GetSide();
                auto strip = channel -> GetStrip();
                if (ContinueRegardingToDataType(det)) continue;
                if (IsPositionSensitiveStrip(channel))
                {
                    double g1 = (calibration==1) ? fC1Parameters[det][0][strip][0] : 1.;
                    double g2 = (calibration==1) ? fC1Parameters[det][0][strip][1] : 1.;
                    auto energy1 = channel -> GetEnergy() * g1;
                    auto energy2 = channel -> GetEnergy2() * g2;
                    if (energy2>0) {
                        auto sum = energy1 + energy2;
                        auto pos = (energy1 - energy2) / sum;
                        for (auto gate=0; gate<fNumGates; ++gate) {
                            if (sum>gatingRange[det][strip][gate][0] && sum<gatingRange[det][strip][gate][1]) {
                                fHistLeftRightGate     [det][side][strip][gate] -> Fill(energy1, energy2);
                                fHistEnergyPositionGate[det][side][strip][gate] -> Fill(pos,sum);
                            }
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // x) Draw examples
    /////////////////////////////////////////////////////////////////////
    auto set = new LKDrawingSet();
    //auto set2 = new LKDrawingSet();
    //set -> CreateGroup("HPJ") -> AddHist(fHistEnergyDetector[0]);
    //set -> CreateGroup("HPO") -> AddHist(fHistEnergyDetector[1]);
    //auto group0 = new LKDrawingGroup("Energy"); set -> AddGroup(group0);
    //for (auto dssGroup : fAllGroupArrayS)
    //{
    //    auto sub = group0 -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
    //    for (auto dss : dssGroup.array)
    //        sub -> AddHist(fHistEnergy[dss.det][dss.side][dss.strip]);
    //}
    auto group1 = set -> CreateGroup("EnergySum");
    auto group2 = set -> CreateGroup("LeftRight");
    auto group3 = set -> CreateGroup("EnergyPosition");
    auto group4 = set -> CreateGroup("LeftRightGate");
    auto group5 = set -> CreateGroup("EnergyPositionGate");
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto sub1 = group1 -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
        auto sub2 = group2 -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
        auto sub3 = group3 -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
        auto sub4 = group4 -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
        auto sub5 = group5 -> CreateSubGroup(Form("%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            sub1 -> AddHist(fHistEnergySum     [dss.det][dss.side][dss.strip]);
            sub2 -> AddHist(fHistLeftRight     [dss.det][dss.side][dss.strip]);
            sub3 -> AddHist(fHistEnergyPosition[dss.det][dss.side][dss.strip]);
            for (auto gate=0; gate<fNumGates; ++gate) {
                sub4 -> AddHist(fHistLeftRightGate     [dss.det][dss.side][dss.strip][gate]);
                sub5 -> AddHist(fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate]);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 7) Write histograms
    /////////////////////////////////////////////////////////////////////
    TString foutName = fHistFileName;
    if (calibration==1)
        foutName = fC1HistFileName;
    auto fileHist = new TFile(foutName,"recreate");
    //for (auto dss : fStripArrayS)
    //{
    //    fHistEnergy[dss.det][dss.side][dss.strip] -> Write();
    //}
    //for (auto dss : fStripArrayR)
    //{
    //    fHistEnergySum     [dss.det][dss.side][dss.strip] -> Write();
    //    fHistLeftRight     [dss.det][dss.side][dss.strip] -> Write();
    //    fHistEnergyPosition[dss.det][dss.side][dss.strip] -> Write();
    //    for (auto gate=0; gate<fNumGates; ++gate) {
    //        fHistLeftRightGate     [dss.det][dss.side][dss.strip][gate] -> Write();
    //        fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate] -> Write();
    //    }
    //}
    //fHistEnergyDetector[0] -> Write();
    //fHistEnergyDetector[1] -> Write();
    //fileHist -> cd();
    //set -> Write("drawings",TObject::kSingleKey);
    //set -> Print();
    set -> Write();//"drawings",TObject::kSingleKey);
    //set2 -> Write();//"drawings",TObject::kSingleKey);

    cout << foutName << endl;
    //(new LKDataViewer(set)) -> Draw();
}
