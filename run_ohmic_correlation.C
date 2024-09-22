#include "si_calibration.h"

void run_ohmic_correlation(bool runViewer=true, int calRun=-1)
{
    bool checkCECorrection = false;

    auto top = new LKDrawingGroup(Form("run_%d%s",fRun,checkCECorrection?"_SC":""));

    MakeRun();
    int ne = fNBinA;
    double e1 = fBinA1;
    double e2 = fBinA2;

    int nx = 100;
    double x1 = fBinX1;
    double x2 = fBinX2;

    e_info << fGateEnergy[0] << endl;
    e_info << fGateEnergy[0]/2. << endl;
    e_info << fGateEnergy[1] << endl;
    e_info << fGateEnergy[1]/2. << endl;

    /////////////////////////////////////////////////////////////////////
    // 1) Get calibration parameter if calibration is 1 (slope correction
    /////////////////////////////////////////////////////////////////////
    TString calName = "";
    calName = "oc_";
    GetC1Parameters(calRun);
    if (checkCECorrection) GetCEParameters(fRun);
    ne = fNBinE;
    e1 = fBinE1;
    //e2 = fBinE2;
    e2 = 6;
    nx = 100;
    x1 = e1;
    x2 = e2;
    int nbinsx = 50;

    /////////////////////////////////////////////////////////////////////
    // 2) Create histograms
    /////////////////////////////////////////////////////////////////////
    for (auto dss : fStripArrayR)
    {
        for (auto gate=0; gate<fNumGates; ++gate)
        {
            for (auto oStrip=0; oStrip<fNumOStrips; ++oStrip)
            {
                TString calName2 = calName + "os" + oStrip + "_";
                fHistLeftRightGateOhmic     [dss.det][dss.side][dss.strip][gate][oStrip] = MakeHist2(calName2+"left",calName2+"right",dss.det,dss.side,dss.strip,gate,ne,e1,e2,ne,e1,e2);
                fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip] = MakeHist2(calName2+"rpos",calName2+"esum", dss.det,dss.side,dss.strip,gate,fNBinX,fBinX1,fBinX2,ne,e1,1.5*e2); // 1
                fHistLRProjGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip][0] = MakeHist1(calName2+"LEnergy","", dss.det,dss.side,dss.strip,gate,nx,x1,x2);
                fHistLRProjGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip][1] = MakeHist1(calName2+"REnergy","", dss.det,dss.side,dss.strip,gate,nx,x1,x2);
            }
            fHistPositionGate[dss.det][dss.side][dss.strip][gate] = MakeHist1(calName+"rpos","", dss.det,dss.side,dss.strip,gate,nbinsx,fBinX1,fBinX2);
        }
    }
    auto histDummy = MakeHist1("dmx","dmy", 0,0,0,0,nx,x1,x2);
    auto gausDummy = new TF1("fit","gaus(0)",x1,x2);

    /////////////////////////////////////////////////////////////////////
    // 3) Fill histogram from stark event tree
    /////////////////////////////////////////////////////////////////////
    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);
    auto numEvents = tree -> GetEntries();
    if (fTest) numEvents = 100000;
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0) cout << "Filling raw histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
        int oStrip = -1;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            if (side==1) {
                oStrip = channel -> GetStrip();
                break;
            }
        }
        if (oStrip<0)
            continue;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            if (side==1)
                continue;
            auto strip = channel -> GetStrip();
            if (ContinueRegardingToDataType(det)) continue;
            if (!IsPositionSensitiveStrip(channel));
            else
            {
                auto energyR = channel -> GetEnergy();
                auto energyL = channel -> GetEnergy2();
                if (energyL<0)
                    continue;
                CalibrateC1(det,side,strip,0,energyL);
                CalibrateC1(det,side,strip,1,energyR);
                if (checkCECorrection && GetNumOhmicStrips(det)>1)
                {
                    CalibrateCE(det,side,strip,0,energyL);
                    CalibrateCE(det,side,strip,1,energyR);
                }
                auto sum = energyR + energyL;
                auto pos = (energyR - energyL) / sum;
                for (auto gate=0; gate<fNumGates; ++gate) {
                    double range1 = fEnergyGateRange[gate][0];
                    double range2 = fEnergyGateRange[gate][1];
                    if (sum>range1 && sum<range2) {
                        fHistLeftRightGateOhmic[det][side][strip][gate][oStrip] -> Fill(energyL,energyR);
                        fHistEnergyPositionGateOhmic[det][side][strip][gate][oStrip] -> Fill(pos,sum);
                        fHistLRProjGateOhmic[det][side][strip][gate][oStrip][0] -> Fill(energyL);
                        fHistLRProjGateOhmic[det][side][strip][gate][oStrip][1] -> Fill(energyR);
                        //fHistPositionGate[det][side][strip][gate] -> Fill(pos,sum);
                        fHistPositionGate[det][side][strip][gate] -> Fill(pos);
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    /*
    auto SetHistColor = [](TH2D* hist, int color)
    {
        auto nx = hist -> GetXaxis() -> GetNbins();
        auto x1 = hist -> GetXaxis() -> GetXmin();
        auto x2 = hist -> GetXaxis() -> GetXmax();
        auto ny = hist -> GetYaxis() -> GetNbins();
        auto y1 = hist -> GetYaxis() -> GetXmin();
        auto y2 = hist -> GetYaxis() -> GetXmax();
        for (auto ix=1; ix<=nx; ++ix) {
            for (auto iy=1; iy<=ny; ++iy) {
                auto value = hist -> GetBinContent(ix,iy);
                if (value>0) {
                    hist -> SetBinContent(ix,iy,color+1);
                    //if      (color==0) hist -> SetBinContent(ix,iy,value);
                    //else if (color==1) hist -> SetBinContent(ix,iy,value*100);
                    //else if (color==2) hist -> SetBinContent(ix,iy,value*10000);
                    //else if (color==3) hist -> SetBinContent(ix,iy,value*1000000);
                }
            }
        }
        hist -> SetMaximum(fNumOStrips+1);
        //hist -> SetMaximum(fNumOStrips*100000000);
    };
    */
    /////////////////////////////////////////////////////////////////////
    double threshold = 10;
    auto FindHistX = [threshold](TH1D* hist0, TH1D* hist1, TH1D* histX)
    {
        histX -> Reset();
        auto nx = hist0 -> GetXaxis() -> GetNbins();
        auto x1 = hist0 -> GetXaxis() -> GetXmin();
        auto x2 = hist0 -> GetXaxis() -> GetXmax();
        for (auto ix=1; ix<=nx; ++ix)
        {
            double value = 0;
            auto value0 = hist0 -> GetBinContent(ix);
            auto value1 = hist1 -> GetBinContent(ix);
            if (value0==0 || value1==0) ;//histX -> SetBinContent(ix,0);
            else if (value0<=value1) value = value0;
            else if (value0> value1) value = value1;
            if (value>threshold) {
                histX -> SetBinContent(ix,value);
            }
        }
    };

    /////////////////////////////////////////////////////////////////////
    // x) find histogram cross section
    /////////////////////////////////////////////////////////////////////
    if (!checkCECorrection) StartWriteCEParameters();
    for (auto dssGroup : fAllGroupArrayR)
    {
        for (auto dss : dssGroup.array)
        {
            fGraphEnergyFit[dss.det][dss.side][dss.strip][0] = new TGraphErrors();
            fGraphEnergyFit[dss.det][dss.side][dss.strip][1] = new TGraphErrors();
            for (auto gate=0; gate<fNumGates; ++gate)
            {
                fGraphLeftRightX[dss.det][dss.side][dss.strip][gate] = new TGraph();
                auto graphLR = fGraphLeftRightX[dss.det][dss.side][dss.strip][gate];
                fGraphEnergyPositionX[dss.det][dss.side][dss.strip][gate] = new TGraph();
                auto graphEP = fGraphEnergyPositionX[dss.det][dss.side][dss.strip][gate];
                graphLR -> SetMarkerStyle(20);
                graphEP -> SetMarkerStyle(20);
                for (auto oStrip=0; oStrip<fNumOStrips-1; ++oStrip)
                {
                    for (auto lr=0; lr<2; ++lr)
                    {
                        auto hist0 = fHistLRProjGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip  ][lr];
                        auto hist1 = fHistLRProjGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip+1][lr];
                        FindHistX(hist0,hist1,histDummy);
                        // tried gaus fit but it doesn't work because shape is not gaussian and number of bin is too small
                        fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][lr][0] = histDummy->GetMaximum();
                        fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][lr][1] = histDummy->GetMean();
                        fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][lr][2] = histDummy->GetStdDev();
                        if (oStrip==1)
                        {
                            double eRatio = ((lr==1) ? double(oStrip+1)/fNumOStrips : 1-double(oStrip+1)/fNumOStrips );
                            auto graphFit = fGraphEnergyFit[dss.det][dss.side][dss.strip][lr];
                            graphFit -> SetPoint(graphFit->GetN(),fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][lr][1],fGateEnergy[gate]*eRatio);
                        }
                    }
                    graphLR -> SetPoint(graphLR->GetN(),fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][0][1],fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][1][1]);
                }
            }
            double pp[2][2];
            for (auto lr=0; lr<2; ++lr)
            {
                fFitEnergy[dss.det][dss.side][dss.strip][lr] = new TF1("f1","pol1",x1,x2);
                auto f1 = fFitEnergy[dss.det][dss.side][dss.strip][lr];
                fGraphEnergyFit[dss.det][dss.side][dss.strip][lr] -> Fit(f1,"RQN0");
                pp[lr][0] = f1 -> GetParameter(0);
                pp[lr][1] = f1 -> GetParameter(1);
                
            }
            if (!checkCECorrection) FillCEParameters(dss.det, dss.side, dss.strip, 2, pp[0][0], pp[0][1], pp[1][0], pp[1][1]);
        }
    }
    if (!checkCECorrection) EndWriteParameters();

    /////////////////////////////////////////////////////////////////////
    // x) Draw examples
    /////////////////////////////////////////////////////////////////////
    int countLR = 0;
    int countFrame;
    auto groupFit = top -> CreateGroup("Fit");
    auto groupEProj = top -> CreateGroup("EnergyProj");
    auto groupPos = top -> CreateGroup("Position");
    auto groupLR = top -> CreateGroup(Form("LRColored%d",countLR));
    auto groupEP = top -> CreateGroup(Form("EPColored%d",countLR));
    for (auto dssGroup : fAllGroupArrayR)
    {
        auto subLR = groupLR -> CreateGroup(Form("LR%d",dssGroup.GetDet()));
        auto subEP = groupEP -> CreateGroup(Form("EP%d",dssGroup.GetDet()));
        if (groupLR -> GetNumAllDrawings()>50) {
            countLR++;
            groupLR = top -> CreateGroup(Form("LRColored%d",countLR));
            groupEP = top -> CreateGroup(Form("EPColored%d",countLR));
        }
        LKDrawingGroup* subFit[2];
        LKDrawingGroup* subXX[2];
        LKDrawingGroup* subEProj[2];
        subFit[0] = groupFit -> CreateGroup(Form("Fit%dL",dssGroup.GetDet()));
        subFit[1] = groupFit -> CreateGroup(Form("Fit%dR",dssGroup.GetDet()));
        subEProj[0] = groupEProj -> CreateGroup(Form("Epj%dL",dssGroup.GetDet()));
        subEProj[1] = groupEProj -> CreateGroup(Form("Epj%dR",dssGroup.GetDet()));
        auto subPos = groupPos -> CreateGroup(Form("Pos%d",dssGroup.GetDet()));
        for (auto dss : dssGroup.array)
        {
            auto drawingEP = subEP -> CreateDrawing();
            auto drawingLR = subLR -> CreateDrawing();
            for (auto gate=0; gate<fNumGates; ++gate)
            {
                auto drawingPos = FitStepHistogram(fHistPositionGate[dss.det][dss.side][dss.strip][gate]);
                subPos -> AddDrawing(drawingPos);
                /*
                {
                    auto histPos = fHistPositionGate[dss.det][dss.side][dss.strip][gate];
                    auto drawingPos = subPos -> CreateDrawing();
                    //auto fslope = new TF1("fslope","( x>[0]&&x<[1] ? ([2]*(x-[0])+[3]) : 0 )",fBinX1,fBinX2);
                    auto hmax = histPos->GetMaximum();
                    auto graphPos = new TGraphErrors();
                    //graphPos -> SetMarkerStyle(24);
                    //graphPos -> SetMarkerSize(0.75);
                    graphPos -> Set(0);
                    double ex = 0.5*(fBinX2-fBinX1)/nbinsx;
                    bool found_bx1 = false;
                    bool found_bx2 = false;
                    double bx1 = 0;
                    double bx2 = 0;
                    for (auto bin=1; bin<=nbinsx; ++bin) {
                        double error0 = 0.005*hmax;
                        double value = histPos -> GetBinContent(bin);
                        double x0 = histPos -> GetBinCenter(bin);
                        graphPos -> SetPoint(bin-1,x0,value);
                        if (value>0)
                            graphPos -> SetPointError(bin-1,ex,hmax*0.1);
                        else
                            graphPos -> SetPointError(bin-1,ex,0);
                        if (value>0.1*hmax && found_bx1==false) {
                            found_bx1 = true;
                            bx1 = x0;
                        }
                        if (value==0 && x0>0 && found_bx1==true && found_bx2==false) {
                            found_bx2 = true;
                            bx2 = x0 - ex*2;
                        }
                    }
                    double slope = -0.8*hmax/(bx2-bx1);
                    double dx = 0.010;
                    TF1* f1Step;
                    if (1) {
                        f1Step = new TF1("f1Step",SiAnaFitStep,fBinX1,fBinX2,5);
                        f1Step -> SetParLimits(0,bx1-0.15,bx1+0.15);
                        f1Step -> SetParLimits(1,bx2-0.15,bx2+0.15);
                        f1Step -> SetParLimits(2,-1.*hmax/(bx2-bx1),-0.5*hmax/(bx2-bx1));
                        f1Step -> SetParLimits(3,hmax*1.2,hmax*0.5);
                        f1Step -> SetParLimits(4,0.001,0.020);
                        f1Step -> SetParameters(bx1,bx2,slope,hmax,dx);
                        f1Step -> SetNpx(500);
                    }
                    else {
                        f1Step = new TF1("f1Step","(x>[0]&&x<[1])*([2]*(x-[0])+[3])",fBinX1,fBinX2);
                        f1Step -> SetParameters(-0.6,0.6,-0.8*hmax/1.2,hmax);
                        f1Step -> SetParLimits(0,-0.9,-0.4);
                        f1Step -> SetParLimits(1,0.4,0.9);
                        f1Step -> SetParLimits(2,-hmax/1.,-0.25*hmax/1.);
                        f1Step -> SetParLimits(3,hmax*1.2,hmax*0.5);
                        f1Step -> SetParameters(-0.6,0.6,-0.8*hmax/1.2,hmax);
                    }
                    graphPos -> Fit(f1Step,"RQN0");
                    drawingPos -> Add(histPos);
                    drawingPos -> Add(graphPos,"samepz");
                    drawingPos -> Add(f1Step,"samel");
                    {
                        auto lg = new TLegend(0.45,0.60,0.95,0.88);
                        lg -> SetBorderSize(0);
                        lg -> SetFillStyle(0);
                        lg -> AddEntry((TObject*)nullptr,Form("x1 = %.2f (%.2f)",f1Step->GetParameter(0),bx1  ),"");
                        lg -> AddEntry((TObject*)nullptr,Form("x2 = %.2f (%.2f)",f1Step->GetParameter(1),bx2  ),"");
                        lg -> AddEntry((TObject*)nullptr,Form("sl = %.2f (%.2f)",f1Step->GetParameter(2),slope),"");
                        lg -> AddEntry((TObject*)nullptr,Form("mx = %.2f (%.2f)",f1Step->GetParameter(3),hmax ),"");
                        lg -> AddEntry((TObject*)nullptr,Form("dx = %.4f (%.4f)",f1Step->GetParameter(4),dx   ),"");
                        drawingPos -> Add(lg,"same");
                    }
                }
                */
                LKDrawing *drawingEProj[2];
                drawingEProj[0] = subEProj[0] -> CreateDrawing();
                drawingEProj[1] = subEProj[1] -> CreateDrawing();
                for (auto oStrip=0; oStrip<fNumOStrips; ++oStrip)
                {
                    {
                        auto histEP = fHistEnergyPositionGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip];
                        SetHistColor(histEP,oStrip,fNumOStrips);
                        //histEP -> SetMaximum(fNumOStrips+1);
                        drawingEP -> Add(histEP);
                    }
                    {
                        auto histLR = fHistLeftRightGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip];
                        SetHistColor(histLR,oStrip,fNumOStrips);
                        //histLR -> SetMaximum(fNumOStrips+1);
                        drawingLR -> Add(histLR);
                    }
                    for (auto lr=0; lr<2; ++lr)
                    {
                        auto histLRP = fHistLRProjGateOhmic[dss.det][dss.side][dss.strip][gate][oStrip][lr];
                        histLRP -> SetLineColor(oStrip+1);
                        drawingEProj[lr] -> Add(histLRP);
                        if (oStrip<fNumOStrips-1)
                            drawingEProj[lr] -> Add(new TMarker(fXParameters[dss.det][dss.side][dss.strip][gate][oStrip][lr][1],10,20));
                    }
                }
                drawingLR -> Add(fGraphLeftRightX[dss.det][dss.side][dss.strip][gate],"samep");
                auto yxline = new TLine(fGateEnergy[gate],0,0,fGateEnergy[gate]);
                drawingLR -> Add(yxline,"samel");
                for (auto lr=0; lr<2; ++lr)
                {
                    auto line = new TLine(0.5*fGateEnergy[gate],0,0.5*fGateEnergy[gate],fHistLRProjGateOhmic[dss.det][dss.side][dss.strip][gate][0][0]->GetMaximum());
                    line -> SetLineColor(kMagenta);
                    //line -> SetLineStyle(2);
                    drawingEProj[lr] -> Add(line);
                }
            }
            for (auto lr=0; lr<2; ++lr)
            {
                auto drawing = subFit[lr] -> CreateDrawing();
                auto graph = fGraphEnergyFit[dss.det][dss.side][dss.strip][lr];
                graph -> SetMarkerStyle(20);
                auto f1 = fFitEnergy[dss.det][dss.side][dss.strip][lr];
                TString title = Form("[%d-%d-%d] %.3f*x + %.3f;Slope calibrated energy;Alpha source energy",dss.det,dss.side,dss.strip,f1->GetParameter(1),f1->GetParameter(0));
                drawing -> Add(new TH2I(Form("hist%d",countFrame++),title,10,0,e2,10,0,e2));
                drawing -> Add(graph,"pl");
                drawing -> Add(f1,"samel");
                auto lg = new TLegend(0.1,0.55,0.55,0.88);
                lg -> SetBorderSize(0);
                lg -> SetFillStyle(0);
                lg -> AddEntry((TObject*)nullptr,Form("p1 = %.3f",fXParameters[dss.det][dss.side][dss.strip][0][1][lr][1]),"");
                lg -> AddEntry((TObject*)nullptr,Form("p2 = %.3f",fXParameters[dss.det][dss.side][dss.strip][1][1][lr][1]),"");
                drawing -> Add(lg,"same");
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 7) Write histograms
    /////////////////////////////////////////////////////////////////////
    if (!checkCECorrection) 
    {
        TString foutName = fCEHistFileName;
        if (fTest)
            foutName.ReplaceAll(".root",".test.root");
        cout << foutName << endl;
        auto fileHist = new TFile(foutName,"recreate");
        top -> Write();
    }

    if (runViewer)
    {
        top -> Draw("viewer");
    }
}
