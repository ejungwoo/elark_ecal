#include "LKLogger.h"

//int fRun = 199;
//int fRun = 303;
//int fRun = 99;
//int fRun = 253;
int fRun = 113;

int fMaxPeaks = 2;

int fDataType = -1;  // 0: 12dE+16E, 1: 12E
int fChooseGate = 1; // 0: Gd, 1: Am
int fEntriesCut = 100; // will be updated as a function of total entries
//int fDrawingExampleDetectors[] = {32,33,34,35,36,37,38,39};
int fDrawingExampleDetectors[] = {4,14,22,34,39};


      int    fNBinA = 200;
const double fBinA1 = 0;
const double fBinA2 = 4000;

      int    fNBinE = 200;
const double fBinE1 = 0;
const double fBinE2 = 10;

const int    fNBinEFix = 400;

      int    fNBinX = 200;
//      int    fNBinX = 120;
const double fBinX1 = -1;
const double fBinX2 = 1;

const int    fNumDetectors = 40;
const int    fNumStrips = 8;
const int    fNumJStrips = 8;
const int    fNumOStrips = 4;
const int    fNumSides = 2;
const int    fNumGates = 2;

const double fExpectedResolution = 0.015;

//////////////////////////////////////////////////////////////////////////////////
TH2D* fHistEnergyDetector[2];
TH2D* fHistEnergyPositionAll;
TH1D* fHistEnergy             [fNumDetectors][2][fNumStrips];
TH1D* fHistEnergySum          [fNumDetectors][2][fNumStrips];
TH2D* fHistLeftRight          [fNumDetectors][2][fNumStrips];
TH2D* fHistEnergyPosition     [fNumDetectors][2][fNumStrips];
TH2D* fHistLeftRightGate      [fNumDetectors][2][fNumStrips][fNumGates];
TH2D* fHistEnergyPositionGate [fNumDetectors][2][fNumStrips][fNumGates];

//////////////////////////////////////////////////////////////////////////////////
TString fRecoFileName;
TString fHistFileName;
TString fC1HistFileName;
TString fC2HistFileName;
TString fDummyFileName;
TString fAllParFileName;
TString fC0ParFileName;
TString fC1ParFileName;
TString fC2ParFileName;

//////////////////////////////////////////////////////////////////////////////////
double fC0Parameters[40][2][8];
double fC1Parameters[40][2][8][8];
double fC2Parameters[40][2][8][3];

//////////////////////////////////////////////////////////////////////////////////
TCanvas *MakeCanvas(TString name, int size=8);

//////////////////////////////////////////////////////////////////////////////////
class strip_info {
    public:
        strip_info(int type_, int det_, int side_, int strip_) { type = type_; det = det_; side = side_; strip = strip_; }
        int type;
        int det;
        int side;
        int strip;
};
class strip_group {
    public:
        strip_group() {}
        vector<strip_info> array;
        int GetDet() {
            if (array.size()==0)
                return -1;
            return array[0].det;
        }
        TCanvas *MakeGroupCanvas(TString name, int npads=-1) {
            TString cvsName = Form("cvs_%s_%d",name.Data(),array.at(0).det);
            if (npads<0) npads = array.size();
            return MakeCanvas(cvsName,npads);
        }
};
vector<strip_info> fStripArrayR;
vector<strip_info> fStripArrayS; /// position sensitive strips
vector<strip_group> fExampleGroupArrayR;
vector<strip_group> fExampleGroupArrayS;
vector<strip_group> fAllGroupArrayR;
vector<strip_group> fAllGroupArrayS;

//////////////////////////////////////////////////////////////////////////////////
SKSiArrayPlane* fStark = nullptr;

bool ContinueRegardingToDataType(int det)
{
    if (fDataType==2)
        return false;
    if ((fDataType==0 && det<12) || (fDataType==1 && det>=12))
        return true;
    return false;
}

bool MakeRun(int run=-1, int type=-1)
{
    if (run<0) {
        if (gSystem -> Getenv("RUN"))
            run = atoi(gSystem -> Getenv("RUN"));
        else
            run = fRun;
    }
    fRun = run;

    fDataType = -1;
    if (type<0) {
        if (gSystem -> Getenv("TYPE"))
            type = atoi(gSystem -> Getenv("TYPE"));
        else
            type = fDataType;
    }
    fDataType = type;

    if (fDataType<0) {
        fDataType = 2;
        //if (fRun==199) fDataType = 1;  // 0: 12dE+16E, 1: 12E
        //if (fRun==303) fDataType = 0;  // 0: 12dE+16E, 1: 12E
        //if (fRun==253) fDataType = 0;  // 0: 12dE+16E, 1: 12E
        //if (fRun==113) fDataType = 2;  // 0: 12dE+16E, 1: 12E
    }

    if (fDataType<0)
        fDataType = 2;

    fStripArrayR.clear();
    fStripArrayS.clear();
    fExampleGroupArrayR.clear();
    fExampleGroupArrayS.clear();
    //fAllGroupArrayR.clear();
    //fAllGroupArrayS.clear();

    if (fRun==303) {
        fNBinA = 500;
        fNBinE = 500;
        fNBinX = 500;
    }

    if (fRun==133) {
        fNBinA = 500;
        fNBinE = 500;
        fNBinX = 500;
    }

    const char* recoDir = "../data_reco";
    const char* dataDir = "data";
    fRecoFileName   = Form("%s/stark_%04d.reco.root",recoDir,run);
    fHistFileName   = Form("%s/stark_%04d.hist.root",dataDir,run);
    fC1HistFileName = Form("%s/stark_%04d.hist_c1.root",dataDir,run);
    fAllParFileName = Form("%s/stark_%04d.par.root",dataDir,run);
    fC0ParFileName  = Form("%s/stark_%04d.c0.root",dataDir,run);
    fC1ParFileName  = Form("%s/stark_%04d.c1.root",dataDir,run);
    fC2HistFileName = Form("%s/stark_%04d.hist_c2.root",dataDir,run);
    fC2ParFileName  = Form("%s/stark_%04d.c2.root",dataDir,run);
    fDummyFileName  = Form("%s/dummy",dataDir);

    auto file = new TFile(fRecoFileName,"read");
    auto tree = (TTree*) file -> Get("event");
    if (tree==nullptr) {
        e_error << fRecoFileName << " " << tree << endl;
        return false;
    }
    fEntriesCut = tree -> GetEntries() / 10000;
    if (fEntriesCut<100) fEntriesCut = 100;
    e_info << "run " << fRun << " " << tree -> GetEntries() << " (" << fEntriesCut << ")" << ", Data-type is " << fDataType << endl;
    file -> Close();

    if (fStark==nullptr)
    {
        fStark = new SKSiArrayPlane();
        fStark -> AddPar("config_stark.mac");
        fStark -> Init();
    }

    auto numDetectors = fStark -> GetNumSiDetectors();
    for (auto det=0; det<numDetectors; ++det) {
        auto detector = fStark -> GetSiDetector(det);
        auto detType = detector -> GetDetType();
        auto numJStrips = detector -> GetNumJunctionStrips();
        auto numOStrips = detector -> GetNumOhmicStrips();
        auto numJDirection = detector -> GetNumJunctionDirection();

        if (ContinueRegardingToDataType(det))
            continue;

        {
            strip_group dss_group1;
            for (auto side=0; side<2; ++side)
            {
                if (side==0&&numJDirection==2)
                    continue;
                auto numStrips = (side==0?numJStrips:numOStrips);
                for (auto strip=0; strip<numStrips; ++strip) {
                    fStripArrayS.push_back(strip_info(detType,det,side,strip));
                    dss_group1.array.push_back(strip_info(detType,det,side,strip));
                }
            }
            bool isExample = false; for (auto ex : fDrawingExampleDetectors) if (det==ex) isExample = true;
            if (isExample) fExampleGroupArrayS.push_back(dss_group1);
            fAllGroupArrayS.push_back(dss_group1);
        }

        if (numJDirection==2) {
            strip_group dss_group2;
            int side = 0;
            for (auto strip=0; strip<numJStrips; ++strip) {
                fStripArrayR.push_back(strip_info(detType,det,side,strip));
                dss_group2.array.push_back(strip_info(detType,det,side,strip));
            }
            bool isExample = false; for (auto ex : fDrawingExampleDetectors) if (det==ex) isExample = true;
            if (isExample) fExampleGroupArrayR.push_back(dss_group2);
            fAllGroupArrayR.push_back(dss_group2);
        }
    }

    return true;
}

TString MakeHistName(TString xname, TString yname, int det, int side, int strip, int gate=-1)
{
    const char* mainName  = yname.IsNull() ? xname.Data() : Form("%s_%s",xname.Data(),yname.Data());
    const char* detName   = ( det  <0 ? "" : Form("_d%d",det   ) );
    const char* sideName  = ( side <0 ? "" : (side==0?"_j":"_o") );
    const char* stripName = ( strip<0 ? "" : Form("_s%d",strip ) );
    const char* gateName  = ( gate <0 ? "" : Form("_g%d",gate  ) );
    TString name = Form("hist%s%s%s%s_%s",detName,sideName,stripName,gateName,mainName);
    //TString name = Form("%s%s%s%s%s",name,detName,sideName,stripName,gateName);
    return name;
}

TString MakeHistTitle(TString xname, TString yname, int det, int side, int strip, int gate=-1)
{
    auto detector = fStark -> GetSiDetector(det);
    const char* detTitle   = ( det  <0 ? "" : Form("detector %d (%s)",det, detector->GetDetTypeName().Data()));
    const char* sideTitle  = ( side <0 ? "" : (side==0?", junction":", ohmic"));
    const char* stripTitle = ( strip<0 ? "" : Form(", strip %d",strip ));
    const char* gateTitle  = ( gate <0 ? "" : Form(", gate %d",gate  ));
    if (TString(detTitle).IsNull()) sideTitle = TString(TString(sideTitle)(2,TString(sideTitle).Sizeof()-3)).Data();
    TString title = Form("%s%s%s%s;%s;%s",detTitle,sideTitle,stripTitle,gateTitle,xname.Data(),yname.Data());
    return title;
}

TH1D* MakeHist1(TString xname, TString yname, int det, int side, int strip, int gate, int nx=0, double x1=0, double x2=0)
{
    if (nx==0) {
        if (xname.Index("c0")==0 || xname.Index("c1")==0 || xname.Index("c2")==0) {
            nx = 200;
            x1 = 0;
            x2 = 10;
        }
        else {
            if      (det<32  && side==1) { nx = 200; x1 = 200; x2 = 1200; }
            else if (det>=32 && side==1) { nx = 200; x1 = 500; x2 = 2200; }
            else if (det==39 && side==0 && (strip==2 || strip==3) ) { nx = 200; x1 = 200; x2 = 2000; }
            else if (det>=32 && side==0) { nx = 200; x1 = 0; x2 = 2800; }
            else                         { nx = 200; x1 = 0; x2 = 4000; }
        }
    }

    auto hist = new TH1D(MakeHistName(xname,yname,det,side,strip,gate),MakeHistTitle(xname,yname,det,side,strip,gate),nx,x1,x2);
    //hist -> SetFillColor(29);
    return hist;
}

TH2D* MakeHist2(TString xname, TString yname, int det, int side, int strip, int gate, int nx, double x1, double x2, int ny, int y1, int y2)
{
    auto hist = new TH2D(MakeHistName(xname,yname,det,side,strip,gate),MakeHistTitle(xname,yname,det,side,strip,gate),nx,x1,x2,ny,y1,y2);
    return hist;
}

TCanvas *MakeCanvas(TString cvsName, int npads)
{
    TCanvas *cvs;
    if (npads==2) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,600,500,0.9);
        cvs -> Divide(1,2,0.002,0.002);
    }
    else if (npads==4) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,550,500,0.7);
        cvs -> Divide(2,2,0.002,0.002);
    }
    else if (npads==8||npads==9) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,550,500,0.95);
        cvs -> Divide(3,3,0.002,0.002);
    }
    else if (npads==20) {
        cvs = LKPainter::GetPainter() -> CanvasResize(cvsName,700,500,1);
        cvs -> Divide(5,4,0.002,0.002);
    }
    else {
        cvs = LKPainter::GetPainter() -> CanvasDefault(cvsName);//,500,500,0.95);
    }
    for (auto i=1; i<=npads; ++i) {
        cvs -> cd(i) -> SetMargin(0.15, 0.12, 0.12, 0.10);
        cvs -> cd(i) -> SetGridx();
        cvs -> cd(i) -> SetGridy();
    }
    cvs -> cd();
    return cvs;
}

bool IsPositionSensitiveStrip(int det, int side=1)
{
    auto detector = fStark -> GetSiDetector(det);
    //detector->Print();
    //lk_debug << "det=" << det << " numjd=" << detector->GetNumJunctionDirection()  <<  " side=" << side << " " << (detector->GetNumJunctionDirection()==2)  <<  " " << (side==0) << endl;
    //if (det==39) cout << det << " " << detector->GetNumJunctionDirection() << endl;
    if (detector->GetNumJunctionDirection()==2 && side==0)
        return true;
    return false;
}

bool IsPositionSensitiveStrip(LKSiChannel* channel)
{
    return IsPositionSensitiveStrip(channel->GetDetID(),channel->GetSide());
}

void GetC0Parameters()
{
    int det, side, strip;
    double entries, a1, m1, s1, a2, m2, s2, g0;
    auto fpar = new TFile(fC0ParFileName,"read");
    auto tree = (TTree*) fpar -> Get("parameters");
    tree -> SetBranchAddress("det"    ,&det    );
    tree -> SetBranchAddress("side"   ,&side   );
    tree -> SetBranchAddress("strip"  ,&strip  );
    tree -> SetBranchAddress("entries",&entries);
    tree -> SetBranchAddress("a1"     ,&a1     );
    tree -> SetBranchAddress("m1"     ,&m1     );
    tree -> SetBranchAddress("s1"     ,&s1     );
    tree -> SetBranchAddress("a2"     ,&a2     );
    tree -> SetBranchAddress("m2"     ,&m2     );
    tree -> SetBranchAddress("s2"     ,&s2     );
    tree -> SetBranchAddress("g0"     ,&g0     );
    auto n = tree -> GetEntries();
    for (auto i=0; i<n; ++i)
    {
        tree -> GetEntry(i);
        fC0Parameters[det][side][strip] = g0;
    }
}

void GetC1Parameters()
{
    int det, side, strip;
    double entries ,x1, y1, x2, y2 ,g1, g2 ,b, slope;
    auto fpar = new TFile(fC1ParFileName,"read");
    auto tree = (TTree*) fpar -> Get("parameters");
    tree -> SetBranchAddress("det"    ,&det    );
    tree -> SetBranchAddress("side"   ,&side   );
    tree -> SetBranchAddress("strip"  ,&strip  );
    tree -> SetBranchAddress("entries",&entries);
    tree -> SetBranchAddress("x1"     ,&x1     );
    tree -> SetBranchAddress("y1"     ,&y1     );
    tree -> SetBranchAddress("x2"     ,&x2     );
    tree -> SetBranchAddress("y2"     ,&y2     );
    tree -> SetBranchAddress("b"      ,&b      );
    tree -> SetBranchAddress("slope"  ,&slope  );
    tree -> SetBranchAddress("g1"     ,&g1     );
    tree -> SetBranchAddress("g2"     ,&g2     );
    auto n = tree -> GetEntries();
    for (auto i=0; i<n; ++i)
    {
        tree -> GetEntry(i);
        fC1Parameters[det][side][strip][0] = g1;
        fC1Parameters[det][side][strip][1] = g2;
        fC1Parameters[det][side][strip][2] = b;
        fC1Parameters[det][side][strip][3] = slope;
        fC1Parameters[det][side][strip][4] = x1;
        fC1Parameters[det][side][strip][5] = y1;
        fC1Parameters[det][side][strip][6] = x2;
        fC1Parameters[det][side][strip][7] = y2;
    }
}

void GetC2Parameters()
{
    int det, side, strip;
    double entries, b0, b1, b2;
    auto fpar = new TFile(fC2ParFileName,"read");
    auto tree = (TTree*) fpar -> Get("parameters");
    tree -> SetBranchAddress("det"    ,&det    );
    tree -> SetBranchAddress("side"   ,&side   );
    tree -> SetBranchAddress("strip"  ,&strip  );
    tree -> SetBranchAddress("entries",&entries);
    tree -> SetBranchAddress("b0"     ,&b0     );
    tree -> SetBranchAddress("b1"     ,&b1     );
    tree -> SetBranchAddress("b2"     ,&b2     );
    auto n = tree -> GetEntries();
    for (auto i=0; i<n; ++i)
    {
        tree -> GetEntry(i);
        fC2Parameters[det][side][strip][0] = b0;
        fC2Parameters[det][side][strip][1] = b1;
        fC2Parameters[det][side][strip][2] = b2;
    }
}

void GetFinalHistograms()
{
    TString calName = "c2_";
    auto fileHist = new TFile(fC2HistFileName,"read");
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = (TH1D*) fileHist ->Get(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) fileHist ->Get(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
        fHistLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) fileHist ->Get(MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) fileHist ->Get(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
    }
    fHistEnergyDetector[0] = (TH2D*) fileHist ->Get(MakeHistName(calName+"esum","det",-1,0,-1,-1));
    fHistEnergyDetector[1] = (TH2D*) fileHist ->Get(MakeHistName(calName+"esum","det",-1,1,-1,-1));

    cout << "fHistEnergy" << endl;
    cout << "fHistEnergySum" << endl;
    cout << "fHistLeftRight" << endl;
    cout << "fHistEnergyPosition" << endl;
    cout << "fHistEnergyDetector" << endl;
    cout << "fHistEnergyDetector" << endl;
}
