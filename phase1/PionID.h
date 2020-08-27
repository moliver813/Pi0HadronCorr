
  #ifndef PIONID_H
  #define PIONID_H

  // --- ROOT system ---
  #include <TFile.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TH3F.h>
  #include <TF1.h>
  #include <TMultiGraph.h>
  #include <TGraphErrors.h>
  #include <TPaveStats.h>
  #include <TSpectrum.h>
  #include <TStopwatch.h>
  #include <TROOT.h>
  #include <TLine.h>
  #include <TCanvas.h>
  #include <TStyle.h>
  #include <TLegend.h>
  #include <TString.h>
  #include <TSystem.h>
  #include <TList.h>
  #include <TLatex.h>
  #include <TColor.h>
  #include <TGaxis.h>
  #include <TStyle.h>
  #include <TEnv.h>
  #include <TMath.h>
  #include <TMinuit.h>
  #include <THnSparse.h>
  #include <TDatime.h>


  #include <Riostream.h>
  #include <TString.h>

  #include <cstdlib>

  #include "fitAlgos.cc"


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - Here are functions that make the plotting of figures nicer- - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void set_plot_style()
  {
      const Int_t NRGBs = 5;
      const Int_t NCont = 99;//max possible?

      //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
      Double_t stops[NRGBs] = { 0.00, 0.25, 0.5, 0.75, 1.00 };
      Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
      Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
      Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
      gStyle->SetNumberContours(NCont);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  class TH1;
  class TH2;
  class TH1F;
  class TH1D;
  class TH2D;
  class TF1;
  class TFile;
  class TList;
  class TBox;
  class TCanvas;
  class TGraphErrors;
  class TLatex;

  class PionID : public TObject {
  public:
    PionID();
    virtual ~PionID() { ; }
    void SetStyle();
    void SetDebugLevel(Int_t input) { fDebugLevel = input; }
    void SetLabel(TString input)                      { sLabel = input; }
    void SetLabel2(TString input)                      { sLabel2 = input; }

    void EnablePerformance(Bool_t input)              { bEnablePerformance = input; }

    void SetPi0CandInputFile(TFile * input)           { fInputFile = input; }
    void SetOutputDirectory(TString input)            { sOutputDir = input; }
    void SetOutputFileName(TString input)             { sOutputFileName = input; }
    void SetListName(TString input)                   { sListName = input; }
    void PrintCanvas(TCanvas * canvas, TString name);

    void DrawWIP(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size);
    void DrawAlicePerf(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size);

    // The Ranges
    void SetEventPlaneBinRange(Int_t input1, Int_t input2) { EventPlaneBinLow = input1; EventPlaneBinHigh = input2; }
    void SetLambdaBinRange(Int_t input1, Int_t input2) { LambdaBinLow = input1; LambdaBinHigh = input2; }
    void SetEnergyBinRange(Int_t input1, Int_t input2) { EnergyBinLow = input1; EnergyBinHigh = input2; }
    void SetAsymBinRange(Int_t input1, Int_t input2)   { AsymBinLow = input1; AsymBinHigh = input2; }
    void SetOpeningAngleBinRange(Int_t input1, Int_t input2) { OpeningAngleBinLow = input1; OpeningAngleBinHigh = input2; }

    // Controls that may actually be used
    void SetPtBinChoice(Int_t input)   { PtBinChoice = input; }
    void SetBkgChoice(Int_t input)     { bkgChoice = input; }
    void SetFitMethod(Int_t input)     { fitMethod = input; }
    void SetFitPeakMethod(Int_t input) { fitPeakMethod = input; }
    void SetFitBkgMethod(Int_t input)  { fitBkgMethod = input; }

    void SetNSigma(Double_t input)      { nSigma = input; }
    void SetNSigmaR(Double_t input)     { nSigmaR = input; }
    void SetFixedMassWindows(Int_t input) { iFixedMassWindows = input; }
    void SetScaleBkg(Bool_t input)      { scaleBkg = input; }
    void SetBkgScaleMin(Double_t input) { bkgScaleMin = input; }
    void SetBkgScaleMax(Double_t input) { bkgScaleMax = input; }
    void SetRemoveMCPi0PS(Bool_t input) { RemoveMCPi0PS = input; }
    void SetRemoveMCEta(Bool_t input)   { RemoveMCEta = input; }

    void SetFitMinX(Double_t input)     { fFitMinX = input; }
    void SetFitMaxX(Double_t input)     { fFitMaxX = input; }
    void SetNRebinMass(Int_t input)     { nRebinMass = input; }
    void SetPtBins(Int_t input)         { nPtBins = input; }
    void SetHighPtCut(Float_t input)    { kHighPtCut = input;}

    void LoadMCPreAnalysis();
    void SetMCPreAnalysisFile(TString input)  {
      sMCPreAnalysisFile = input;
      if (input.Length() > 0) bUseMCPreAnalysis = true;
    }

    void LoadThetaModelParameters();
    void SetEnableThetaModel(Bool_t input)    { bEnableThetaModel = input; }
    void SetUseThetaLookUpTable(Bool_t input) { bUseThetaLookUpTable = input; }
    void SetThetaModelRootFile(TString input) { sThetaModelRootFile = input; }
    void SetThetaModelChoice(Int_t input)     { iThetaModelParamChoice = input; }
    void SetThetaModelTrigger(Int_t input)    { iThetaModelTrigger = input; }
    void SetThetaModelCent(Int_t input)       { 
      iThetaModelCent = input; 
    //  if(input==1) kMagicScale=1.5; // Doing so much more with this function
    }
    
    void SetEnablePSScaleMethod(Bool_t input) { bEnablePSScaleMethod = input; }
    void SetEnablePSDirectMethod(Bool_t input) { bEnablePSDirectMethod = input; }
    // Setting initial values for signal model
    void SetWPower(Double_t input)      { fWPower = input; }
    void SetWYield(Double_t input)      { fWYield = input; }
    void SetWPt0(Double_t input)      { fWPt0 = input; }
    void SetWMass(Double_t input)       { fWMass = input; }
    void SetWSigma(Double_t input)       { fWSigma = input; }

    void PrintSettings();
    void Run();

  protected:
    Bool_t LoadHistograms();
    void OpeningAngleAnalysis();   // obsolete part of analysis
    void SetPtBins();              // Sets up the Pt Bins. Can generate optimized bins
    void MakeBasicPlots();
    void DrawMassPlots();
    void DrawResultGraphs();
    void DrawDemoPlots();
    void ProduceDeltaPsiPlots();
    void MeasureVn();
    void DoProjections();           // Makes the projections from the THnSparse
    void Pi0MassAnalysis();

    bool bUseMCPreAnalysis = true;
    TString sMCPreAnalysisFile = "/home/moliver/cern/gammaHadron/wrk/phase1/output/T38/FitScan/CentN/Fit_6/SecondAnalysis.root"; // if blank, no file used for MC fit parameters
//    TString sMCPreAnalysisFile = "/home/alidock/cern/gammaHadron/wrk/phase1/output/T38/CentN/EP-1/SecondAnalysis.root"; // if blank, no file used for MC fit parameters
    vector<vector<TF1 *>> fMCPreAnalysisPi0Fits = {}; // [Cent Bin] [Pt Bin]
    int iMCPreAnalysis_Cent    = 0; // Which centrality's MC Pre-Analysis to use?

    bool bEnableThetaModel = true;
    bool bUseThetaLookUpTable = true;
    // Will there need to be differences for different centralities?
    TString sThetaModelRootFile = "/home/moliver/cern/gammaHadron/wrk/phase1/AngleEffect/T40/GA/CentN_HadCorr/AngleAnalysis.root"; //FIXME
     //TString sThetaModelRootFile = "/home/moliver/cern/gammaHadron/wrk/phase1/AngleEffect/T40/GA/CentN_HadCorr/AngleAnalysis.root";
    int  iThetaModelParamChoice = 0; // Which parameter set to use

    int iNPtBinsThetaModel = 6; // How many pt bins to load
    int iThetaModelTrigger = 0; // 0 for MB, 1 for GA, 2 for MC
    int iThetaModelCent    = 0; // Which centrality's angle analysis to use? // Also used by fixed pi0 mass windows // And by the Vn analysis

    vector<vector<TGraphErrors *>> fThetaMassPrimeGraphs = {};
    vector<vector<TGraphErrors *>> fThetaLambdaGraphs = {};
    vector<TGraphErrors *> fMassPrimePar1Graphs = {};
    vector<TGraphErrors *> fLambdaPar1Graphs = {};


    // void LoadThetaModelParameters(); // Implement if parameters stored in histogram
    void GetThetaModelParameters(double fPt, double fThetaC, double &ThetaModelLambda, double &ThetaModelMPrime);

    void FitMCTruthPi0();

    // For MC, make plots showing potential PSCorrection with Same-Pos Pi0 pairs
    void DrawPosSwapMCSub();

    // Advanced Background Modelling 
    void InitializeWSignal();


    // Extending the method outlined by E+W
    bool bEnablePSScaleMethod = true;
    void DrawPSScaleCorrectionPlots(); // Experimental scaling correction to Position Swap method
    void BuildPSScaleDistributions();  // Use U,V and W to build E and P
    double IntegrateScaleModification(TH2 * fDist, Double_t lMass, Double_t lPt); 

    // A more direct method involving direct mass pt maps
    bool bEnablePSDirectMethod = true;
    void DrawPSDirectCorrectionPlots(); 
    void BuildPSDirectDistributions();


    std::vector <double> FindEqualBins(TH1D * intHist, int nBins = 10);

    void CopyTF1Details(TF1 * target, TF1 * source);
    void SetTH1Histo(TH1 *hist,TString xTitle,TString yTitle,Bool_t big);
    void AddFunctionToHist(TH1 *Hist, TF1 * func);
    TH1D * IntegralHist(TH1D * hist);

    void SaveResults();


    TString sLabel  = "";
    TString sLabel2 = "";


    Int_t fDebugLevel;                      ///<For Debugging Purposes
    Int_t fCanvasWidth = 600;
    Int_t fCanvasHeight = 800;
    Int_t kCanvasWidthThin = 300;
    Int_t kCanvasHeightThin = 400;

    Int_t kCanvasWidthWide = 600;
    Int_t kCanvasHeightWide = 400;

    TFile * fInputFile;
    TString sOutputDir;
    TString sOutputFileName;
    TString sListName = "";

    static const Int_t kNCentBins = 4;

    Int_t EventPlaneAxis = 7; 
    Int_t EventPlaneBinLow = 0;
    Int_t EventPlaneBinHigh = 2; 
    Int_t LambdaBinLow = 1;
    Int_t LambdaBinHigh = 5;
    Int_t EnergyBinLow = 1;
    Int_t EnergyBinHigh = -2; 
    Int_t AsymBinLow = 1;
    Int_t AsymBinHigh = -1; 
    Int_t OpeningAngleBinLow = 1;
    Int_t OpeningAngleBinHigh = -1; 
    Int_t bkgChoice = 6;
    Int_t fitMethod = 0;
    Int_t fitPeakMethod = 0;
    Int_t fitBkgMethod = 0;
    
    std::vector<TString> sPeakNames = {"Gaussian","ExpDecay/Gaussian","Brent-Wigner","Crystal Ball (Left Side)","Crystal Ball (Right Side)","ExpDecay/Gaussian (Right Side)","ExpDecay/Gaussian (Both Sides)","Voigt Profile"};
    std::vector<TString> sPeakFormula = {};

/*
FitPeakMethod: 6
#             0 Gaussian
#             1 ExpDecay/Gaussian
#             2 Brent-Wigner
#             3 Crystal Ball (Left Side)
#             4 Crystal Ball (Right Side)
#             5 ExpDecay/Gaussian (Right Side)
#             6 ExpDecay/Gaussian (Both Sides)
#             7 Voigt Profile
*/


    // Upper cut variables after interpreting -1 or -2 entries
    Int_t nOpeningAngleBinHigh = OpeningAngleBinHigh;
    Int_t nLambdaBinHigh = LambdaBinHigh;
    Int_t nEnergyBinHigh = EnergyBinHigh;


    std::vector <Double_t> equalBinHeights;

    Bool_t bDrawMassWindowLines = true;

    Bool_t scaleBkg = true; 

    Int_t iFixedMassWindows = 0; // (not implemented yet) Whether to use already determined mass window cuts. 0 uses freshly determined windows
    Double_t nSigma = 2;       // Size of Mass Window
    Double_t nSigmaR = -1;      // Size of Right side of mass window. if < 0, use symmetric window

    Int_t nRebinMass = 5;      // How much to rebin the mass plots // 5 is good
    // 4

    Int_t nRebinDeltaPsi = 6;  // Rebin the 60 bin DeltaPsi histograms

    Bool_t bApplyOpeningAngleCorrection = false; // An older correction attempt
    Bool_t bPtOpAngleCorr = false; // Whether to make angle correction in separate pt bins

    Bool_t bUsingClusPairRot = true; // this replaces pos swap

    //Int_t nPtBins = 8;         // How many pt bins
    Int_t nPtBins = 5;         // How many pt bins
    Int_t PtBinChoice = 1;
    Float_t kHighPtCut = 30.;  // Don't plot above here

    Int_t kUsedPi0TriggerPtBins = 5; // How many Pt bins do we actually use

    Double_t fOpeningAngleCut = 0; // This is set during LoadHistogram stage based on yaml

    Float_t kEtaDrawMin = 0.4;
    Float_t kEtaThreshold = 5.; //5.//11.; // when to start using eta peak
    Float_t kEtaDrawMax = 0.7;

    // Graphical Things
    string kMonthList[12] = {"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec."};
    // char *kMonthList[12] = {"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec."};
    Int_t kInvarMassColor = 1;
    Int_t kInvarMassStyle = kOpenSquare;
    Int_t kBkgSubColor = kRed+1;
    Int_t kBkgSubStyle = kFullCircle;
    Int_t kCombBkgColor   = kGray + 3;
    Int_t kCombBkgMarkerStyle = kFullSquare;

    Int_t kFitNpx = 256;

    float kMagicScale = 1.1;
    float kMagicNegativeScale = 0.03;

    Int_t kMassWindowLineStyle = 2;
    Int_t kMassWindowLineWidth = 2;
    Int_t kMassWindowLineColor = kOrange+8;

    Int_t kSBMassWindowLineStyle = 4;
    Int_t kSBMassWindowLineWidth = 2;
    Int_t kSBMassWindowLineColor = 8;

    Int_t kPi0BkgColor = kGreen+1; // This isn't drawn normally
    Int_t kPi0FitColor = kAzure-6; //kAzure+8

    Int_t kEtaPeakLineColor = kAzure+8; //kAzure-6
    Int_t kEtaPeakLineStyle = 1;
    Int_t kTotalFitColor = kOrange+10;
    Int_t kTotalBkgColor = kGray+1;

    Int_t kTotalFitLineWidth = 3;
    Int_t kTotalBkgLineWidth = 3;
    Int_t kTotalBkgLineStyle = 2;

    Int_t kCombBkgLineWidth = 3;

    Int_t kUnmodPi0FitColor = kOrange+10;
    Int_t kUnmodPi0FitLineStyle = 2;
    Int_t kUnmodPi0FitLineWidth = 3;

    Int_t kMCYieldTotalColor = kRed; //kRed+3;
    Int_t kYieldTotalRatioMarkerSize = 1;

    Int_t kRecMCYieldRatioColor = kBlack;
    Int_t kRecMCYieldRatioMarkerStyle = kFullCircle;

    //  std::vector<TString> sMCIdNames  = {"NoMatch","1MCPart","Pi02Gamma","Pi0Dalitz","Eta2Gamma","Eta3Pi0","EtaPi0Dalitz","EtaGammaDalitz","GammaPC","OtherShared"};
    std::vector<Int_t> kMCIdColor = {kRed+3,kSpring-1,kBlue,kAzure-2,kRed-4,kRed-4,kRed-4,kRed-4,kOrange+7,kMagenta};
    std::vector<Int_t> kMCIdStyle = {kOpenSquare,22,kFullCircle,kOpenCircle,kFullCircle,23,kOpenTriangleUp,kOpenCircle,kFullStar,kOpenStar};

    Int_t kMCNoPeakColor = kRed+2;
    Int_t kMCNoEtaColor = kOrange-7;
    Int_t kMCRotBkgPi0Color = kAzure - 3;
    Int_t kMCRotBkgEtaColor = kRed-4;

    Int_t kMCNoMatchColor  = kOrange+10; //kRed+3;
    Int_t kMCPi0MatchColor = kBlue; //kViolet - 3; 
    Int_t kMCEtaMatchColor = kAzure-2;
    Int_t kMCGammaPCMatchColor= kOrange+10;
    Int_t kMCSinglePartMatchColor= kSpring-1;
    Int_t kMCSharedEtaAncMatchColor= kGray+3;
    Int_t kMCSharedAncMatchColor= kMagenta;//kOrange+10;

    Int_t kMCNoPeakMarkerStyle = kFullSquare;
    Int_t kMCNoEtaMarkerStyle = kFullCircle;
    Int_t kMCRotBkgPi0MarkerStyle = kOpenTriangleUp;
    Int_t kMCRotBkgEtaMarkerStyle = kFullTriangleDown;

    Int_t kMCNoMatchMarkerStyle   = kOpenSquare;
    Int_t kMCPi0MatchMarkerStyle  = 22;
    Int_t kMCEtaMatchMarkerStyle  = kOpenCircle;
    Int_t kMCGammaPCMatchMarkerStyle = kFullTriangleDown;
    Int_t kMCSinglePartMatchMarkerStyle = kOpenTriangleUp;
    Int_t kMCSharedEtaAncMatchMarkerStyle = kFullDiamond;
    Int_t kMCSharedAncMatchMarkerStyle = kOpenStar;



    Int_t kPSPi0EnergyColor = kBlue-3;
    Int_t kPSPi0PosColor = kAzure+10;
    Int_t kPSEtaEnergyColor = kOrange+1;
    Int_t kPSEtaPosColor = kOrange+7;

    Int_t kRotBkgMinusMCPi0PosColor = kMagenta+1;
    Int_t kRotBkgMinusMCPi0PosStyle = kOpenSquare;

    Int_t kRotBkgMinusMCAllColor = kMagenta-9;
    Int_t kRotBkgMinusMCAllStyle = kOpenCircle;


    Bool_t havePatchStatus = false; // Whether or not the tree has the Patch Candidate Status axis. (added July 31, 2018)
    Int_t iPatchStatusAxis = 8; // may be 7

    Bool_t restrictToPatch = false;  // temporary switch to use only patch candidates
    Bool_t restrictPatchOnlyInME = false; // Only require GA Patch status in Mixed Events
    Int_t fSelectTriggerType = 0; // 1 for MB, 2 for GA
    // 0 will use the last Pi0 TList in the file

    Bool_t drawMCInfo = false; // Whether to draw MC info.  Could go into a different plot entirely.

    Bool_t RemoveMCPi0PS = false; // Whether to try removing the PosSwapped Pi0 peak (in MC)
    Bool_t RemoveMCEta = false; // Whether to remove MC Eta from Same Event and from Pos Swap (if applicable)

    Bool_t bZoomInvarMass = true; // Switch for zooming the invariant mass

    Double_t bkgScaleMin = 0.7; // Where to normalize 
    Double_t bkgScaleMax = 1.0; 
    Double_t fitDrawRangeMin = 0.05;
    Double_t fitDrawRangeMax = 0.75;

    // These don't really do anything anymore
    Int_t nThetaBins = 0;
    Int_t nHighLambdaBin = 6; // Lambda min for uncorrelated real clusters 
    int nHighLambdaBinMax = 6; // max is 6
    bool useLambdaBkgOpenAngle = false;
    Int_t nLowPtBinMin = 4; // Pt range for finding the angular correction factors
    Int_t nLowPtBinMax = 4;
    Int_t nLowPtEnergyBinMin = 4; // min energy for the low pt angular correction //2 ?
    Bool_t usePtBkgOpenAngle = false; // Set to false to use the entire pt range for angular corrections

    Bool_t drawFits = true;
    Bool_t bDrawUnmod = true;
    Bool_t drawResBkg = true; //Whether to draw the residual bkg fit
    Bool_t drawBkg = true; // Whether to draw the ME/Rot background

    Double_t kMCPi0WindowLow = 0.05; // Obsolete window for MC pi0 peak
    Double_t kMCPi0WindowHigh = 0.3;

    Double_t fMCRescaling = 1.0; // How much to rescale Other Shared Ancestor (Jet) MC component
    // Nothing done if 1.0

    Bool_t bEnablePerformance = true;

    // These will be set automatically based on the input file
    Bool_t haveRotBkg = false;
    Bool_t haveMEBkg  = false;
    Bool_t havePSBkg  = false;
    Bool_t haveMCStatus = false;

    // These two numbers are updated as needed so we don't store the +
    Int_t nSkipPoints = 0; // number of pt bins with zero content
    Int_t iFirstRealBin = -1; // index of first bin with content

    // These are set based on other input settings
    bool fitBkgSub = false; // Whether to subtract the background from the histograms first.  Otherwise, they will be ignored or included as a term in the fit
    int bkgType = 3;

    // MC Label Information
    Int_t iMCAxis = 7;
    Int_t nMCId = 10; // This does not include 11-14, the labels for rotation
    std::vector<TString> sMCIdNames  = {"NoMatch","1MCPart","Pi02Gamma","Pi0Dalitz","Eta2Gamma","Eta3Pi0","EtaPi0Dalitz","EtaGammaDalitz","GammaPC","OtherShared"};
    std::vector<TString> sMCIdTitles = {"No Match","1 MC Part.","#pi^{0}#rightarrow#gamma#gamma","#pi^{0}#rightarrow#gammae^{+}e^{-}","#eta#rightarrow#gamma#gamma","#eta#rightarrow3#pi^{0}","#eta#rightarrow#pi^{0}#pi^{+}#pi^{-},#pi^{0}#gamma#gamma","#eta#rightarrow#gamma#pi^{+}#pi^{-},#gammae^{+}e^{-}","#gamma#rightarrowe^{+}e^{-}","Other Shared Ancestor"};

    // Important Shared Histograms

    THnSparse * Pi0Cands = 0;
    THnSparse * ClusterProp = 0;


    TH1F * fClusEnergy = 0;
    TH2F * fHistEvsPt = 0;
  
    // The Main Two
    TH2F * fInvarMasspT = 0;
    TH2F * fInvarMasspTRotBkg = 0;

    TH2F * fInvarMasspTRotBkgAngleScaled = 0;
    TH2F * fInvarMassPtRotBkgGAPatch = 0;
    TH2F * fInvarMassPtRotBkgNoPatch = 0;
    TH2F * fInvarMasspTRotBkgMCPi0 = 0;
    TH2F * fInvarMasspTRotBkgMCEta = 0;
    TH2F * fInvarMasspTRotBkgMCPi0EnergyPair = 0;
    TH2F * fInvarMasspTRotBkgMCEtaEnergyPair = 0;
    TH2F * fInvarMasspTRotBkgMCPi0PosPair = 0;
    TH2F * fInvarMasspTRotBkgMCEtaPosPair = 0;


    TH2F * fMatchDeltaPhiTrackPt = 0;
    TH2F * fMatchDeltaEtaTrackPt = 0;
    TH2F * fMatchCondDeltaPhiTrackPt = 0;
    TH2F * fMatchCondDeltaEtaTrackPt = 0;

    // Subtracting MCPi0Pos 
    TH2F * fInvarMasspTRotBkgMinusMCPi0Pos = 0;
    TH2F * fInvarMasspTRotBkgMinusMCAll = 0;

    std::vector<TH1D *> hInvarMasspTBin;
    // Residual to fit
    std::vector<TH1D *> hInvarMasspTBinResid;
    // Bkg Method Histograms
    std::vector<TH1D *> hInvarMassBkgSubpTBin;
    std::vector<TH1D *> hInvarMasspTBinRotBkg;
    std::vector<TH1D *> hInvarMassPtBinRotSub;
    // MC histograms
    std::vector<std::vector<TH1D *>> hInvarMassPtBinMCId;
    std::vector<TH1D *> hInvarMassPtBinMCNoPeak;
    std::vector<TH1D *> hInvarMassPtBinMCNoEta;
    // bkg MC info
    // The old Method (combined)
    std::vector<TH1D *> hInvarMasspTBinRotBkgMCPi0;
    std::vector<TH1D *> hInvarMasspTBinRotBkgMCEta;

    // The New Method
    std::vector<TH1D *> hInvarMasspTBinRotBkgMCPi0EnergyPair;
    std::vector<TH1D *> hInvarMasspTBinRotBkgMCEtaEnergyPair;

    std::vector<TH1D *> hInvarMasspTBinRotBkgMCPi0PosPair;
    std::vector<TH1D *> hInvarMasspTBinRotBkgMCEtaPosPair;

    std::vector<TH1D *> hInvarMasspTBinRotBkgMinusMCPi0PosPair;
    std::vector<TH1D *> hInvarMasspTBinRotBkgMinusMCAll;

    std::vector<TH1D *> hInvarMasspTBinRotBkgMinusMCPi0PosPairRatio; // Ratio over true MC backgruond

    // MC 2D Histos
    std::vector<TH2F *> fInvarMassPtMCId = {};
    std::vector<TH2F *> fInvarMassPtRotBkgMCId = {};

    TH2F * fInvarMassPtMCNoPeak = 0;
    TH2F * fInvarMassPtMCNoEta = 0;  // No eta contribution
    TH2F * fInvarMassPtMCNoEtaBkg = 0;  // Only eta contribution from eta->2g
    TH2F * fInvarMasspTNoCuts = 0;
    TH1D * hPairPt = 0;
    TH1D * hPairPtRotBkg = 0;
    TH1D * hPairOpeningAngle = 0;
    TH1D * hPairOpeningAngleRotBkg = 0;

    // For pt dependent opening angle
    TH2F * hPairPtOpAngle = 0;
    TH2F * hPairPtOpAngleBkg = 0;

    TH1D * hPairLambdaBkgOpeningAngle = 0;
    TH1D * hPairPtBkgOpeningAngle = 0;
    TH2F * hPtAngleDetail, * hPtAngleDetailRotBkg;
    TH2F * hAngleEnergyCut, * hAngleEnergyCutRotBkg;


    // Delta Psi Histograms
    // Event Plane
    TH2F * hPtEPAnglePionAcc=0;       // Pi0 Candidate
    TH2F * hPtEPAngleMCPion=0;        // MC Pi0
    TH2F * hPtEPAngleTrueRecMCPion=0; // MC Pi0 reconstructed as Pi0 Candidate
    TH3F * hHistTrackPsiEPPtCent=0;   // Tracks
    // Reaction Plane (Angle 0 in Data)
    TH2F * hPtRPAnglePionAcc=0;
    TH2F * hPtRPAngleMCPion=0;
    TH2F * hPtRPAngleTrueRecMCPion=0;
    TH3F * hHistTrackPsiRPPtCent=0;

    // Projections of the above
    std::vector<TH1F *> hPtEPAnglePionAcc_Proj;
    std::vector<TH1F *> hPtEPAngleMCPion_Proj;
    std::vector<TH1F *> hPtEPAngleTrueRecMCPion_Proj;

    std::vector<TH1F *> hPtRPAnglePionAcc_Proj;
    std::vector<TH1F *> hPtRPAngleMCPion_Proj;
    std::vector<TH1F *> hPtRPAngleTrueRecMCPion_Proj;

    TH2F * hHistTrackPsiEPPt=0;   // Tracks
    TH2F * hHistTrackPsiRPPt=0;
    std::vector<TH1F *> hPtEPAngleTrack_Proj;

    // Vn Measurements
    // track pt info, so far only used for studying v_n
    const int kNTrackPtBins = 11;
    std::vector <double> fTrackPtBins = {0.15, 0.25, 0.5,   1, 1.5, 2,   3, 4, 5,   6, 8, 10};


    TGraphErrors * gTrigger_Bv = 0;
    TGraphErrors * gTrigger_V2 = 0;
    TGraphErrors * gTrigger_V4 = 0;
    TGraphErrors * gTrigger_V6 = 0;

    TGraphErrors * gTrack_Bv = 0;
    TGraphErrors * gTrack_V2 = 0;
    TGraphErrors * gTrack_V4 = 0;
    TGraphErrors * gTrack_V6 = 0;

    // Functions and Data Points
    std::vector <double> Pi0PtBins;

    std::vector<TF1  *> fPi0Fit;
    std::vector<TF1  *> fUnmodPi0Fit;
    std::vector<TF1  *> fPi0Peak;
    std::vector<TF1  *> fPi0Bkg;
    std::vector<TF1  *> fEtaPeak;

    std::vector<TH1D *> hTotalBkg; //add background function 
    std::vector<TH1D *> hTotalFit;

    std::vector<TF1  *> fMCPi0Fit;


    std::vector<double> ptPointsForTGraph = {};
    std::vector<double> ptErrorsForTGraph = {};

    std::vector<double> Pi0MassCutLow = {};
    std::vector<double> Pi0MassCutHigh = {};

    std::vector<double> SBMassCut0 = {};
    std::vector<double> SBMassCut1 = {};
    std::vector<double> SBMassCut2 = {};
    std::vector<double> SBMassCut3 = {};
    std::vector<double> SBMassCut4 = {};


    std::vector<double> Pi0MassArr = {};
    std::vector<double> Pi0MassArrUn = {};
    std::vector<double> Pi0SigmaArr = {};
    std::vector<double> Pi0SigmaArrUn = {};

    // Parameter Yields
    std::vector<double> Pi0YieldArr = {};
    std::vector<double> Pi0YieldArrUn = {};
    std::vector<double> Pi0NormYieldArr = {};
    std::vector<double> Pi0NormYieldArrUn = {};
    // Integral Yields
    std::vector<double> Pi0IntegralArr = {};
    std::vector<double> Pi0IntegralArrUn = {};
    std::vector<double> Pi0IntegralTotalArr = {}; // Total = Signal + Background (or data without subtraction)
    std::vector<double> Pi0IntegralTotalArrUn = {};
    std::vector<double> Pi0NormIntegralArr = {};
    std::vector<double> Pi0NormIntegralArrUn = {};

    // MC Integral Yields (directly using MC info)
    std::vector<double> Pi0MCIntegralArr = {};
    std::vector<double> Pi0MCIntegralArrUn = {};
    std::vector<double> Pi0MCNormIntegralArr = {};
    std::vector<double> Pi0MCNormIntegralArrUn = {};

    // Scaling used for combinatorial background (from specific range normalization or from fit)
    std::vector<double> Pi0BkgScaleArr = {};
    std::vector<double> Pi0BkgScaleArrUn = {};

    // Background Estimation
    std::vector<double> Pi0BkgArr = {};  // How to get uncertainty?
    std::vector<double> Pi0BkgArrUn = {};

    // MC Background
    std::vector<double> MCBkgArr = {};
    std::vector<double> MCBkgArrUn = {};
    std::vector<double> MCPi0PeakSigRatioArr = {};
    std::vector<double> MCPi0PeakSigRatioArrUn = {};
    std::vector<double> MCPi0YieldBkgRatioArr = {};
    std::vector<double> MCPi0YieldBkgRatioArrUn = {};
    std::vector<double> MCPi0YieldTotalRatioArr = {};
    std::vector<double> MCPi0YieldTotalRatioArrUn = {};

    std::vector<double> RecMCPi0YieldRatioArr = {}; // reconstructed (rec. yield / MC yield) ratio
    std::vector<double> RecMCPi0YieldRatioArrUn = {};

    std::vector<double> Pi0PeakSigRatioArr = {};
    std::vector<double> Pi0PeakSigRatioArrUn = {};
    std::vector<double> Pi0YieldBkgRatioArr = {};
    std::vector<double> Pi0YieldBkgRatioArrUn = {};
    std::vector<double> Pi0YieldTotalRatioArr = {};
    std::vector<double> Pi0YieldTotalRatioArrUn = {};

    std::vector<double> Pi0ChiSquareArr = {};
    std::vector<double> MCPi0ChiSquareArr = {};

    // Output Objects

    TH1D * hOpeningAngleCorrection = 0;
    TH2F * h2DOpeningAngleCorr = 0; // now with pt dependence!!


    TGraphErrors * Pi0Yield = 0;
    TGraphErrors * Pi0Spectrum = 0;
    TGraphErrors * Pi0IntYield = 0;
    TGraphErrors * Pi0IntTotal = 0;
    TGraphErrors * Pi0IntSpectrum = 0;

    TGraphErrors * Pi0MCIntYield = 0;
    TGraphErrors * Pi0MCIntSpectrum = 0;

    TGraphErrors * Pi0Bkg = 0;
    TGraphErrors * Pi0PeakSigRatio = 0;
    TGraphErrors * Pi0YieldBkgRatio = 0;
    TGraphErrors * Pi0YieldTotalRatio = 0;

    TGraphErrors * MCBkg = 0;
    TGraphErrors * MCPi0PeakSigRatio = 0;
    TGraphErrors * MCPi0YieldBkgRatio = 0;
    TGraphErrors * MCPi0YieldTotalRatio = 0;
    TGraphErrors * RecMCPi0YieldRatio = 0;

    TGraphErrors * Pi0Mass = 0;
    TGraphErrors * Pi0Sigma = 0;
    TGraphErrors * Pi0ChiSquare = 0;
    TGraphErrors * MCPi0ChiSquare = 0;

    TF1 * pi0MassFit = 0;
    TF1 * pi0SigmaFit = 0; 

/*
    THnSparse * fUDist = 0;
    THnSparse * fUTildeDist = 0;
    THnSparse * fVDist = 0;
    THnSparse * fVTildeDist = 0;
*/
    THnSparse * fUScaleMatrix = 0;
    THnSparse * fVScaleMatrix = 0;
    THnSparse * fPSMassPtMap = 0;
    THnSparse * fESMassPtMap = 0;

    TH2F * fU2DCut = 0;
    TH2F * fV2DCut = 0;
    TH2F * fPSInitialMap = 0;
    TH2F * fPSFinalMap = 0;

    // A histogram of the predicted transformed peak
    TH2F * fPSFinalPeak = 0;
    // Array of predicted transformed peaks
    std::vector<TH1D *> hInvarMasspTBinPSCorr;
    // Array of post subtraction rotbkg histograms
    std::vector<TH1D *> hInvarMasspTRotBkgPSSub;

    // Model for the input signal
    TF2 * W2DSignalModel = 0;
    TH2 * W2DSignalHist  = 0;
    // Initial Signal Parameters
    Double_t fWPower = 4;
    Double_t fWYield = 1e3;
    Double_t fWPt0    = 4.5;
    // These may be made pT dependent
    Double_t fWMass  = 0.135;
    Double_t fWSigma = 0.016;

    // Mass Dependent Features


    Double_t fWD_Mass = 9.042345e+00;
    Double_t fWE_Mass = 1.387551e-01;
    Double_t fWM1_Mass = 1.038669e-02;
    Double_t fWM2_Mass = 2.903288e-03;
  
    Double_t fWD_Sigma = 9.042345e+00;
    Double_t fWE_Sigma = 1.232047e-02;
    Double_t fWM1_Sigma = -3.715172e-03;
    Double_t fWM2_Sigma = 4.273227e-04;


    // First version
    /*
    Double_t fWD_Mass = 10.816;
    Double_t fWE_Mass = 0.141;
    Double_t fWM1_Mass = 0.003;
    Double_t fWM2_Mass = 0.002;
  
    Double_t fWD_Sigma = 12.;
    Double_t fWE_Sigma = 0.014;
    Double_t fWM1_Sigma = 0.001;
    Double_t fWM2_Sigma = -0.001;
    */



    TH2F * fPSEnergyPairDist = 0;
    TH2F * fPSPosPairDist    = 0;
    
    Bool_t bPSCorrLogMod = true;



    // Fixed Mass Windows:
    // These are the values for GA Data on 20191202
    // Lambda Range: [0.10 - 0.70]     Bins: 1 7
    // Energy Range: [2.00 - 100.00]     Bins: 1 3
    // Asym Range:   [0.00 - 0.70]     Bins: 1 1
    // OpeningAngle Range:   [0.017 - 3.142]     Bins: 3 12
    // Track Cluster Correction: Subtraction

    // note: last three bins are basically fake
    Double_t fPi0MassFixedValue_3[kNCentBins][9] = {
      { 0.129265, 0.130740, 0.138696, 0.147779, 0.160138, 0.177054, 0.177054, 0.177054, 0.177054},
      { 0.131305, 0.131078, 0.137113, 0.143671, 0.153884, 0.171041, 0.171041, 0.171041, 0.171041},
      { 0.130332, 0.129157, 0.134172, 0.140263, 0.151399, 0.162242, 0.162242, 0.162242, 0.162242},
      { 0.127595, 0.128214, 0.132793, 0.137396, 0.148063, 0.167322, 0.167322, 0.167322, 0.167322}
    };
    Double_t fPi0SigmaFixedValue_3[kNCentBins][9] = {
      { 0.019000, 0.019000, 0.017793, 0.017149, 0.017567, 0.011105, 0.011105, 0.011105, 0.011105},
      { 0.019000, 0.019000, 0.013467, 0.012001, 0.014112, 0.013643, 0.011105, 0.011105, 0.011105},
      { 0.019000, 0.016856, 0.010454, 0.010690, 0.011508, 0.011226, 0.011226, 0.011226, 0.011226},
      { 0.019000, 0.003000, 0.009186, 0.011264, 0.014117, 0.012515, 0.012515, 0.012515, 0.012515}
    };









  private:
    PionID(const PionID&); // not imp.
    ClassDef(PionID,1);
};



#endif
