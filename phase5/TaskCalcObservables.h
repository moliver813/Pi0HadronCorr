#ifndef TASKCALCOBSERVABLES_H
#define TASKCALCOBSERVABLES_H

#include <Riostream.h>
#include <TString.h>
using namespace std;
#include <vector>


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
class TLatex;


class TaskCalcObservables : public TObject {

public:

  TaskCalcObservables();

  virtual ~TaskCalcObservables() { ; }

  void SetStyle();
  void LoadHistograms();
  void LoadSystematics();
  void LoadModels();
  void CalculateSystematics();
  void InitArrays();

  void SetPtBin(Int_t input)              { iPtBin       = input; }
  void SetObservable(Int_t input)         { iObservable  = input; }

  void SetCentralInputFile(TFile * inputFile) { fInputFileCentral = inputFile; }

  void AddSystematicComparison(TFile * inputFile, TString fLabel) {
    printf("Adding systematic file for %s\n",fLabel.Data());
    //fSystematicsNames.push_back(fLabel);
    TString sFormattedLabel = fLabel;
    sFormattedLabel.ReplaceAll("_"," ");
    fSystematicsNames.push_back(sFormattedLabel);
    fSystematicsFiles.push_back(inputFile);
  }

  void AddModelComparison(TFile * inputFile, TString fLabel, TString fTitle) {
    printf("Adding model file %s\n",fLabel.Data());
    fModelNames.push_back(fLabel);
    TString sFormattedTitle = fTitle;
    sFormattedTitle.ReplaceAll("_"," ");
    fModelTitles.push_back(sFormattedTitle);
    //fModelTitles.push_back(fTitle);
    fModelFiles.push_back(inputFile);
  }

  void AddSecondaryAxis(TString fLabel) {
    fInputFilesSecondary.push_back({});
    fAxesSecondary.push_back(fLabel);
  }

  // Could also add abilitiy to use string instead of iAxis
  void AddSecondaryFile(TFile * inputFile, int iAxis, TString fLabel = "") {
    //fInputFilesSecondary.push_back(inputFile);
    //fLabelsSecondary.push_back(fLabel);
    if (iAxis > (int) fInputFilesSecondary.size()) {
      fprintf(stderr,"Error, this secondary axis doesn't exist yet\n");
      return;
    }
    fInputFilesSecondary[iAxis].push_back(inputFile);
    fLabelsSecondary[iAxis].push_back(fLabel);
  }

 
  void SetPlotOptions(TString input)      { fPlotOptions = input; }
  void SetOutputDir(TString input)        { fOutputDir   = input; }
  void SetOutputFile(TFile * outputFile)  { fOutputFile  = outputFile; }

  TLegend * DrawAliceLegend(TObject *obj, Float_t x, Float_t y, Float_t x_size, Float_t y_size);


  void SetMCGenMode(Int_t input = 1)      { fIsMCGenMode = input; }
  Int_t GetMCGenMode()                    { return fIsMCGenMode;  }
  void SetCentralityBin(Int_t input)      { iCentBin     = input; }
  void SetRPFMethod(Int_t input)          { iRPFMethod   = input; }


  void SetSigmaRangeNS(double input)      { fSigmaRangeNS = input; }
  void SetSigmaRangeAS(double input)      { fSigmaRangeAS = input; }

  void SetSigmaRangeNSBinChange(int input) { iSigmaNSRangeBinChange = input; }
  void SetSigmaRangeASBinChange(int input) { iSigmaASRangeBinChange = input; }

  void SetNSkipPoints(Int_t input)        { nSkipPoints  = input; }
  int  GetNSkipPoints()                   { return nSkipPoints; }

  void Debug(Int_t input);
  void SetDebugLevel(Int_t input)         { fDebugLevel = input; }

  TString GetLabel() { return sLabel; }
  TString GetLabel2() { return sLabel2; }
  TString GetMyTitle() { return sTitle; }
  TString GetMyTitle2() { return sTitle2; }

  void SetLabel(TString input) { sLabel = input; }
  void SetLabel2(TString input) { sLabel2 = input; }
  void SetMyTitle(TString input) { sTitle = input; }
  void SetMYTitle2(TString input) { sTitle2 = input; }

  void SaveOutput();
  void Run();

protected:

  TGraphErrors * CalculateSystematicIndiv(TGraphErrors * graph, vector<TGraphErrors *> graph_sys_errors);


  TH1D * ProduceSystematicFromHists(vector<TH1D*> fHistArray, TH1D * fInputHist);
  TGraphErrors * ProduceSystematicFromGraphs(vector<TGraphErrors*> fGraphArray, TGraphErrors * fInputGraph);


  void ProduceSystematicUncertPlots();
  void ProduceSystematicUncertPlotIndiv(TGraphErrors * fCentral, vector<TGraphErrors *> fSysErrors, TGraphErrors * fTotalSysErrors, float legendX, float legendY);

  TString sLabel  = "";
  TString sLabel2 = "";
  TString sTitle  = "";
  TString sTitle2 = "";


  // Constants
  static const Int_t kNEPBins=3;

  static const Int_t kGammaNBINS=9;  //9    ///< Number of 2D histograms for Gamma energy
  static const Int_t kZtNBINS=7;            ///< Number of 2D histograms for Zt of g-h pair
  static const Int_t kXiNBINS=8;            ///< Number of 2D histograms for Xi of g-h pair
  static const Int_t kNoHPtBins=8;          ///< Bins in hadron pT

  static const Int_t kCentBINS= 4;

  const Int_t kFitLineColor = kViolet-5;
  const Float_t kOmniMarkerSize = 0.5;

  const Int_t kNonSelectColor = kGray;
  const Int_t kSelectColor = kBlack;

  // General Model Color/Style
  vector<Int_t> kModelColors = {kRed,4,kSpring-1,6,kRed+2,kOrange+10, 6, 8, 9 , 42};
  vector<Int_t> kModelStyles = {kFullSquare,kOpenSquare,kFullCircle,kOpenCircle,kFullStar,kOpenStar,kOpenCross,kFullCross};
  vector<Int_t> kModelLineStyles = {kSolid,kDashed,kDotted,kDashDotted,kSolid,kDashed,kDotted,kDashDotted,kSolid,kDashed,kDotted,kDashDotted};

  // Event Plane Bin Styles
  // in, mid, out, incl
  vector<vector<Int_t>> kModelColorsEP = {{kAzure+1,kAzure+0,kCyan+2,kViolet+10},
  //Int_t kModelColorsEP[4][4] = {{kAzure+1,kAzure+0,kCyan+2,kViolet+10},
    {kGreen-4,kTeal+3,kSpring-6,kSpring+9},
    {kRed-7,kPink-3,kOrange+8,kOrange+2},
    {kGray+3,kGray+2,kGray+1,kGray}
    };


  // Systematic Error Source color/sizes
  vector<Int_t> kSysErrorColor = {kRed-7,kAzure+1,kPink-3,kAzure+0,kOrange+8,kCyan+2,kOrange+2,kViolet+10};
  vector<Int_t> kSysErrorStyle = {kFullSquare,kOpenSquare,kFullCircle,kOpenCircle,kFullStar,kOpenStar,kOpenCross,kFullCross};

  void DrawDirectComparisons();

  void CalculateResults();                   ///< Loop over obs bins
  void CalculateResultsObsBinEPBin(int iObsBin, int iEPBin, TCanvas * canv);  ///< Calculate yields, widths for this bin in iObsBin, EP bin

  void CalculateRPFErrorYieldsRms();

  void CreateRatioAndDifferencesGraphs();

  // Calculate Ratios
  void CalculateResultsRatiosObsBin(int iObsBin, TCanvas * canv);

  void CalculateRPFErrorRatios();

  // Calculate Differences
  void CalculateResultsDifferencesObsBin(int iObsBin, TCanvas * canv);

  void CleanResults(); ///< Remove nSkipPoints, where the analysis is ineffective

  void AddTrackingUncertainties();
  //void AddTrackingUncertaintiesIndiv(TGraphErrors * graph, vector<TGraphErrors *> fErrorsArray, double value);
  void AddTrackingUncertaintiesIndiv(TGraphErrors * graph, vector<TGraphErrors *> * fErrorsArray, double value);
  TGraphErrors * BuildConstantGraph(TGraphErrors * graph, double value);

  void CombineStatisticalErrors();

  void CombineStatisticalErrorsIndiv(TGraphErrors * fRawStatErrors, TGraphErrors * fRPFErrors);

  void SetGraphStyles(); ///< Set the colors, markers on the graphs

 
  void DrawResults(); ///< Draw the plots

  //void DrawObservable(vector<TGraphErrors *> fObsGraphs, vector<TGraphErrors *> fObsRPFErrors, vector<TGraphErrors *> fObsSysErrors = {});
  // Draws an observable with stat, RPF, and Sys errors, and then with models
  void DrawObservable(vector<TGraphErrors *> fObsGraphs, vector<TGraphErrors *> fObsRPFErrors, vector<TGraphErrors *> fObsSysErrors = {}, vector<vector<TGraphErrors *>> fModels = {});

  void DrawFinalResults(); ///< Draw the final plots

  void DrawRatio();


  // Draw Observable or Ratio
  void DrawFinalObservable(TGraphErrors * fObsGraph, TGraphErrors * fObsGraphSysErrors, vector<TGraphErrors*> fModels, TCanvas * cFinal, TString sName, TString sTitle);

//  CalculateYield(TH1D * hist, bool bAwaySide, double &Error);


  int fYieldCanvasWidth = 900;
  int fYieldCanvasHeight = 600;  

  //TF1 * FitAndCalculateSigmaFromHisto(TH1D * hist, double &Sigma, double &SigmaErr);
  TF1 * FitAndCalculateSigmaFromHisto(TH1D * hist, double *Sigma, double *SigmaErr, int iSide);


private:

  Int_t fIsMCGenMode;                       ///< 0 = data, 1 = mcGen mode

  Int_t fDebugLevel;                       ///< For Debugging Purposes

  Bool_t bFitSigmaSlices=true;            ///< True = use FitSlices for fitting the sigmas of the variants. False = project and fit each slice as a TH1


  Float_t fGlobalTrackingUncertainty = 0.04; ///< Global tracking efficiency uncertainty
  Float_t fTrackingEventPlaneUncertainty = 0.0042; ///< Uncertainty on how much the tracking efficiency may change between event plane angles.



  TFile *fInputFileCentral;               ///< File for central (systematic) file

  // This method may be obsolete; the variations will be done and compared with sysCompare
  vector<vector<TFile *>> fInputFilesSecondary;   ///< Files for systematically changing things
    // 1st index is the axis (v3 change, purity change, etc)
    // 2nd index is just within that list

  vector<TString> fAxesSecondary;     ///< Labels for each "axis" of the variations
  vector<vector<TString>> fLabelsSecondary;       ///< Corresponding labels
  // Could make this have more dimensions, for being more differential

  vector<TString> fSystematicsNames = {};
  vector<TFile *> fSystematicsFiles = {};

  vector<TString> fModelNames = {};
  vector<TString> fModelTitles = {};
  vector<TFile *> fModelFiles = {};

  TString fOutputDir;                       ///< Output directory to save the plots
  
  Int_t iPtBin = 4;                         ///< Which Pt Bin is used (for fixed V_n^t values
  // update use iPtBin also for restricting pt trigger distributions

  Double_t PtBins[6] = {5,7,9,11,14,17}; // might be able to get this information from a pt trigger histograms
  Int_t iObservable = 1;                    ///< Which observable (0=Pt,1=Zt) is being iterated over

  Int_t nRebinDPhi = 1;                     ///< Rebinning in Delta Phi

  Int_t iCentBin = 0;                           ///< Which centrality bin (0-10,10-30,30-50,50-80)

  Int_t iRPFMethod = 0;                         ///< Which RPF method to use. 0 is my c++, 1 is Raymond's python implementation

  TString fPlotOptions="COLZ";              ///< Style for plotting 3D histograms.  Default is colz

  const bool fPreliminary = false;          ///< Switch on in case we get performance approval

//  TFile *fInputCentral;                     ///< File with central analysis

  TFile *fOutputFile;                       ///< File where output objects are to be stored

  Int_t fObservable;                        ///< Observable for the current analysis
  TString fObservableName;                  ///< Name of the current observable (for plot labels)

  Int_t nObsBins;                           ///< How many bins we have of the current observable. Determined at run time.

  vector<Double_t> fObsBins;                ///< Bin Edges for whichever observable we are looking at

  TString fEPBinTitles[kNEPBins+1] = {"In-Plane","Mid-Plane","Out-of-Plane","All EP Angles"};
 // TString fEPBinTitles[kNEPBins+1] = {"In-Plane","Mid-Plane","Out-of-Plane","All EP Angles"};
  TString fEPBinTitlesShort[kNEPBins+1] = {"EP0","EP1","EP2","Incl"};

  Float_t kAliceLegendWidth=0.225;
  Float_t kAliceLegendHeight=0.2;

  //Int_t kEPColorList[4] = {kBlack, kBlue-4, kGreen-3, kRed+1};
  //Int_t kEPMarkerList[4] = {kOpenSquare, kFullSquare, 39, kFullDiamond};
  //Int_t kEPColorList[4] = {kBlue-4, kGreen-3, kRed+1,kBlack};
  Int_t kEPColorList[4] = {kBlue-4, kGreen+2, kRed+1,kBlack};
  Int_t kEPMarkerList[4] = {kFullSquare, 39, kFullDiamond, kOpenSquare};

  Float_t kOutInRPFErrorAlpha = 0.4;
  Float_t kMidInRPFErrorAlpha = 0.4;
  Float_t kOutInErrorAlpha = 0.4;
  Float_t kMidInErrorAlpha = 0.4;

  // Red / Blue
  Int_t kOutInColor = kViolet-1;
  // Green / Blue
  Int_t kMidInColor = kEPColorList[1];

  Int_t kOutInRPFErrorColor = kOutInColor;
  Int_t kMidInRPFErrorColor = kMidInColor;
  Int_t kOutInErrorColor = kOutInColor+5;
  Int_t kMidInErrorColor = kMidInColor+5;

  Int_t kOutOverInMarker = kFullSquare;
  Int_t kMidOverInMarker = kFullCircle;

  


  // Red / Blue, Green / Blue
//  int kOutInErrorColor = kOutInColor;
//  int kMidInErrorColor = kMidInColor;


  // Sys Error Fills
  int kOutInSysFillStyle = 3001;
  int kMidInSysFillStyle = 3002;
  int kOutInNSSysFillStyle = 3001;
  int kMidInNSSysFillStyle = 3002;



  // RPF Error Fills
  int kOutInFillStyle = 3245;
  int kMidInFillStyle = 3254;
  int kOutInNSFillStyle = 3245;
  int kMidInNSFillStyle = 3254;

  // RPF errors
  Int_t kEPRPFFillColorList[4] = {kBlue-4, kGreen+2, kRed+1,kGray};
  //Int_t kEPRPFFillColorList[4] = {kBlue-4, kGreen-3, kRed+1,kBlack};
  //Int_t kEPRPFFillStyleList[4] = {3245,3254,3245,3254};
  Int_t kEPRPFFillStyleList[4] = {3004,3005,3004,3005};

  Int_t kEPSysFillColorList[4] = {kBlue-7, kGreen+4, kRed+3,kBlack};
  Int_t kEPSysFillStyleList[4] = {3002,3002,3002,3002};
  //Int_t kEPSysFillStyleList[4] = {3004,3005,3004,3005};


  Int_t nSkipPoints = 0;  // How many of the first points to skip

//  Double_t fRMSRange = 1.047; // Range in DeltaPhi to calulate truncated RMS

  // May want to vary these for systematic uncertainty

  Double_t fRmsRangeNS = TMath::Pi()/3.;
  Double_t fRmsRangeAS = TMath::Pi()/3.;

  Double_t fSigmaRangeNS = TMath::Pi()/3.;
  Double_t fSigmaRangeAS = TMath::Pi()/3.;

  Double_t fYieldRangeNS = TMath::Pi()/3.;
  Double_t fYieldRangeAS = TMath::Pi()/3.;

  // variables for varying the ranges by a number of bins in delta phi
  // 
  int iSigmaNSRangeBinChange  = 0;
  int iSigmaASRangeBinChange  = 0;


  // The histograms to integrate

  // Delta Eta Projections
  vector<TH1D *> fNearSideSubDEtaFinalAll;    ///< SB Sub Nearside Delta Eta (all EP)
  vector<vector<TH1D *>> fNearSideSubDEtaFinalEP;    ///< SB Sub Nearside Delta Eta (EP Bins)

  // First index is observable bin, 2nd is the event plane bin
  vector<vector<TH1D *>> fFullDPhiProj_Sub;  ///< Full projections in DPhi after subtracting flow and rescaling for number of triggers in the EP bin
  vector<vector<TH1D *>> fNearEtaDPhiProj_Sub;  ///< Full projections in DPhi after subtracting flow and rescaling for number of triggers in the EP bin
  vector<vector<TH1D *>> fFarEtaDPhiProj_Sub;  ///< Full projections in DPhi after subtracting flow and rescaling for number of triggers in the EP bin

  // TH2s with variations of the RPF 
  // Delta Eta Projections
  // I don't apply RPF to these?
  //vector<TH2D *> fNearSideSubDEtaFinalAll_RPFVar;  
  //vector<vector<TH2D *>> fNearSideSubDEtaFinalEP_RPFVar;

  // First index is observable bin, 2nd is the event plane bin
  vector<vector<TH2F *>> fFullDPhiProj_Sub_RPFVar;
  vector<vector<TH2F *>> fNearEtaDPhiProj_Sub_RPFVar;
  vector<vector<TH2F *>> fFarEtaDPhiProj_Sub_RPFVar;


  // Items for applying errors from reaction plane fit
  int nRPFVariants = 0; 
 
  // Primary observable graphs

  TGraphErrors * fNSYieldsInc; ///< Near-Side Yields in all EP
  vector<TGraphErrors *> fNSYieldsEP; ///< Near-Side Yields in the EP bins
  TGraphErrors * fASYieldsInc; ///< Away-Side Yields in all EP
  vector<TGraphErrors *> fASYieldsEP; ///< Away-Side Yields in the EP bins

  // TGraphErrors with the error bars as the calculated uncertainty from RPF variations
  TGraphErrors * fNSYieldsInc_RPFError; ///< Near-Side Yields in all EP
  vector<TGraphErrors *> fNSYieldsEP_RPFError; ///< Near-Side Yields in the EP bins
  TGraphErrors * fASYieldsInc_RPFError; ///< Away-Side Yields in all EP
  vector<TGraphErrors *> fASYieldsEP_RPFError; ///< Away-Side Yields in the EP bins

  // TGraphErrors with the error bars as the calculated uncertainty from systematic variations
  TGraphErrors * fNSYieldsInc_SysError = 0; ///< Near-Side Yields in all EP
  vector<TGraphErrors *> fNSYieldsEP_SysError={}; ///< Near-Side Yields in the EP bins
  TGraphErrors * fASYieldsInc_SysError = 0; ///< Away-Side Yields in all EP
  vector<TGraphErrors *> fASYieldsEP_SysError={}; ///< Away-Side Yields in the EP bins

  vector<TGraphErrors *> fNSYieldsInc_SysErrorBySource= {}; ///< Near-Side Yields in all EP
  vector<vector<TGraphErrors *>> fNSYieldsEP_SysErrorBySource= {}; ///< Near-Side Yields in the EP bins
  vector<TGraphErrors *> fASYieldsInc_SysErrorBySource= {}; ///< Away-Side Yields in all EP
  vector<vector<TGraphErrors *>> fASYieldsEP_SysErrorBySource= {}; ///< Away-Side Yields in the EP bins

  vector<TGraphErrors *> fNSYieldsInc_Models = {}; ///< Near-Side Yields in all EP
  vector<vector<TGraphErrors *>> fNSYieldsEP_Models = {}; ///< Near-Side Yields in the EP bins
  vector<TGraphErrors *> fASYieldsInc_Models = {}; ///< Away-Side Yields in all EP
  vector<vector<TGraphErrors *>> fASYieldsEP_Models = {}; ///< Away-Side Yields in the EP bins

  // RPF variations used to calculate uncertainty
 // vector<TGraphErrors *> fNSYieldsInc_RPFVariants; ///< Near-Side Yields in all EP
 // vector<vector<TGraphErrors *>> fNSYieldsEP_RPFVariants; ///< Near-Side Yields in the EP bins
 // vector<TGraphErrors *> fASYieldsInc_RPFVariants; ///< Away-Side Yields in all EP
 // vector<vector<TGraphErrors *>> fASYieldsEP_RPFVariants; ///< Away-Side Yields in the EP bins


  TGraphErrors * fNSRmsInc; ///< Near-Side RMS in all EP
  vector<TGraphErrors *> fNSRmsEP; ///< Near-Side RMS in the EP bins
  TGraphErrors * fASRmsInc; ///< Away-Side RMS in all EP
  vector<TGraphErrors *> fASRmsEP; ///< Away-Side RMS in the EP bins

  // TGraphErrors with the error bars as the calculated uncertainty from RPF variations
  TGraphErrors * fNSRmsInc_RPFError; ///< Near-Side Rms in all EP
  vector<TGraphErrors *> fNSRmsEP_RPFError; ///< Near-Side Rms in the EP bins
  TGraphErrors * fASRmsInc_RPFError; ///< Away-Side Rms in all EP
  vector<TGraphErrors *> fASRmsEP_RPFError; ///< Away-Side Rms in the EP bins

  // calculated uncertainty from systematic variations
  TGraphErrors * fNSRmsInc_SysError = 0; ///< Near-Side Rms in all EP
  vector<TGraphErrors *> fNSRmsEP_SysError; ///< Near-Side Rms in the EP bins
  TGraphErrors * fASRmsInc_SysError = 0; ///< Away-Side Rms in all EP
  vector<TGraphErrors *> fASRmsEP_SysError; ///< Away-Side Rms in the EP bins


  vector<TGraphErrors *> fNSRmsInc_SysErrorBySource= {}; ///< Near-Side Rms in all EP
  vector<vector<TGraphErrors *>> fNSRmsEP_SysErrorBySource= {}; ///< Near-Side Rms in the EP bins
  vector<TGraphErrors *> fASRmsInc_SysErrorBySource= {}; ///< Away-Side Rms in all EP
  vector<vector<TGraphErrors *>> fASRmsEP_SysErrorBySource= {}; ///< Away-Side Rms in the EP bins


  vector<TGraphErrors *> fNSRmsInc_Models = {}; ///< Near-Side Rms in all EP
  vector<vector<TGraphErrors *>> fNSRmsEP_Models = {}; ///< Near-Side Rms in the EP bins
  vector<TGraphErrors *> fASRmsInc_Models = {}; ///< Away-Side Rms in all EP
  vector<vector<TGraphErrors *>> fASRmsEP_Models = {}; ///< Away-Side Rms in the EP bins

  // RPF variations used to calculate uncertainty in RMS
  //vector<TGraphErrors *> fNSRmsInc_RPFVariants; ///< Near-Side Rms in all EP
 // vector<vector<TGraphErrors *>> fNSRmsEP_RPFVariants; ///< Near-Side Rms in the EP bins
  //vector<TGraphErrors *> fASRmsInc_RPFVariants; ///< Away-Side Rms in all EP
  //vector<vector<TGraphErrors *>> fASRmsEP_RPFVariants; ///< Away-Side Rms in the EP bins

  // Sigma (width via fit parameter)


  TGraphErrors * fNSSigmasInc; ///< Near-Side Sigmas in all EP
  vector<TGraphErrors *> fNSSigmasEP; ///< Near-Side Sigmas in the EP bins
  TGraphErrors * fASSigmasInc; ///< Away-Side Sigmas in all EP
  vector<TGraphErrors *> fASSigmasEP; ///< Away-Side Sigmas in the EP bins

  // TGraphErrors with the error bars as the calculated uncertainty from RPF variations
  TGraphErrors * fNSSigmasInc_RPFError; ///< Near-Side Sigmas in all EP
  vector<TGraphErrors *> fNSSigmasEP_RPFError; ///< Near-Side Sigmas in the EP bins
  TGraphErrors * fASSigmasInc_RPFError; ///< Away-Side Sigmas in all EP
  vector<TGraphErrors *> fASSigmasEP_RPFError; ///< Away-Side Sigmas in the EP bins

  // TGraphErrors with the error bars as the calculated uncertainty from systematic variations
  TGraphErrors * fNSSigmasInc_SysError = 0; ///< Near-Side Sigmas in all EP
  vector<TGraphErrors *> fNSSigmasEP_SysError={}; ///< Near-Side Sigmas in the EP bins
  TGraphErrors * fASSigmasInc_SysError = 0; ///< Away-Side Sigmas in all EP
  vector<TGraphErrors *> fASSigmasEP_SysError={}; ///< Away-Side Sigmas in the EP bins

  vector<TGraphErrors *> fNSSigmasInc_SysErrorBySource= {}; ///< Near-Side Sigmas in all EP
  vector<vector<TGraphErrors *>> fNSSigmasEP_SysErrorBySource= {}; ///< Near-Side Sigmas in the EP bins
  vector<TGraphErrors *> fASSigmasInc_SysErrorBySource= {}; ///< Away-Side Sigmas in all EP
  vector<vector<TGraphErrors *>> fASSigmasEP_SysErrorBySource= {}; ///< Away-Side Sigmas in the EP bins

  vector<TGraphErrors *> fNSSigmasInc_Models = {}; ///< Near-Side Sigmas in all EP
  vector<vector<TGraphErrors *>> fNSSigmasEP_Models = {}; ///< Near-Side Sigmas in the EP bins
  vector<TGraphErrors *> fASSigmasInc_Models = {}; ///< Away-Side Sigmas in all EP
  vector<vector<TGraphErrors *>> fASSigmasEP_Models = {}; ///< Away-Side Sigmas in the EP bins




  // Ratios Graphs
  TGraphErrors * OutOverIn_AS;
  TGraphErrors * OutOverIn_NS;
  TGraphErrors * MidOverIn_AS;
  TGraphErrors * MidOverIn_NS;

  // Variations of RPF parameters
  vector<TGraphErrors *> OutOverIn_AS_RPFVariants;
  vector<TGraphErrors *> OutOverIn_NS_RPFVariants;
  vector<TGraphErrors *> MidOverIn_AS_RPFVariants;
  vector<TGraphErrors *> MidOverIn_NS_RPFVariants;

  vector<TGraphErrors *> RmsOutOverIn_AS_RPFVariants;
  vector<TGraphErrors *> RmsOutOverIn_NS_RPFVariants;
  vector<TGraphErrors *> RmsMidOverIn_AS_RPFVariants;
  vector<TGraphErrors *> RmsMidOverIn_NS_RPFVariants;

  vector<TGraphErrors *> SigmasOutOverIn_AS_RPFVariants;
  vector<TGraphErrors *> SigmasOutOverIn_NS_RPFVariants;
  vector<TGraphErrors *> SigmasMidOverIn_AS_RPFVariants;
  vector<TGraphErrors *> SigmasMidOverIn_NS_RPFVariants;

  // Ratios with RPF errors as the calculated uncertainty from all variations
  TGraphErrors * OutOverIn_AS_RPFError;
  TGraphErrors * OutOverIn_NS_RPFError;
  TGraphErrors * MidOverIn_AS_RPFError;
  TGraphErrors * MidOverIn_NS_RPFError;

  TGraphErrors * RmsOutOverIn_AS_RPFError;
  TGraphErrors * RmsOutOverIn_NS_RPFError;
  TGraphErrors * RmsMidOverIn_AS_RPFError;
  TGraphErrors * RmsMidOverIn_NS_RPFError;

  TGraphErrors * SigmasOutOverIn_AS_RPFError;
  TGraphErrors * SigmasOutOverIn_NS_RPFError;
  TGraphErrors * SigmasMidOverIn_AS_RPFError;
  TGraphErrors * SigmasMidOverIn_NS_RPFError;

  // Variations of Systematics
//  vector<vector<TGraphErrors *>> OutOverIn_AS_SysVariants; ///< First axis is the systematic type, 2nd is variant ID
//  vector<vector<TGraphErrors *>> OutOverIn_NS_SysVariants;
//  vector<vector<TGraphErrors *>> MidOverIn_AS_SysVariants;
//  vector<vector<TGraphErrors *>> MidOverIn_NS_SysVariants;


  // Systematic Uncertainty TGraphErrors from sysCompare output
  // The x,y are the central values in those. The error bars are the sys error
  vector<TGraphErrors *> OutOverIn_AS_SysErrorBySource; ///< Axis is the systematic type
  vector<TGraphErrors *> OutOverIn_NS_SysErrorBySource;
  vector<TGraphErrors *> MidOverIn_AS_SysErrorBySource;
  vector<TGraphErrors *> MidOverIn_NS_SysErrorBySource;

  vector<TGraphErrors *> RmsOutOverIn_AS_SysErrorBySource;
  vector<TGraphErrors *> RmsOutOverIn_NS_SysErrorBySource;
  vector<TGraphErrors *> RmsMidOverIn_AS_SysErrorBySource;
  vector<TGraphErrors *> RmsMidOverIn_NS_SysErrorBySource;

  vector<TGraphErrors *> SigmasOutOverIn_AS_SysErrorBySource;
  vector<TGraphErrors *> SigmasOutOverIn_NS_SysErrorBySource;
  vector<TGraphErrors *> SigmasMidOverIn_AS_SysErrorBySource;
  vector<TGraphErrors *> SigmasMidOverIn_NS_SysErrorBySource;

  /*
  vector<TGraphErrors *> OutOverIn_AS_SysErrorBySource; ///< Axis is the systematic type
  vector<TGraphErrors *> OutOverIn_NS_SysErrorBySource;
  vector<TGraphErrors *> MidOverIn_AS_SysErrorBySource;
  vector<TGraphErrors *> MidOverIn_NS_SysErrorBySource;
*/
  // TGraphErrors of these observables, where the errors are the systematic erros from
  // the sources added in quadrature.

  TGraphErrors * OutOverIn_AS_SysError;
  TGraphErrors * OutOverIn_NS_SysError;
  TGraphErrors * MidOverIn_AS_SysError;
  TGraphErrors * MidOverIn_NS_SysError;

  TGraphErrors * RmsOutOverIn_AS_SysError;
  TGraphErrors * RmsOutOverIn_NS_SysError;
  TGraphErrors * RmsMidOverIn_AS_SysError;
  TGraphErrors * RmsMidOverIn_NS_SysError;

  TGraphErrors * SigmasOutOverIn_AS_SysError;
  TGraphErrors * SigmasOutOverIn_NS_SysError;
  TGraphErrors * SigmasMidOverIn_AS_SysError;
  TGraphErrors * SigmasMidOverIn_NS_SysError;


  // Systematic Uncertainty TGraphErrors from model output
  vector<TGraphErrors *> OutOverIn_AS_Models; ///< Axis is the model
  vector<TGraphErrors *> OutOverIn_NS_Models;
  vector<TGraphErrors *> MidOverIn_AS_Models;
  vector<TGraphErrors *> MidOverIn_NS_Models;
  vector<TGraphErrors *> RmsOutOverIn_AS_Models;
  vector<TGraphErrors *> RmsOutOverIn_NS_Models;
  vector<TGraphErrors *> RmsMidOverIn_AS_Models;
  vector<TGraphErrors *> RmsMidOverIn_NS_Models;
  vector<TGraphErrors *> SigmasOutOverIn_AS_Models;
  vector<TGraphErrors *> SigmasOutOverIn_NS_Models;
  vector<TGraphErrors *> SigmasMidOverIn_AS_Models;
  vector<TGraphErrors *> SigmasMidOverIn_NS_Models;



  TGraphErrors * RmsOutOverIn_AS;
  TGraphErrors * RmsOutOverIn_NS;
  TGraphErrors * RmsMidOverIn_AS;
  TGraphErrors * RmsMidOverIn_NS;

  TGraphErrors * SigmasOutOverIn_AS;
  TGraphErrors * SigmasOutOverIn_NS;
  TGraphErrors * SigmasMidOverIn_AS;
  TGraphErrors * SigmasMidOverIn_NS;

  // Differences Graphs
  TGraphErrors * OutMinusIn_AS;
  TGraphErrors * OutMinusIn_NS;
  TGraphErrors * MidMinusIn_AS;
  TGraphErrors * MidMinusIn_NS;







  // Idea: first index is 0 for the central, n for any comparison, error, etc?
  // New Idea: use these arrays to store calculations with different bins
  // Yields
  vector<TGraphErrors *> fNSYieldsInc_RPFVariants = {};  ///< Near-Side Yields in the EP inclusive region
  vector<vector<TGraphErrors *>> fNSYieldsEP_RPFVariants = {}; ///< Near-Side Yields in each EP bin.
  vector<TGraphErrors *> fASYieldsInc_RPFVariants = {};  ///< Away-Side Yields in the EP inclusive region
  vector<vector<TGraphErrors *>> fASYieldsEP_RPFVariants = {}; ///< Away-Side Yields in each EP bin.

  // Widths: Truncated RMS
  vector<TGraphErrors *> fNSRmsInc_RPFVariants = {};  ///< Near-Side Rms in the EP inclusive region
  vector<vector<TGraphErrors *>> fNSRmsEP_RPFVariants = {}; ///< Near-Side Rms in each EP bin.
  vector<TGraphErrors *> fASRmsInc_RPFVariants = {};  ///< Away-Side Rms in the EP inclusive region
  vector<vector<TGraphErrors *>> fASRmsEP_RPFVariants = {}; ///< Away-Side Rms in each EP bin.
  
  vector<TGraphErrors *> fNSSigmasInc_RPFVariants = {};  ///< Near-Side Sigmas in the EP inclusive region
  vector<vector<TGraphErrors *>> fNSSigmasEP_RPFVariants = {}; ///< Near-Side Sigmas in each EP bin.
  vector<TGraphErrors *> fASSigmasInc_RPFVariants = {};  ///< Away-Side Sigmas in the EP inclusive region
  vector<vector<TGraphErrors *>> fASSigmasEP_RPFVariants = {}; ///< Away-Side Sigmas in each EP bin.

  // Comparison plots.

};


#endif

