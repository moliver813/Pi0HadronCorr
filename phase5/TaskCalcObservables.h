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
  void InitArrays();

  void SetPtBin(Int_t input)              { iPtBin       = input; }
  void SetObservable(Int_t input)         { iObservable  = input; }

  void SetCentralInputFile(TFile * inputFile) { fInputFileCentral = inputFile; }


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

  void SetMCGenMode(Int_t input = 1)      { fIsMCGenMode = input; }
  Int_t GetMCGenMode()                    { return fIsMCGenMode;  }
  void SetCentralityBin(Int_t input)      { iCentBin     = input; }
  void SetRPFMethod(Int_t input)          { iRPFMethod   = input; }

  void SetNSkipPoints(Int_t input)        { nSkipPoints  = input; }
  int  GetNSkipPoints()                   { return nSkipPoints; }

  void Debug(Int_t input);
  void SetDebugLevel(Int_t input)         { fDebugLevel = input; }

  void SaveOutput();
  void Run();

protected:
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

  void SetGraphStyles(); ///< Set the colors, markers on the graphs

  void DrawObservable(vector<TGraphErrors *> fObsGraphs, vector<TGraphErrors *> fObsRPFErrors, vector<TGraphErrors *> fObsSysErrors = {});


  void DrawResults(); ///< Draw the plots

//  CalculateYield(TH1D * hist, bool bAwaySide, double &Error);


  int fYieldCanvasWidth = 900;
  int fYieldCanvasHeight = 600;  


private:

  Int_t fIsMCGenMode;                       ///< 0 = data, 1 = mcGen mode

  Int_t fDebugLevel;                       ///< For Debugging Purposes

  TFile *fInputFileCentral;               ///< File for central (systematic) file
  vector<vector<TFile *>> fInputFilesSecondary;   ///< Files for systematically changing things
    // 1st index is the axis (v3 change, purity change, etc)
    // 2nd index is just within that list

  vector<TString> fAxesSecondary;     ///< Labels for each "axis" of the variations
  vector<vector<TString>> fLabelsSecondary;       ///< Corresponding labels
  // Could make this have more dimensions, for being more differential

  TString fOutputDir;                       ///< Output directory to save the plots

  
  Int_t iPtBin = 4;                         ///< Which Pt Bin is used (for fixed V_n^t values
  // update use iPtBin also for restricting pt trigger distributions

  Double_t PtBins[6] = {5,7,9,11,14,17}; // might be able to get this information from a pt trigger histograms
  Int_t iObservable = 1;                    ///< Which observable (0=Pt,1=Zt) is being iterated over

  Int_t nRebinDPhi = 1;                     ///< Rebinning in Delta Phi

  Int_t iCentBin = 0;                           ///< Which centrality bin (0-10,10-30,30-50,50-80)

  Int_t iRPFMethod = 0;                         ///< Which RPF method to use. 0 is my c++, 1 is Raymond's python implementation

  TString fPlotOptions="COLZ";              ///< Style for plotting 3D histograms.  Default is colz

//  TFile *fInputCentral;                     ///< File with central analysis

  TFile *fOutputFile;                       ///< File where output objects are to be stored

  Int_t fObservable;                        ///< Observable for the current analysis
  TString fObservableName;                  ///< Name of the current observable (for plot labels)

  Int_t nObsBins;                           ///< How many bins we have of the current observable. Determined at run time.

  vector<Double_t> fObsBins;                ///< Bin Edges for whichever observable we are looking at

  TString fEPBinTitles[kNEPBins+1] = {"In-Plane","Mid-Plane","Out-of-Plane","All EP Angles"};
 // TString fEPBinTitles[kNEPBins+1] = {"In-Plane","Mid-Plane","Out-of-Plane","All EP Angles"};
  TString fEPBinTitlesShort[kNEPBins+1] = {"EP0","EP1","EP2","Incl"};

  //Int_t kEPColorList[4] = {kBlack, kBlue-4, kGreen-3, kRed+1};
  //Int_t kEPMarkerList[4] = {kOpenSquare, kFullSquare, 39, kFullDiamond};
  Int_t kEPColorList[4] = {kBlue-4, kGreen-3, kRed+1,kBlack};
  Int_t kEPMarkerList[4] = {kFullSquare, 39, kFullDiamond, kOpenSquare};



  // Red / Blue
  Int_t kOutInColor = kViolet-1;
  // Green / Blue
  Int_t kMidInColor = kEPColorList[1];

  Int_t kOutInErrorColor = kOutInColor;
  Int_t kMidInErrorColor = kMidInColor;

  Int_t kOutOverInMarker = kFullSquare;
  Int_t kMidOverInMarker = kFullCircle;

  // Red / Blue, Green / Blue
//  int kOutInErrorColor = kOutInColor;
//  int kMidInErrorColor = kMidInColor;
  
  int kOutInFillStyle = 3245;
  int kMidInFillStyle = 3254;
  int kOutInNSFillStyle = 3245;
  int kMidInNSFillStyle = 3254;

  // RPF errors
  Int_t kEPRPFFillColorList[4] = {kBlue-4, kGreen-3, kRed+1,kBlack};
  Int_t kEPRPFFillStyleList[4] = {3245,3254,3245,3254};

  Int_t nSkipPoints = 0;  // How many of the first points to skip

//  Double_t fRMSRange = 1.047; // Range in DeltaPhi to calulate truncated RMS

  // May want to vary these for systematic uncertainty

  Double_t fRmsRangeNS = TMath::Pi()/3.;
  Double_t fRmsRangeAS = TMath::Pi()/3.;

  Double_t fYieldRangeNS = TMath::Pi()/3.;
  Double_t fYieldRangeAS = TMath::Pi()/3.;

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

  // RPF variations used to calculate uncertainty in RMS
  //vector<TGraphErrors *> fNSRmsInc_RPFVariants; ///< Near-Side Rms in all EP
 // vector<vector<TGraphErrors *>> fNSRmsEP_RPFVariants; ///< Near-Side Rms in the EP bins
  //vector<TGraphErrors *> fASRmsInc_RPFVariants; ///< Away-Side Rms in all EP
  //vector<vector<TGraphErrors *>> fASRmsEP_RPFVariants; ///< Away-Side Rms in the EP bins




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

  // Ratios with RPF errors as the calculated uncertainty from all variations
  TGraphErrors * OutOverIn_AS_RPFError;
  TGraphErrors * OutOverIn_NS_RPFError;
  TGraphErrors * MidOverIn_AS_RPFError;
  TGraphErrors * MidOverIn_NS_RPFError;
  TGraphErrors * RmsOutOverIn_AS_RPFError;
  TGraphErrors * RmsOutOverIn_NS_RPFError;
  TGraphErrors * RmsMidOverIn_AS_RPFError;
  TGraphErrors * RmsMidOverIn_NS_RPFError;

  // Variations of Systematics
  vector<vector<TGraphErrors *>> OutOverIn_AS_SysVariants; ///< First axis is the systematic type, 2nd is variant ID
  vector<vector<TGraphErrors *>> OutOverIn_NS_SysVariants;
  vector<vector<TGraphErrors *>> MidOverIn_AS_SysVariants;
  vector<vector<TGraphErrors *>> MidOverIn_NS_SysVariants;

  // Ratios with Systematic errors as the calculated uncertainty from all variations
  vector<TGraphErrors *> OutOverIn_AS_SysError; ///< Axis is the systematic type
  vector<TGraphErrors *> OutOverIn_NS_SysError;
  vector<TGraphErrors *> MidOverIn_AS_SysError;
  vector<TGraphErrors *> MidOverIn_NS_SysError;


  TGraphErrors * RMSOutOverIn_AS;
  TGraphErrors * RMSOutOverIn_NS;
  TGraphErrors * RMSMidOverIn_AS;
  TGraphErrors * RMSMidOverIn_NS;

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
  // Yields
  vector<TGraphErrors *> fNSRmsInc_RPFVariants = {};  ///< Near-Side Rms in the EP inclusive region
  vector<vector<TGraphErrors *>> fNSRmsEP_RPFVariants = {}; ///< Near-Side Rms in each EP bin.
  vector<TGraphErrors *> fASRmsInc_RPFVariants = {};  ///< Away-Side Rms in the EP inclusive region
  vector<vector<TGraphErrors *>> fASRmsEP_RPFVariants = {}; ///< Away-Side Rms in each EP bin.
  

  // Comparison plots.

};


#endif

