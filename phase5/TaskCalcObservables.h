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
  void AddSecondaryFile(TFile * inputFile, TString fLabel = "") {
  fInputFilesSecondary.push_back(inputFile);
  fLabelsSecondary.push_back(fLabel);
}

 
  void SetPlotOptions(TString input)      { fPlotOptions = input; }
  void SetOutputDir(TString input)        { fOutputDir   = input; }
  void SetOutputFile(TFile * outputFile)  { fOutputFile  = outputFile; }

  void SetMCGenMode(Int_t input = 1)      { fIsMCGenMode = input; }
  Int_t GetMCGenMode()                    { return fIsMCGenMode;  }
  void SetCentralityBin(Int_t input)      { iCentBin     = input; }
  void SetRPFMethod(Int_t input)          { iRPFMethod   = input; }

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

  void CalculateResults();                   ///< Loop over obs bins
  void CalculateResultsObsBinEPBin(int iObsBin, int iEPBin, TCanvas * canv);  ///< Calculate yields, widths for this bin in iObsBin, EP bin


  void CreateRatioAndDifferencesGraphs();

  // Calculate Ratios
  void CalculateResultsRatiosObsBin(int iObsBin, TCanvas * canv);

  // Calculate Differences
  void CalculateResultsDifferencesObsBin(int iObsBin, TCanvas * canv);


  void SetGraphStyles(); ///< Set the colors, markers on the graphs
  void DrawResults(); ///< Draw the plots

//  CalculateYield(TH1D * hist, bool bAwaySide, double &Error);


  int fYieldCanvasWidth = 900;
  int fYieldCanvasHeight = 600;  


private:

  Int_t fIsMCGenMode;                       ///< 0 = data, 1 = mcGen mode

  Int_t fDebugLevel;                       ///< For Debugging Purposes

  TFile *fInputFileCentral;               ///< File for central (systematic) file
  vector<TFile *> fInputFilesSecondary;   ///< Files for systematically changing things
  vector<TString> fLabelsSecondary;       ///< Corresponding labels
  // Could make this have more dimensions, for being more differential

  TString fOutputDir;                       ///< Output directory to save the plots

  
  Int_t iPtBin = 4;                         ///< Which Pt Bin is used (for fixed V_n^t values
  // update use iPtBin also for restricting pt trigger distributions

  Double_t PtBins[6] = {5,7,9,11,14,17}; // might be able to get this information from a pt trigger histograms
  Int_t iObservable = 1;                    ///< Which observable (0=Pt,1=Zt) is being iterated over



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

  Int_t kEPColorList[4] = {kBlack, kBlue-4, kGreen-3, kRed+1};
  Int_t kEPMarkerList[4] = {kOpenSquare, kFullSquare, 39, kFullDiamond};


//  Double_t fRMSRange = 1.047; // Range in DeltaPhi to calulate truncated RMS

  // May want to vary these for systematic uncertainty

  Double_t fRmsRangeNS = TMath::Pi()/3.;
  Double_t fRmsRangeAS = TMath::Pi()/3.;

  Double_t fYieldRangeNS = TMath::Pi()/3.;
  Double_t fYieldRangeAS = TMath::Pi()/3.;

  // The histograms to integrate

  // First index is observable bin, 2nd is the event plane bin
  vector<vector<TH1D *>> fFullDPhiProj_Sub;  ///< Full projections in DPhi after subtracting flow and rescaling for number of triggers in the EP bin
  vector<vector<TH1D *>> fNearEtaDPhiProj_Sub;  ///< Full projections in DPhi after subtracting flow and rescaling for number of triggers in the EP bin
  vector<vector<TH1D *>> fFarEtaDPhiProj_Sub;  ///< Full projections in DPhi after subtracting flow and rescaling for number of triggers in the EP bin

  // Items for applying errors from reaction plane fit

  
  // Primary observable graphs

  TGraphErrors * fNSYieldsInc; ///< Near-Side Yields in all EP
  vector<TGraphErrors *> fNSYieldsEP; ///< Near-Side Yields in the EP bins
  TGraphErrors * fASYieldsInc; ///< Away-Side Yields in all EP
  vector<TGraphErrors *> fASYieldsEP; ///< Away-Side Yields in the EP bins

  TGraphErrors * fNSRmsInc; ///< Near-Side RMS in all EP
  vector<TGraphErrors *> fNSRmsEP; ///< Near-Side RMS in the EP bins
  TGraphErrors * fASRmsInc; ///< Away-Side RMS in all EP
  vector<TGraphErrors *> fASRmsEP; ///< Away-Side RMS in the EP bins



  // Ratios Graphs
  TGraphErrors * OutOverIn_AS;
  TGraphErrors * OutOverIn_NS;
  TGraphErrors * MidOverIn_AS;
  TGraphErrors * MidOverIn_NS;

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
  // Yields
  vector<TGraphErrors *> fNSYieldsInc_Array;  ///< Near-Side Yields in the EP inclusive region
  vector<vector<TGraphErrors *>> fNSYieldsEP_Array; ///< Near-Side Yields in each EP bin.
  vector<TGraphErrors *> fASYieldsInc_Array;  ///< Away-Side Yields in the EP inclusive region
  vector<vector<TGraphErrors *>> fASYieldsEP_Array; ///< Away-Side Yields in each EP bin.

  // Widths: Truncated RMS

  

  // Comparison plots.


};


#endif

