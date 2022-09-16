#ifndef TASKSIDEBAND_H
#define TASKSIDEBAND_H

// --- ROOT system ---
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
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

#include <vector>
using namespace std;
#include <Riostream.h>
#include <TString.h>

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

class TaskSideband : public TObject {
public:
	TaskSideband();
	virtual ~TaskSideband() { ; }


	void Run();

	void Debug(Int_t input);

	void SetDebugLevel(Int_t input)         { fDebugLevel = input; }
  void SetMCMode(Int_t input)             { iMCMode = input; }

//  void SetSidebandMode(Int_t input)       { fSidebandMode = input; if (input==1) {kNSB=3; fNSBFit=3;}}
  void SetSidebandMode(Int_t input);
	void SetPlotOptions(TString input)      { fPlotOptions = input; }
	void SetOutputDir(TString input)        { fOutputDir = input; }
	void SetOutputFileName(TString input)   { fOutputFileName = input; }

  // This uses the SetSidebandMode so that the order of calling SetBackgroundSelection vs SetSideband Mode doesn't matter
	void SetBackgroundSelection(Int_t input)                      { fBackgroundSelection = input; SetSidebandMode(fSidebandMode);  }
  void SetSidebandFitMask(Int_t input);
  void SetScalingFitFunction(Int_t input)                       { fScalingFitFunction = input; }
	void SetPtBin(Int_t input)                                    { iPtBin = input; fGlobalMinPt=fPtBins[input-1]; fGlobalMaxPt=fPtBins[input]; }

	void SetPi0CorrInputFile(TFile * inputFile)                   { fPi0CorrFile = inputFile; }
	void SetSidebandInputFile(Int_t iSideband, TFile * inputFile) { fSidebandFile[iSideband] = inputFile; }
	void SetPi0PurityFile(TFile * inputFile)                      { fPi0PurityFile = inputFile; }
  void SetPurityChoice(Int_t input)                             { iPurityChoice = input; }
  void SetFixedPurity(Float_t input)                            { fFixedPurity = input; }
  void SetUseMCPurity(Int_t input)                              { fUseMCPurity = input; }
  void SetEPBin(Int_t input)                                    { fEPBin = input; }
  void SetCentralityBin(Int_t input)                            { fCent = input; }

  void SetLabel(TString input)                      { sLabel = input; }
  void SetLabel2(TString input)                      { sLabel2 = input; }

protected:
	void SetStyle();
	void PrintCanvas(TCanvas * canvas, TString name);
	void LoadPurity();
	void LoadHistograms();
	void InitArrays();
	void MassAnalysis();
	void ProduceSidebandFigure();
	void ProduceSidebandComparison();
	void SimpleNormFit();
	void ProduceBackground();
	TH1D * MergeAndScaleBkg(Int_t index, Int_t iType);
	TH1D * MergeAndScaleDEtaBkg(Int_t index, Int_t iType);
  void PlotBkgAndSignal();
  void Subtract();

  void ProcessFlow();
//  void CalculateVnPreSub();
//  void CalculateVnPostSub();
  void CalculateVnSub(int iPostSub); // iPostSub=0 -> pre subtraction, 1 -> postsubtraction
  void CalculateV3Sub(int iPostSub); // iPostSub=0 -> pre subtraction, 1 -> postsubtraction

  void RunSidebandExtrapolator(TH1D * fPredBkg, vector<TH1D *> fSBHists, vector<double> fSBMasses, double fLocalPi0Mass);

	void SaveResults();
  void DrawWIP(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size);
  void DrawAlicePerf(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size);

	Int_t fDebugLevel;                      ///<For Debugging Purposes

  Int_t fSidebandMode = 0;                ///< 0 -> Sidebands 3,4,5,6; 1 -> Sidebands 1,2,3

  TString sLabel  = "";
  TString sLabel2 = "";
  Bool_t bEnablePerformance = false;

  Bool_t bEnableArbUnits = false;            ///< Whether to use arbitrary units on axes for sharing
  Bool_t bNoYLabel = false;                  ///< Whether to remove Y/Z axis labels for performance purposes


  Int_t fCanvasWidth = 1400;
  Int_t fCanvasHeight = 900;

  Bool_t fSidebandFitMask[4] = {1,1,1,1};   ///< Which sidebands are included in the fit
  Bool_t fSidebandMask[4] = {1,1,1,1};      ///< Which sidebands are summed together for the subtraction

  Float_t fGlobalMinPt = 3;
  Float_t fGlobalMaxPt = 5; // debug values

  Int_t iMCMode = 0;                        ///< What to do with MC information
                                            // 0 -> include all (equiv to data)
                                            // 1 -> Background only
                                            // 2 -> True Pi0s
                                            // 3 -> True Etas

  Int_t fObservable;                        ///< Observable for the current analysis
  TString fTriggerName;                     ///< Name of the current trigger. e.g. "Pi0"
  TString fObservableName;                  ///< Name of the current observable (for plot labels)
  Int_t nObsBins;                           ///< How many bins we have of the current observable. Determined at run time.
  vector<Double_t> fObsBins;                ///< Bin Edges for whichever observable we are looking at

 // Double_t                    fArray_G_Bins[10];         ///< 10=kNoGammaBins+1
 // Double_t fPtBins[10];                 ///< Bin Edges for the Pt Bins

  Double_t fPtBins[6] = {5,7,9,11,14,17}; // might be able to get this information from a pt trigger histograms


  Bool_t bEnableV6 = false;               ///< Whether to include V6 in trigger flow fits

	Int_t iPtBin;                       ///< PtBin of analysis. Not used if fObs = 0 (Pt is the observable)
	Int_t fBackgroundSelection;             ///< 0: all 4 Sidebands combined, 1: SB3+SB4, 2: SB5+SB6, 3: SB4+SB5
  Int_t fScalingFitFunction;              ///< Which function to use in fit. 0: constant 1: linear 2: quadratic

  vector< Double_t> fPurityArray = {};    ///< Array storing the purity values to be stored
  vector< Double_t> fPurityArray_Err = {};///< Array storing the purity error values to be stored. Assumed to be statistical errors for now
  Int_t iPurityChoice = 1;                ///< Which purity value to use. 0 = purity=0, 1 = standard purity from graph, 2 = 1, 3 = standard - error, 4 = standard + error
      // 5 = standard - 2 * error, 6 = standard + 2 * error

  float fFixedPurity = -1;
  Int_t fUseMCPurity = 0;                 ///< Whether to use the MC purity from the phase1 file (fUseMCPurity = 1) or the MC purity from phase2 (2)
  // int in case we want to add an alternate MC purity determination (based on data from phase 2 files)

	TFile * fPi0CorrFile;                   ///< File with Pi0 candidate - hadron correlations
	TFile * fSidebandFile[4];               ///< Files with SB correlations
	TFile * fPi0PurityFile;                 ///< File containing TGraphErrors with S/(S+B) (Purity) 

  int iExtrapolatorMode = 1;              ///<
                                          ///< 0 -> old scalar method
                                          ///< 1 -> bin-by-bin extrapolation (linear)

	TString fPlotOptions;                   ///< Style for 2D plots
  TString fOutputDir;                     ///< Output directory to save the plots
  TString fOutputFileName;                ///< Name of Root file to which to save results 


  // Histogram useful for counting events
  TH1F * fHistEventHash = 0;

  // May need one of these for each input file
  Bool_t bNeedToRenormalize = false;      ///< Set if the correlations need to be renormalized by 1 / num triggers
                                          ///< This is needed in the newest mode, where phase2 saves them unnormalized (for the purpose of merging)
  vector<Bool_t> bNeedToRenormalizeSB = {};


  TH1D  * VariableInfo = 0;               ///< Histogram tracking some variable information

  TH1D  * fTriggerPt;                     ///< Distribution of trigger Pts (will be needed later)
  TH1D  * fTriggerPtWithinEPBin;          ///< Distribution of trigger within selected EP range

  // fTrigger Pts for Sidebands
  vector<TH1D *> fTriggerPtSB = {};

  // Track information
  TH1D * fTrackPtProjectionSE = 0;        ///< Histogram of track pT made from projecting Corr THnSparse for Same Events (this is a biased distribution)
  TH1D * fTrackPtProjectionME = 0;        ///< Histogram of track pT made from projecting Corr THnSparse for Mixed Events (this should be an unbiased distribution)

  // These should come from the Pi0 file, and should be the same for all of them
  TH1D * fTrackPtFromTrackPsi = 0;        ///< Histogram of track pT made from projecting the TrackPsiEPPtCent TH3 
  TH3F *hHistTrackPsiEPPtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
  TH3F *hHistTrackPsiEP3PtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
  TH3F *hHistTrackPsiEP4PtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)



  // Pi0_Cand and Sidebands vs the event plane.
  // Second order event plane
  std::vector<TH1D *> hPtEPAnglePionAcc_Proj_Pion;
  std::vector<std::vector<TH1D *>> hPtEPAnglePionAcc_Proj_SB; // First index is the pt bin, second index is the sideband index
  std::vector<TH1D *> hPtEPAngleBackgroundEstimate;
  std::vector<TH1D *> hPtEPAnglePionAcc_PionPostSub;

  // Third order
  std::vector<TH1D *> hPtEP3AnglePionAcc_Proj_Pion;
  std::vector<std::vector<TH1D *>> hPtEP3AnglePionAcc_Proj_SB; // First index is the pt bin, second index is the sideband index
  std::vector<TH1D *> hPtEP3AngleBackgroundEstimate;
  std::vector<TH1D *> hPtEP3AnglePionAcc_PionPostSub;

  // Fourth order
  std::vector<TH1D *> hPtEP4AnglePionAcc_Proj_Pion;
  std::vector<std::vector<TH1D *>> hPtEP4AnglePionAcc_Proj_SB; // First index is the pt bin, second index is the sideband index
  std::vector<TH1D *> hPtEP4AngleBackgroundEstimate;
  std::vector<TH1D *> hPtEP4AnglePionAcc_PionPostSub;



  TGraphErrors * gTriggerFlowPreSub_V2 = 0;
  TGraphErrors * gTriggerFlowPreSub_V3 = 0;
  TGraphErrors * gTriggerFlowPreSub_V4 = 0;
  TGraphErrors * gTriggerFlowPreSub_V6 = 0;

  std::vector<TGraphErrors *> gTriggerFlowSidebands_V2;
  std::vector<TGraphErrors *> gTriggerFlowSidebands_V3;
  std::vector<TGraphErrors *> gTriggerFlowSidebands_V4;
  std::vector<TGraphErrors *> gTriggerFlowSidebands_V6;


  // These may now be done in phase 4
  TGraphErrors * gTriggerFlowPostSub_V2 = 0;
  TGraphErrors * gTriggerFlowPostSub_V4 = 0;
  TGraphErrors * gTriggerFlowPostSub_V6 = 0;

  TGraphErrors * gTriggerFlowPostSub_V3 = 0;


	vector<TH1D *> fMassPtBinPi0;           ///< Histograms of the mass distribution
  vector<TH1D *> fMassPtBinAll;           ///< Histogram of all mass (vector in pt bins)
	vector<vector<TH1D *>> fMassPtBinSB;           ///< Histograms of the mass distribution (first index is pt, 2nd is SB index)

	vector<Double_t> fMeanMassPi0Val;       ///< Mean Mass values for the Pi0 Candidates
	vector<Double_t> fMeanMassPi0Val_Un;		///< Uncertainties in above
	vector<vector<Double_t>> fMeanMassSBVal;     ///< Mean Mass values for the Pi0 Candidates
	vector<vector<Double_t>> fMeanMassSBVal_Un;  ///< Uncertainties in above

	vector<TGraphErrors *> fMassVsIntegral_Full;         ///< Graph of dPhi integral vs mass.  Full DPhi projections
	vector<TF1 *> fMassEffectFit;                        ///< Fits of the Sideband mass effect 

	vector<TH1D *> fFullDPhiPi0;            ///< Full DPhi Projections for Pi0 Candidates
	vector<vector<TH1D *>> fFullDPhiSB;     ///< Full DPhi Projections for the Sidebands (first index is obs, 2nd is SB index)
	vector<TH1D *> fNearEtaDPhiPi0;            ///< Near Eta DPhi Projections for Pi0 Candidates
	vector<vector<TH1D *>> fNearEtaDPhiSB;     ///< Near Eta DPhi Projections for the Sidebands (first index is obs, 2nd is SB index)
	vector<TH1D *> fFarEtaDPhiPi0;            ///< Far Eta DPhi Projections for Pi0 Candidates
	vector<vector<TH1D *>> fFarEtaDPhiSB;     ///< Far Eta DPhi Projections for the Sidebands (first index is obs, 2nd is SB index)

  // Experiment with DEta subtractions
	vector<TH1D *> fFullDEtaPi0;            ///< Full DPhi Projections for Pi0 Candidates
	vector<vector<TH1D *>> fFullDEtaSB;     ///< Full DPhi Projections for the Sidebands (first index is obs, 2nd is SB index)

	vector<TH1D *> fNearSideDEtaPi0;            ///< Full DPhi Projections for Pi0 Candidates
	vector<TH1D *> fAwaySideDEtaPi0;            ///< Full DPhi Projections for Pi0 Candidates
	vector<vector<TH1D *>> fNearSideDEtaSB;     ///< Full DPhi Projections for the Sidebands (first index is obs, 2nd is SB index)
	vector<vector<TH1D *>> fAwaySideDEtaSB;     ///< Full DPhi Projections for the Sidebands (first index is obs, 2nd is SB index)




  vector<TH1D *> fNearSideSubDEtaPi0;  ///< Nearside DEta minus scaled Awayside (pi0 cands)
  vector<vector<TH1D *>> fNearSideSubDEtaSB;  ///< Nearside DEta minus scaled Awayside (sidebands)


	TGraphErrors * Pi0YieldTotalRatio;      ///< Graph with S/(S+B) from purity input file
	TGraphErrors * MCPi0YieldTotalRatio=0;  ///< Graph with MC True S/(S+B) from purity input file (if phase1 file has MC)
  TH1D * MCPhase2Pi0YieldTotalRatio=0;    ///< Histogram with MC True S/(S+B) from phase 2 (valid if in MC)

  TH1D * fMCTriggerDistPi0 = 0;           ///< Histogram of MC true status of triggers from Pi0-h Corr. file
  vector<TH1D *> fMCTriggerDistSBs;       ///< Array of Hists of MC True status of triggers from SB-h Corr. files

  // V_n Information from Purity Phase1 file
  TGraphErrors * gTrigger_Bv = 0;
  TGraphErrors * gTrigger_V2 = 0;
  TGraphErrors * gTrigger_V4 = 0;
  TGraphErrors * gTrigger_V6 = 0;

  TGraphErrors * gTrack_Bv = 0;
  TGraphErrors * gTrack_V2 = 0;
  TGraphErrors * gTrack_V4 = 0;
  TGraphErrors * gTrack_V6 = 0;


	vector<TH1D *> fFullPredBkgDPhi;        ///< Predicted background.  Scaled by mass effect and (1-purity)	
	vector<TH1D *> fNearEtaPredBkgDPhi;     ///< Predicted background.  Scaled by mass effect and (1-purity)	
	vector<TH1D *> fFarEtaPredBkgDPhi;      ///< Predicted background.  Scaled by mass effect and (1-purity)	
  vector<TH1D *> fFullPredBkgDPhi_Unscaled; ///< Predicted background, scaled to awayside for performance plots
  vector<TH1D *> fFullPredBkgDPhi_ArbScaled; ///< Predicted background, scaled arbtrarily to show shape

  vector<TH1D *> fFullDPhiFinal;          ///< Final subtracted correlations  
  vector<TH1D *> fNearEtaDPhiFinal;          ///< Final subtracted correlations  
  vector<TH1D *> fFarEtaDPhiFinal;          ///< Final subtracted correlations  

  // Delta Eta Correlations
	vector<TH1D *> fNearSideSubPredBkgDEta;      ///< Predicted background.  Scaled by mass effect and (1-purity)	
	vector<TH1D *> fAwaySideSubPredBkgDEta;      ///< Predicted background.  Scaled by mass effect and (1-purity)	

  vector<TH1D *> fNearSideSubDEtaFinal;          ///< Final subtracted correlations  
  vector<TH1D *> fAwaySideSubDEtaFinal;          ///< Final subtracted correlations  




  // Constants
	//static const Int_t kNSB=4;              ///< Total Number of Sidebands (a constant for now)
	Int_t kNSB=4;                           ///< Total Number of Sidebands (a constant no longer)
	Int_t fNSBFit=4;                        ///< Number of Sidebands to be used for mass scaling fit
	Int_t fNSB=4;                           ///< Number of Sidebands to be used for final sum

  Int_t fEPBin = -1;                      ///< Which event plane bin are we looking at
  Int_t fCent = 2;                        ///< Which centrality bin were looking at

  static const Int_t kGammaNBINS=9;  //9    ///< Number of 2D histograms for Gamma energy
  static const Int_t kZtNBINS=7;            ///< Number of 2D histograms for Zt of g-h pair
  static const Int_t kXiNBINS=8;            ///< Number of 2D histograms for Xi of g-h pair

  static const Int_t kNPtBins=5;



  static const Int_t kRawPi0Color=kBlack;
  static const Int_t kScaledSBColor=kRed+1;

private:
	TaskSideband(const TaskSideband&); // not imp.
	ClassDef(TaskSideband,1);
};


#endif
