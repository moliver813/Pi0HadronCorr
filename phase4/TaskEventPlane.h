#ifndef TASKEVENTPLANE_H
#define TASKEVENTPLANE_H

#include <Riostream.h>
#include <TString.h>
#include "TaskEventPlaneMathTools.cxx"
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

class RPF_Functor;
class RPF_Functor_Single;



/*CColor=9
CMarkerStyle=22

PyColor=8
PyMarkerStyle=26

PyRPDepColor=46
PyRPDepMarkerStyle=25
*/


class TaskEventPlane : public TObject {


public:
	TaskEventPlane();

	virtual ~TaskEventPlane()  { ; }
 
	void SetStyle();
	void LoadHistograms();
  void FitFlow();
	void InitArrays();
	void DrawRawOmniPlots();
	void DrawFitOmniPlots();
	void DrawRescaleOmniPlots();
  void PrelimCalculation();
  void PrelimCalculation_Step(Int_t iV);
	void DrawOmniPlots_Type(vector<TH1D *> fHists, TString fLabel, vector<TF1 *> fFits = {});
//	void DrawOmniPlots_Type(vector<TH1D *> fHists, TString fLabel, Bool_t bIncludeFit = false);
	void DoRPFThing();
	void DoRPFThing_Step(vector<TH1D *> fHists, TString fLabel,Int_t iObsBin, Double_t fV2T_Fixed);
  // RPF functions (written in TaskEventPlaneMathTools.cxx

//  Double_t RPF_BR(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0);
//  Double_t RPF_VR_2(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0);
//  Double_t RPF_VR_4(Double_t fPhi_S, Double_t fC_S, Double_t vt_2, Double_t vt_4 = 0);
//  Double_t RPFFunction_Single(Double_t * x, Double_t * par);
//  Double_t RPFFunction(Double_t * x, Double_t * par );

//  TF1 * FitRPF(TH1D * fHist, RPF_Functor *fFit, TString fName, Double_t fV2T_Fixed);
//  void RPF_Prefit(TF1 * fit,TH1D * fHist, RPF_Functor * funct);
//  TH1D * MergeEvtPlanesForRPF(vector<TH1D *> fHists, TString fName);


  void FormatPythonRPFs();
  void CompareParameters();
	void ProcessFitParams();
	void PlotFitParams();
	void SubtractBackground();
  void DrawOmniSandwichPlots();
	void DrawOmniSandwichPlots_Step(Int_t iV, Int_t iObsBin);
  void Rescale();
  void RescaleRegion(Int_t iV, Int_t iObsBin, Int_t iRegion); // 0 for full, 1 for near, 2 for far
  void SaveOutput();
	void Run_Part1();
  void Run_Part2();


	void Debug(Int_t input);
	void SetDebugLevel(Int_t input)         { fDebugLevel = input; }

  void SetMCGenMode(Int_t input = 1)      { fIsMCGenMode = input; }
  Int_t GetMCGenMode()                    { return fIsMCGenMode;  }
  void SetPtBin(Int_t input)              { iPtBin       = input; }
  void SetObservable(Int_t input)         { iObservable  = input; }
	void SetRPFMode(Int_t input)            { iRPFMode     = input; }
	void SetPlotOptions(TString input)      { fPlotOptions = input; }
	void SetOutputDir(TString input)        { fOutputDir   = input; }
  void SetOutputFile(TFile * outputFile)  { fOutputFile  = outputFile; }
  TFile * GetOutputFile()  { return fOutputFile; }
	void SetSavePlots(Bool_t input)         { fSavePlots   = input; }
  void SetEPRSet(Int_t input)             { iEPRSet      = input; }
  void SetCentralityBin(Int_t input)      { iCentBin     = input; }
  void SetFixV2TMode(Bool_t input)        { bFixV2T      = input; }
  void SetFixV3To0(Bool_t input)          { bFixV3To0    = input; }

	void SetIntermediateInputFile(Int_t evtPlaneIndex, TFile * inputFile) {
		if (evtPlaneIndex == -1) { fInputFileAll = inputFile; return; }
		if (evtPlaneIndex >= 0 && evtPlaneIndex <= 2)  { fInputFileEvt[evtPlaneIndex] = inputFile; }
	}

  int GetObservable()                     { return fObservable; }
  TString GetObservableName()             { return fObservableName; }
  int GetNObsBins()                       { return nObsBins;    }
  int GetCentBin()                        { return iCentBin;    }
  int GetEPRSet()                         { return iEPRSet;     }

  double GetMCRescaleFactor()             { return fMCRescaleFactor; }

  int GetFlowTermModeTrigger()                   { return iFlowTermModeTrigger; }
  void SetFlowTermModeTrigger(int input)       { iFlowTermModeTrigger = input; }
  int GetFlowTermModeAssoc()                   { return iFlowTermModeAssoc; }
  void SetFlowTermModeAssoc(int input)       { iFlowTermModeAssoc = input; }

  int GetFlowV5Mode()                     { return iFlowV5Mode; }
  void SetFlowV5Mode(int input)           { iFlowV5Mode = input; }
  int GetFlowV6TMode()                     { return iFlowV6TMode; }
  void SetFlowV6TMode(int input)           { iFlowV6TMode = input; }
  int GetFlowV6AMode()                     { return iFlowV6AMode; }
  void SetFlowV6AMode(int input)           { iFlowV6AMode = input; }

  double GetNumTriggers()                    { return fAllTriggerPt->Integral(); }
  double GetNumTriggersEP(int iEP)        { return fEPBinTriggerPt[iEP]->Integral(); }
  TH1D * GetFullDPhiProjAll(int i)         { return fFullDPhiProjAll[i]; }
  TH1D * GetFullDPhiProjEP(int i, int iEP)  { return fFullDPhiProj[i][iEP]; }
  TH1D * GetFarEtaDPhiProjAll(int i)         { return fFarEtaDPhiProjAll[i]; }
  TH1D * GetFarEtaDPhiProjEP(int i, int iEP)  { return fFarEtaDPhiProj[i][iEP]; }
  TH1D * GetNearEtaDPhiProjAll(int i)         { return fNearEtaDPhiProjAll[i]; }
  TH1D * GetNearEtaDPhiProjEP(int i, int iEP)  { return fNearEtaDPhiProj[i][iEP]; }

  TGraphErrors * GetParamGraph(int i)     { return fParGraphs[i]; }
  TGraphErrors * GetChiSquareGraph()      { return fChiSqGraph; }

  TGraphErrors * GetTriggerBv()           { return gTrigger_Bv; }
  TGraphErrors * GetTriggerV2()           { return gTrigger_V2; }
  TGraphErrors * GetTriggerV4()           { return gTrigger_V4; }
  TGraphErrors * GetTriggerV6()           { return gTrigger_V6; }
  TGraphErrors * GetTrackBv()             { return gTrack_Bv; }
  TGraphErrors * GetTrackV2()             { return gTrack_V2; }
  TGraphErrors * GetTrackV4()             { return gTrack_V4; }
  TGraphErrors * GetTrackV6()             { return gTrack_V6; }

  void InputPyBkgChiSqGraph(TGraphErrors * input)        { fPyBkgChiSqGraph=input;   }
  void InputPyBkgParamGraph(int i, TGraphErrors * input) { fPyBkgParGraphs[i]=input; }
  void InputPyRPSChiSqGraph(TGraphErrors * input)        { fPyRPSChiSqGraph=input;   }
  void InputPyRPSParamGraph(int i, TGraphErrors * input) { fPyRPSParGraphs[i]=input; }


  static const Int_t kCMethodColor = 9;
  static const Int_t kCMethodMarkerStyle = 22;

  static const Int_t kPyBkgColor = 8;
  static const Int_t kPyBkgMarkerStyle = 26;

  static const Int_t kPyRPDepColor = 46;
  static const Int_t kPyRPDepMarkerStyle = 25;

  static const Int_t kFlowFitColor = kOrange+10;
  static const Int_t kFlowFitMarkerStyle = kOpenCircle;
  static const Int_t kFlowFillStyle = 3004;

  static const Bool_t kEnableGridX = 1;
  static const Bool_t kEnableGridY = 1;


protected:
	// Constants
	static const Int_t kNEPBins=3;

  static const Int_t kGammaNBINS=9;  //9    ///< Number of 2D histograms for Gamma energy
  static const Int_t kZtNBINS=7;            ///< Number of 2D histograms for Zt of g-h pair
  static const Int_t kXiNBINS=8;            ///< Number of 2D histograms for Xi of g-h pair
  static const Int_t kNoHPtBins=8;          ///< Bins in hadron pT

// following 3 probable not important at this stage 
  static const Int_t kNvertBins=10;          ///< Z-Vertex bins in which the ME are mixed
  static const Int_t kNvertBins_alt=8;      ///< Instead of the usual 1cm binning we will do a 4cm binning at the edge and a 2cm binning in the center
  static const Int_t kCentBINS= 4;    

  const Int_t kFitLineColor = kViolet-5;
  const Float_t kOmniMarkerSize = 0.5;

  // Variables that change in Evt Plane Bins
//  const Double_t kPhi_S[3] = {0., 0.25, 0.5}; // To be multiplied by pi
    // what about 3/4?
//  const Double_t kC_S[3]   = {1./6.,1./12.,1./6.};

  // Preset : Cent 10-30%
  //Double_t fEPRes[6] = {0.0,0.885,0.605,0.245,0.0,0.1};
  Double_t fEPRes[6] = {0,-1,-1,-1,-1,0};


  void SetEPTitles(vector<vector<TH1D *>> HistArray);

  Int_t fIsMCGenMode;                       ///< 0 = data, 1 = mcGen mode

	Int_t fDebugLevel;                       ///< For Debugging Purposes

	TFile *fInputFileAll;                     ///< File with all event plane windows
	TFile *fInputFileEvt[3];                  ///< Files In, Mid, And Out plane windows

	TFile *fOutputFile;                       ///< File where output objects are to be stored

  TString fDPhiHistName = "SBSub";         ///< Name at beginning of full, nearEta, farEta files
                                           // Should be SBSub for post Sideband Subtraction hists
	// Options

	vector<string> fPlaneLabels;              ///< General Labels for each EP bin.  Last is for All Angles
  TString fPlotOptions;                     ///< Style for plotting 2D histograms.  Default is lego2
  TString fOutputDir;                       ///< Output directory to save the plots
  Bool_t fSavePlots;                        ///< Whether to save plots

  Int_t iPtBin = 4;                         ///< Which Pt Bin is used (for fixed V_n^t values
  // update use iPtBin also for restricting pt trigger distributions

  Double_t PtBins[6] = {5,7,9,11,14,17}; // might be able to get this information from a pt trigger histograms
  Int_t iObservable = 1;                    ///< Which observable (0=Pt,1=Zt) is being iterated over

  Double_t fMCRescaleFactor = 1e7;          ///< Rescale for MC to give realistic scale

	Int_t iRPFMode;                           ///< Which mode to use RPF method in #FIXME write them down somewhere

  Int_t nRPFMethods = 1;                   ///< How many versions of the RPF to use? 1 is just the C++ implementation, 3 includes both the Bkg and RPSignal fit versions

  Int_t iEPRSet = 0;                            ///< Which set of event plane resolutions to use. Default to MB data measurement
  Int_t iCentBin = 0;                           ///< Which centrality bin for EPR? (0-10,10-30,30-50,50-80)
	
  Bool_t bFixV2T = 0;                       ///< Whether to fix the V2T to the value found in the first z_t bin.
  Bool_t bFixV3To0 = 0;                     ///< Whether to fix the V3AV3T to 0

  // Flow analysis info
  const int kNTrackPtBins = 14;
  std::vector <double> fTrackPtBins = {0.15, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 17};



  Int_t iFlowTermModeTrigger = 0;           ///< How to apply the phase 1 flow information for trigger
                                            // 0: none
                                            // 1: fix V2T, V4T
                                            // 2: Set V2T, V4T ranges based on 1 sigma 

  Int_t iFlowTermModeAssoc = 0;             ///< How to apply the phase 1 flow information for associated
                                            // 0: none
                                            // 1: Fix V2A, V4A
                                            // 2: Fix V2A, V4A ranges based on 1 sigma

  Int_t iFlowV5Mode = 0;                    ///< How to deal with with v5tv5a term
                                            // 0 : none
                                            // 1 : enable

  Int_t iFlowV6TMode = 0;                   ///< How to deal with with v6t term
                                            // 0 : none
                                            // 1 : enable, free
                                            // 2 : enable, fix to measured EP value
                                            // 3 : enable, fix near measured EP value

  Int_t iFlowV6AMode = 0;                   ///< How to deal with with v6a term
                                            // Same definitions as above



	Int_t fObservable;                        ///< Observable for the current analysis
	TString fObservableName;                  ///< Name of the current observable (for plot labels)

	Int_t nObsBins;                           ///< How many bins we have of the current observable. Determined at run time.

	vector<Double_t> fObsBins;                ///< Bin Edges for whichever observable we are looking at

  TString fEPBinTitles[kNEPBins+1] = {"In-Plane","Mid-Plane","Out-of-Plane","All EP Angles"};

  TString fParNames[9] = {"B","V2T","V2A","V3","V4T","V4A","V5","V6T","V6A"};
  TString fParTitles[9] = {"B","v_{2}^{t}","v_{2}^{a}","v_{3}^{t}v_{3}^{a}","v_{4}^{t}","v_{4}^{a}","v_{5}^{t}v_{5}^{a}","v_{6}^{t}","v_{6}^{a}"};

  Int_t kEPColorList[4] = {kBlack, kBlue-4, kGreen-3, kRed+1};

  // Trigger Information
  TH1D *                 fAllTriggerPt;             ///< Distribution of trigger Pt in all EP bins
  vector<TH1D *>         fEPBinTriggerPt;           ///< Distribution of trigger Pt in each EP bin

  // Track information
  TH1D * fTrackPtProjectionSE = 0;        ///< Histogram of track pT made from projecting Corr THnSparse for Same Events (this is a biased distribution)
  TH1D * fTrackPtProjectionME = 0;        ///< Histogram of track pT made from projecting Corr THnSparse for Mixed Events (this should be an unbiased distribution)

  Int_t kUsedPi0TriggerPtBins = 5; // How many Pt bins do we actually use

  TH2F * hHistTrackPsiEPPt=0;   // Tracks
  TH2F * hHistTrackPsiRPPt=0;
  std::vector<TH1F *> hPtEPAngleTrack_Proj;

  TH1D * fTrackPtFromTrackPsi = 0;        ///< Histogram of track pT made from projecting the TrackPsiEPPtCent TH3 


  // Vn Information from fits in Phase 1 or to be recalculated
  // To be recalculated here from the raw pi0 candidate vs event plane
  TGraphErrors * gTrigger_Bv = 0;
  TGraphErrors * gTrigger_V2 = 0;
  TGraphErrors * gTrigger_V4 = 0;
  TGraphErrors * gTrigger_V6 = 0;

  // The trigger vs event plane prior to sideband subtraction
  TGraphErrors * gTrigger_Bv_Presub = 0;
  TGraphErrors * gTrigger_V2_Presub = 0;
  TGraphErrors * gTrigger_V4_Presub = 0;
  TGraphErrors * gTrigger_V6_Presub = 0;

  TGraphErrors * gTrack_Bv = 0;
  TGraphErrors * gTrack_V2 = 0;
  TGraphErrors * gTrack_V4 = 0;
  TGraphErrors * gTrack_V6 = 0;

	// Histograms
	vector<TH1D *>         fFullDPhiProjAll;          ///< Full Projections in DPhi.
	vector<vector<TH1D *>> fFullDPhiProj;             ///< Full Projections in DPhi.  First index is observable bin, second is event plane bin
  TH3F *hHistTrackPsiEPPtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
 // TH3F *hHistTrackPsiEPPtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
	vector<TH1D *>         fNearEtaDPhiProjAll;       ///< Near Eta (Signal) Projections in DPhi.
	vector<vector<TH1D *>> fNearEtaDPhiProj;          ///< Near Eta (Signal) Projections in DPhi.  First index is observable bin, second is event plane bin
	vector<TH1D *>         fFarEtaDPhiProjAll;       ///< Far Eta (Background) Projections in DPhi.
	vector<vector<TH1D *>> fFarEtaDPhiProj;          ///< Far Eta (Background) Projections in DPhi.  First index is observable bin, second is event plane bin


  // RPF Fit Residuals

  // FIXME add additional dimension of RPF Method

	// Histograms after background subtraction
	/*vector<TH1D *>         fFullDPhiProjAll_Sub;          ///< Full Projections in DPhi.
	vector<vector<TH1D *>> fFullDPhiProj_Sub;             ///< Full Projections in DPhi.  First index is observable bin, second is event plane bin
	vector<TH1D *>         fNearEtaDPhiProjAll_Sub;       ///< Near Eta (Signal) Projections in DPhi.
	vector<vector<TH1D *>> fNearEtaDPhiProj_Sub;          ///< Near Eta (Signal) Projections in DPhi.  First index is observable bin, second is event plane bin
	vector<TH1D *>         fFarEtaDPhiProjAll_Sub;       ///< Far Eta (Background) Projections in DPhi.
	vector<vector<TH1D *>> fFarEtaDPhiProj_Sub;          ///< Far Eta (Background) Projections in DPhi.  First index is observable bin, second is event plane bin*/
	vector<vector<TH1D *>>         fFullDPhiProjAll_Sub;          ///< Full Projections in DPhi.
	vector<vector<vector<TH1D *>>> fFullDPhiProj_Sub;             ///< Full Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<TH1D *>>         fNearEtaDPhiProjAll_Sub;       ///< Near Eta (Signal) Projections in DPhi.
	vector<vector<vector<TH1D *>>> fNearEtaDPhiProj_Sub;          ///< Near Eta (Signal) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<TH1D *>>         fFarEtaDPhiProjAll_Sub;        ///< Far Eta (Background) Projections in DPhi.
	vector<vector<vector<TH1D *>>> fFarEtaDPhiProj_Sub;           ///< Far Eta (Background) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin

  // After background sub. and trigger rescaling
	vector<vector<vector<TH1D *>>> fFullDPhiProj_Rescale;             ///< Full Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<vector<TH1D *>>> fNearEtaDPhiProj_Rescale;          ///< Near Eta (Signal) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<vector<TH1D *>>> fFarEtaDPhiProj_Rescale;          ///< Far Eta (Background) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin


  // Fits and Maths
//	vector<TF1 *> fRPFFits;                   ///< Array of the fits
//	vector<vector<TF1 *>> fRPFFits_Indiv;     ///< Array of resultant fits (1st index is Obs, 2nd index is EP bin)

//  vector<TH1D *> fRPF_Residuals = {};                ///< Array of (Data - Fit) / Fit [iObsBin] // the original fit
//  vector<vector<TH1D *>> fRPF_Residuals_Indiv = {};  ///< Array of (Data - Fit) / Fit [iObsBin,iEPBin]

  vector<vector<TF1 *>> fRPFFits;                   ///< Array of the fits
	vector<vector<vector<TF1 *>>> fRPFFits_Indiv;     ///< Array of resultant fits (1st index is RPF Method, 2nd index is Obs, 3rd index is EP bin)

  // this one still needs to be updated for the outer loop, though I'm not sure if it is acutally used
  vector<vector<TH1D *>> fRPF_Residuals = {};                ///< Array of (Data - Fit) / Fit [iRPFMethod][iObsBin] // the original fit

  vector<vector<vector<TH1D *>>> fRPF_Residuals_Indiv = {};  ///< Array of (Data - Fit) / Fit [iRPFMethod][iObsBin,iEPBin]


//	vector<TString> fParNames;                ///< Names of the parameters
  // Not used yet
	vector<vector<Double_t>> fParameters;     ///< Parameters from the fit.  First dim is Obs
	vector<vector<Double_t>> fParErrors;      ///< Errors from ^^^

  // C++ Implementation Results
  TGraphErrors * fChiSqGraph=0;             // TGraphErrors of Chi Sq over NDF
  vector<TGraphErrors *> fParGraphs;        // Parameter TGraphs
  
  // Python Bkg Fit Results
  TGraphErrors * fPyBkgChiSqGraph=0;             // TGraphErrors of Chi Sq over NDF
  vector<TGraphErrors *> fPyBkgParGraphs;        // Parameter TGraphs

  // Python Reaction Plane Signal Fit Results
  TGraphErrors * fPyRPSChiSqGraph=0;             // TGraphErrors of Chi Sq over NDF
  vector<TGraphErrors *> fPyRPSParGraphs;        // Parameter TGraphs


 // vector<TGraphErrors *> fArrayChiSqGraph; // Array of Chi Squared / nDOF TGraph Errors

  vector<TGraphErrors *> fPrelimASYieldsInc_Array;         ///< Preliminary calculations of yields in the inclusive region
  vector<TGraphErrors *> fPrelimNSYieldsInc_Array;         ///< Preliminary calculations of yields in the inclusive region
  vector<vector<TGraphErrors *>> fPrelimASYieldsEP_Array;  ///< Preliminary calculations of the yields in each EP region
  vector<vector<TGraphErrors *>> fPrelimNSYieldsEP_Array;  ///< Preliminary calculations of the yields in each EP region

  vector<TGraphErrors *> fPrelimNSYieldsOutOverIn_Array; ///< Preliminary out-of-plane yield / in-plane yield
  vector<TGraphErrors *> fPrelimASYieldsOutOverIn_Array; ///< Preliminary out-of-plane yield / in-plane yield


private:
	TaskEventPlane              (const TaskEventPlane&); // not implemented
	ClassDef(TaskEventPlane, 1) // Class to analyze intermediate results of PlotGHCorrelations2


};

#endif
