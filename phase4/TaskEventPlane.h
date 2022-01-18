#ifndef TASKEVENTPLANE_H
#define TASKEVENTPLANE_H

#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
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
  void ProcessMCGenFlow();
  void FitFlow();
  void LoadPublishedFlow();
	void InitArrays();
	void DrawRawOmniPlots();
	void DrawFitOmniPlots();
	void DrawRescaleOmniPlots();
  void PrelimCalculation();
  void PrelimCalculation_Step(Int_t iV);
	void DrawOmniPlots_Type(vector<TH1D *> fHists, TString fLabel, vector<TF1 *> fFits = {});
//	void DrawOmniPlots_Type(vector<TH1D *> fHists, TString fLabel, Bool_t bIncludeFit = false);
	void DoRPFThing();
	void DoRPFThing_Step(vector<TH1D *> fHists, TString fLabel,Int_t iObsBin, double fV2T_Fixed, double fV4T_Fixed);
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
  void AnalyzeFlowParameters();
  void PrepareParameterArrays();
  void ProduceVariants();
  TH2F * SubtractVariantFits(TH1D * fHist, int iVersion, int iObs, int iEPBin);
	void ProcessFitParams();
	void PlotFitParams();
  TH1D * BuildOverSubQAHist(TH1D * fHist);
	void SubtractBackground();
  void DrawOmniSandwichPlots();
	void DrawOmniSandwichPlots_Step(Int_t iV, Int_t iObsBin);
  void Rescale();
  void RescaleRegion(Int_t iV, Int_t iObsBin, Int_t iRegion); // 0 for full, 1 for near, 2 for far
  double GetFlowVNAFromObsBin(int N, int iObsBin);
  double GetFlowVNAeFromObsBin(int N, int iObsBin);
  void SaveOutput();
	void Run_Part1();
  void Run_Part2();


	void Debug(Int_t input);
	void SetDebugLevel(Int_t input)         { fDebugLevel = input; }

  void SetOverallMode(Int_t input)        { iOverallMode = input; }
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
  Int_t GetFixV2TMode()                   { return iFixV2T;       }
  void SetFixV2TMode(Int_t input)         { iFixV2T      = input; }
  void SetFixV3To0(Bool_t input)          { bFixV3To0    = input; }
  void SetDisableFlow(Int_t input)        { iDisableFlow = input; }
  Int_t GetDisableFlow()                  { return iDisableFlow; }

  Int_t GetEnableDeltaEtaScaling()        { return iEnableDeltaEtaScaling; }
  void SetEnableDeltaEtaScaling(Int_t input) { iEnableDeltaEtaScaling = input; }

  Int_t GetNRebin()                       { return nRebin;        }
  void SetNRebin(Int_t input)             { nRebin      = input;  }

	void SetIntermediateInputFile(Int_t evtPlaneIndex, TFile * inputFile) {
		if (evtPlaneIndex == -1) { fInputFileAll = inputFile; return; }
		if (evtPlaneIndex >= 0 && evtPlaneIndex <= 2)  { fInputFileEvt[evtPlaneIndex] = inputFile; }
	}

  int GetObservable()                     { return fObservable; }
  TString GetObservableName()             { return fObservableName; }
  int GetNObsBins()                       { return nObsBins;    }
  int GetCentBin()                        { return iCentBin;    }
  int GetEPRSet()                         { return iEPRSet;     }

  int GetNumVariants()                    { return iNumVariants; }
  void SetNumVariants(int input)          { iNumVariants=input; }


  RPF_Functor * GetRPFFunctor(int iObsBin)  { return fRPFFunctors[iObsBin]; }

  void SetFlowFinderMode(int input)       { iFlowFinderMode=input; }
  int GetFlowFinderMode()                 { return iFlowFinderMode; }

  double GetMCRescaleFactor()             { return fMCRescaleFactor; }

  int GetFlowTermModeTrigger()                   { return iFlowTermModeTrigger; }
  void SetFlowTermModeTrigger(int input)       { iFlowTermModeTrigger = input; }
  int GetFlowTermModeAssoc()                   { return iFlowTermModeAssoc; }
  void SetFlowTermModeAssoc(int input)       { iFlowTermModeAssoc = input; }
  int GetFlowSource()                     { return iFlowSource; }
  void SetFlowSource(int input)           { iFlowSource = input; }

  int GetFlowV3Mode()                    { return iFlowV3Mode; }
  void SetFlowV3Mode(int input)            { iFlowV3Mode = input; }

  int GetFlowV1Mode()                     { return iFlowV1Mode; }
  void SetFlowV1Mode(int input)           { iFlowV1Mode = input; }
  int GetFlowV5Mode()                     { return iFlowV5Mode; }
  void SetFlowV5Mode(int input)           { iFlowV5Mode = input; }
  int GetFlowV6TMode()                     { return iFlowV6TMode; }
  void SetFlowV6TMode(int input)           { iFlowV6TMode = input; }
  int GetFlowV6AMode()                     { return iFlowV6AMode; }
  void SetFlowV6AMode(int input)           { iFlowV6AMode = input; }

  void SetV3CalcChoice(float input)         { fV3CalcChoice = input; }
  float GetV3CalcChoice()                  { return fV3CalcChoice; }

  int GetFixV4Threshold()                  { return iFixV4Threshold; }
  void SetFixV4Threshold(int input)        { iFixV4Threshold = input; }

  int GetNegativeVnMode()                   { return iNegativeVnMode; }
  void SetNegativeVnMode(int input)         { iNegativeVnMode = input; }

  double GetGlobalV1Max()                  { return fV1_AbsMax; }
  void SetGlobalV1Max(double input)        { fV1_AbsMax = input; }

  double GetGlobalV2TMax()                  { return fV2T_AbsMax; }
  void SetGlobalV2TMax(double input)        { fV2T_AbsMax = input; }
  double GetGlobalV2AMax()                  { return fV2A_AbsMax; }
  void SetGlobalV2AMax(double input)        { fV2A_AbsMax = input; }

  double GetGlobalV3Max()                  { return fV3_AbsMax; }
  void SetGlobalV3Max(double input)        { fV3_AbsMax = input; }

  double GetGlobalV4TMax()                  { return fV4T_AbsMax; }
  void SetGlobalV4TMax(double input)        { fV4T_AbsMax = input; }
  double GetGlobalV4AMax()                  { return fV4A_AbsMax; }
  void SetGlobalV4AMax(double input)        { fV4A_AbsMax = input; }

  double GetGlobalV5Max()                  { return fV5_AbsMax; }
  void SetGlobalV5Max(double input)        { fV5_AbsMax = input; }

  double GetGlobalV6TMax()                  { return fV6T_AbsMax; }
  void SetGlobalV6TMax(double input)        { fV6T_AbsMax = input; }
  double GetGlobalV6AMax()                  { return fV6A_AbsMax; }
  void SetGlobalV6AMax(double input)        { fV6A_AbsMax = input; }


  // Global Minimum Values. These involve the Negative Vn Mode and abs max
  double GetGlobalV1Min() { if (iNegativeVnMode==3) return 0; else return -fV1_AbsMax; }

  double GetGlobalV2TMin() { return 0;}
  double GetGlobalV2AMin() { return 0;}

  double GetGlobalV3Min() { if (iNegativeVnMode==1 || iNegativeVnMode==3) return 0; else return -fV3_AbsMax; }

  double GetGlobalV4TMin() { return 0;}
  double GetGlobalV4AMin() { return 0;}

  double GetGlobalV5Min() { if (iNegativeVnMode==1 || iNegativeVnMode==3) return 0; else return -fV5_AbsMax; }

  double GetGlobalV6TMin() { return 0;}
  double GetGlobalV6AMin() { return 0;}

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

  TGraphErrors * GetTriggerV3()           { return gTrigger_V3; }

  TGraphErrors * GetTrackBv()             { return gTrack_Bv; }
  TGraphErrors * GetTrackV2()             { return gTrack_V2; }
  TGraphErrors * GetTrackV4()             { return gTrack_V4; }
  TGraphErrors * GetTrackV6()             { return gTrack_V6; }

  TGraphErrors * GetTrackV3()             { return gTrack_V3_EP3; }

  int GetNTotalParMethod(int input)      { return fNTotalParametersMethod[input]; }

  // Methods for adding to the NumFreePar, FreeMaskArray, ParNamesArray, ParMuArray, ParSigmaArray
  void UpdateNumFreePar(int iVersion, int iObsBin, int nFreePar);
  void UpdateFreeMaskArray(int iVersion, int iObsBin, int iParIndex, bool isFree);
  void UpdateParNamesArray(int iVersion, int iObsBin, int iParIndex, TString sParName);
  void UpdateParMuArray(int iVersion, int iObsBin, int iParIndex, double fParMu);
  void UpdateParSigmaArray(int iVersion, int iObsBin, int iParIndex, double fParSigma);

  // This method works if we create the TMatrixDSym objecs in runRPF.py
  void InputPyCovMatrix(int i, TMatrixDSym matrix) { fCovMatrices[i].push_back(matrix); }
  void CreateCovarianceMatrix(int iVersion, int iObsBin, int nPar);
  void UpdateCovarianceMatrix(int iVersion, int iObsBin, int iPar, int jPar, double fCovValue);



  void InputPyBkgChiSqGraph(TGraphErrors * input)        { fPyBkgChiSqGraph=input;   }
  void InputPyBkgParamGraph(int i, TGraphErrors * input) { fPyBkgParGraphs[i]=input; }
  void InputPyRPSChiSqGraph(TGraphErrors * input)        { fPyRPSChiSqGraph=input;   }
  void InputPyRPSParamGraph(int i, TGraphErrors * input) { fPyRPSParGraphs[i]=input; }


  static const Int_t kRPFColor = kViolet-5;

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


  TString GetLabel() { return sLabel; }
  TString GetLabel2() { return sLabel2; }

  void SetLabel(TString input) { sLabel = input; }
  void SetLabel2(TString input) { sLabel2 = input; }

  TLegend * DrawGeneralInfo(TCanvas * canv, double xMin, double xMax, double yMin, double yMax);

protected:

  TString sLabel  = "";
  TString sLabel2 = "";

  bool bAllowNegativeVn = false;

  int iNegativeVnMode = 2;
  // Mode 0 -> No limitations.
  // Mode 1 -> Only V1 can be negative
  // Mode 2 -> No Even negative Vn
  // Mode 3 -> No Negative Vn

  // B parameters are fractions of the average value of the nearside far eta histogram
  double fB_Init = 1.0;
  double fB_GlobalMin = 0.3;
  double fB_GlobalMax = 1.5;

  double fV1_AbsMax = 0.2; // this is v1*v1
  double fV2T_AbsMax = 0.2;
  double fV2A_AbsMax = 0.3;
  double fV3_AbsMax = 0.05; // this is v3*v3
  // note that sqrt(.1) = 0.316
  double fV4T_AbsMax = 0.1;
  double fV4A_AbsMax = 0.3;

  double fV5_AbsMax = 0.05; // this is v5*v5
  double fV6T_AbsMax = 0.1;
  double fV6A_AbsMax = 0.1;




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
  Double_t fEPRes[6]  = {0,1,1,1,1,0};
  Double_t fEPRes_Err[6]  = {0,0,0,0,0,0};
  Double_t fEP3Res[6] = {0,-1,-1,-1,-1,0};
  Double_t fEP4Res[6] = {0,-1,-1,-1,-1,0};


  void SetEPTitles(vector<vector<TH1D *>> HistArray);

  Int_t fIsMCGenMode;                       ///< 0 = data, 1 = mcGen mode

  Int_t iOverallMode = 0;                  ///< 0=RPF,1=ZYAM,2=FarEtaAve

	Int_t fDebugLevel;                        ///< For Debugging Purposes

	TFile *fInputFileAll;                     ///< File with all event plane windows
	TFile *fInputFileEvt[3];                  ///< Files In, Mid, And Out plane windows

	TFile *fOutputFile;                       ///< File where output objects are to be stored

  TGraphErrors * fEP2RGraph = 0; // EPR with respect to 2nd order EP
  TGraphErrors * fEP3RGraph = 0;

  TString sFlowGraphPath = "/home/moliver/cern/gammaHadron/wrk/FlowMeasurements/FlowGraphs.root";
  TGraph2DErrors * fV2Graph2D = 0;          ///< charged v2 vs pt and cent from ALICE published (1804.02944)
  TGraph2DErrors * fV3Graph2D = 0;          ///< charged v3 vs pt and cent
  TGraph2DErrors * fV4Graph2D = 0;          ///< charged v4 vs pt and cent
  // 2D Plots whose points are centered on point+err
  TGraph2DErrors * fV2Graph2DErrUp = 0;
  TGraph2DErrors * fV3Graph2DErrUp = 0;
  TGraph2DErrors * fV4Graph2DErrUp = 0;
  // 2D Plots whose points are centered on point-err
  TGraph2DErrors * fV2Graph2DErrDown = 0;
  TGraph2DErrors * fV3Graph2DErrDown = 0;
  TGraph2DErrors * fV4Graph2DErrDown = 0;

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

  Int_t iNumVariants = 500;                  ///< Number of variations of the RPF parameters

	Int_t iRPFMode;                           ///< Which mode to use RPF method in #FIXME write them down somewhere

  Int_t nRPFMethods = 1;                   ///< How many versions of the RPF to use? 1 is just the C++ implementation, 3 includes both the Bkg and RPSignal fit versions

  Int_t iEPRSet = 0;                            ///< Which set of event plane resolutions to use. Default to MB data measurement
  Int_t iCentBin = 0;                           ///< Which centrality bin for EPR? (0-10,10-30,30-50,50-80)
	
  Double_t fCentArray[5] = {0.,10.,30.,50.,80.};

  Int_t nRebin = 0;                         ///< How much to rebin the delta phi histograms


  Int_t iDisableFlow = 0;                  ///< Set to 1 to disable all Vn

  Int_t iFixV2T = 0;                       ///< Whether to fix the V2T and V4T to the value found in the first iFixV2T bins;
  Bool_t bFixV3To0 = 0;                     ///< Whether to fix the V3AV3T to 0

  Int_t iFlowFinderMode = 0;                ///< How to get VN estimates
  // See DoRPFThing_Step() for the most up to date definitions of the modes
  

  void RunDeltaEtaScaling();

  // Flow analysis info
  // First version (higher pt resolution)
  //const int kNTrackPtBins = 14;
  //std::vector <double> fTrackPtBins = {0.15, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 17};
  // Second version (close to correlation track pt bins)

  const int kNTrackPtBins = 8;
  std::vector <double> fTrackPtBins ={0.2,0.4,0.8,1.5,2.5,4,7,11,17};


  Int_t iFlowTermModeTrigger = 0;           ///< How to apply the phase 1 flow information for trigger
                                            // 0: none
                                            // 1: fix V2T, V4T
                                            // 2: Set V2T, V4T ranges based on 1 sigma 

  Int_t iFlowTermModeAssoc = 0;             ///< How to apply the phase 1 flow information for associated
                                            // 0: none
                                            // 1: Fix V2A, V4A
                                            // 2: Fix V2A, V4A ranges based on 1 sigma
  // FIXME not fixing V4A now


  Int_t iFlowSource = 0;                    ///< Where to get flow graphs from
                                            // 0: This analysis
                                            // 1: ALICE published results


  Int_t iFlowV1Mode = 0;                    ///< Whether to include a free V1 parameter. Only applicable in python (RPF2)
                                            // 0: none (default)
                                            // 1: free

  Int_t iFlowV3Mode = 0;                    ///< How to deal with V3 term (in addiont to FixV30 variable)
                                            // 0: free (default)
                                            // 1: Fixed to Calculated v3*v3 value
                                            // 2: Bound by Calculated v3*v3 range

  Int_t iFlowV5Mode = 0;                    ///< How to deal with with v5tv5a term
                                            // 0 : none (default)
                                            // 1 : free

  Int_t iFlowV6TMode = 0;                   ///< How to deal with with v6t term
                                            // 0 : none (default)
                                            // 1 : enable, free
                                            // 2 : enable, fix to measured EP value
                                            // 3 : enable, fix near measured EP value

  Int_t iFlowV6AMode = 0;                   ///< How to deal with with v6a term
                                            // Same definitions as above


  Float_t fV3CalcChoice = 0.0;              // Choice of varying the calculated v3 value. Use value of best estimate + iV3CalcChoice * V3V3Error

  Int_t iFixV4Threshold = 30;                ///< For iObsBin >= iFixV4Threshold,
                                            // the V4a and V4t are fixed to 0, unless already fixed to a value


	Int_t fObservable;                        ///< Observable for the current analysis
	TString fObservableName;                  ///< Name of the current observable (for plot labels)

	Int_t nObsBins;                           ///< How many bins we have of the current observable. Determined at run time.

	vector<Double_t> fObsBins;                ///< Bin Edges for whichever observable we are looking at

  TString fEPBinTitles[kNEPBins+1] = {"In-Plane","Mid-Plane","Out-of-Plane","All EP Angles"};

//  TString fParNames[9] = {"B","V2T","V2A","V3","V4T","V4A","V5","V6T","V6A"};
//  TString fParTitles[9] = {"B","v_{2}^{t}","v_{2}^{a}","v_{3}^{t}v_{3}^{a}","v_{4}^{t}","v_{4}^{a}","v_{5}^{t}v_{5}^{a}","v_{6}^{t}","v_{6}^{a}"};
  TString fParNames[10] = {"B","V1","V2T","V2A","V3","V4T","V4A","V5","V6T","V6A"};
  TString fParTitles[10] = {"B","v_{1}^{t}v_{1}^{a}","v_{2}^{t}","v_{2}^{a}","v_{3}^{t}v_{3}^{a}","v_{4}^{t}","v_{4}^{a}","v_{5}^{t}v_{5}^{a}","v_{6}^{t}","v_{6}^{a}"};

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
  std::vector<TH1F *> hPtEP3AngleTrack_Proj;
  std::vector<TH1F *> hPtEP4AngleTrack_Proj;

  TH1D * fTrackPtFromTrackPsi = 0;        ///< Histogram of track pT made from projecting the TrackPsiEPPtCent TH3 

  // Trigger and Track VN from this analysis

  // Vn Information from fits in Phase 1 or to be recalculated
  // To be recalculated here from the raw pi0 candidate vs event plane
  TGraphErrors * gTrigger_Bv = 0;
  TGraphErrors * gTrigger_V2 = 0;
  TGraphErrors * gTrigger_V4 = 0;
  TGraphErrors * gTrigger_V6 = 0;

  TGraphErrors * gTrigger_V3 = 0;

  // The trigger vs event plane prior to sideband subtraction
  TGraphErrors * gTrigger_Bv_Presub = 0;
  TGraphErrors * gTrigger_V2_Presub = 0;
  TGraphErrors * gTrigger_V4_Presub = 0;
  TGraphErrors * gTrigger_V6_Presub = 0;

  TGraphErrors * gTrigger_V3_Presub = 0;

  TGraphErrors * gTrack_Bv = 0;
  TGraphErrors * gTrack_V2 = 0;
  TGraphErrors * gTrack_V4 = 0;
  TGraphErrors * gTrack_V6 = 0;

  TGraphErrors * gTrack_Bv_EP3 = 0;
  TGraphErrors * gTrack_V3_EP3 = 0;

  TGraphErrors * gTrack_Bv_EP4 = 0;
  TGraphErrors * gTrack_V4_EP4 = 0;




  // Additional Flow Histograms from MCGen Toy and Flow Studies


  TH1F * hEP2CosDeltaPsiIncl = 0;
  TH1F * hEP2CosDeltaPsiToyOnly = 0;
  TH1F * hEP3CosDeltaPsiIncl = 0;
  TH1F * hEP3CosDeltaPsiToyOnly = 0;
  TH1F * hEP4CosDeltaPsiIncl = 0;
  TH1F * hEP4CosDeltaPsiToyOnly = 0;

  double fMCEP2Incl = -1;
  double fMCEP2ToyOnly = -1;
  double fMCEP3Incl = -1;
  double fMCEP3ToyOnly = -1;
  double fMCEP4Incl = -1;
  double fMCEP4ToyOnly = -1;



  TH1F * hToyV2EP = 0;
  TH1F * hToyV3EP = 0;
  TH1F * hToyV4EP = 0;
  TH1F * hToyV2RP = 0;
  TH1F * hToyV3RP = 0;
  TH1F * hToyV4RP = 0;

  TH1F * hInclusiveV2EP = 0;
  TH1F * hInclusiveV3EP = 0;
  TH1F * hInclusiveV4EP = 0;
  TH1F * hInclusiveV2RP = 0;
  TH1F * hInclusiveV3RP = 0;
  TH1F * hInclusiveV4RP = 0;

  TH1F * hToyTriggerV2EP = 0;
  TH1F * hToyTriggerV3EP = 0;
  TH1F * hToyTriggerV4EP = 0;
  TH1F * hToyTriggerV2RP = 0;
  TH1F * hToyTriggerV3RP = 0;
  TH1F * hToyTriggerV4RP = 0;

  TH1F * hInclusiveTriggerV2EP = 0;
  TH1F * hInclusiveTriggerV3EP = 0;
  TH1F * hInclusiveTriggerV4EP = 0;
  TH1F * hInclusiveTriggerV2RP = 0;
  TH1F * hInclusiveTriggerV3RP = 0;
  TH1F * hInclusiveTriggerV4RP = 0;


  // Calculated V2Trigger graphs for comparison to RPF parameter
  TGraphErrors * fGraphFlowToyV2TriggerEP = 0;
  TGraphErrors * fGraphFlowToyV2TriggerRP = 0;
  TGraphErrors * fGraphFlowInclusiveV2TriggerEP = 0;
  TGraphErrors * fGraphFlowInclusiveV2TriggerRP = 0;


  // Track Vn from ALICE published inclusive charged particle Vn
  TGraphErrors * gAliTrack_V2 = 0;
  TGraphErrors * gAliTrack_V3 = 0;
  TGraphErrors * gAliTrack_V4 = 0;

  // Calculated V3TV3A from ALICE published inclusive charged particle Vn
  TGraphErrors * gAliTrack_CalcV3TV3A = 0;

  // Calculated V3TV3A from my analysis
  TGraphErrors * gCalcV3TV3A = 0;

  // Histograms from toy model V3
  TH1F * hMCGenToyV3V3EP = 0; // Calculatied with Toy V3 vs EP3 using Toy V3 as Triggers
  TH1F * hMCGenToyV3V3RP = 0; // Calculatied with Toy V3 vs RP3

  TH1F * hMCGenInclusiveTriggerV3InclusiveV3EP = 0; // Inclusive trigger and toy
  TH1F * hMCGenInclusiveTriggerV3InclusiveV3RP = 0; // Inclusive trigger and toy


	// Histograms

  // Delta Eta Projections

  vector<TH1D *> fNearSideSubDEtaFinalAll;    ///< SB Sub Nearside Delta Eta (all EP)
  vector<vector<TH1D *>> fNearSideSubDEtaFinalEP;    ///< SB Sub Nearside Delta Eta (EP Bins)


  // Delta Phi Projections

	vector<TH1D *>         fFullDPhiProjAll;          ///< Full Projections in DPhi.
	vector<vector<TH1D *>> fFullDPhiProj;             ///< Full Projections in DPhi.  First index is observable bin, second is event plane bin
  TH3F *hHistTrackPsiEPPtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
  TH3F *hHistTrackPsiEP3PtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
  TH3F *hHistTrackPsiEP4PtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
 // TH3F *hHistTrackPsiEPPtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
	vector<TH1D *>         fNearEtaDPhiProjAll;       ///< Near Eta (Signal) Projections in DPhi.
	vector<vector<TH1D *>> fNearEtaDPhiProj;          ///< Near Eta (Signal) Projections in DPhi.  First index is observable bin, second is event plane bin
	vector<TH1D *>         fFarEtaDPhiProjAll;       ///< Far Eta (Background) Projections in DPhi.
	vector<vector<TH1D *>> fFarEtaDPhiProj;          ///< Far Eta (Background) Projections in DPhi.  First index is observable bin, second is event plane bin

  double fEtaScaleRange = TMath::Pi();  ///< range used for calculating awayside scales.


  // Scale information for far-eta vs full-eta and near-eta
  int iEnableDeltaEtaScaling = 0; ///< Control for delta eta scaling
  // 0 is off, 1 is central scaling, 3 and 4 are down and up 1 sigma
  // 7 is pedestal subtraction (only valid for JEWEL keep Recoils)

  // RPF subtraction for the Full eta or Near eta, which ever is not used here
  int iScaleToFullEta = 0; // 0 for scaling FarEta to match Near Eta, 1 for scaling FarEta to match Full Eta
  // The Scale of the RPF function needs to be adjusted, as the Full Eta (used for the awayside) and the near eta (used for the nearside) have different
  // scales or pedestals that need to be applied.


  // First axis is ObsBin, Second is iEP Index [all,in,mid,out] or [in,mid,out,all]
  vector<vector<double>> fScaleFullOverFar={};       ///< Scale of awayside in full delta eta vs near delta eta
  vector<vector<double>> fScaleFullOverFarErr={};       ///< 
  // uncertainty?
  vector<vector<double>> fScaleNearOverFar={};
  vector<vector<double>> fScaleNearOverFarErr={};      ///<

  // For Pedestal subtraction
  vector<vector<double>> fScaleFullMinusFar={};       ///< Integral of awayside in full delta eta vs near delta eta
  vector<vector<double>> fScaleFullMinusFarErr={};    ///< 

  vector<vector<double>> fScaleNearMinusFar={};
  vector<vector<double>> fScaleNearMinusFarErr={};      ///<

  //  Scales actually used
  vector<vector<double>> fUsedDEtaScalingScales = {}; ///<
  vector<vector<double>> fUsedDEtaScalingScalesErr = {}; ///<



  // could store scales in a histogram -> easy to analyze and to track errors.



  // RPF Fit Residuals
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

  // Variants
	vector<vector<TH2F *>>         fFullDPhiProjAll_Sub_Variants;          ///< Full Projections in DPhi.
	vector<vector<vector<TH2F *>>> fFullDPhiProj_Sub_Variants;             ///< Full Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<TH2F *>>         fNearEtaDPhiProjAll_Sub_Variants;       ///< Near Eta (Signal) Projections in DPhi.
	vector<vector<vector<TH2F *>>> fNearEtaDPhiProj_Sub_Variants;          ///< Near Eta (Signal) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<TH2F *>>         fFarEtaDPhiProjAll_Sub_Variants;        ///< Far Eta (Background) Projections in DPhi.
	vector<vector<vector<TH2F *>>> fFarEtaDPhiProj_Sub_Variants;           ///< Far Eta (Background) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin


  // Oversubtraction QA Histograms
	vector<vector<TH1D *>>         fFullDPhiProjAll_OverSubQA={};          ///< Full Projections in DPhi.
	vector<vector<vector<TH1D *>>> fFullDPhiProj_OverSubQA={};             ///< Full Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<TH1D *>>         fNearEtaDPhiProjAll_OverSubQA={};       ///< Near Eta (Signal) Projections in DPhi.
	vector<vector<vector<TH1D *>>> fNearEtaDPhiProj_OverSubQA={};          ///< Near Eta (Signal) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<TH1D *>>         fFarEtaDPhiProjAll_OverSubQA={};        ///< Far Eta (Background) Projections in DPhi.
	vector<vector<vector<TH1D *>>> fFarEtaDPhiProj_OverSubQA={};           ///< Far Eta (Background) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin


  // After background sub. and trigger rescaling
	vector<vector<vector<TH1D *>>> fFullDPhiProj_Rescale;             ///< Full Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<vector<TH1D *>>> fNearEtaDPhiProj_Rescale;          ///< Near Eta (Signal) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<vector<TH1D *>>> fFarEtaDPhiProj_Rescale;          ///< Far Eta (Background) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin

  // Variants
	vector<vector<vector<TH2F *>>> fFullDPhiProj_Rescale_Variants;             ///< Full Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<vector<TH2F *>>> fNearEtaDPhiProj_Rescale_Variants;          ///< Near Eta (Signal) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin
	vector<vector<vector<TH2F *>>> fFarEtaDPhiProj_Rescale_Variants;          ///< Far Eta (Background) Projections in DPhi.  First index is RPF version, 2nd is observable bin, 3rd is event plane bin




  // Fits and Maths
//	vector<TF1 *> fRPFFits;                   ///< Array of the fits
//	vector<vector<TF1 *>> fRPFFits_Indiv;     ///< Array of resultant fits (1st index is Obs, 2nd index is EP bin)

//  vector<TH1D *> fRPF_Residuals = {};                ///< Array of (Data - Fit) / Fit [iObsBin] // the original fit
//  vector<vector<TH1D *>> fRPF_Residuals_Indiv = {};  ///< Array of (Data - Fit) / Fit [iObsBin,iEPBin]

  vector<RPF_Functor *> fRPFFunctors = {};            ///< Array of functors

  vector<vector<TF1 *>> fRPFFits;                   ///< Array of the fits
	vector<vector<vector<TF1 *>>> fRPFFits_Indiv;     ///< Array of resultant fits (1st index is RPF Method, 2nd index is Obs, 3rd index is EP bin)

  // Old approach array of different variants:
  vector<vector<vector<TF1 *>>> fRPFFits_Variants={};    ///< Array of global fits, final index is the variant id
	vector<vector<vector<vector<TF1 *>>>> fRPFFits_Indiv_Variants={};     ///< Array of resultant fits (1st index is RPF Method, 2nd index is Obs, 3rd index is EP bin, 4th index is the variant ID)

  // New approach: array of one variant TF1 per function, with arrays of parameters to cycle through
  //vector<vector<TF1 *>> fRPFFits_Variant;                   ///< Array of the fits
	vector<vector<TF1 *>> fRPFFits_Indiv_Variant;     ///< Array of resultant fits (1st index is RPF Method, 2nd index is Obs, 3rd index is EP bin)
  //vector<vector<vector<vector<vector<double>>>>> fRPFFits_Indiv_Parameters_Variants = {}; ///< Array of parameters for each variant
  // removing dimension of EPBin
  vector<vector<vector<vector<double>>>> fRPFFits_Parameters_Variants = {}; ///< Array of parameters for each variant

  // Vectors of parameters information
  // Since different obs bins can have different free parameters, need to track closely;
  // For all, the first index is the method, the second is the obsbin
  // Third is the parameter index
  // Initializing with two empty sub arrays. More can be added for Incl, RPDep, Redux RPF options
  vector<vector<int>> fNumFreeParameters = {{},{}};
  vector<vector<vector<bool>>> fParFreeMaskArray = {{},{}};
  // These arrays should only be filled with free parameters
  vector<vector<vector<TString>>> fParNamesArray = {{},{}};
  vector<vector<vector<double>>> fParMuArray {{},{}};
  vector<vector<vector<double>>> fParSigmaArray {{},{}};

  vector<vector<TTree *>> fParameterTreeArray;  ///< Array of trees storing the sets of parameter variants. First bin is RPF method, 2nd bin is Obs Bin

  // Vectors of the covariance matrices
  vector<vector<TMatrixDSym>> fCovMatrices;         ///< Array of fit covariance matrices. First index is the RPF method, second is the Obs Bin
  // Vectors of the parameter correlation matrices
  vector<vector<TMatrixDSym>> fCorMatrices;         ///< Array of fit correlation matrices. First index is the RPF method, second is the Obs Bin


  // this one still needs to be updated for the outer loop, though I'm not sure if it is actually used
  vector<vector<TH1D *>> fRPF_Residuals = {};                ///< Array of (Data - Fit) / Fit [iRPFMethod][iObsBin] // the original fit

  vector<vector<vector<TH1D *>>> fRPF_Residuals_Indiv = {};  ///< Array of (Data - Fit) / Fit [iRPFMethod][iObsBin,iEPBin]


  vector<int> fNTotalParametersMethod = {10,8}; ///< Array of the total parameters in the fits for each mode. It is higher for Method1 because it can have v5, v6a, v6t parameters


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
