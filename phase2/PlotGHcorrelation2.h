#ifndef PLOTGHCORRELATION2_H
#define PLOTGHCORRELATION2_H

#include <Riostream.h>
#include <TString.h>
//#include <TObject.h>

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

class PlotGHcorrelation2 : public TObject {

public:
	PlotGHcorrelation2();
	PlotGHcorrelation2(Int_t observable,Int_t centrality, Int_t eventPlane);

	virtual ~PlotGHcorrelation2()  { ; }

	void Run();

	//..setters for the analysis
    void SetInputTrainFile(TFile *inputFileSE, TFile *inputFileME,TString inputNameSE,TString inputNameME) { fInputFileSE= inputFileSE; fInputFileME = inputFileME; fWagonNameSE = inputNameSE; fWagonNameME = inputNameME;}
    //..load external root file with already projected histograms
    void SetInputRootFile(TFile *inputFile) { fRootFile    = inputFile;}
    void SetInputRootFileME(TFile *inputFile) { fRootFileME    = inputFile;}
    void SetGammaPi0(Int_t input)           { fGammaOrPi0  = input;  
                                             //if (input!=0) fTriggerName = "#pi^{0}";
                                             //else fTriggerName = "#gamma"; 
                                             fTriggerName = "trig";}
    void SetMCMode(Int_t input)             { iMCMode = input; }
    void SetPlotOptions(TString input)      { fPlotOptions = input; }
    void SetSavePlots(Bool_t input)         { fSavePlots   = input;}
    void SetMergeMEplots(Bool_t input)      { fMergeMEplots   = input;}
    void SetPtBinRange(Int_t input1, Int_t input2) { fPtMinBin=input1; fPtMaxBin = input2;}
    void SetPlotAdvanced(Bool_t input)      { fPlotAdvancedSB = input;}
    void SetDrawVtzBins(Bool_t input)       { fPlotVtzBins = input;}
    void SetPlot2DHistos(Bool_t input)      { fPlot2DHistos = input;}
		void SetNSigma(Double_t input)          { fNSigma = input;}

    // epsilon = 0.0001 to avoid dumb binning effects
    void SetMinDEtaSignalRange(Double_t input) { fMinDEtaSignalRange = input + 0.0001; }
    void SetMaxDEtaSignalRange(Double_t input) { fMaxDEtaSignalRange = input - 0.0001; }

    void SetMEDetaRangeForNorm(Double_t input) { fMEDEtaRangeForNorm = input;}
    void SetMENormSideBins(Int_t input)     { nMENormSideBins = input;}
    void SetRebinMEForNorm(Int_t input)     { nRebinMEForNorm = input;}
    void SetUseFindLastGoodBin(Bool_t input){ fUseFindLastGoodBin = input;}
    void SetPlotMEcomparision(Bool_t input) { fPlotMoreMEstrategies = input;}
    void SetWait(Bool_t input)              { fWait              = input;}
    void SetOutputDir(TString input)        { fOutputDir = input; }
		void SetOutputLabel(TString label)      { fLabel = label; }

//    void SetNormMode(Int_t input)          { fNormMode       = input;}
protected:

	void InitArrays();
	void ReadjustVariables();
	void LoadHistograms();
	TH2D* Get2DHistoFromFile(TFile* RootFile, TString name);
	TH1D* Get1DHistoFromFile(Bool_t smaMix, TString subListName, TString name);
	void DrawWIP(TH1 *Histo,Float_t x, Float_t y, Float_t x_size, Float_t y_size);
	void SetTH1Histo(TH1 *Histo,TString Xtitle,TString Ytitle,Bool_t big=0);
	void SetTH2Histo(TH2 *Histo,TString Xtitle,TString Ytitle,TString zTitle="",Int_t rows=2);
	void ZoomYRange(TH1D *Histo,Double_t border=0.1,Double_t Range1=-1,Double_t Range2=-1);
	void DrawEtaPhi2D(TH2 *Histo, TString sZTitle="");
	void SetPlotStyle();
	TLatex* PlotTopLegend(const char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0);
	void PlotVerLineRange(Double_t x_val, Double_t yLow, Double_t yHigh, Int_t Line_Col);
	void PlotVerLine3(Double_t x_val,TH1* Histo, Double_t y_fac, Int_t Line_Col);
	void DrawAliceInfoBox(TObject* Histo);
	void DrawAlicePerf(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size);

	void DrawSEplots();
	void AnalyzePionAccRej();
	void MergeMEplots();
	void NormalizeSEsignal();
	void NormalizeMEsignal();
	void ScaleMEbackground(TH2D* Histo,Double_t lowRange,Double_t highRange,TCanvas* Can,Int_t CanvasPad);
	void Plot2DHistograms();
	void PlotCorrected2DHistograms();
	void PlotCompare2DHistoCorrection();

	void DetermineEtaWidths(TH2D* corrHistoSE[]);
	void DeterminePhiWidths();
	void DrawWidths();
	void FitGaussAndDraw(Int_t bin,TH1D* corrProjHistoSE,TF1* Func1, TF1* Func2,TF1* Func3,Bool_t EtaPhi);
	static Double_t PolyTwoGaussFitFunc(Double_t* x_val, Double_t* par);

	void FitEtaSides(TH2D* corrHistoSE[],Double_t sigma1,TCanvas* can,TCanvas* can2,TCanvas* can3);
	void FitEtaSide(TH2D* corrHistoSE,Double_t width,Double_t Sigma1,TCanvas* Can,TCanvas* Can2,TCanvas* Can3,Int_t CanvasPad);
	void FindLastGoodBin(TH2D* corrHistoSE,Int_t startSearch,Int_t direction,Int_t *lastGoodBin);
	void UpdateEtaCanvas(Int_t bin,Double_t edgeLow,Double_t edgeSigL,Double_t edgeSigR,Double_t edgeHigh);

  void ProduceDeltaPsiPlots();
  void MeasureVn();
  void DPhiQA();

	//..last functions of part1
	void SaveIntermediateResult(Int_t stage);
	void Wait();

// Zero Out a selected region of a histogram
  void zeroRegion(TH1D * hist, int minBin, int maxBin) {
    for (int i = minBin; i <= maxBin; i++) {
      hist->SetBinContent(i,0);
      hist->SetBinError(i,0);
    }
  }

	//..This belongs to part 2
	void DrawSupComponentsJoel(TF1 *func, Int_t harmonic);
	void DrawSupComponentsEli(TF1 *func, Int_t harmonic);
	TF1* allFitFuncVn(TString name,Int_t VnTerms,Double_t allPhiMin,Double_t allPhiMax);
	void SetJoelsParLimits(TF1 *func, TH1 *histoToFit,Double_t par_V10[]);
	static Double_t JoelsVnFunctionValue(const double * x, const double * p);
	static Double_t FlowFunction(Double_t* x_val, Double_t* par);


	TFile *fInputFileSE;                      ///< Input file for the SameEvent histograms
	TFile *fInputFileME;                      ///< Input file for the MixedEvent histograms
	TString fWagonNameSE;                     ///< Name of the wagon or subwagon in which the SE histograms are saved
	TString fWagonNameME;                     ///< Name of the wagon or subwagon in which the ME histograms are saved
	TString fPlotOptions;                     ///< Style for plotting 2D histograms.  Default is lego2
	TString fOutputDir;                       ///< Output directory to save the plots
	TFile *fRootFile;                         ///< Input root file with already projected histograms
	TFile *fRootFileME;                       ///< Input file for mixed event histograms (same file as SE used if NULL)
	Bool_t fSavePlots;                        ///< Whether to save plots
	Bool_t fMergeMEplots;                     ///< switch to Merge ME plots for different bins, and QA that process
  Bool_t fProducePi0AnalysisPlots = 0;      ///< Whether to make the pi0 acceptded vs rejected plots
//  Int_t  fNormMode;                         ///< 0 to normalize after ME corr, 1 to normalize SE before ME Corr
	Bool_t fPlotAdvancedSB;                   ///< switch to plot more advanced histograms
  Bool_t fPerformDPhiQA;                    ///< Switch for additional QA of corrected DPhi Correlations
//	Bool_t fPlot2MEstrategies;                ///< switch on a detailed comparision of diffent ME correction strategies
	Bool_t fUseFindLastGoodBin;               ///< Switch for using the FindLastGoodBin Method
	Bool_t fPlotMoreMEstrategies;             ///< switch on a detailed comparision of diffent ME correction strategies
	Int_t fObservable;                        ///< observable for the analysis 0=Ga bins, 1=zt bins, 2=xi bins...
	TString fObservableName;                  ///< Name of the current observable (for plot labels)
	Int_t fCent;                              ///< centrality selection for the analysis -1=all
	Int_t fEventPlane;                        ///< event plane selection for the analysis -1=all
	Bool_t fUseHistogramFile;                 ///< Identifies the analysis to run with an external root file with alredy projected histograms
	TString fLabel;                           ///< Label for creating external root file [fLabel_Observable%d_....root]
  // FIXME merge these 2 variables
	const bool fPerformance = false;          ///< Switch on in case we get performance approval again
  Bool_t bEnablePerformance = false;
  Bool_t bEnableArbUnits = false;           ///< Whether to use arbitrary units on axes for sharing
  Bool_t bNoYLabel = false;                  ///< Whether to remove Y/Z axis labels for performance purposes
	Bool_t fWait;                             ///< This enables you to look at canvases interactively

  Int_t iMCMode = 0;                        ///< What to do with MC information
                                            // 0 -> include all (equiv to data)
                                            // 1 -> Background only
                                            // 2 -> True Pi0s
                                            // 3 -> True Etas

  static const Int_t kProjFullStyle=kFullSquare;
  static const Int_t kProjFullColor=kRed+1;//kOrange+10;

	static const Int_t kGammaNBINS=9;  //9    ///< Number of 2D histograms for Gamma energy
	static const Int_t kZtNBINS=7;            ///< Number of 2D histograms for Zt of g-h pair
	static const Int_t kXiNBINS=8;            ///< Number of 2D histograms for Xi of g-h pair
  static const Int_t kNoHPtBins=8;          ///< Bins in hadron pT

	static const Int_t kNvertBins=10;          ///< Z-Vertex bins in which the ME are mixed

	static const Int_t kNvertBins_alt=8;      ///< Instead of the usual 1cm binning we will do a 4cm binning at the edge and a 2cm binning in the center
	static const Int_t kCentBINS= 4;          ///< Number of centrality bins in the analysis
	Int_t fmaxBins;                           ///< Number of bins in the analysis specified by fObservable
	Int_t fminZvtx;                           ///< minimum range for z-Vertex, if different from -10cm
	Int_t fmaxZvtx;                           ///< maximum range for z-Vertex, if different from +10cm

	//..binning in E_g, zT, and xi
	Double_t fZtStep;                          ///<
	Double_t fXiStep;                          ///<

  Bool_t fPlotVtzBins = 0;                   ///< Whether or not to draw the 2d raw SE and ME in vtz bins
  Bool_t fPlot2DHistos =  1;                 ///< Whether or not to make the big 2D plots

	Double_t fNSigma;                          ///< Number of sigma of NS peak to use to defined NearEta,FarEta peak
  Double_t fMEDEtaRangeForNorm = 0.1;        ///< Range in delta eta to project ME for normalization determination
  Double_t fMinDEtaSignalRange = 0.4 + 0.001;      ///< Minimum range in Delta eta for the signal dominated region (epsilon = 0.001 needed cause of bin rounding)
//  Double_t fMinDEtaSignalRange = 0.8 + 0.001;      ///< Minimum range in Delta eta for the signal dominated region (epsilon = 0.001 needed cause of bin rounding)
  Double_t fMaxDEtaSignalRange = 0.8 - 0.001;      ///< Maximum range in Delta eta for the signal dominated region (epsilon = 0.001 needed cause of bin rounding)
  Double_t fMaxDeltaEtaRange = 1.35-0.001; // epsilon = 0.001 to avoid bin effects

  Double_t fMaxDeltaEtaPlotRange = 1.2;      ///< How far to plot the 2D plots

  Int_t nRebinMEForNorm = 3;                 ///< Rebin number for mixed event normalization scaling
  Int_t nMENormSideBins = 1;                 ///< Number of side bins included in Mixed Event normalization

	Int_t fPtMinBin,fPtMaxBin;                 ///< Bin range for Pt (if using Zt or Xi)

	Double_t fArray_G_Bins[10];                ///< 10=kGammaNBINS+1
	Double_t fArray_ZT_Bins[8];                ///< 8=kZtNBINS+1
	Double_t fArray_XI_Bins[9];                ///< 9=kXiNBINS+1
  Double_t fArray_HPT_Bins[9];               ///< 9=kNoHPtBins+1


	Double_t fscaleFactorSBtoNS[10];           ///< scales the SB eta range to the NS eta range
	Double_t fArray_cent_Bins[5];              ///< 5=kCentBINS+1
	Double_t fArray_zVtx_Bins[kNvertBins+1];   ///< orig. zVtx binning
	Double_t fArray_zVtx_BinsAlt[kNvertBins_alt+1];///< 9=kNvertBins_alt+1 alternative zVtx binning

  static const Int_t kRebinDEtaThreshold=-1;  ///< Rebin the lowest bins due to bad statistics

	//..delta phi integration ranges
	static const Int_t kNDeltaPhiBins=8;       ///<

	Double_t fDeltaPhiBins[8];                 ///<  8=kNDeltaPhiBins  <={180,125,180,24,156,24,132,24};
	const Int_t fDoublePhiBins=16;             ///<

	Color_t fColorSceme[6];                    ///< saved color palett for plotting various histograms
	TBox *fBoxes[16];                          ///< 16= fDoublePhiBins

	//..2D corr histograms
	TH2D **fDetaDphi_SE[10];                   ///< 2D array of TH2 SE in bins of Eg/Zt/Xi and z-vertex.
	TH2D **fDetaDphi_ME[10];                   ///< 2D array of TH2 ME in bins of Eg/Zt/Xi and z-vertex.
	TH2D *fDetaDphi_ME_alt1[10];               ///< add all z-vertex bins - keep individual fObservable bins
	TH2D *fDetaDphi_ME_alt2[10][8];            ///< add all fobservable bins - keep zvertex bins
	TH2D *fDetaDphi_ME_alt3[10];               ///< add all fobservable bins - keep zvertex bins [kNvertBins]
	TH2D *fDetaDphi_SE_alt1[10];               ///< add all z-vertex bins - keep individual fObservable bins
	TH2D *fDetaDphi_SE_alt2[10][8];            ///< add all fobservable bins - keep zvertex bins
	TH2D *fDetaDphi_SE_alt3[10][10];           ///< add all fobservable bins - keep zvertex bins [fObserv][kNvertBins]
	TH2D *fsumCorrSE[10];                      ///< 1D array of corrected TH2 one histogram per bin Eg/Zt/Xi
	TH2D *fsumCorrSE_alt1[10];                 ///< 1D array of corrected TH2 one histogram per bin Eg/Zt/Xi
	TH2D *fsumCorrSE_alt2[10];                 ///< 1D array of corrected TH2 one histogram per bin Eg/Zt/Xi
	TH2D *fsumCorrSE_alt3[10];                 ///< 1D array of corrected TH2 one histogram per bin Eg/Zt/Xi


	// TH1 arrays for projections to be saved to Intermediate results
	TH1D * fDeta_Proj[10];                     ///< 1D array of projections of deta in NS
	TH1D * fDeta_ProjSub[10];                  ///< 1D array of projections of deta after AS subtraction
	TH1D * fsumCorrSE_ProjFull[10];            ///< 1D array of projections of the full dphi projection
	TH1D * fsumCorrSE_NearEta[10];             ///< 1D array of projections of the near eta (signal) region
	TH1D * fsumCorrSE_FarEta[10];              ///< 1D array of projections of the far eta (background) region

  TH1D * fFFTsumCorrSE_ProjFull[10];         ///< 1D array of FFTs of corrected full dphi projection

	// Trigger Information
	Int_t fGammaOrPi0;                         ///< 0 for Gammas, 1 for Pi0s 2 p0 SB1, 3 Pi0 SB2.  Only needed for normalization right now.
	TString fTriggerName;
	TH1D *fTrigger_SE[10];                     ///< 1D array of z-vertex histograms for each trigger class in SE.

	TH1D *fTriggerPt;                          ///< Histogram of trigger (pT of Pi0 or E of Gamma)
  TH1D *fTriggerPtWithinEPBin;               ///< Histogram of trigger pt within the chosen EvtPlane bin

  // Track information
  TH1D * fTrackPtProjectionSE = 0;           ///< Histogram of track pT made from projecting Corr THnSparse for Same Events (this is a biased distribution)
  TH1D * fTrackPtProjectionME = 0;           ///< Histogram of track pT made from projecting Corr THnSparse for Mixed Events (this should be an unbiased distribution)

  TH1D * fTrackPtFromTrackPsi = 0;           ///< Histogram of track pT made from projecting the TrackPsiEPPtCent TH3 

  TH1D * fPhase2Purity = 0;                  ///< Purity made from Triggerhist ThnSparse (only in MC)
  TH1D * fMCTriggerDist = 0;                 ///< Distribution of MC Truth Status of Triggers

	TH2F *fMassPtPionAcc;                      ///< Histogram of Mass vs Pt for Accepted Pi0ns
	TH2F *fMassPtPionRej;                      ///< Histogram of Mass vs Pt for Rejected Pi0ns
	TH3F *fMassPtCentPionAcc;                      ///< Histogram of Mass vs Pt vs Centrality for Accepted Pi0ns
	TH3F *fMassPtCentPionRej;                      ///< Histogram of Mass vs Pt vs Centrality for Rejected Pi0ns

	TH1D *fMassPtPionAccProj[kGammaNBINS];     ///< Histograms storing the mass dist. in each pt bin
	TH1D *fMassPtPionRejProj[kGammaNBINS];     ///< Histograms storing the rejected mass dist. in each pt bin

  TH2F *hPtEPAnglePionAcc = 0;               ///< Accepted Pi0 Candidates vs event plane (n=2)
  // may have to get final PtEPAnglePionAcc from Pi0Cand trains (that were broken up into centrality)

  TH3F *hHistTrackPsiEPPtCent = 0;           ///< Accepted Tracks vs event plane (broken down by centrality)
  std::vector<TH1F *> hPtEPAnglePionAcc_Proj;



	//..new histograms from the analysis
	TH1 *fYield_VS_Eg[4];                      ///<
	TH1 *fYield_VS_Zt[4];                      ///<
	TH1 *fYield_VS_Xi[4];                      ///<

	TH1D* fEtaWidth;                           ///<
	TH1D* fPhiWidth;                           ///<


  // Debug Tree for ME Scaling
  //
  TTree * fMEScaleTree = 0;                   ///< Tree for investigating ME normalization
  Double_t fMEValueAtZero   = 0;
  Double_t fMEValueAtPi     = 0;
  Double_t fMEScaleValue    = 0;
  Double_t fMEMean1DValue   = 0;
  Double_t fMEValueAt1DMax  = 0;
  Int_t fMEIndexOfNormRegion = 0; ///< Bin of normalization region
  // From 2D Normalized Plot
  Double_t fMEValueAtOrigin = 0; 
  Double_t fMEMean2DValue   = 0;
  Double_t fMEValueAt2DMax  = 0;

	TCanvas *fCanBinCheck;                     ///<
	TCanvas *fCanNormCheck;                    ///<
	TCanvas *fMEPlots2DGamma[10];              ///< [kNvertBins]
	TCanvas *fMEPlots1DGamma[10];              ///< [kNvertBins]
	TCanvas *fMEPlots2DXi[10];                 ///< [kNvertBins]
	TCanvas *fMEPlots2DZt[10];                 ///< [kNvertBins]
//	TCanvas *fRaw_Plots_2D_GammaV1[20];        ///<
//	TCanvas *fPlots_2D_Sum;                    ///< Canvas with corrected 2D distribution (summed over z-vertexes)
//	TCanvas *fPlots_2D_Corr_alt;               ///< Canvas with alternative way of correcting SE histogram
//	TCanvas *fPlots_2D_MEcompareV1;            ///< Canvas comparing the two ways of correcting histograms (errors)
//	TCanvas *fPlots_2D_MEcompareV2;            ///< Canvas comparing the two ways of correcting histograms (yield)

//	TCanvas *fCorrected_Plots_2D_Gamma[20];    ///<
	TCanvas *fRaw_Plots_2D[10];                ///< [kNvertBins]
	TCanvas *fRaw_Plots_2D_alt2[10];                ///< [kNvertBins]
	TCanvas *fRaw_Plots_2D_alt3[10];                ///< [kNvertBins]
	TCanvas *fPlots_2D_ME_alt1;                ///< Canvas with alternative way of correcting SE histogram
	TCanvas *fPlots_2D_ME_alt2[10];            ///< Canvas with alternative way of correcting SE histogram
	TCanvas *fPlots_2D_ME_alt3;                ///< Canvas with alternative way of correcting SE histogram
	TCanvas *fPlots_2D_Corr[20];               ///< This is the actual corrected spectrum, still in several zvertex bins
	TCanvas *fPlots_2D_Corr_alt1;              ///< Canvas with alternative way of correcting SE histogram
	TCanvas *fPlots_2D_Corr_alt2[10];          ///< Canvas with alternative way of correcting SE histogram
	TCanvas *fPlots_2D_Corr_alt3[10];          ///< Canvas with alternative way of correcting SE histogram
	TCanvas *fPlots_1D_MEcompare;              ///< Canvas comparing the two ways of correcting histograms (errors)
	TCanvas *fPlots_1D_MEcompareEta;           ///< Canvas comparing the two ways of correcting histograms (errors)
	TCanvas *fPlots_1D_MEcompareRatio;         ///<
	TCanvas *fPlots_1D_MEcompareErrors;        ///<

	TCanvas *fPlots_2D_CorrSum;                ///< Canvas with corrected 2D distribution (summed over z-vertexes)
	TCanvas *fPlots_2D_CorrSum_alt1;           ///< Canvas with corrected 2D distribution (summed over z-vertexes)
	TCanvas *fPlots_2D_CorrSum_alt2;           ///< Canvas with corrected 2D distribution (summed over z-vertexes)
	TCanvas *fPlots_2D_CorrSum_alt3;           ///< Canvas with corrected 2D distribution (summed over z-vertexes)
	TCanvas *fCanCorr1D;                       ///< Projection on the Delta Eta axis to estimate the with of the NS peak
	TCanvas *fCanCorr1D_Sub;                   ///< Subtracted Projection on the Delta Eta axis to estimate the with of the NS peak
	TCanvas *fCanvWidth;                       ///<
	TCanvas *fCanProj;                         ///<
	TCanvas *fCanProj2;                        ///<
	TCanvas *fCanProj3;                        ///<
	TCanvas *fCanProjFull;                     ///<
	TCanvas *fMEplotsG;    				      ///<
	TCanvas *fSEplotsG;    				      ///<

	const double fDetaLimit = 1.0;             ///< Range used to find min/max for 2D plots.
//	const bool fPerformance = true;           ///< Switch on in case we get performance approval again

private:
	PlotGHcorrelation2           (const PlotGHcorrelation2&);            // not implemented
	PlotGHcorrelation2 &operator=(const PlotGHcorrelation2&); // not implemented

	ClassDef(PlotGHcorrelation2, 1) // Class to analyse gamma hadron correlations
};
#endif



//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - Here are functions that get the historgams from the files - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/*TH1D* Get1DHistoFromFile(TFile* RootFile,TString inputListName,TString SubListName,TString Name)
{
	TList* IntermediatList;
	TList* FinalList;

	if(SubListName=="")
	{
		FinalList    =(TList*)RootFile->Get(inputListName);
	}
	else
	{
		IntermediatList=(TList*)RootFile       ->Get(inputListName);
		FinalList      =(TList*)IntermediatList->FindObject(SubListName);
	}
	TH1D* Histo   =(TH1D*)FinalList->FindObject(Name);

	return Histo;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TH2D* Get2DHistoFromFile(TFile* RootFile,TString inputListName,TString SubListName,TString Name)
{
	TList* IntermediatList;
	TList* FinalList;

	if(SubListName=="")
	{
		FinalList    =(TList*)RootFile->Get(inputListName);
	}
	else
	{
		IntermediatList=(TList*)RootFile       ->Get(inputListName);
		FinalList      =(TList*)IntermediatList->FindObject(SubListName);
	}
	TH2D* Histo   =(TH2D*)FinalList->FindObject(Name);

	return Histo;
}

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
void SetTH1Histo(TH1 *Histo,TString Xtitel,TString Ytitel,Bool_t big=0)
{
	Histo->SetStats(0);
	Histo->SetTitle("");
	if(big==0)	Histo->GetYaxis()->SetTitleOffset(1.4);
	if(big==0)	Histo->GetXaxis()->SetTitleOffset(1.4);
	if(big==1)	Histo->GetYaxis()->SetTitleOffset(0.8);
	if(big==1)	Histo->GetXaxis()->SetTitleOffset(1.0);
	//if(big==1)	Histo->GetYaxis()->SetLabelOffset(0.015);
	//if(big==1)	Histo->GetXaxis()->SetLabelOffset(0.015);
	if(big==0)	Histo->GetXaxis()->SetLabelSize(0.05);
	if(big==0)  Histo->GetYaxis()->SetLabelSize(0.05);
	if(big==1)	Histo->GetXaxis()->SetLabelSize(0.07);
	if(big==1)  Histo->GetYaxis()->SetLabelSize(0.07);
	if(big==0)	Histo->GetXaxis()->SetTitleSize(0.045);
	if(big==0)	Histo->GetYaxis()->SetTitleSize(0.045);
	if(big==1)  Histo->GetXaxis()->SetTitleSize(0.08);
	if(big==1)	Histo->GetYaxis()->SetTitleSize(0.08);
	Histo->GetXaxis()->CenterTitle();
	Histo->GetYaxis()->CenterTitle();
	Histo->GetXaxis()->SetNdivisions(505);
	Histo->GetYaxis()->SetNdivisions(505);
	//make nice font
    Histo->GetXaxis()->SetLabelFont(42);
    Histo->GetYaxis()->SetLabelFont(42);
    Histo->GetXaxis()->SetTitleFont(42);
    Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitel!="")Histo->GetXaxis()->SetTitle(Xtitel);
	if(Ytitel!="")Histo->GetYaxis()->SetTitle(Ytitel);

	Histo->SetLineColor(1);
	Histo->SetMarkerColor(1);
	Histo->SetMarkerStyle(20);
	Histo->SetMarkerSize(0.5);

}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Double_t ScaleMEbackground(TH2D* Histo,Double_t lowRange,Double_t highRange,TCanvas* Can,Int_t CanvasPad)
{
//	TF1* LinFit = new TF1("LinFit","","[0]",-50,50,1);
	TF1* LinFit = new TF1("pol0","pol0",-50,50);

	Can->cd(CanvasPad+1);
	//..project to the x-axis. But only around y=0.
	TString ProjectionName;
	ProjectionName= Histo->GetName();
	ProjectionName+="_projX_range";
	TH1D *PprojX=Histo->ProjectionX((char)ProjectionName,Histo->GetYaxis()->FindBin(lowRange),Histo->GetYaxis()->FindBin(highRange));
	SetTH1Histo(PprojX,"#Delta #phi","dN^{#gamma-h}/dN^{#gamma}",1);
	PprojX->DrawCopy("E");
	PprojX->Fit("pol0","Q","",-50,50);//Q = quiet mode, no printout

	//cout<<"param: "<<LinFit->GetParameter(0)<<endl;
	//..determine/etabin yield (count bins over which it was integrated)
	Int_t nBins= Histo->GetYaxis()->FindBin(highRange)-Histo->GetYaxis()->FindBin(lowRange);
	Double_t Scalef = LinFit->GetParameter(0)/(1.0*nBins);
	TString TopLegendText=Form("fit value: %0.2f/%i",LinFit->GetParameter(0),nBins);
	plotTopLegend((char)TopLegendText,0.2,0.25,0.07,kGreen-2);
	Histo->Scale(1/Scalef);
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0)
{
	// coordinates in NDC!
	// plots the string label in position x and y in NDC coordinates
	// size is the text size
	// color is the text color

	if(x<0||y<0)
	{   // defaults
		x=gPad->GetLeftMargin()*1.15;
		y=(1-gPad->GetTopMargin())*1.04;
	}
	TLatex* text=new TLatex(x,y,label);
	text->SetTextSize(size);
	text->SetNDC();
	text->SetTextColor(color);
	text->SetTextAngle(angle);
	text->Draw();
	return text;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void SetTH2Histo(TH2 *Histo,TString Xtitel,TString Ytitel)
{
	Histo->SetStats(0);
	Histo->SetTitle("");
	Histo->GetYaxis()->SetTitleOffset(1.7);
	Histo->GetXaxis()->SetTitleOffset(1.7);
	Histo->GetXaxis()->SetLabelSize(0.05);
	Histo->GetYaxis()->SetLabelSize(0.05);
	Histo->GetXaxis()->SetTitleSize(0.045);
	Histo->GetYaxis()->SetTitleSize(0.045);
	Histo->GetXaxis()->CenterTitle();
	Histo->GetYaxis()->CenterTitle();
	Histo->GetXaxis()->SetNdivisions(505);
	Histo->GetYaxis()->SetNdivisions(505);
	//make nice font
    Histo->GetXaxis()->SetLabelFont(42);
    Histo->GetYaxis()->SetLabelFont(42);
    Histo->GetXaxis()->SetTitleFont(42);
    Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitel!="")Histo->GetXaxis()->SetTitle(Xtitel);
	if(Ytitel!="")Histo->GetYaxis()->SetTitle(Ytitel);
	Histo->SetLineColorAlpha(kBlue+2,0.095);
//	Histo->GetYaxis()->SetRangeUser(-1.3,1.3);

}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TH1D* FitEtaSides(TH2* Histo,Double_t Array[],Int_t Sigmas,TCanvas* Can,Int_t CanvasPad)
{
	TString TopLegendText;
    TString ProjectionName;
    TF1* BackgroundFunction = new TF1("FlowFunction",FlowFunction,-100,300,4);
    Int_t vN=3;
    TF1 *allFit             = allFitFuncVn("JoelsFlowFunction",vN,-100,300);

	//..select the eta range in which you want to project your signal
	Int_t SignalCenter    = Histo->GetYaxis()->FindBin(0);
	Int_t SignalEdgeLow   = Histo->GetYaxis()->FindBin(0-Array[0]*Sigmas);
	Int_t SignalEdgeHigh  = Histo->GetYaxis()->FindBin(0+Array[0]*Sigmas);
	Int_t SignalEdgeLowT  = Histo->GetYaxis()->FindBin(0-Array[0]*(Sigmas+1));//..increase sigma ragne by 1 to check signal left over
	Int_t SignalEdgeHighT = Histo->GetYaxis()->FindBin(0+Array[0]*(Sigmas+1));//..increase sigma ragne by 1 to check signal left over
	Int_t LowestBin       = Histo->GetYaxis()->FindBin(-1.4);
	Int_t HighestBin      = Histo->GetYaxis()->FindBin(1.4);

	//..check that the mean+-sigma is not larger or smaller than the histogram range
	if(SignalEdgeLow<LowestBin || SignalEdgeHigh>HighestBin)
	{
		cout<<"Error: Problem detected!"<<endl;
		cout<<"In Histo: "<<Histo->GetName()<<endl;
		cout<<"Signal range is reaching outside the histogram boundaries - please correct"<<endl;
		cout<<"bins lowes"<<LowestBin<<", edge "<<SignalEdgeLow<<", center"<<SignalCenter <<", uppedge"<< SignalEdgeHigh<<", highestbin "<< HighestBin<<endl;
		SignalEdgeLow=LowestBin;
		SignalEdgeHigh=HighestBin;
	}
	//..determine a scale factor to fit NS to away side
	//..this is necessary since the flow function is only fit
	//..to the side band but needs to be scaled to the whole eta range
	Double_t scaleFactorSBtoNS;// = ScaleSBtoNS(Histo,LowestBin,SignalEdgeLow,SignalEdgeHigh,HighestBin);

	Int_t nBinsSB = (SignalEdgeLow-LowestBin)+(HighestBin-SignalEdgeHigh);
	Int_t nBinsNS = SignalEdgeHigh-SignalEdgeLow;
	scaleFactorSBtoNS=1.0*nBinsSB/nBinsNS;
    cout<<"^-^"<<"bins SB: "<<nBinsSB<<", bins NS: "<<nBinsNS<<", scaling: "<<scaleFactorSBtoNS<<endl;

    ProjectionName=Histo->GetName();
    ProjectionName+="PprojXSig";
	TH1D *PprojXSig  =Histo->ProjectionX(ProjectionName,SignalEdgeLow,SignalEdgeHigh);
    ProjectionName=Histo->GetName();
    ProjectionName+="PprojXSide1";
	TH1D *PprojXSide1=Histo->ProjectionX(ProjectionName,LowestBin,SignalEdgeLow);
    ProjectionName=Histo->GetName();
    ProjectionName+="PprojXSide2";
    TH1D *PprojXSide2=Histo->ProjectionX(ProjectionName,SignalEdgeHigh,HighestBin);
    //.. Test part
    ProjectionName+="PprojXSide1_Test";
	TH1D *PprojXSide1_T=Histo->ProjectionX(ProjectionName,LowestBin,SignalEdgeLowT);
    ProjectionName=Histo->GetName();
    ProjectionName+="PprojXSide2_Test";
    TH1D *PprojXSide2_T=Histo->ProjectionX(ProjectionName,SignalEdgeHighT,HighestBin);


	//..perform now a projection of the 2Dhistogam in a given range x-sigmas outside the jet region
	Can->cd(CanvasPad*2+1);
	//..project to plot the delta phi distribution in a certain Delta Eta intervall
	SetTH1Histo(PprojXSig,"","",1);
	ZoomYRange(PprojXSig);
	PprojXSig->DrawCopy("E");
	TopLegendText=Form("Projection in range:");
	plotTopLegend((char)TopLegendText,0.53,0.85,0.07);
	TopLegendText=Form("%0.2f < #eta < %0.2f",0-Array[0]*Sigmas,0+Array[0]*Sigmas);
	plotTopLegend((char)TopLegendText,0.53,0.77,0.07);

	//..large eta
	Can->cd(CanvasPad*2+2);
	PprojXSide1  ->Add(PprojXSide2);
	PprojXSide1_T->Add(PprojXSide2_T);
	SetTH1Histo(PprojXSide1,"","",1);
	ZoomYRange(PprojXSide1);
	PprojXSide1->DrawCopy("E");
	//..scale the test histo to the awway side yield
//	Double_t Integral1=PprojXSide1  ->Integral(PprojXSide1->FindBin(160),PprojXSide1->FindBin(200));
//	Double_t Integral2=PprojXSide1_T->Integral(PprojXSide1_T->FindBin(160),PprojXSide1_T->FindBin(200));
//	PprojXSide1_T->SetLineColor(6);
//	PprojXSide1_T->SetMarkerColor(6);
//	Double_t ratio=(Double_t)Integral1/(Double_t)Integral2;
//	PprojXSide1_T->Scale(ratio);
//	Double_t Integral3=PprojXSide1_T->Integral(PprojXSide1_T->FindBin(160),PprojXSide1_T->FindBin(200));
//	PprojXSide1_T->DrawCopy("same E");


    //fit with the flow function
	BackgroundFunction->SetParNames("B", "comb v1", "comb v2", "comb v3");

	for(Int_t g = 0; g < 4; g++)
	{
		BackgroundFunction->ReleaseParameter(g);
		BackgroundFunction->SetParameter(g,0);
		BackgroundFunction->SetParError(g,0.0);
	}
	//flat line
//	BackgroundFunction->FixParameter(1,0.0);
//	BackgroundFunction->FixParameter(2,0.0);
//	BackgroundFunction->FixParameter(3,0.0);

	BackgroundFunction->SetParLimits(1,0.0,100); //..do not allow negative values
	BackgroundFunction->SetParLimits(2,0.0,1);   //..do not allow negative values, v2<v1 (v2 strength fraction of v1 strength)
	BackgroundFunction->SetParLimits(3,0.0,1);   //..do not allow negative values, v3<v2 (v3 strength fraction of v2 strength)

	// arbirtary - set for you, I'm just copying and generalizing stuff here
	//double par_V10[15] = {450, 0.8e-1, 5e-2, 1e-4, 1e-2, 1e-2,1e-4, 1e-3, 1e-3, 1e-5, 1e-3, 1e-3, 1e-5, 1e-3, 1e-3};
	double par_V10[15] = {450, 0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0};
	//.. Fit Range
	Double_t limits = 90; //..pi
	//limits=1.25*180.0/TMath::Pi();  //..joels paper range 1 (71.6deg too small)
	limits=1.57*180.0/TMath::Pi();  //..joels paper range 2 (90deg)

	SetJoelsParLimits(allFit,PprojXSide1,par_V10);

	TString funcName    = BackgroundFunction->GetName();
	TString JoelfuncName= allFit->GetName();
	TFitResultPtr r  = PprojXSide1->Fit(funcName,"","",-limits,limits);//Q = quiet mode, no printout
	TFitResultPtr r2 = PprojXSide1->Fit(JoelfuncName,"S0","",-limits,limits);

	allFit->SetLineColor(kRed-9);
	allFit->SetRange(-100,300);
	allFit->DrawCopy("same");
	BackgroundFunction->SetLineStyle(9);//
	BackgroundFunction->SetLineColor(kCyan-5);//
	BackgroundFunction->SetRange(-100,300);
	BackgroundFunction->DrawCopy("same");

	TopLegendText=Form("Projection in range:");
	plotTopLegend((char)TopLegendText,0.2,0.85,0.07);
	TopLegendText=Form("%0.2f < #eta < %0.2f",-1.5,0-Array[0]*Sigmas);
	plotTopLegend((char)TopLegendText,0.2,0.77,0.07);
	TopLegendText=Form(" %0.2f < #eta < %0.2f",0+Array[0]*Sigmas,1.5);
	plotTopLegend((char)TopLegendText,0.2,0.69,0.07);

	//..Return the background subtracted correlation spectrum
	ProjectionName = PprojXSide1->GetName();
	ProjectionName+="_backgroundSubtracted";
	BackgroundFunction->SetLineColor(kWhite);
	TH1D* BackgroundSubtraction = (TH1D*)PprojXSide1->Clone(ProjectionName);
	BackgroundSubtraction->Add(BackgroundFunction,-1);
//	BackgroundSubtraction->GetFunction(funcName)->Delete();//..to not draw the fit function at any point later

    return BackgroundSubtraction;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void DetermineWidths(TH2* Histo,Double_t Array[],TCanvas* Can,Int_t CanvasPad)
{
	//cout<<"inside of DetermineWidths()"<<endl;
	TF1* GaussFunc  = new TF1("GaussFunc",PolyTwoGaussFitFunc,-100,300,11);
	TF1* GaussFunc1 = new TF1("GaussFunc1",PolyTwoGaussFitFunc,-100,300,11);
	TF1* GaussFunc2 = new TF1("GaussFunc2",PolyTwoGaussFitFunc,-100,300,11);
	TString ProjectionName;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Can->cd(CanvasPad*2+1);
	//project to the y-axis. But only the near side.
	ProjectionName= Histo->GetName();
	ProjectionName+="_projY";
	TH1D *projY=Histo->ProjectionY((char)ProjectionName,0,Histo->GetYaxis()->FindBin(180));
	SetTH1Histo(projY,"","dN^{#gamma-h}/dN^{#gamma}",1);
	ZoomYRange(projY);
	projY->DrawCopy("E");
	//..fit, draw and save widths
	FitGaussAndDraw(projY,GaussFunc,GaussFunc1,GaussFunc2,0,Array);

	Can->cd(CanvasPad*2+2);
	TH1D *projX=Histo->ProjectionX();
	SetTH1Histo(projX,"","dN^{#gamma-h}/dN^{#gamma}",1);
	ZoomYRange(projX);
	projX->DrawCopy("E");
	FitGaussAndDraw(projX,GaussFunc,GaussFunc1,GaussFunc2,1,Array);

}

void ZoomYRange(TH1 *Histo,Double_t border=0.1)
{
	//Double_t min=Histo->GetBinContent(Histo->GetMinimumBin());
	Double_t min=Histo->GetMinimum(0);
	Double_t max=Histo->GetBinContent(Histo->GetMaximumBin());
	//cout<<"zoom into range: "<<min<<"-"<<max<<", "<<min*(1-border)<<"-"<<max*(1+border)<<endl;
	Histo->GetYaxis()->SetRangeUser(min*(1-border),max*(1+2*border));
}


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void FitGaussAndDraw(TH1D* Hist,TF1* Func1, TF1* Func2,TF1* Func3,Bool_t EtaPhi,Double_t Array[])
{
	Double_t Width;

	for(Int_t g = 0; g < 11; g++)
	{
		Func1->ReleaseParameter(g);
		Func1->SetParameter(g,0);
		Func1->SetParError(g,0.0);
	}
	Func1->SetParameter(1,0); //..start value for mean1
	Func1->SetParameter(4,0); //..start value for mean2


	if(EtaPhi==0)//..Eta case
	{
		//- - - - - - - - - - - - - - - -
		//..Flat background
		//..for eta -> estimate the background level on eta=-0.7&+0.7
		Double_t backgroundLevel=Hist->GetBinContent(Hist->FindBin(-0.7));
	    backgroundLevel+=Hist->GetBinContent(Hist->FindBin(-0.7));
	    backgroundLevel*=0.5;
	    //cout<<"Background level="<<backgroundLevel<<endl;
	    Func1->SetParameter(6,backgroundLevel);
	    Func1->SetParLimits(6,backgroundLevel*0.9,backgroundLevel*1.1);  //..allow a variation of +-10%

	    //- - - - - - - - - - - - - - - -
		//..big, narrow gaussian
	    //..etimate amplitude by 0 heigth and background level
		Double_t amplEst=Hist->GetBinContent(Hist->FindBin(0));
		amplEst-=backgroundLevel;
		Func1->SetParameter(0,amplEst);    //..amplitude
		Func1->SetParameter(2,0.05);       //..width
		Func1->SetParLimits(0,amplEst*0.9,amplEst*1.1);
		Func1->SetParLimits(1,-0.1,0.1);  //..mean limits
		Func1->SetParLimits(2,0.05,0.5);  //..width limits

	    //- - - - - - - - - - - - - - - -
		//..small, wide gaussian
		Func1->SetParameter(3,0.05);       //..amplitude
		Func1->SetParameter(5,1.1);        //..width
		Func1->SetParLimits(3,0.05,0.5);   //..amplitude limits 5-50% of the main peak (ampl2 = param0*param3)
		Func1->SetParLimits(4,-0.1,0.1);   //..mean limits
		Func1->SetParLimits(5,1.05,3.0);   //..width limits 105%-300% of the big one (width2 = param2*param5)

	    //- - - - - - - - - - - - - - - -
		//..plot range for eta projection
	    Func1->SetRange(-1.5,1.5);
	    Func2->SetRange(-1.5,1.5);
	    Func3->SetRange(-1.5,1.5);
	}

	if(EtaPhi==1)//..phi case
	{
	    //- - - - - - - - - - - - - - - -
		//..Flat background
		//..for phi -> estimate the background level on phi=-85&+85
		Double_t backgroundLevel=Hist->GetBinContent(Hist->FindBin(-85));
	    backgroundLevel+=Hist->GetBinContent(Hist->FindBin(+85));
	    backgroundLevel*=0.5;
	    //cout<<"Background level="<<backgroundLevel<<endl;
	    Func1->SetParameter(6,backgroundLevel);
	    Func1->SetParLimits(6,backgroundLevel*0.9,backgroundLevel*1.1);  //..allow a variation of +-10%

	    //- - - - - - - - - - - - - - - -
		//..big, narrow gaussian
	    //..etimate amplitude by 0 heigth and background level
		Double_t amplEst=Hist->GetBinContent(Hist->FindBin(0));
		amplEst-=backgroundLevel;
		Func1->SetParameter(0,amplEst);   //..amplitude
		Func1->SetParameter(2,15);        //..width
		Func1->SetParLimits(0,amplEst*0.9,amplEst*1.1);
		Func1->SetParLimits(1,-0.1,0.1);  //..mean limits
		Func1->SetParLimits(2,5,30);      //..width limits

	    //- - - - - - - - - - - - - - - -
		//..small, wide gaussian
		Func1->SetParameter(3,0.05);       //..amplitude
		Func1->SetParameter(5,1.1);        //..width
		Func1->SetParLimits(3,0.05,0.5);   //..amplitude limits 5-50% of the main peak (ampl2 = param0*param3)
		Func1->SetParLimits(4,-0.1,0.1);   //..mean limits
		Func1->SetParLimits(5,1.05,3.0);   //..width limits 105%-300% of the big one (width2 = param2*param5)

	    //- - - - - - - - - - - - - - - -
		//..plot range for eta projection
	    Func1->SetRange(-90,90);
	    Func2->SetRange(-90,90);
	    Func3->SetRange(-90,90);
	}

    //.. CAREFUL YOU CAN ALSO SET THE SECOND GAUSS TO 0.
	Func1->FixParameter(3,0);

	Func1->SetLineColor(15);
	Func1->SetLineStyle(2);

	TString Name= Func1->GetName();
	if(EtaPhi==0)Hist->Fit(Name,"Q","",-1,1);//Q = quiet mode, no printout
	if(EtaPhi==1)Hist->Fit(Name,"Q","",-70,70);//Q = quiet mode, no printout

	//..width that is used to define an eta range not contaminated by the near side peak
	//..you can define the width in multiple ways.
	//..THINK ABOUT THAT!
	//Width   =GaussFunc->GetParameter(5)*GaussFunc->GetParameter(2);//bigger width in delta phi
	Width   =GaussFunc->GetParameter(2); //bigger width in delta phi
	Double_t ErrWidth=GaussFunc->GetParError(2);
	if(EtaPhi==0)//..eta case
	{
		Array[0]=Width;
		Array[1]=ErrWidth;
	}
	if(EtaPhi==1)//..phi case
	{
		Array[2]=Width;
		Array[3]=ErrWidth;
	}
	//  cout<<"Width gauss2: "<<GaussFunc->GetParameter(5)<<" times of width 1"<<endl;

	for(Int_t g = 0; g < 11; g++)
	{
		Func2->SetParameter(g,Func1->GetParameter(g));
		Func3->SetParameter(g,Func1->GetParameter(g));
	}
	//..small, wide gaussian
	//..due to the fact that param 0 and param 3 are proportional
	//..we do a little hack here. Setting param0 to 0 is neseccary
	//..to see only the small wiede gaussian. If we set param0 to 0
	//..however, param3 will become 0 by defualt. We can however
	//..set param0 to a negligibly small value x and multiply param3
	//..by the inverse of x. (normally param3 is in the range 0-1, but we omit this for this specific case)
	Double_t Shrinkage=0.00001;
	//Func2
	Func2->SetParameter(0,Shrinkage);
	Func2->SetParameter(3,1.0*Func1->GetParameter(3)*Func1->GetParameter(0)/Shrinkage);
	Func2->SetLineColor(kPink-9);
	if(GaussFunc->GetParameter(3)!=0)Func2 ->DrawCopy("same"); //..only when the small-broad gaussian is not set to 0

    //..big, narrow gaussian (green)
	Func3->SetParameter(3,0);
	Func3->SetLineColor(kGreen-2);
	Func3 ->DrawCopy("same");

	//..flat background
	Func1->SetParameter(0,0);
	Func1->SetParameter(3,0);
	Func1 ->DrawCopy("same");

	//..plot 3 sigma range (enough to split off near side peak)
	PlotVerLine3(Width*3,Hist,0.8,17);
	PlotVerLine3(-(Width*3),Hist,0.8,17);

	//..Draw legend
	TString TopLegendText;
	Double_t x_Pos;
	if(EtaPhi==0)x_Pos=0.2;//Eta case
	if(EtaPhi==1)x_Pos=0.6;//Phi case
	TopLegendText=Form("#mu_{1}: %0.3f",Func1->GetParameter(1));
	plotTopLegend((char)TopLegendText,x_Pos,0.85,0.07,kGreen-2);
	TopLegendText=Form("#sigma_{1}: %0.2f",Func1->GetParameter(2));
	plotTopLegend((char)TopLegendText,x_Pos,0.78,0.07,kGreen-2);

	if(GaussFunc->GetParameter(3)!=0)
	{
		TopLegendText=Form("#mu_{2}: %0.3f",Func1->GetParameter(4));
		plotTopLegend((char)TopLegendText,x_Pos,0.71,0.07,kPink-9);
		TopLegendText=Form("#sigma_{2}: %0.2f",Func1->GetParameter(5)*Func1->GetParameter(2));
		plotTopLegend((char)TopLegendText,x_Pos,0.64,0.07,kPink-9);
	}
}
Double_t PolyTwoGaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5, par6, par7, par8, par9, par10,CommomMean;
    par0  = par[0]; //amplitude gauss 1
    par1  = par[1]; //mean gauss 1
    par2  = par[2]; //width gauss 1
    //the second gaus is smaller in amplitude and larger in width
    par3  = par[3]*par[0]; //amplitude gauss 2 (parameter3 ranges from 0-1)
    par4  = par[4]; //mean gauss 2
    par5  = par[5]*par[2]; //width gauss 2 (parameter5 is larger than 1)
    par6  = par[6]; //a
    par7  = par[7]; //b x^1
    par8  = par[8]; //c x^2
    par9  = par[9]; //d x^3
    par10 = par[10];//e x^4
    x = x_val[0];

    //Do that so that the mean of the two gaussians are the same
    CommomMean=(par1+par4)*0.5;

 //   cout<<"current p0: "<<par0<<", current p3: "<<par3<<endl;

    y = par0*(TMath::Gaus(x,CommomMean,par2,0))+par3*(TMath::Gaus(x,CommomMean,par5,0))+par6;//+(par6+par7*x+par8*x*x+par9*x*x*x+par10*x*x*x*x);
    return y;
}
void PlotVerLine3(Double_t x_val,TH1* Histo, Double_t y_fac, Int_t Line_Col)
{
	Double_t min=Histo->GetMinimum(0);
	Double_t max=Histo->GetBinContent(Histo->GetMaximumBin());
	Double_t maxlocal=Histo->GetBinContent(Histo->FindBin(x_val));
    if(y_fac*max<maxlocal*1.1)max=maxlocal*1.1*(1/y_fac);

    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x_val);
    Zero_line -> SetX2(x_val);
    Zero_line -> SetY1(min);
    Zero_line -> SetY2(y_fac*max);
    //cout << "x_val = " << x_val << ", Bin = " << Histo->FindBin(x_val) << ", Y2 = " << Histo->GetBinContent(Histo->FindBin(x_val)) << endl;
    Zero_line -> SetLineWidth(2);
    Zero_line -> SetLineStyle(2);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TF1 * allFitFuncVn(char *name,Int_t VnTerms,Double_t allPhiMin,Double_t allPhiMax)
{
  if(VnTerms == 10) TF1* f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 16);
  if(VnTerms == 7)  TF1* f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 11);
  if(VnTerms == 5)  TF1* f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 8);
  if(VnTerms == 4)  TF1* f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 7);
  if(VnTerms == 3)  TF1* f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 5);

  return f1;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Double_t JoelsVnFunctionValue(const double * x, const double * p)
{
  Int_t VnTerms=3; //!!!!!!!!!!!!! do not hard code
  double B  = p[0];
  double v1       = p[1];
  double v2jet    = p[2];
  double v2assoc  = p[3];
  double v3 = p[4];
  double v4jet    = p[5];
  double v4assoc  = p[6];
  double v5 = p[7];
  double v6jet    = p[8];
  double v6assoc  = p[9];
  double v7 = p[10];
  double v8jet    = p[11];
  double v8assoc  = p[12];
  double v9 = p[13];
  double v10jet   = p[14];
  double v10assoc = p[15];

  //- - - - - - - - - - - - - - - - -
  //..set not needed terms to 0
  if(VnTerms<=7)
  {
	  v8jet    = 0.;
	  v8assoc  = 0.;
	  v9       = 0.;
	  v10jet   = 0.;
	  v10assoc = 0.;
  }
  if(VnTerms<=5)
  {
	  v6jet   = 0.;
	  v6assoc = 0.;
	  v7      = 0.;
  }
  if(VnTerms<=4)
  {
	  v5 = 0.;
  }
  if(VnTerms<=3)
  {
	  v4jet   = 0.;
	  v4assoc = 0.;
  }
  //- - - - - - - - - - - - - - - - -

  double result;
  double phi = x[0]*TMath::Pi()/180.0;  //transform from deg to rad


  // changed the following line, because normalizing by triggers not # of events
  //..why is there no v1?   -> 2.0*v1*TMath::Cos(1*phi);
  Float_t Part1  = 60.*TMath::Pi();
  Float_t Part2  = 2.0*v2jet*v2assoc*TMath::Cos(2.0*phi);
  Float_t Part4  = 2.0*v4jet*v4assoc*TMath::Cos(4.0*phi);
  Float_t Part6  = 2.0*v6jet*v6assoc*TMath::Cos(6.0*phi);
  Float_t Part8  = 2.0*v8jet*v8assoc*TMath::Cos(8.0*phi);
  Float_t Part10 = 2.0*v10jet*v10assoc*TMath::Cos(10.*phi);

  Float_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Float_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  Float_t Part7 = 2.0*v7*TMath::Cos(7*phi);
  Float_t Part9 = 2.0*v9*TMath::Cos(9*phi);

  // final func
  //Joel original result = (B/(60.*TMath::Pi()))*Part1*(1 + Part2 + Part3 + Part4 + Part5 + Part6 + Part7 + Part8 + Part9 + Part10);
  //..Eliane adding v1 into the function
  Part1  = 2.0*v1*TMath::Cos(phi);
  result = B*(1 + Part1 + Part2 + Part3 + Part4 + Part5 + Part6 + Part7 + Part8 + Part9 + Part10);

  return result;
}
//-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
Double_t FlowFunction(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3;
    par0  = par[0]; //Background level
    par1  = par[1]; //
    par2  = par[2]; //0-1
    par3  = par[3]; //0-1

    x = x_val[0]*TMath::Pi()/180.0;  //transform from deg to rad

    y = par0*(1+2*par1*TMath::Cos(x)+2*par1*par2*TMath::Cos(2*x)+2*par2*par3*TMath::Cos(3*x));
    return y;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void SetJoelsParLimits(TF1 *func, TH1 *histoToFit,Double_t par_V10[])
{
  par_V10[0] = (histoToFit->GetMinimum()+histoToFit->GetMaximum());
  func->SetParameters(&par_V10[0]);
  func->SetParNames("B","v1", "v2jet", "v2assoc", "v3", "v4jet", "v4assoc", "v5", "v6jet", "v6assoc", "v7");
  func->SetParName(11, "v8jet");
  func->SetParName(12, "v8assoc");
  func->SetParName(13, "v9");
  func->SetParName(14, "v10jet");
  func->SetParName(15, "v10assoc");

  // - completely arbirtary - set whats right for you
  func->SetParLimits(0,1e-6,histoToFit->GetMaximum());
  func->SetParLimits(1,1e-6,10);
  func->SetParLimits(2,1e-6,10);
  func->SetParLimits(3,1e-6,10);
  func->SetParLimits(4,1e-6,10);
  func->SetParLimits(5,1e-6,10);
  func->SetParLimits(6,1e-6,10);
  func->SetParLimits(7,1e-4,10);
  func->SetParLimits(8,1e-4,10);
  func->SetParLimits(9,1e-6, 10);
  func->SetParLimits(10,1e-4,10);
  func->SetParLimits(11,1e-4,10);
  func->SetParLimits(12,1e-6,10);
  func->SetParLimits(13,1e-4,10);
  func->SetParLimits(14,1e-4,10);
  func->SetParLimits(15,1e-4,10);
}
 */
/*



   |
   |
   |
  \  /
   \/
 to be updated


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - Here are functions that make the plotting of figures nicer- - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Draw_Legend(TH1* File1Histo,TString File1HistoName,TH1* File2Histo,TString File2HistoName)
{
	TLegend *leg1=new TLegend(0.35,0.65,0.80,0.88);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->AddEntry(File1Histo,File1HistoName,"L");
	leg1->AddEntry(File2Histo,File2HistoName,"L");
	leg1->Draw();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Draw_Legend_V2(TH1* File1Histo,TString HistoName1,TH1* File2Histo,TString HistoName2,TH1* File3Histo,TString HistoName3,TH1* File4Histo,TString HistoName4)
{
	TLegend *leg1=new TLegend(0.45,0.70,0.70,0.88);
	leg1->SetBorderSize(0);
	leg1->SetTextSize(0.025);
	leg1->SetFillColor(0);
	leg1->AddEntry(File1Histo,HistoName1,"LP");
	leg1->AddEntry(File2Histo,HistoName2,"LP");
	leg1->AddEntry(File3Histo,HistoName3,"LP");
	leg1->AddEntry(File4Histo,HistoName4,"LP");
	leg1->Draw();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Draw_Boxes(TH1* Histo,TBox *Box[],const Int_t Bins,Color_t colorSceme[])
{
	Double_t min=Histo->GetMinimum(-5000);
	Double_t max=Histo->GetBinContent(Histo->GetMaximumBin());

	for(Int_t i=0;i<Bins/2;i++)
	{
		Box[i*2]->SetY1(max); //set maximum
		Box[i*2]->SetY2(min); //set minimum
		Box[i*2+1]->SetY1(max); //set maximum
		Box[i*2+1]->SetY2(min); //set minimum
		if(i==0)
		{
			Box[i*2]->SetY1(max*1.1); //set maximum
			Box[i*2+1]->SetY1(max*1.1); //set maximum
			Box[i*2]->SetY2(min*0.9); //set maximum
			Box[i*2+1]->SetY2(min*0.9); //set maximum
			Box[i*2]  ->SetFillColorAlpha(colorSceme[i],1);
			Box[i*2+1]->SetFillColorAlpha(colorSceme[i],1);
		}
		else
		{
			Box[i*2]  ->SetFillColorAlpha(colorSceme[i],0.9);
			Box[i*2+1]->SetFillColorAlpha(colorSceme[i],0.9);
		}
		Box[i*2]  ->Draw("");
		Box[i*2+1]->Draw("");
	}
}
//-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-




//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void PlotHorLine(Double_t x1_val, Double_t x2_val, Double_t y_val, Int_t Line_Col)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y_val);
    Zero_line -> SetY2(y_val);
    Zero_line -> SetLineWidth(2);
    Zero_line -> SetLineStyle(2);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void PlotVerLine(Double_t x_val, Double_t y_val_low, TH1* Histo, Double_t y_fac, Int_t Line_Col)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x_val);
    Zero_line -> SetX2(x_val);
    Zero_line -> SetY1(y_val_low);
    Zero_line -> SetY2(y_fac*Histo->GetBinContent(Histo->FindBin(x_val)));
    //cout << "x_val = " << x_val << ", Bin = " << Histo->FindBin(x_val) << ", Y2 = " << Histo->GetBinContent(Histo->FindBin(x_val)) << endl;
    Zero_line -> SetLineWidth(2);
    Zero_line -> SetLineStyle(1);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void PlotVerLine2(Double_t x_val, Double_t y_val_low, Double_t y_val_high, Int_t Line_Col)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x_val);
    Zero_line -> SetX2(x_val);
    Zero_line -> SetY1(y_val_low);
    Zero_line -> SetY2(y_val_high);
    //cout << "x_val = " << x_val << ", Bin = " << Histo->FindBin(x_val) << ", Y2 = " << Histo->GetBinContent(Histo->FindBin(x_val)) << endl;
    Zero_line -> SetLineWidth(2);
    Zero_line -> SetLineStyle(1);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - Here are functions that take the histograms and modify them - - - - - - - - - - - - -
//- - - - - - - In order to obtain new prepresenations with clearer information - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Double_t ScaleMEbackground2(TH2D* Histo,Double_t lowRange,Double_t highRange,TCanvas* Can,Int_t CanvasPad)
{
	TF1* LinFit = new TF1("pol0","pol0",-0.2,0.2);

	Can->cd(CanvasPad+1);
	//..project to the x-axis. But only around y=0.
	TString ProjectionName;
	ProjectionName= Histo->GetName();
	ProjectionName+="_projY_range";
	TH1D *PprojY=Histo->ProjectionY((char)ProjectionName,Histo->GetXaxis()->FindBin(lowRange),Histo->GetXaxis()->FindBin(highRange));
	SetTH1Histo(PprojY,"#Delta #eta","dN^{#gamma-h}/dN^{#gamma}",1);
	PprojY->GetXaxis()->SetRangeUser(-1,1);
	PprojY->DrawCopy("E");
	PprojY->Fit("pol0","Q","",-0.2,0.2);//Q = quiet mode, no printout

	cout<<"param: "<<LinFit->GetParameter(0)<<endl;
	//..determine/etabin yield (count bins over which it was integrated)
	Int_t nBins= Histo->GetXaxis()->FindBin(highRange)-Histo->GetXaxis()->FindBin(lowRange);
	TString TopLegendText=Form("fit value: %0.2f/%i",LinFit->GetParameter(0),nBins);
	plotTopLegend((char)TopLegendText,0.2,0.25,0.07,kGreen-2);
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Double_t ScaleSBtoNS(TH2D* Histo,)
{

}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - Here are the high end functions that take the histograms, - - - - - - - - - - - - - -
//- - - - - - - integrate them and plot the extracted information - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void ExtractInfo_Y_vs_Egamma(TH1D* Histo,Int_t Bin,TH1D* SummaryHisto[])
{
	Double_t BinWidth = Histo->GetBinWidth(1);
	for(Int_t i=0;i<4;i++)
	{
		//							        		   __
 		//integrate both sides (1,2) of the peak 1_/  \_2
		//Define limits below 180
		Int_t BIN_Start=Histo->FindBin(180-3*BinWidth*(i));
		Int_t BIN_end  =Histo->FindBin(180-3*BinWidth*(i-1));
		//summ bin contents:
		Double_t Yield=0;
        Double_t Error=0;

        //integrate the "full" range (currently 55¡-180¡)
        if(i==0)BIN_Start=Histo->FindBin(55);
        	if(i==0)BIN_end  =Histo->FindBin(180);

		for(Int_t j=BIN_Start;j<BIN_end;j++)
		{
		  Yield+= Histo->GetBinContent(j);
		  Error+= pow(Histo->GetBinError(j),2);
		  //cout<<"1) Adding bin nr. "<<j<<", with content: "<<Yield<<endl;
		}

		//Define limits above 180
		BIN_Start=Histo->FindBin(180+3*BinWidth*(i-1));
		BIN_end  =Histo->FindBin(180+3*BinWidth*(i));
		if(180+3*BinWidth*(i) >270) cout<<"ExtractInfo_Y_vs_Egamma::Error : End of histogram. Please correct!!"<<endl;

        //integrate the "full" range (currently 55¡-180¡ =125¡)
        if(i==0)BIN_Start=Histo->FindBin(180);
        	if(i==0)BIN_end  =Histo->FindBin(270);   //

        	for(Int_t j=BIN_Start;j<BIN_end;j++)
		{
		  Yield+= Histo->GetBinContent(j);
		  Error+= pow(Histo->GetBinError(j),2);
		  //cout<<"2) Adding bin nr. "<<j<<", with content: "<<Yield<<endl;
		}

        	//for large ranges integrate also the very beginning of the histogram!
        	if(i==0)
        	{
        		BIN_Start=Histo->FindBin(90);
        		BIN_end  =Histo->FindBin(90+35);
        		for(Int_t j=BIN_Start;j<BIN_end;j++)
        		{
        			Yield+= Histo->GetBinContent(j);
        			Error+= pow(Histo->GetBinError(j),2);
        		}
        	}

		if(Error<0)cout<<"ExtractInfo_Y_vs_Egamma:: strange error"<<endl;
		Error=sqrt(Error);
		SummaryHisto[i]->SetBinContent(Bin+1,Yield/SummaryHisto[i]->GetBinWidth(Bin+1));
		SummaryHisto[i]->SetBinError(Bin+1,Error/SummaryHisto[i]->GetBinWidth(Bin+1));
	}
}
 */

