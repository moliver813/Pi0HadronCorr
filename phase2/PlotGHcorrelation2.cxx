#if defined(__CINT__)
#define _SYS_TYPES_H_
#endif
#define DTR TMath::DegToRad()
// --- ROOT system ---
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPaveText.h>
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
#include <stdio.h>
#include <fstream>
#include <iostream>

#include "PlotGHcorrelation2.h"


using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;
/// \cond CLASSIMP
ClassImp(PlotGHcorrelation2);
/// \endcond

///
/// Default constructor
///
//________________________________________________________________________
PlotGHcorrelation2::PlotGHcorrelation2():
TObject(),
fInputFileSE(),fInputFileME(),fWagonNameSE(),fWagonNameME(),fRootFile(),fRootFileME(0),fUseHistogramFile(0),fObservable(-1),fObservableName(),
fPtMinBin(1),fPtMaxBin(-1),fmaxBins(0),fPlotAdvancedSB(0),fUseFindLastGoodBin(1),fWait(0),
//fPtMinBin(1),fPtMaxBin(-1),fmaxBins(0),fMergeMEplots(0),fPlotAdvancedSB(0),fUseFindLastGoodBin(1),fWait(0),
fDetaDphi_SE(),fDetaDphi_ME(),fDeta_Proj(),fDeta_ProjSub(),fsumCorrSE(),fsumCorrSE_alt1(),fsumCorrSE_alt2(),fsumCorrSE_alt3(),fsumCorrSE_ProjFull(),fsumCorrSE_NearEta(),fsumCorrSE_FarEta(),
fTrigger_SE(),fTriggerPt(0), fMassPtPionAcc(0), fMassPtPionRej(0), fMassPtCentPionAcc(0), fMassPtCentPionRej(0),
fMassPtPionAccProj(), fMassPtPionRejProj()

{
	fObservable=0;
	fCent      =-1;
	fEventPlane=-1;

	fOutputDir   = "output";
	fLabel       = "Proj";

	fTriggerName = "#gamma";
	fPlotOptions = "LEGO2";

//	fNSigma = 2.5; 
//	fNSigma = 3.5; 
//	fNSigma = 4.; 

	InitArrays();
}
///
/// constructor
///
//________________________________________________________________________
PlotGHcorrelation2::PlotGHcorrelation2(Int_t observable,Int_t centrality, Int_t eventPlane):
TObject(),
fInputFileSE(0),fInputFileME(0),fWagonNameSE(),fWagonNameME(),fRootFile(),fRootFileME(0),fUseHistogramFile(0), fObservable(-1), fObservableName(),//fMergeMEplots(),
fPtMinBin(1),fPtMaxBin(-1),fmaxBins(0), fPlotAdvancedSB(0),fPerformDPhiQA(1),fUseFindLastGoodBin(1), fWait(0),
//fPtMinBin(1),fPtMaxBin(-1),fmaxBins(0),fMergeMEplots(0), fPlotAdvancedSB(0),fUseFindLastGoodBin(1), fWait(0),
fDetaDphi_SE(),fDetaDphi_ME(),fDeta_Proj(),fDeta_ProjSub(),fsumCorrSE(),fsumCorrSE_alt1(),fsumCorrSE_alt2(),fsumCorrSE_alt3(),fsumCorrSE_ProjFull(),fsumCorrSE_NearEta(),fsumCorrSE_FarEta(),
fTrigger_SE(),fTriggerPt(0),fMassPtPionAcc(0), fMassPtPionRej(0), fMassPtCentPionAcc(0), fMassPtCentPionRej(0),
fMassPtPionAccProj(), fMassPtPionRejProj()
{
	//...set parameter values
	fObservable    =observable;
	fCent          =centrality;
	fEventPlane    =eventPlane;

	fOutputDir   = "output";
	fLabel       = "Proj";

	fTriggerName = "#gamma";
	//fPlotOptions = "LEGO2";
	fPlotOptions = "COLZ";

//	fNSigma = 2.5; 
//	fNSigma = 3; 
//	fNSigma = 4.; 

  // Preparing ME Scale Tree 
  fMEScaleTree = new TTree("MEScaleTree","ME Normalization Debug Tree");
  fMEScaleTree->SetDirectory(0);
  // 1D Branches
  fMEScaleTree->Branch("ValueAtZero",&fMEValueAtZero);
  fMEScaleTree->Branch("ValueAtPi",&fMEValueAtPi);
  fMEScaleTree->Branch("Mean1DValue",&fMEMean1DValue);
  fMEScaleTree->Branch("ValueAt1DMax",&fMEValueAt1DMax);
  fMEScaleTree->Branch("IndexOfNormRegion",&fMEIndexOfNormRegion);
  fMEScaleTree->Branch("ScaleValue",&fMEScaleValue);
  // 2D 
  fMEScaleTree->Branch("ValueAtOrigin",&fMEValueAtOrigin);
//  fMEScaleTree->Branch("Mean2DValue",&fMEMean2DValue);
  fMEScaleTree->Branch("ValueAt2DMax",&fMEValueAt2DMax);


  //Double_t fMEValueAtZero   = 0;
  //Double_t fMEValueAtPi     = 0;
  //Double_t fMEValueAt1DMax  = 0;
  //Int_t fMEIndexOfNormRegion = 0; ///< Bin of normalization region
  // From 2D Normalized Plot
  //   Double_t fMEValueAtOrigin = 0;
  //     Double_t fMEMeanValue     = 0;
  //       Double_t fMEValueAt2DMax  = 0;



	InitArrays();
}
///
/// Initializes values of arrays
///
//________________________________________________________________________
void PlotGHcorrelation2::InitArrays()
{
	fZtStep =1.0/(7-1.0);
	fXiStep =2.5/(8-1.0);

	//..assign valued to initialized arrays
  Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,20,23,30,60};
//	Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,22,30,60,90};
	//	Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,22};
	//	Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9};
	Double_t array_ZT_BinsValue[kZtNBINS+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
	Double_t array_XI_BinsValue[kXiNBINS+1]   ={-100,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};
 // Double_t fArray_HPT_BinsValue[kNoHPtBins+1]   ={0.15,0.4,0.8,1.45,2.5,4.2,6.95,11.4,18.6};
  Double_t fArray_HPT_BinsValue[kNoHPtBins+1]  ={0.2,0.4,0.8,1.5,2.5,4,7,11,17};


	//..assign values to z-vertex array
	//	Double_t array_zVtx_BinsValueAlt[kNvertBins_alt+1] ={-10,-6,-4,-2,0,2,4,6,10};
	Double_t array_zVtx_BinsValue[kNvertBins+1]        ={-10,-8,-6,-4,-2,0,2,4,6,8,10};
	Double_t array_zVtx_BinsValueAlt[kNvertBins_alt+1] ={-10,-6,-4,-2,0,2,4,6,10};

	//..assign values to Centrality array
	Double_t array_Cent_BinsValue[kCentBINS+1] = {0.0,10.0,30.0,50.0,90.0};
  // FIXME is it 90 or 80 here?

	//..assign values to phi binning
	//Double_t deltaPhiBinsValue[kNDeltaPhiBins]={180,125,180,24,156,24,132,24}; // degrees
	Double_t deltaPhiBinsValue[kNDeltaPhiBins]={180*DTR,125*DTR,180*DTR,24*DTR,156*DTR,24*DTR,132*DTR,24*DTR}; // degrees
	//Double_t deltaPhiBinsValue[kNDeltaPhiBins]={TMath::Pi(),25.*TMath::Pi()/36,TMath::Pi(),24,156,24,132,24}; // radians

	//..color sceme for plots
	Color_t fColorScemeValue[6]={kSpring+10,kCyan-10,kCyan-6,kOrange+7,kOrange-3,kOrange};
	//Double_t fColorScemeValue[6]={kGreen-9,kBlue-9,kRed-9,kOrange+7,kOrange-3,kOrange};

	memcpy (fArray_G_Bins, array_G_BinsValue, sizeof (fArray_G_Bins));
	memcpy (fArray_ZT_Bins, array_ZT_BinsValue, sizeof (fArray_ZT_Bins));
	memcpy (fArray_XI_Bins, array_XI_BinsValue, sizeof (fArray_XI_Bins));
  memcpy (fArray_HPT_Bins, fArray_HPT_BinsValue, sizeof (fArray_HPT_Bins));

	memcpy (fArray_zVtx_BinsAlt, array_zVtx_BinsValueAlt, sizeof (array_zVtx_BinsValueAlt));
	memcpy (fArray_zVtx_Bins, array_zVtx_BinsValue, sizeof (array_zVtx_BinsValue));
	memcpy (fArray_cent_Bins, array_Cent_BinsValue, sizeof (array_Cent_BinsValue));
	memcpy (fDeltaPhiBins, deltaPhiBinsValue, sizeof (fDeltaPhiBins));
	memcpy (fColorSceme, fColorScemeValue, sizeof (fColorSceme));

	if(fObservable==0)fmaxBins=6;
//	if(fObservable==0)fmaxBins=kGammaNBINS;
	if(fObservable==1)fmaxBins=kZtNBINS;
//	if(fObservable==2)fmaxBins=kXiNBINS;
	if(fObservable==2)fmaxBins=kNoHPtBins;

	fminZvtx=0;
	fmaxZvtx=10;

}


/**
  * Sets the global style
  */
void PlotGHcorrelation2::SetStyle() {

	gStyle->SetOptStat(0);
	//gStyle->SetPalette(53);  //standard is 1
  gStyle->SetPalette(kTemperatureMap);
	gStyle->SetCanvasColor(10);
//	TGaxis::SetMaxDigits(4);  //..ELI I don't remember why I wanted 4 changed to 2
	TGaxis::SetMaxDigits(3);

  // Margin Settings
//	gStyle->SetPadTopMargin(0.07);//0.05
//	gStyle->SetPadBottomMargin(0.18);//0.15
////	gStyle->SetPadRightMargin(0.045);
//	gStyle->SetPadRightMargin(0.08);
//	gStyle->SetPadLeftMargin(0.21);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.07);
  //gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);


	gStyle->SetFrameFillColor(10);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(5.0,"X");
	gStyle->SetTitleSize(5.0,"Y");
	TGaxis::SetMaxDigits(2);
	gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas
	gStyle->SetOptTitle(0);


}

///
/// Main run function
///
//________________________________________________________________________
void PlotGHcorrelation2::Run()
{
	//..security check for fPtMaxBin & fPtMaxBin
	if(fPtMaxBin>kGammaNBINS)
	{
		fPtMaxBin=kGammaNBINS;
		cout<<"---- Please fix: had to reset fPtMaxBin!"<<endl;
	}
	if(fPtMinBin<1)
	{
		fPtMinBin=1;
		cout<<"---- Please fix: had to reset fPtMinBin!"<<endl;
	}
	if(fmaxBins==0) cout<<"Something is wrong here!"<<endl;

	if(fCent>=kCentBINS)
	{
		fCent=kCentBINS-1;
		cout<<"---- Please fix: had to reset fCent!"<<endl;
	}

	if (fObservable==0)      fObservableName = "#it{p}_{T}";
	else if (fObservable==1) fObservableName = "z_{T}";
	else if (fObservable==2) fObservableName = "#it{p}_{T}^{assoc}";
//	else if (fObservable==2) fObservableName = "#xi";

	if (fGammaOrPi0==0)     fTriggerName = "#gamma";
	else if (fGammaOrPi0==1)fTriggerName = "#pi^{0}";
	else if (fGammaOrPi0==2)fTriggerName = "#pi^{0}-SB1";
	else if (fGammaOrPi0==3)fTriggerName = "#pi^{0}-SB2";
	else if (fGammaOrPi0==4)fTriggerName = "#pi^{0}-SB3";
	else if (fGammaOrPi0==5)fTriggerName = "#pi^{0}-SB4";
	else if (fGammaOrPi0==6)fTriggerName = "#pi^{0}-SB5";
	else if (fGammaOrPi0==7)fTriggerName = "#pi^{0}-SB6";

  // FIXME  For poster performance
  fTriggerName = "trig";

  SetStyle();

	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..This avoids that all the pdf/png printing status is posted on the screen

	//cout<<"in Run()"<<endl;
	if(!fInputFileSE && !fInputFileME)fUseHistogramFile=1;
	if(!fUseHistogramFile && !(fInputFileSE && fInputFileME)) {
		cout<<"Error: Missing one of the two train input files"<<endl;
		return;
	}
	if(!fInputFileSE && !fInputFileME && !fRootFile) {
		cout<<"Error, input file needs to be set"<<endl;
		return;
	}

	//..if you have an input root file reajust the input variables
	//..according to the file name
	ReadjustVariables();
	cout<<"o Analysis is run for:"<<endl;
	if(fObservable==0) cout<<"  o "<<fmaxBins<<" bins in GA"<<endl;
	if(fObservable==1) cout<<"  o "<<fmaxBins<<" bins in zT"<<endl;
	if(fObservable==2) cout<<"  o "<<fmaxBins<<" bins in pT^{h}"<<endl;
//	if(fObservable==2) cout<<"  o "<<fmaxBins<<" bins in xi"<<endl;
	if(fCent==-1)      cout<<"  o all centralities"<<endl;
	if(fCent>=0)       cout<<"  o centrality bin "<<fCent<<", this is "<<fArray_cent_Bins[fCent]<<"-"<<fArray_cent_Bins[fCent+1]<<"%"<<endl;
	if(fEventPlane==-1)cout<<"  o all event plane orientations "<<endl;
	if(fEventPlane==0) cout<<"  o gammas in plane "<<endl;
	if(fEventPlane==1) cout<<"  o gammas mid plane "<<endl;
	if(fEventPlane==2) cout<<"  o gammas out of plane "<<endl;




	//	SetPlotStyle();
/* // Might be taken care of in python script now
	//..Create first the general output folder, in case you run this analysis for the very first time
	//..Should be expanded later by the extention of whether this is Cent/peri or EvtPlane etc...
	if(fGammaOrPi0==0)fOutputDir="./output/GammaHadron";
	if(fGammaOrPi0==1)fOutputDir="./output/Pi0Hadron";
	if(fGammaOrPi0==2)fOutputDir="./output/Pi0SB1";
	if(fGammaOrPi0==3)fOutputDir="./output/Pi0SB2";
	gSystem->mkdir("./output");
	gSystem->mkdir(Form("%s",fOutputDir.Data()));
	gSystem->mkdir(Form("%s/CFiles",fOutputDir.Data()));
*/
	//..Load the histograms for the analyis
	LoadHistograms();
	cout<<"o Loaded histograms"<<endl;
	//not used currently 	DrawSEplots();


  ProduceTriggerPhiEtaPlots();

  fprintf(stderr,"Processing event plane histograms\n");
  ProduceDeltaPsiPlots();
  MeasureVn();

	// For Pi0 Analysis, Analyze Accepted/Rejected Pion Pt vs Mass
//	if(fGammaOrPi0>=1 && fProducePi0AnalysisPlots)
	if(fGammaOrPi0>=1)
	{
		cout<<"Analyzing Accepted Vs Rejected Pion Candidates"<<endl;
		AnalyzePionAccRej();
    cout<<"Done."<<endl;
	}


	if(fPlotMoreMEstrategies==1)
	{
		MergeMEplots();
	} 
  cout<<"Starting Normalization..."<<endl;
//	//..Normalize g-h pairs to No of triggers
//	if (fNormMode==1) NormalizeSEsignal();
//	cout<<"o Normalized SE"<<endl;
	//..Normalize the ME so that it is 1 at its plateau
	NormalizeMEsignal();
	cout<<"o Normalized ME and SE"<<endl;

  

  if(fPlot2DHistos) {
    Plot2DHistograms();
    cout<<"o Plotted raw ME and SE"<<endl;
  }

	PlotCorrected2DHistograms(); 
	cout<<"o Plotted corrected and normalized SE correlations"<<endl;

	if(fPlotMoreMEstrategies==1)
	{
		//..plot and compre the two different Mixed Event strategies
		//..How do they differ in yield?
		//..How do they differ in errors?
		PlotCompare2DHistoCorrection();
		//FIXME !!!to be done!!
		//for absolute differences and their error size (double check sumw2 stuff!!)
	}


	//-.-.-.-.-.-.Checkout the projections to the Delta Phi axis.-.-.-.-.-.-.-.-.-.
	//-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	if(fPlotAdvancedSB==1)
	{
		DetermineEtaWidths(fsumCorrSE);
		/// think about thatDetermine PhiWidths();
		cout<<"o Determined widths of NS peak"<<endl;

		DrawWidths();
		cout<<"o Plotted width of NS peak vs bin"<<endl;

		FitEtaSides(fsumCorrSE,fNSigma,fCanProj,fCanProj2,fCanProj3);
//		FitEtaSides(fsumCorrSE,2.5,fCanProj,fCanProj2,fCanProj3);
		cout<<"o Fitted eta side bands"<<endl;

		// Saving Plots
		fCanCorr1D->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fCanCorr1D->GetName()));
		fCanCorr1D->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanCorr1D->GetName()));
//		fCanCorr1D->Print(TString::Format("%s/%s.eps",fOutputDir.Data(),fCanCorr1D->GetName()));
		fCanCorr1D->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanCorr1D->GetName()));

		fCanCorr1D_Sub->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fCanCorr1D_Sub->GetName()));
		fCanCorr1D_Sub->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanCorr1D_Sub->GetName()));
//		fCanCorr1D_Sub->Print(TString::Format("%s/%s.eps",fOutputDir.Data(),fCanCorr1D_Sub->GetName()));
		fCanCorr1D_Sub->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanCorr1D_Sub->GetName()));

		fCanProj->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fCanProj->GetName()));
		fCanProj->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanProj->GetName()));
//		fCanProj->Print(TString::Format("%s/%s.eps",fOutputDir.Data(),fCanProj->GetName()));
		fCanProj->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanProj->GetName()));

		fCanProj2->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanProj2->GetName()));
		fCanProj2->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanProj2->GetName()));

		fCanProj3->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanProj3->GetName()));
		fCanProj3->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanProj3->GetName()));

		fCanProjFull->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fCanProjFull->GetName()));
		fCanProjFull->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanProjFull->GetName()));
//		fCanProjFull->Print(TString::Format("%s/%s.eps",fOutputDir.Data(),fCanProjFull->GetName()));
		fCanProjFull->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanProjFull->GetName()));

		fCanvWidth->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fCanvWidth->GetName()));
		fCanvWidth->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanvWidth->GetName()));
		fCanvWidth->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fCanvWidth->GetName()));
	}

  if (fPerformDPhiQA) {
    DPhiQA();
  }

	//..Save the width of the eta peak and the projections to an intermediate file
	//..this is done for each event plane orientation
	//..in step 2 the 3 different orientations are loaded and fitted together
	// FIXME removing fUseHistogramCheck (SaveIntermediate produces outputfile)
//	if(fUseHistogramFile==0)SaveIntermediateResult(2);
  if (bDenormalize) DenormalizeHists();
	SaveIntermediateResult(2);
	Wait();
	cout<<"o End of Program"<<endl;
}


///
///
///
//________________________________________________________________________
void PlotGHcorrelation2::DenormalizeHists() {
    cout<<"Undoing per-trigger normalization on output 1D histograms so they can be merged"<<endl;

    // Denorm scale
    double fDenormScale = 1;

    if (fObservable==0) {

      fprintf(stderr,"Error: using Observable = Trigger Pt, but I haven't coded the denormalization for this\n");
     // iIntegralPtBin = i+1;
     // fScale = fTriggerPt->Integral(iIntegralPtBin,iIntegralPtBin);
     // fScale = fTriggerPt->Integral(fPtMinBin,fPtMaxBin);
     // printf(" DEBUG: fMassPtPionAccProj_%d gives %f.\n",iIntegralPtBin-1,fMassPtPionAccProj[iIntegralPtBin-1]->Integral());
    } else { // 1 or 2
      double fScale = fTriggerPt->Integral();
      fDenormScale = fScale;
      //if (fScale != 0) fDenormScale = 1./fScale;
    }
    printf("Denorm Scale = %e\n",fDenormScale);
		//..Saving Delta eta projections
		for(Int_t i=0;i<fmaxBins;i++)
		{
			if (fDeta_Proj[i]) {
				fDeta_Proj[i]->Scale(fDenormScale);
			}
		}
		for(Int_t i=0;i<fmaxBins;i++)
		{
			if (fDeta_ProjSub[i]) {
				fDeta_ProjSub[i]->Scale(fDenormScale);
			}
		}
		for(Int_t i=0;i<fmaxBins;i++)
    {
      if (fDeta_AwaySide[i]) {
        fDeta_AwaySide[i]->Scale(fDenormScale);
      }
    }
		for (Int_t i=0;i<fmaxBins;i++)
		{
			if (fsumCorrSE_ProjFull[i]) {
				fsumCorrSE_ProjFull[i]->Scale(fDenormScale);
			}
		}

		for (Int_t i=0;i<fmaxBins;i++)
		{
			if (fsumCorrSE_NearEta[i]) {
				fsumCorrSE_NearEta[i]->Scale(fDenormScale);
			}
		}

		// Saving Far Eta histograms
		for (Int_t i=0;i<fmaxBins;i++)
		{
			if (fsumCorrSE_FarEta[i]) {
				fsumCorrSE_FarEta[i]->Scale(fDenormScale);
			}
		}



}

///
/// Readjust InputVariable read from the root file
///
//________________________________________________________________________
void PlotGHcorrelation2::ReadjustVariables()
{
	if(fRootFile)
	{
		cout<<"o Retrive settings from root file: "<<fRootFile->GetName()<<endl;
		TH1D* variables =(TH1D*)fRootFile->Get("VariableInfo");
		if(!variables)cout<<"Error: could not retrieve settings!"<<endl;
		fObservable=variables->GetBinContent(1);
		fCent      =variables->GetBinContent(2);
		fEventPlane=variables->GetBinContent(3);
	}
}
///
/// Allows to plot the results on the screen
///
//________________________________________________________________________
void PlotGHcorrelation2::Wait()
{
	TCanvas* z1= new TCanvas("z1","z1",0);
	if(fWait==1)
	{
		z1->Update();
		z1->WaitPrimitive();
	}
}
///
/// load the histograms from the provided files
///
//________________________________________________________________________
void PlotGHcorrelation2::LoadHistograms()
{
	for(Int_t j=0; j<fmaxBins;j++)  //-->kGammaNBINS+1??
	{
		fDetaDphi_SE[j]    = new TH2D*[kNvertBins+1];
		fDetaDphi_ME[j]    = new TH2D*[kNvertBins+1];
	}
	TString histName="";

	TStopwatch* watch = new TStopwatch();
	watch->Start();
	cout<<"o Start loading histograms"<<endl;
	double pi = TMath::Pi();



	//..get histograms from the THnSparse
	THnSparseF* corrVsParamSE = 0;
	THnSparseF* corrVsParamME = 0;
	THnSparseF* triggerHistSE = 0;
	bool haveTriggerHist = 0;
	if(fUseHistogramFile==0)
	{
		cout<<"  o Load SE ThnS from Train output: "<<fInputFileSE->GetName()<<endl;
		cout<<"  o Load ME ThnS from Train output: "<<fInputFileME->GetName()<<endl;

		if(!fInputFileSE)cout<<"Error: No SE input file with name: "<<fInputFileSE->GetName()<<" found "<<endl;
		if(!fInputFileME)cout<<"Error: No ME input file with name: "<<fInputFileME->GetName()<<" found "<<endl;
		TList* FinalListSE         =(TList*)fInputFileSE->Get(fWagonNameSE);
		TList* FinalListME         =(TList*)fInputFileME->Get(fWagonNameME);
		if(!FinalListSE)cout<<"Error: No SE wagon with name: "<<fWagonNameSE<<" found "<<endl;
		if(!FinalListME)cout<<"Error: No ME wagon with name: "<<fWagonNameME<<" found "<<endl;
		corrVsParamSE  =(THnSparseF*)FinalListSE->FindObject("CorrVsManyThings");
		corrVsParamME  =(THnSparseF*)FinalListME->FindObject("CorrVsManyThings");
		triggerHistSE  =(THnSparseF*)FinalListSE->FindObject("TriggerHist");
		cout<<"  o Project SE histograms from the THnSparse: "<<corrVsParamSE->GetName()<<endl;
		cout<<"  o Project ME histograms from the THnSparse: "<<corrVsParamME->GetName()<<endl;

    fHistEventHash = (TH1F *) FinalListSE->FindObject("HistEventHash");

		haveTriggerHist = (triggerHistSE != 0);
		if (haveTriggerHist)
		{
			cout<<"  o Found TriggerHist THnSparses"<<endl;
		}
		// Checking if gamma or Pi0 analysis
		//fGammaOrPi0 = (0 != FinalListSE->FindObject("fMassPtPionAcc"));
		if (fGammaOrPi0)
		{
			cout<<"  o Using Pi0 Triggers"<<endl;
			fMassPtPionAcc = (TH2F *) FinalListSE->FindObject("fMassPtPionAcc");
			fMassPtPionRej = (TH2F *) FinalListSE->FindObject("fMassPtPionRej");
		    //cout<<"fMassPtPionAcc "<<fMassPtPionAcc->GetName()<<endl;
		    //cout<<"fMassPtPionRej "<<fMassPtPionRej->GetName()<<endl;
			fMassPtPionAcc->SetDirectory(0);
			fMassPtPionRej->SetDirectory(0);
		    //cout<<"fMassPtPionAcc Dir. "<<fMassPtPionAcc->GetDirectory()->GetName()<<endl;
		    //cout<<"fMassPtPionRej Dir. "<<fMassPtPionRej->GetDirectory()->GetName()<<endl;
			fMassPtCentPionAcc = (TH3F *) FinalListSE->FindObject("fMassPtCentPionAcc");
			fMassPtCentPionRej = (TH3F *) FinalListSE->FindObject("fMassPtCentPionRej");
			fMassPtCentPionAcc->SetDirectory(0);
			fMassPtCentPionRej->SetDirectory(0);

      fEtaPhiPionAcc = (TH2F *) FinalListSE->FindObject("fEtaPhiPionAcc");
      fEtaPhiPionAcc->SetDirectory(0);



      // Saving some Vs Event Plane histograms that will be useful down the line
      hPtEPAnglePionAcc = (TH2F *) FinalListSE->FindObject("PtEPAnglePionAcc");
      hPtEPAnglePionAcc->SetDirectory(0);

      fPtEPAnglePionAccCent = (TH3F *) FinalListSE->FindObject("PtEPAnglePionAccCent");
      if (!fPtEPAnglePionAccCent) {
        fprintf(stderr,"Could not find fPtEPAnglePionAccCent\n");
      } else {
        fPtEPAnglePionAccCent->SetDirectory(0);
      }

      fPtEP3AnglePionAccCent = (TH3F *) FinalListSE->FindObject("PtEP3AnglePionAccCent");
      if (!fPtEP3AnglePionAccCent) {
        fprintf(stderr,"Could not find fPtEP3AnglePionAccCent\n");
      } else {
        fPtEP3AnglePionAccCent->SetDirectory(0);
      }
      fPtEP4AnglePionAccCent = (TH3F *) FinalListSE->FindObject("PtEP4AnglePionAccCent");
      if (!fPtEP4AnglePionAccCent) {
        fprintf(stderr,"Could not find fPtEP4AnglePionAccCent\n");
      } else {
        fPtEP4AnglePionAccCent->SetDirectory(0);
      }

      hHistTrackPsiEPPtCent = (TH3F *) FinalListSE->FindObject("fHistTrackPsiEPPtCent");
      hHistTrackPsiEPPtCent->SetDirectory(0);
      hHistTrackPsiEP3PtCent = (TH3F *) FinalListSE->FindObject("fHistTrackPsiEP3PtCent");
      if (hHistTrackPsiEP3PtCent) hHistTrackPsiEP3PtCent->SetDirectory(0);
      hHistTrackPsiEP4PtCent = (TH3F *) FinalListSE->FindObject("fHistTrackPsiEP4PtCent");
      if (hHistTrackPsiEP4PtCent) hHistTrackPsiEP4PtCent->SetDirectory(0);

		}
		else
		{
			cout <<"  o Using Gamma Triggers"<<endl;
		}

		//..Set the THN sparse ranges for the projections
		if (fObservable != 0)
		{   //..For Xi or Zt, set pt range to ptMinBin, ptMaxBins
			if (fPtMaxBin < 0) fPtMaxBin = kGammaNBINS;
			printf("  o Setting Pt bins: %i  -  %i\n",fPtMinBin,fPtMaxBin);
			printf("  o Setting Pt Range: %f  -  %f GeV/c \n",fArray_G_Bins[fPtMinBin-1],fArray_G_Bins[fPtMaxBin]);

			triggerHistSE->GetAxis(0)->SetRange(fPtMinBin,fPtMaxBin);
        printf("  o triggerTHn has range %f - %f\n",triggerHistSE->GetAxis(0)->GetBinLowEdge(fPtMinBin),triggerHistSE->GetAxis(0)->GetBinUpEdge(fPtMaxBin));

			corrVsParamSE->GetAxis(2)->SetRange(fPtMinBin,fPtMaxBin);
        printf("  o SameEvent Corr has range %f - %f\n",corrVsParamSE->GetAxis(2)->GetBinLowEdge(fPtMinBin),corrVsParamSE->GetAxis(2)->GetBinUpEdge(fPtMaxBin));

			corrVsParamME->GetAxis(2)->SetRange(fPtMinBin,fPtMaxBin);
        printf("  o MixedEvent Corr has range %f - %f\n",corrVsParamME->GetAxis(2)->GetBinLowEdge(fPtMinBin),corrVsParamME->GetAxis(2)->GetBinUpEdge(fPtMaxBin));

		}
	}

  // Old location of fTriggerPt and fTriggerPtWithinEPBin

/*	if(haveTriggerHist)
		fTriggerPt = (TH1D*)triggerHistSE->Projection(0); // Pt of trigger (after limiting pt range)
		fTriggerPt->SetDirectory(0);
		fTriggerPt->SetName("fTriggerPt"); // this will be used for normalizing later
		fTriggerPt->SetTitle("Trigger #it{p}_{T} (All Event Plane Angles);#it{p}_{T} (GeV/c)");
    // FIXME check if this has the right centrality limits
		
    // Now get the triggers within the specific event plane bin (used in phase 4)
    triggerHistSE->GetAxis(2)->SetRange(fEventPlane+1, fEventPlane+1);
    fTriggerPtWithinEPBin = (TH1D*)triggerHistSE->Projection(0);
    fTriggerPtWithinEPBin->SetDirectory(0);
    fTriggerPtWithinEPBin->SetName("fTriggerPtWithinEPBin");
    fTriggerPtWithinEPBin->SetTitle(Form("Trigger #it{p}_{T} in EP Bin %d;#it{p}_{T} (GeV/c)",fEventPlane+1));
  */
	//..get already projected histograms from the root file
	if(fUseHistogramFile==1)
	{
		cout<<"  o Load histograms from a root file: "<<fRootFile->GetName()<<endl;
		if (fRootFileME) {
			cout<<"  o Load ME histograms from separate root file: "<<fRootFileME->GetName()<<endl;
		}


    fHistEventHash = (TH1F *) fRootFile->Get("HistEventHash");


    fPhase2Purity = (TH1D *) fRootFile->Get("fPhase2Purity");
    if (fPhase2Purity) {
      fPhase2Purity->SetDirectory(0);
    } else {
			fprintf(stderr,"Warning: fPhase2Purity Missing\n");
    }
    fMCTriggerDist = (TH1D *) fRootFile->Get("fMCTriggerDist");
    if (fMCTriggerDist) {
      fMCTriggerDist->SetDirectory(0);
    } else {
      fprintf(stderr,"Warning: fMCTriggerDist Missing\n");
    }

		fTriggerPt = (TH1D*) fRootFile->Get("fTriggerPt");
		if (fTriggerPt) {
			fTriggerPt->SetDirectory(0);
		} else {
			fprintf(stderr,"Warning; fTriggerPt Missing\n");
		}
		fTriggerPtWithinEPBin = (TH1D*) fRootFile->Get("fTriggerPtWithinEPBin");
		if (fTriggerPtWithinEPBin) {
			fTriggerPtWithinEPBin->SetDirectory(0);
		} else {
			fprintf(stderr,"Warning; fTriggerPtWithinEPBin Missing\n");
		}


    // Getting these from the projections file
    // Saving some Vs Event Plane histograms that will be useful down the line
    hPtEPAnglePionAcc = (TH2F *) fRootFile->Get("PtEPAnglePionAcc");
    if (!hPtEPAnglePionAcc) {
      fprintf(stderr,"Could not find hPtEPAnglePionAcc\n");
    }
    hPtEPAnglePionAcc->SetDirectory(0);

    // PtEPAnglePionAccCent
    fPtEPAnglePionAccCent = (TH3F *) fRootFile->Get("PtEPAnglePionAccCent");
    if (!fPtEPAnglePionAccCent) {
      fprintf(stderr,"Could not find fPtEPAnglePionAccCent\n");
    } else {
      fPtEPAnglePionAccCent->SetDirectory(0);
    }
    // PtEP3AnglePionAccCent
    fPtEP3AnglePionAccCent = (TH3F *) fRootFile->Get("PtEP3AnglePionAccCent");
    if (!fPtEP3AnglePionAccCent) {
      fprintf(stderr,"Could not find fPtEP3AnglePionAccCent\n");
    } else {
      fPtEP3AnglePionAccCent->SetDirectory(0);
    }

    // PtEP4AnglePionAccCent
    fPtEP4AnglePionAccCent = (TH3F *) fRootFile->Get("PtEP4AnglePionAccCent");
    if (!fPtEP4AnglePionAccCent) {
      fprintf(stderr,"Could not find fPtEP4AnglePionAccCent\n");
    } else {
      fPtEP4AnglePionAccCent->SetDirectory(0);
    }




    // Track flow histograms
    hHistTrackPsiEPPtCent = (TH3F *) fRootFile->Get("fHistTrackPsiEPPtCent");
    if (!hHistTrackPsiEPPtCent) {
      fprintf(stderr,"Could not find hHistTrackPsiEPPtCent\n");
    }
    hHistTrackPsiEPPtCent->SetDirectory(0);
    hHistTrackPsiEP3PtCent = (TH3F *) fRootFile->Get("fHistTrackPsiEP3PtCent");
    if (!hHistTrackPsiEP3PtCent) {
      fprintf(stderr,"Could not find hHistTrackPsiEP3PtCent\n");
    } else {
      hHistTrackPsiEP3PtCent->SetDirectory(0);
    }
    hHistTrackPsiEP4PtCent = (TH3F *) fRootFile->Get("fHistTrackPsiEP4PtCent");
    if (!hHistTrackPsiEP4PtCent) {
      fprintf(stderr,"Could not find hHistTrackPsiEP4PtCent\n");
    } else {
      hHistTrackPsiEP4PtCent->SetDirectory(0);
    }

    //fMassPtCentPionAcc->SetDirectory(0);
    //fMassPtCentPionRej->SetDirectory(0);

    // Get the trackpT histograms 
    fTrackPtProjectionSE = (TH1D *) fRootFile->Get("TrackPtProjectionSE");
    fTrackPtProjectionSE->SetDirectory(0);
    if (!fTrackPtProjectionSE) {
      fprintf(stderr,"Could not find TrackPtProjectionSE\n");
    } else {
      printf("Got object %s from RootFile %s, calling it fTrackPtProjectionSE\n",fTrackPtProjectionSE->GetName(),fRootFile->GetName());
    }
    fTrackPtProjectionME = (TH1D *) fRootFile->Get("TrackPtProjectionME");
    fTrackPtProjectionME->SetDirectory(0);

    fTrackPtFromTrackPsi = (TH1D *) fRootFile->Get("TrackPtFromTrackPsi");
    fTrackPtFromTrackPsi->SetDirectory(0);

	} else { // Projecting from THnSparses
    printf("Projecting Histograms from THnSparses\n");
    //..Set centrality axis
    if(fCent>=0)
    {
      corrVsParamSE->GetAxis(7)->SetRange(fCent+1, fCent+1);
      corrVsParamME->GetAxis(7)->SetRange(fCent+1, fCent+1);
      if (haveTriggerHist)
      {
        triggerHistSE->GetAxis(3)->SetRange(fCent+1, fCent+1);
      }
    }
    // Note: these are currently subject to the trigger pt bin limits
    // this may be interesting for SE, but useless for ME
    fTrackPtProjectionSE = (TH1D *) corrVsParamSE->Projection(4);
    fTrackPtProjectionSE->SetName("TrackPtProjectionSE");
    fTrackPtProjectionME = (TH1D *) corrVsParamME->Projection(4);
    fTrackPtProjectionME->SetName("TrackPtProjectionME");

    //
    if (hHistTrackPsiEPPtCent != 0) {
      // Set Centrality limits (z-axis)
      hHistTrackPsiEPPtCent->GetZaxis()->SetRange(fCent+1,fCent+1);
      

      fTrackPtFromTrackPsi = (TH1D *) hHistTrackPsiEPPtCent->ProjectionY("TrackPtFromTrackPsi");
      // Reset Centrality limits
      hHistTrackPsiEPPtCent->GetZaxis()->SetRange(1,4); // ncentbins
    }


    if(haveTriggerHist)
    {
      Int_t iTrigMCAxis = 4;
      Int_t iTrigMCTruePi0Axis = 3; // (2-1)
      // Create Purity Histogram
      // Note that this is ignoring event plane for now.
      // This will be 0 in data
      fInclusiveTriggerPt = (TH1D*)triggerHistSE->Projection(0); // Pt of trigger (after limiting pt range)
      fInclusiveTriggerPt->SetName("fInclusiveTriggerPt"); // Just used for purity calculation.

      if (triggerHistSE->GetNdimensions() > iTrigMCAxis) { // current check for MC mode
      //if (fMCTriggerDist != 0) { // current check for MC mode
        // Create MC True distributions
        fMCTriggerDist = (TH1D *) triggerHistSE->Projection(iTrigMCAxis);
        fMCTriggerDist->SetName("fMCTriggerDist");

        triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(iTrigMCTruePi0Axis,iTrigMCTruePi0Axis);
        fPhase2Purity = (TH1D*)triggerHistSE->Projection(0); // Pt of trigger (after limiting pt range)
        fPhase2Purity->SetDirectory(0);
        fPhase2Purity->SetName("fPhase2Purity"); // this will be used for normalizing later
        fPhase2Purity->Divide(fInclusiveTriggerPt);

        // Move the switch for the MCMode here??
        // iMCMode
        int maxbin = triggerHistSE->GetAxis(iTrigMCAxis)->GetNbins(); // for resetting
        switch (iMCMode) {
          default:
          case 0:
            triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(0,0);// Resets the bins to even include underflow
            //triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(1,maxbin);
            printf("DEBUG: Setting Trigger THn range to %d %d (Resetting to include all).\n (All)",0,0);
            //printf("DEBUG: Setting Trigger THn range to %d %d.\n (All)",1,maxbin);
            break;
          case 1: // Background Only
            triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(1,1);
            printf("DEBUG: Setting Trigger THn range to %d %d.\n (True Background)",1,1);
            break;
          case 2: // True Pi0s
            triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(3,3);
            printf("DEBUG: Setting Trigger THn range to %d %d (True Pi0s).\n",3,3);
            break;
          case 3: // True Eta(2 gamma)
            triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(2,2);
            printf("DEBUG: Setting Trigger THn range to %d %d (True Etas).\n",2,2);
            break;
        }
      }

      fTriggerPt = (TH1D*)triggerHistSE->Projection(0); // Pt of trigger (after limiting pt range)
      fTriggerPt->SetDirectory(0);
      fTriggerPt->SetName("fTriggerPt"); // this will be used for normalizing later
      fTriggerPt->SetTitle("Trigger #it{p}_{T} (All Event Plane Angles);#it{p}_{T} (GeV/c)");
      
      // Now get the triggers within the specific event plane bin (used in phase 4)
			if(fEventPlane>=0) {
        triggerHistSE->GetAxis(2)->SetRange(fEventPlane+1, fEventPlane+1);
      }
      fTriggerPtWithinEPBin = (TH1D*)triggerHistSE->Projection(0);
      fTriggerPtWithinEPBin->SetDirectory(0);
      fTriggerPtWithinEPBin->SetName("fTriggerPtWithinEPBin");
      fTriggerPtWithinEPBin->SetTitle(Form("Trigger #it{p}_{T} in EP Bin %d;#it{p}_{T} (GeV/c)",fEventPlane+1));
			if(fEventPlane>=0) {
        // Resetting this axis.
        triggerHistSE->GetAxis(2)->SetRange(0,0); // This resets the axis
      }
    }

  }
  printf("Finished prepping histograms\n");


	for(Int_t i=0;i<fmaxBins;i++)
	{
		cout<<endl;
		cout<<"- Bin No. "<<i<<" -"<<endl;
		//..axes: 0-phi, 1-eta, 2-Et gamma, 3-zT, 4-Xi, 5-z-vertex, 6-event plane, 7-centr.
		if(fUseHistogramFile==0)
		{
			//..Set observable axis
			corrVsParamSE->GetAxis(fObservable+2)->SetRange(i+1, i+1); //..Set Et bin, zT bin, or pTA bin
			corrVsParamME->GetAxis(fObservable+2)->SetRange(i+1, i+1); //..Set Et bin, zT bin, or pTA bin

			if (haveTriggerHist && (fObservable == 0))
			{   // Only need to restrict trigger range for pt/Et obs.
				// Restricting pt range bin by bin
				triggerHistSE->GetAxis(0)->SetRange(i+1,i+1);
        printf("Setting trigger THn axis %d to range %d %d.\n",0,i+1,i+1);
        // FIXME reenable this for vertex vs pt histos to work
			}

			//..Set centrality axis
/*			if(fCent>=0)
			{
				corrVsParamSE->GetAxis(7)->SetRange(fCent+1, fCent+1);
				corrVsParamME->GetAxis(7)->SetRange(fCent+1, fCent+1);
				if (haveTriggerHist)
				{
					triggerHistSE->GetAxis(3)->SetRange(fCent+1, fCent+1);
				}
			}*/
			//..Set Event plane axis
			if(fEventPlane>=0)
			{
				corrVsParamSE->GetAxis(6)->SetRange(fEventPlane+1, fEventPlane+1);
        // Don't Use Event plane for mixed events.
				//if (bSplitMEbyEventPlane) corrVsParamME->GetAxis(6)->SetRange(fEventPlane+1, fEventPlane+1);



				if (haveTriggerHist)
				{
					// FIXME if we want the normalization in EP bins to be for all triggers,
					// as opposed to for triggers in the given event plane, don't do this
          // August 16, 2019: reenabled this
					//triggerHistSE->GetAxis(2)->SetRange(fEventPlane+1, fEventPlane+1);
          // August 29, 2019: redisabled this, as it is incorrect for what goes into the RPF
          // the per trigger in each bin must be done after the RPF fit.
          // Oct 3, 2019: this now done later, to produce a second trigger histogram
				}
			}
      // Extracting trigger pt distributions
			// Extracting trigger z-vertex histograms
			// FIXME make ability to get from histogram file
			// FIXME what about limiting the zVertex range?
			//       that will effectivley throw out events and thus also triggers

      // This line of code is setting the event plane angle bin to all angles.
      triggerHistSE->GetAxis(2)->SetRange(1, triggerHistSE->GetAxis(2)->GetNbins());
      printf("Setting the TriggerTHn to event plane range 1, %d\n",triggerHistSE->GetAxis(2)->GetNbins());
			if(haveTriggerHist)
			{

				fTrigger_SE[i] = (TH1D*)triggerHistSE->Projection(1);
				fTrigger_SE[i]->SetName(TString::Format("fTrigger_SE_proj_1_bin_%d",i));
				cout<<"Trigger integral: "<<fTrigger_SE[i]->Integral()<<endl;
        // no longer using fTrigger_SE for the trigger counting
			}
/*
      if(haveTriggerHist)
      {// move this outside the i loop
        fTriggerPt = (TH1D*)triggerHistSE->Projection(0); // Pt of trigger (after limiting pt range)
        fTriggerPt->SetDirectory(0);
        fTriggerPt->SetName("fTriggerPt"); // this will be used for normalizing later
        fTriggerPt->SetTitle("Trigger #it{p}_{T} (All Event Plane Angles);#it{p}_{T} (GeV/c)");
        
        // Now get the triggers within the specific event plane bin (used in phase 4)
        triggerHistSE->GetAxis(2)->SetRange(fEventPlane+1, fEventPlane+1); // FIXME check that this doesn't need to 
        // be undone
        fTriggerPtWithinEPBin = (TH1D*)triggerHistSE->Projection(0);
        fTriggerPtWithinEPBin->SetDirectory(0);
        fTriggerPtWithinEPBin->SetName("fTriggerPtWithinEPBin");
        fTriggerPtWithinEPBin->SetTitle(Form("Trigger #it{p}_{T} in EP Bin %d;#it{p}_{T} (GeV/c)",fEventPlane+1));
      }*/
		}

		for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
		{
			cout<<". "<<flush;
			if(j==9)cout<<endl;
			if(fUseHistogramFile==0)
			{
				TString xTitle,yTitle;
				if(fObservable==0)yTitle=Form("#Delta#eta^{%s-h} %0.1f<#it{p}_{T}<%0.1f",fTriggerName.Data(),fArray_G_Bins[i],fArray_G_Bins[i+1]);
				if(fObservable==1)yTitle=Form("#Delta#eta^{%s-h} %0.1f<z_{T}<%0.1f",fTriggerName.Data(),fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]);
//				if(fObservable==2)yTitle=Form("#Delta#eta^{%s-h} %0.1f<#xi<%0.1f",fTriggerName.Data(),fArray_XI_Bins[i],fArray_XI_Bins[i+1]);
				if(fObservable==2)yTitle=Form("#Delta#eta^{%s-h} %0.1f<#it{p}_{T}^{assoc}<%0.1f",fTriggerName.Data(),fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]);
				if(fObservable==0)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<#it{p}_{T}<%0.1f",fTriggerName.Data(),fArray_G_Bins[i],fArray_G_Bins[i+1]);
				if(fObservable==1)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<z_{T}<%0.1f",fTriggerName.Data(),fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]);
//				if(fObservable==2)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<#xi<%0.1f",fTriggerName.Data(),fArray_XI_Bins[i],fArray_XI_Bins[i+1]);
				if(fObservable==2)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<#it{p}_{T}^{assoc}<%0.1f",fTriggerName.Data(),fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]);


        Int_t iCorrMCAxis = 8;
        Int_t iTrigMCAxis = 4;
        // Setting MC Information FIXME
        // iMCMode
        switch (iMCMode) {
          default:
          case 0:
            break;
          case 1: // Background Only
            corrVsParamSE->GetAxis(iCorrMCAxis)->SetRange(1,1);
 //           triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(1,1);
            break;
          case 2: // True Pi0s
            corrVsParamSE->GetAxis(iCorrMCAxis)->SetRange(3,3);
   //         triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(3,3);
            break;
          case 3: // True Eta(2 gamma)
            corrVsParamSE->GetAxis(iCorrMCAxis)->SetRange(2,2);
     //       triggerHistSE->GetAxis(iTrigMCAxis)->SetRange(2,2);
            break;
        }


				corrVsParamSE->GetAxis(5)->SetRange(j+1, j+1); //..bins start at 1
				fDetaDphi_SE[i][j] = corrVsParamSE->Projection(1,0);
				fDetaDphi_SE[i][j] ->SetName(Form("fHistDEtaDPhiG%d_Id1_Vert%d",i,j));
				fDetaDphi_SE[i][j]->GetXaxis()->SetTitle(xTitle);
				fDetaDphi_SE[i][j]->GetYaxis()->SetTitle(yTitle);
				fDetaDphi_SE[i][j] ->Sumw2();

        // FIXME Around here we can choose to not restrict pt range for mixed events
				corrVsParamME->GetAxis(5)->SetRange(j+1, j+1); //..bins start at 1
				fDetaDphi_ME[i][j] = corrVsParamME->Projection(1,0);
				fDetaDphi_ME[i][j] ->SetName(Form("fHistDEtaDPhiG%d_Id0_Vert%d",i,j));
				fDetaDphi_ME[i][j] ->Sumw2();
			}
			if(fUseHistogramFile==1)
			{
				fDetaDphi_SE[i][j] =(TH2D*)fRootFile->Get(Form("fHistDEtaDPhiG%d_Id1_Vert%d",i,j));
				if (fRootFileME) fDetaDphi_ME[i][j] =(TH2D*)fRootFileME->Get(Form("fHistDEtaDPhiG%d_Id0_Vert%d",i,j));
				else fDetaDphi_ME[i][j] =(TH2D*)fRootFile->Get(Form("fHistDEtaDPhiG%d_Id0_Vert%d",i,j));

				if(!fDetaDphi_SE[i][j])cout<<"did not find SE histogram: "<<j<<" (zVtx) "<<i<<"(Bin)"<<endl;
				if(!fDetaDphi_ME[i][j])cout<<"did not find ME histogram: "<<j<<" (zVtx) "<<i<<"(Bin)"<<endl;
				TString xTitle,yTitle;
				if(fObservable==0)yTitle=Form("#Delta#eta^{%s-h} %0.1f<#it{p}_{T}<%0.1f",fTriggerName.Data(),fArray_G_Bins[i],fArray_G_Bins[i+1]);
				if(fObservable==1)yTitle=Form("#Delta#eta^{%s-h} %0.1f<z_{T}<%0.1f",fTriggerName.Data(),fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]);
//				if(fObservable==2)yTitle=Form("#Delta#eta^{%s-h} %0.1f<#xi<%0.1f",fTriggerName.Data(),fArray_XI_Bins[i],fArray_XI_Bins[i+1]);
				if(fObservable==2)yTitle=Form("#Delta#eta^{%s-h} %0.1f<#it{p}_{T}^{assoc}<%0.1f",fTriggerName.Data(),fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]);
				if(fObservable==0)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<#it{p}_{T}<%0.1f",fTriggerName.Data(),fArray_G_Bins[i],fArray_G_Bins[i+1]);
				if(fObservable==1)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<z_{T}<%0.1f",fTriggerName.Data(),fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]);
//				if(fObservable==2)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<#xi<%0.1f",fTriggerName.Data(),fArray_XI_Bins[i],fArray_XI_Bins[i+1]);
				if(fObservable==2)xTitle=Form("#Delta#varphi^{%s-h} %0.1f<#it{p}_{T}^{assoc}<%0.1f",fTriggerName.Data(),fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]);

				fDetaDphi_SE[i][j]->GetXaxis()->SetTitle(xTitle);
				fDetaDphi_SE[i][j]->GetYaxis()->SetTitle(yTitle);

				fDetaDphi_SE[i][j] ->SetDirectory(0);
				fDetaDphi_ME[i][j] ->SetDirectory(0);
			}
		}
		//..get properties of the 2D histogram and then empty the entries
		//..this is the histo which is used for adding the different z-vertex bins
//		fsumCorrSE[i] = (TH2D*)fDetaDphi_SE[i][0]->Clone(Form("%s_Clone",fDetaDphi_SE[i][0]->GetName()));
		fsumCorrSE[i] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("%s_Clone",fDetaDphi_SE[i][fminZvtx]->GetName()));
		fsumCorrSE[i]->Reset();
		fsumCorrSE_alt1[i] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("%s_Clone_alt1",fDetaDphi_SE[i][fminZvtx]->GetName()));
		fsumCorrSE_alt1[i]->Reset();
		fsumCorrSE_alt2[i] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("%s_Clone_alt2",fDetaDphi_SE[i][fminZvtx]->GetName()));
		fsumCorrSE_alt2[i]->Reset();
		fsumCorrSE_alt3[i] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("%s_Clone_alt3",fDetaDphi_SE[i][fminZvtx]->GetName()));
		fsumCorrSE_alt3[i]->Reset();

		fsumCorrSE[i]      ->SetDirectory(0);
		fsumCorrSE_alt1[i] ->SetDirectory(0);
		fsumCorrSE_alt2[i] ->SetDirectory(0);
		fsumCorrSE_alt3[i] ->SetDirectory(0);
	}
	cout<<endl;

	// Find trigger pT distribution
	if (fUseHistogramFile==1)
	{
		//Eli Don't understand difference to "A"
		//..Get these from the observable root file
		if (fGammaOrPi0)
		{
			fMassPtPionAcc = (TH2F *) fRootFile->Get("fMassPtPionAcc");
			fMassPtPionRej = (TH2F *) fRootFile->Get("fMassPtPionRej");
			if (!fMassPtPionAcc) cout<<"Problem loading fMassPtPionAcc"<<endl;
			if (!fMassPtPionRej) cout<<"Problem loading fMassPtPionRej"<<endl;
			fMassPtPionAcc ->SetDirectory(0);
			fMassPtPionRej ->SetDirectory(0);
			// Now in centrality bins:
			fMassPtCentPionAcc = (TH3F *) fRootFile->Get("fMassPtCentPionAcc");
			fMassPtCentPionRej = (TH3F *) fRootFile->Get("fMassPtCentPionRej");
			if (!fMassPtCentPionAcc) cout<<"Problem loading fMassPtCentPionAcc"<<endl;
			if (!fMassPtCentPionRej) cout<<"Problem loading fMassPtCentPionRej"<<endl;
			if (fMassPtCentPionAcc) fMassPtCentPionAcc->SetDirectory(0);
			if (fMassPtCentPionAcc) fMassPtCentPionRej->SetDirectory(0);

      fEtaPhiPionAcc = (TH2F *) fRootFile->Get("fEtaPhiPionAcc");
      if (!fEtaPhiPionAcc) cout<<"Problem loading fEtaPhiPionAcc"<<endl;
      if (fEtaPhiPionAcc) fEtaPhiPionAcc->SetDirectory(0);

		}

		for(Int_t i=0;i<fmaxBins;i++)
		{
			//Eli don#'t understand difference to "B"
			fTrigger_SE[i] = (TH1D*)fRootFile->Get(Form("fTrigger_SE_proj_1_bin_%d",i));
			fTrigger_SE[i] ->SetDirectory(0);
			if (fTrigger_SE[i]) haveTriggerHist = 1;
		}

		//..close the root file
		if(fRootFile)
		{
			fRootFile->Close();
			cout<<"o Input Rootfile closed"<<endl;
		}
	}

	//..test an alternative way of using the ME background
	if(fPlotMoreMEstrategies==1)
	{
		for(Int_t i=0;i<fmaxBins;i++)
		{
			fDetaDphi_ME_alt1[i] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("ME_alt1_ObservableBin%i",i));
			fDetaDphi_ME_alt1[i] ->Reset();
			fDetaDphi_SE_alt1[i] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("SE_alt1_ObservableBin%i",i));
			fDetaDphi_SE_alt1[i] ->Reset();

			for(Int_t j=0;j<kNvertBins_alt;j++)
			{
				fDetaDphi_ME_alt2[i][j] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("ME_alt2_ObservableBin%i_Vtx%i",i,j));
				fDetaDphi_ME_alt2[i][j] ->Reset();
				fDetaDphi_SE_alt2[i][j] = (TH2D*)fDetaDphi_SE[i][fminZvtx]->Clone(Form("SE_alt2_ObservableBin%i_Vtx%i",i,j));
				fDetaDphi_SE_alt2[i][j] ->Reset();
			}
			for(Int_t j=0;j<fmaxZvtx;j++)
			{
				if(i==0)fDetaDphi_ME_alt3[j] = (TH2D*)fDetaDphi_SE[0][fminZvtx]->Clone(Form("ME_alt3_AnyObservableBin_Vtx%i",j));
				if(i==0)fDetaDphi_ME_alt3[j] ->Reset();
				fDetaDphi_SE_alt3[i][j] = (TH2D*)fDetaDphi_SE[i][j]->Clone(Form("SE_alt3_AnyObservableBin%i_Vtx%i",i,j));
				//fDetaDphi_SE_alt3[i][j] ->Reset(); //..Do not reset keep the original binning from fObservables and z-Vtx bins
			}
		}
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Histograms created new from the analysis information
/*	for(Int_t i=0;i<4;i++)
	{
		if(fObservable==0)fYield_VS_Eg[i] = new TH1D(Form("Yield_VS_EGamma_Angle%0d",i),Form("Yield_VS_EGamma_Angle%0d",i), kGammaNBINS, fArray_G_Bins);
		if(fObservable==1)fYield_VS_Eg[i] = new TH1D(Form("Yield_VS_Zt_Angle%0d",i),Form("Yield_VS_Zt_Angle%0d",i), kZtNBINS, fArray_ZT_Bins);
//		if(fObservable==2)fYield_VS_Eg[i] = new TH1D(Form("Yield_VS_Xi_Angle%0d",i),Form("Yield_VS_Xi_Angle%0d",i), kXiNBINS, fArray_XI_Bins);
		if(fObservable==2)fYield_VS_Eg[i] = new TH1D(Form("Yield_VS_HPt_Angle%0d",i),Form("Yield_VS_HPt_Angle%0d",i), kNoHPtBins, fArray_HPT_Bins);
	}*/

	if(fObservable==0)
	{
		fEtaWidth = new TH1D("EtaWidth","EtaWidth",kGammaNBINS, fArray_G_Bins);
		fPhiWidth = new TH1D("PhiWidth","PhiWidth",kGammaNBINS, fArray_G_Bins);
	}
	if(fObservable==1)
	{
		fEtaWidth = new TH1D("EtaWidth","EtaWidth",kZtNBINS, fArray_ZT_Bins);
		fPhiWidth = new TH1D("PhiWidth","PhiWidth",kZtNBINS, fArray_ZT_Bins);
	}
	if(fObservable==2)
	{
	//	fEtaWidth = new TH1D("EtaWidth","EtaWidth",kXiNBINS, fArray_XI_Bins);
	//	fPhiWidth = new TH1D("PhiWidth","PhiWidth",kXiNBINS, fArray_XI_Bins);
		fEtaWidth = new TH1D("EtaWidth","EtaWidth",kNoHPtBins, fArray_HPT_Bins);
		fPhiWidth = new TH1D("PhiWidth","PhiWidth",kNoHPtBins, fArray_HPT_Bins);
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Create boxes for the histograms
	for(Int_t i=0;i<kNDeltaPhiBins/2;i++)
	{
		fBoxes[i*2]   =new TBox(fDeltaPhiBins[i*2]-fDeltaPhiBins[i*2+1],0, fDeltaPhiBins[i*2],100);

		if(fDeltaPhiBins[i*2]==pi)
			//		if(fDeltaPhiBins[i*2]==180)
		{
			fBoxes[i*2+1] =new TBox(fDeltaPhiBins[i*2],0,fDeltaPhiBins[i*2]+fDeltaPhiBins[i*2+1],100);
		}
		else
		{
			Double_t Delta=TMath::Pi()-fDeltaPhiBins[i*2];
			fBoxes[i*2+1] =new TBox(TMath::Pi()+Delta,0,pi+Delta+fDeltaPhiBins[i*2+1],100);
			//			Double_t Delta=180-fDeltaPhiBins[i*2];
			//			fBoxes[i*2+1] =new TBox(180+Delta,0,180+Delta+fDeltaPhiBins[i*2+1],100);
		}
	}
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Create Canvases
	for(Int_t j=0; j<fmaxBins; j++)
	{
		fMEPlots1DGamma[j]       = new TCanvas(Form("ME_1D_Plots_gBin%0d",j),Form("A) 1D ME norm. Bin%0d",j),1750,1400);
		fMEPlots2DGamma[j]       = new TCanvas(Form("ME_2D_Plots_gBin%0d",j),Form("B) 2D ME norm. Bin%0d",j),1750,1400);
		fRaw_Plots_2D[j]         = new TCanvas(Form("Raw_2D_Plots1Bin%0d",j),Form("C) 2D SE raw Bin%0d",j),1750,1400);
		fRaw_Plots_2D_alt2[j]    = new TCanvas(Form("Raw_2D_Plots1Bin%0d_alt2",j),Form("C) 2D SE raw Bin%0d alt2",j),1750,1400);
		fRaw_Plots_2D_alt3[j]    = new TCanvas(Form("Raw_2D_Plots1Bin%0d_alt3",j),Form("C) 2D SE raw Bin%0d alt3",j),1750,1400);
		fPlots_2D_Corr[j]        = new TCanvas(Form("Corrected_2D_Plots%0d",j),Form("D) 2D Correct. multiple Z-vtx for Bin %0d",j),1750,1400);
	}
	fPlots_2D_CorrSum       = new TCanvas("Corr_2D_Plots","E) 2D Corr Sum",1400,900);
	fPlots_2D_CorrSum_alt1  = new TCanvas("Corr_2D_Plots_alt1","E) 2D Corr Sum",1400,900);
	fPlots_2D_CorrSum_alt2  = new TCanvas("Corr_2D_Plots_alt2","E) 2D Corr Sum",1400,900);
	fPlots_2D_CorrSum_alt3  = new TCanvas("Corr_2D_Plots_alt3","E) 2D Corr Sum",1400,900);
	fCanBinCheck   = new TCanvas("Can_BinCheck","Check for the binning",800,800);
	fCanNormCheck  = new TCanvas("Can_NormCheck","Check for the Normalization of triggers",800,800);

	if(fPlotMoreMEstrategies==1)
	{
		fPlots_2D_ME_alt1    = new TCanvas("ME_2D_PlotsAlt1","E) 2D ME Sum (one zVert bin)",1400,900);
		fPlots_2D_Corr_alt1  = new TCanvas("Corr_2D_PlotsAlt1","E) 2D Corr Sum (one zVert bin)",1400,900);
		fPlots_2D_ME_alt3    = new TCanvas(Form("ME_2D_PlotsAlt3"),Form("E) 2D ME Sum (summed over triggers) all Bins"),1400,900);
		for(Int_t j=0; j<fmaxBins; j++)
		{
			fPlots_2D_ME_alt2[j]  = new TCanvas(Form("ME_2D_PlotsAlt2_%0d",j),Form("E) 2D ME Sum (reduced zVert bins) Bin%0d",j),1400,900);
			fPlots_2D_Corr_alt2[j]= new TCanvas(Form("Corr_2D_PlotsAlt2_%0d",j),Form("E) 2D Corr Sum (reduced zVert bins) Bin%0d",j),1400,900);
			fPlots_2D_Corr_alt3[j]= new TCanvas(Form("Corr_2D_PlotsAlt3_%0d",j),Form("E) 2D Corr Sum (summed over triggers) Bin%0d",j),1400,900);
		}
		fPlots_1D_MEcompare      = new TCanvas("CompareCorr_1D_PlotsAlt","E) 1D Comparision of many methods",1400,900);
		fPlots_1D_MEcompareEta   = new TCanvas("CompareCorr_1D_PlotsAltEta","E) 1D Comparision of many methods DeltaEta",1400,900);
		fPlots_1D_MEcompareRatio = new TCanvas("CompareCorr_1D_PlotsAltRatio","E) 1D Ratio comparision of many methods",1400,900);
		fPlots_1D_MEcompareErrors= new TCanvas("CompareCorr_1D_PlotsAltErrors","E) 1D Error comparision of many methods",1400,900);
	}

	if(fPlotAdvancedSB==1)
	{
		fCanCorr1D    = new TCanvas("Corr1D_Plots","F) Corrected 1D Plots 1) Full",1400,900); //1400
		fCanCorr1D_Sub= new TCanvas("Corr1D_Plots_Sub","F) Corrected 1D Plots 1) UE subtracted",1400,900);
		fCanProj      = new TCanvas("Corr1D_Plots_SB","F) Corrected 1D Plots 2) SB",1400,900);
		fCanProj2     = new TCanvas("Corr1D_Plots_SB2","F) Corrected 1D Plots 2) SB compare",1400,900);
		fCanProj3     = new TCanvas("Corr1D_Plots_SB3","F) Corrected 1D Plots 2) SB compare",1400,900);
		fCanProjFull  = new TCanvas("CanProjFull","F) Corrected 1D Plots 3) All",1400,900);
		fCanvWidth    = new TCanvas("Canv_Width","F) Width of NS peak vs bin",1400,700);
		fCanProjFull->Divide(3,3,0.001,0.0012);
	//	fCanProjFull->Divide(3,2,0.001,0.0012);
		//fCanProjFull->Divide(3,3,0.001,0.0012);
	}
	watch->Stop();
	cout<<"o Finished loading all histograms in "<<watch->RealTime()/60<<" min"<<endl;

///* //Why was this commented out?
	//..If there is no input root file, save the loaded histograms to a file
	if(!fRootFile)
	{
		TString fileName;
		//fileName = Form("./Observable%i_Histograms.root",fObservable);
		fileName = Form("Projections/%s_Observable%i_Cent%i_EvtPlane%i_Histograms.root",fLabel.Data(),fObservable,fCent,fEventPlane);
		fRootFile = new TFile(fileName,"recreate");
	}
	if (fUseHistogramFile==0) {
		for(Int_t i=0;i<fmaxBins;i++)
		{
			for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
			//for(Int_t j=0;j<kNvertBins;j++)
			{
				fRootFile->WriteObject(fDetaDphi_SE[i][j],fDetaDphi_SE[i][j]->GetName());
				fRootFile->WriteObject(fDetaDphi_ME[i][j],fDetaDphi_ME[i][j]->GetName());
			}
		}
    if (fHistEventHash) fRootFile->WriteObject(fHistEventHash,fHistEventHash->GetName());

		TH1D* VariableInfo = 0;
		if(fUseHistogramFile==0)
		{
			VariableInfo = new TH1D("VariableInfo","VariableInfo",3,0,3);
			VariableInfo->SetBinContent(1,fObservable);
			VariableInfo->SetBinContent(2,fCent);
			VariableInfo->SetBinContent(3,fEventPlane);
		}

		fRootFile->WriteObject(VariableInfo,VariableInfo->GetName()); //..pass variable info into root file

    if (fTrackPtProjectionSE != 0) fRootFile->WriteObject(fTrackPtProjectionSE,fTrackPtProjectionSE->GetName());
    else fprintf(stderr,"TrackPtProjectionSE is not being saved\n");
    if (fTrackPtProjectionME != 0) fRootFile->WriteObject(fTrackPtProjectionME,fTrackPtProjectionME->GetName());
    if (fTrackPtFromTrackPsi != 0) fRootFile->WriteObject(fTrackPtFromTrackPsi,fTrackPtFromTrackPsi->GetName());

    if (hPtEPAnglePionAcc != 0) fRootFile->WriteObject(hPtEPAnglePionAcc,hPtEPAnglePionAcc->GetName());
    if (fPtEPAnglePionAccCent != 0) fRootFile->WriteObject(fPtEPAnglePionAccCent,fPtEPAnglePionAccCent->GetName());
    if (fPtEP3AnglePionAccCent != 0) fRootFile->WriteObject(fPtEP3AnglePionAccCent,fPtEP3AnglePionAccCent->GetName());
    if (fPtEP4AnglePionAccCent != 0) fRootFile->WriteObject(fPtEP4AnglePionAccCent,fPtEP4AnglePionAccCent->GetName());

    if (hHistTrackPsiEPPtCent != 0) fRootFile->WriteObject(hHistTrackPsiEPPtCent,hHistTrackPsiEPPtCent->GetName());
    if (hHistTrackPsiEP3PtCent != 0) fRootFile->WriteObject(hHistTrackPsiEP3PtCent,hHistTrackPsiEP3PtCent->GetName());
    if (hHistTrackPsiEP4PtCent != 0) fRootFile->WriteObject(hHistTrackPsiEP4PtCent,hHistTrackPsiEP4PtCent->GetName());


		if (haveTriggerHist) {
			if (fPhase2Purity) fRootFile->WriteObject(fPhase2Purity,fPhase2Purity->GetName());
      if (fMCTriggerDist) fRootFile->WriteObject(fMCTriggerDist,fMCTriggerDist->GetName());
			fRootFile->WriteObject(fTriggerPt,fTriggerPt->GetName());
      fRootFile->WriteObject(fTriggerPtWithinEPBin,fTriggerPtWithinEPBin->GetName());
			for(Int_t i=0;i<fmaxBins;i++) {
				//fRootFile->Add(fTrigger_SE[i]);
				fRootFile->WriteObject(fTrigger_SE[i],fTrigger_SE[i]->GetName());
			}
		}

		if (fMassPtPionAcc) fRootFile->WriteObject(fMassPtPionAcc,fMassPtPionAcc->GetName());
		if (fMassPtPionRej) fRootFile->WriteObject(fMassPtPionRej,fMassPtPionRej->GetName());
		if (fMassPtCentPionAcc) fRootFile->WriteObject(fMassPtCentPionAcc,fMassPtCentPionAcc->GetName());
		if (fMassPtCentPionRej) fRootFile->WriteObject(fMassPtCentPionRej,fMassPtCentPionRej->GetName());
		// Axis Reference Histograms
	//	    if (fGammaPt2) fRootFile->Add(fGammaPt2);
  //  else fprintf(stderr,"Missing fGammaPt2\n");
  //  if (fGammahXi) fRootFile->Add(fGammahXi);
  //  else fprintf(stderr,"Missing fGammahXi\n");
  //  if (fGammahZt) fRootFile->Add(fGammahZt);
  //  else fprintf(stderr,"Missing fGammahZt\n");
		 
		printf("MHO: Writing fRootFile ...\n");
			fRootFile->Write();
		printf("MHO: Closing fRootFile ...\n");
			fRootFile->Close();
		printf("MHO: done with fRootFile\n");
	}
//*/
	// FIXME MHO trying removing fUseHistogramFile check
	SaveIntermediateResult(1);
//	if(fUseHistogramFile==0)SaveIntermediateResult(1);
}
///
/// load the 2D histograms from input file
///
//________________________________________________________________________
TH2D* PlotGHcorrelation2::Get2DHistoFromFile(TFile* RootFile, TString name)
{
	TH2D* Histo   =(TH2D*)RootFile->FindObject(name);
	return Histo;
}
///
/// load the 1D histograms from input file
///
//________________________________________________________________________
TH1D* PlotGHcorrelation2::Get1DHistoFromFile(Bool_t smaMix,TString subListName,TString name)
{
	TFile* RootFile= nullptr;
	TString WagonName;
	if(smaMix==0)RootFile=fInputFileME; //..mixed event
	if(smaMix==1)RootFile=fInputFileSE;
	if(smaMix==0)WagonName=fWagonNameME; //..mixed event
	if(smaMix==1)WagonName=fWagonNameSE;

	TList* IntermediatList;
	TList* FinalList;

	if(subListName=="")
	{
		FinalList    =(TList*)RootFile->Get(WagonName);
	}
	else
	{
		IntermediatList=(TList*)RootFile       ->Get(WagonName);
		FinalList      =(TList*)IntermediatList->FindObject(subListName);
	}
	TH1D* Histo   =(TH1D*)FinalList->FindObject(name);

	return Histo;
}

void PlotGHcorrelation2::ProduceTriggerPhiEtaPlots() {
  if (!fEtaPhiPionAcc) return;

  TCanvas * cPhiEta = new TCanvas("PhiEta","PhiEta");


  fEtaPhiPionAcc->Draw("COLZ");
  fEtaPhiPionAcc->GetXaxis()->SetMaxDigits(5);

  cPhiEta->SetRightMargin(0.1);

  cPhiEta->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),"EtaPhiAcceptPion"));
  cPhiEta->Print(TString::Format("%s/%s.png",fOutputDir.Data(),"EtaPhiAcceptPion"));
  cPhiEta->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),"EtaPhiAcceptPion"));

  // I want these in output -> compare different triggers
  fEtaPionAcc = (TH1F *) fEtaPhiPionAcc->ProjectionX();
  fPhiPionAcc = (TH1F *) fEtaPhiPionAcc->ProjectionY();

  // FIXME add the projections to output files

}





///
/// Saves an intermediate result for a later simultaneous fitting
///
//________________________________________________________________________
void PlotGHcorrelation2::SaveIntermediateResult(Int_t stage)
{
  // Lazily updating iFixedDEtaCutIndex
  //if (fFixedDEtaCut == -1) iFixedDEtaCutIndex = 0;
  //if (fFixedDEtaCut == 0.8) iFixedDEtaCutIndex = 1;
  //if (fFixedDEtaCut == 0.7) iFixedDEtaCutIndex = 2;

	TString fileName;
	cout<<"o Save intermediate result for: obs: "<<fObservable<<" cent: "<<fCent<<" EvtPlane"<<fEventPlane<<endl;
	cout<<"o Saving stage: "<<stage<<endl;
	fileName = Form("output/IntermediateResult_%s_Observable%i_Cent%i_DEta%i_EvtPlane%i.root",fLabel.Data(),fObservable,fCent,iFixedDEtaCutIndex,fEventPlane);
	//fileName = Form("output/IntermediateResult_%s_Observable%i_Cent%i_EvtPlane%i.root",fLabel.Data(),fObservable,fCent,fEventPlane);
//	if(fGammaOrPi0==0)fileName = Form("output/IntermediateResult_GH_Observable%i_Cent%i_EvtPlane%i.root",fObservable,fCent,fEventPlane);
//	if(fGammaOrPi0==1)fileName = Form("output/IntermediateResult_PiH_Observable%i_Cent%i_EvtPlane%i.root",fObservable,fCent,fEventPlane);
//	if(fGammaOrPi0==2)fileName = Form("output/IntermediateResult_PiHSB1_Observable%i_Cent%i_EvtPlane%i.root",fObservable,fCent,fEventPlane);
//	if(fGammaOrPi0==3)fileName = Form("output/IntermediateResult_PiHSB2_Observable%i_Cent%i_EvtPlane%i.root",fObservable,fCent,fEventPlane);
	cout<<"o File Name: "<<fileName<<endl;
	TFile* outputRootFile=nullptr;


	if(stage==1)
	{
		for(Int_t i=0;i<fmaxBins;i++)
		{
			for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
			//for(Int_t j=0;j<kNvertBins;j++)
			{
				fDetaDphi_SE[i][j]->SetDirectory(0); //FIXME now it needs to be manually deleted
				fDetaDphi_ME[i][j]->SetDirectory(0);
			}
		}
  //  fTrackPtProjectionSE->SetDirectory(0);
 //   fTrackPtProjectionME->SetDirectory(0);
	//	fMassPtPionAcc->SetDirectory(0); //FIXME now it needs to be manually deleted
//		fMassPtPionRej->SetDirectory(0); //FIXME now it needs to be manually deleted
    printf("Debug: This part of the intermediate saving file is running.\n");

		outputRootFile = new TFile(fileName,"recreate");

//    outputRootFile->Add(fTrackPtProjectionSE);




		TH1D* VariableInfo = 0;
		VariableInfo = new TH1D("VariableInfo","VariableInfo",5,0,5);
		VariableInfo->SetBinContent(1,fObservable);
		VariableInfo->SetBinContent(2,fCent);
		VariableInfo->SetBinContent(3,fEventPlane);
		VariableInfo->SetBinContent(4,fPtMinBin);
		VariableInfo->SetBinContent(5,fPtMaxBin);
		VariableInfo->Write();

		for(Int_t i=0;i<fmaxBins;i++)
		{
			for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
			//for(Int_t j=0;j<kNvertBins;j++)
			{
				//fDetaDphi_SE[i][j]->SetDirectory(0); //FIXME now it needs to be manually deleted
				//fDetaDphi_ME[i][j]->SetDirectory(0);
				//outputRootFile->WriteObject(fDetaDphi_SE[i][j],fDetaDphi_SE[i][j]->GetName());
				//outputRootFile->WriteObject(fDetaDphi_ME[i][j],fDetaDphi_ME[i][j]->GetName());
				
				fDetaDphi_SE[i][j]->Write();
				fDetaDphi_ME[i][j]->Write();
			}
		}
    cout << "About to try saving some more intermediate (stage 1) histograms" << endl;
/*    if (fTrackPtProjectionSE) {
      printf("Saving out histogram %s (%s), which I expect to be fTrackPtProjectionSE\n",fTrackPtProjectionSE->GetName(),fTrackPtProjectionSE->GetTitle());
      fTrackPtProjectionSE->Write();
    }
    else fprintf(stderr,"TrackPtProjectionSE is not being saved\n");
*/

    if (fHistEventHash) fHistEventHash->Write();

		if (fPhase2Purity) fPhase2Purity->Write();
    if (fMCTriggerDist) fMCTriggerDist->Write();
		if (fTriggerPt) {
      printf("  Writing fTriggerPt = %s (%s)\n",fTriggerPt->GetName(),fTriggerPt->GetTitle());
      fTriggerPt->Write();
    }
    if (fTriggerPtWithinEPBin) fTriggerPtWithinEPBin->Write();
		if (fTrigger_SE[0])
		{
			for(Int_t i=0;i<fmaxBins;i++)
			{
				fTrigger_SE[i]->Write();
			}
		}

 
    if (fTrackPtProjectionSE) {
      printf("Saving out histogram %s (%s), which I expect to be fTrackPtProjectionSE\n",fTrackPtProjectionSE->GetName(),fTrackPtProjectionSE->GetTitle());
      fTrackPtProjectionSE->Write();
    }
    else fprintf(stderr,"TrackPtProjectionSE is not being saved\n");

    if (fTrackPtProjectionME) fTrackPtProjectionME->Write();
    if (fTrackPtFromTrackPsi) fTrackPtFromTrackPsi->Write();

		if (fMassPtPionAcc) fMassPtPionAcc->Write();
		if (fMassPtPionRej) fMassPtPionRej->Write();
		if (fMassPtCentPionAcc) fMassPtCentPionAcc->Write();
		if (fMassPtCentPionRej) fMassPtCentPionRej->Write();

    if (fEtaPhiPionAcc) fEtaPhiPionAcc->Write();


	}
	if(stage==2)
	{
		outputRootFile = new TFile(fileName,"UPDATE");



    printf("Trying to write out EP projections\n");
    printf("   list has size %d\n",(int)hPtEPAnglePionAcc_Proj.size());
    for (int i = 0; i < (int) hPtEPAnglePionAcc_Proj.size(); i++) hPtEPAnglePionAcc_Proj[i]->Write();

    printf("Adding the trigger EP histograms for the centrality selection\n");
    for (int i = 0; i < (int)  fPtEPAnglePionAccCent_Proj.size(); i++)  fPtEPAnglePionAccCent_Proj[i]->Write();
    for (int i = 0; i < (int) fPtEP3AnglePionAccCent_Proj.size(); i++) fPtEP3AnglePionAccCent_Proj[i]->Write();
    for (int i = 0; i < (int) fPtEP4AnglePionAccCent_Proj.size(); i++) fPtEP4AnglePionAccCent_Proj[i]->Write();

    printf("Trying to write out trigger EP th2\n");
    if (hPtEPAnglePionAcc) hPtEPAnglePionAcc->Write();

    if (fPtEPAnglePionAccCent) fPtEPAnglePionAccCent->Write();
    if (fPtEP3AnglePionAccCent) fPtEP3AnglePionAccCent->Write();
    if (fPtEP4AnglePionAccCent) fPtEP4AnglePionAccCent->Write();


    if (hHistTrackPsiEPPtCent) hHistTrackPsiEPPtCent->Write();
    if (hHistTrackPsiEP3PtCent) hHistTrackPsiEP3PtCent->Write();
    if (hHistTrackPsiEP4PtCent) hHistTrackPsiEP4PtCent->Write();
/*
    cout << "Debug I am once again trying to save fTrackPtProjectionSE, now in stage 2" << endl;
    if (fTrackPtProjectionSE) {
      printf("Saving out histogram %s (%s), which I expect to be fTrackPtProjectionSE\n",fTrackPtProjectionSE->GetName(),fTrackPtProjectionSE->GetTitle());
      fTrackPtProjectionSE->Write();
    }
    else fprintf(stderr,"TrackPtProjectionSE is not being saved\n");
*/
    TH1F * fRenorm = new TH1F("Renorm","Renorm",5,0,6);
    if (bDenormalize) {
      fRenorm->Write();
    }



		// FIXME trying out saving projections of raw SE, normalized ME
		for(Int_t i=0;i<fmaxBins;i++)
		{
			if(fDetaDphi_ME_alt1[i]) {
				TH1F * fDPhi_ME = (TH1F *) fDetaDphi_ME_alt1[i]->ProjectionX();
				fDPhi_ME->SetName(Form("DPhi_ME_%d",i));
				fDPhi_ME->SetMarkerStyle(kFullSquare);
				fDPhi_ME->SetFillColor(1);
				fDPhi_ME->SetLineStyle(1);
				fDPhi_ME->Write();
				TH1F * fDEta_ME = (TH1F *) fDetaDphi_ME_alt1[i]->ProjectionY();
				fDEta_ME->SetName(Form("DEta_ME_%d",i));
				fDEta_ME->SetMarkerStyle(kFullSquare);
				fDEta_ME->SetFillColor(1);
				fDEta_ME->SetLineStyle(1);
				fDEta_ME->Write();
			} 
		}

		//..Saving Delta eta projections
		for(Int_t i=0;i<fmaxBins;i++)
		{
			if (fDeta_Proj[i]) {
				fDeta_Proj[i]->SetTitle(Form("%s Bin %d",fObservableName.Data(),i));
				fDeta_Proj[i]->Write();
			}
		}
		for(Int_t i=0;i<fmaxBins;i++)
		{
			if (fDeta_ProjSub[i]) {
				fDeta_ProjSub[i]->SetTitle(Form("%s Bin %d",fObservableName.Data(),i));
				fDeta_ProjSub[i]->Write();
			}
		}
		for(Int_t i=0;i<fmaxBins;i++)
    {
      if (fDeta_AwaySide[i]) {
        fDeta_AwaySide[i]->SetTitle(Form("%s Bin %d",fObservableName.Data(),i));
        fDeta_AwaySide[i]->Write();
      }
    }


		//..eta width
		outputRootFile->WriteObject(fEtaWidth,fEtaWidth->GetName());
    //outputRootFile->WriteObject(fMEScaleTree,"MEScaleTree");

    // Trying to add the TTree
    printf("Trying to add the tree\n");
    outputRootFile->Add(fMEScaleTree);

		for(Int_t i=0;i<fmaxBins;i++)
		{
			//..2D corrected SE distribution
			if (fsumCorrSE[i]) outputRootFile->WriteObject(fsumCorrSE[i],fsumCorrSE[i]->GetName());
			if (fsumCorrSE_alt1[i]) outputRootFile->WriteObject(fsumCorrSE_alt1[i],fsumCorrSE_alt1[i]->GetName());
			if (fsumCorrSE_alt2[i]) outputRootFile->WriteObject(fsumCorrSE_alt2[i],fsumCorrSE_alt2[i]->GetName());
			if (fsumCorrSE_alt3[i]) outputRootFile->WriteObject(fsumCorrSE_alt3[i],fsumCorrSE_alt3[i]->GetName());
		}

		// Saving Full Projections, and adding final titles
		for (Int_t i=0;i<fmaxBins;i++)
		{
			if (fsumCorrSE_ProjFull[i]) {
//				fsumCorrSE_ProjFull[i]->SetTitle(Form("%s Bin %d",fObservableName.Data(),i));
				if(fObservable==0) fsumCorrSE_ProjFull[i]->SetTitle(Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]));
				if(fObservable==1) fsumCorrSE_ProjFull[i]->SetTitle(Form("%0.1f < z_{T}^{%s} < %0.1f",fArray_ZT_Bins[i],fTriggerName.Data(),fArray_ZT_Bins[i+1]));
//				if(fObservable==2) fsumCorrSE_ProjFull[i]->SetTitle(Form("%0.1f < #xi^{%s} < %0.1f",fArray_XI_Bins[i],fTriggerName.Data(),fArray_XI_Bins[i+1]));
				if(fObservable==2) fsumCorrSE_ProjFull[i]->SetTitle(Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]));
				fsumCorrSE_ProjFull[i]->Write();
			}
		}
		// Note: consider storing eta ranges in title	
	
		// Saving Signal Region (Near Eta) histograms
		for (Int_t i=0;i<fmaxBins;i++)
		{
			if (fsumCorrSE_NearEta[i]) {
	//			fsumCorrSE_NearEta[i]->SetTitle(Form("%s Bin %d",fObservableName.Data(),i));
				if(fObservable==0) fsumCorrSE_NearEta[i]->SetTitle(Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]));
				if(fObservable==1) fsumCorrSE_NearEta[i]->SetTitle(Form("%0.1f < z_{T}^{%s} < %0.1f",fArray_ZT_Bins[i],fTriggerName.Data(),fArray_ZT_Bins[i+1]));
//				if(fObservable==2) fsumCorrSE_NearEta[i]->SetTitle(Form("%0.1f < #xi^{%s} < %0.1f",fArray_XI_Bins[i],fTriggerName.Data(),fArray_XI_Bins[i+1]));
				if(fObservable==2) fsumCorrSE_NearEta[i]->SetTitle(Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]));
				fsumCorrSE_NearEta[i]->Write();
			}
		}

		// Saving Far Eta histograms
		for (Int_t i=0;i<fmaxBins;i++)
		{
			if (fsumCorrSE_FarEta[i]) {
	//			fsumCorrSE_FarEta[i]->SetTitle(Form("%s Bin %d",fObservableName.Data(),i));
				if(fObservable==0) fsumCorrSE_FarEta[i]->SetTitle(Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]));
				if(fObservable==1) fsumCorrSE_FarEta[i]->SetTitle(Form("%0.1f < z_{T}^{%s} < %0.1f",fArray_ZT_Bins[i],fTriggerName.Data(),fArray_ZT_Bins[i+1]));
//				if(fObservable==2) fsumCorrSE_FarEta[i]->SetTitle(Form("%0.1f < #xi^{%s} < %0.1f",fArray_XI_Bins[i],fTriggerName.Data(),fArray_XI_Bins[i+1]));
				if(fObservable==2) fsumCorrSE_FarEta[i]->SetTitle(Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]));
				fsumCorrSE_FarEta[i]->Write();
			}
		}

	//fMassPtPionAccProj[kGammaNBINS];
    printf("Trying to save out the Accepted Trigger Projections\n");
		if (fObservable == 0) {
      printf("Adding all projections\n");
			for (Int_t i = 0; i < kGammaNBINS; i++) {
				if (fMassPtPionAccProj[i]) fMassPtPionAccProj[i]->Write();
			}
			for (Int_t i = 0; i < kGammaNBINS; i++) {
				if (fMassPtPionRejProj[i]) fMassPtPionRejProj[i]->Write();
			}
		} else {
      printf("Adding projections within the range %d - %d\n",fPtMinBin,fPtMaxBin);
			for (Int_t i = fPtMinBin; i <= fPtMaxBin; i++) {
				if (fMassPtPionAccProj[i-1]) fMassPtPionAccProj[i-1]->Write();
        else printf("Projection for pt bin %d not found\n",i-1);
			}
			for (Int_t i = fPtMinBin; i <= fPtMaxBin; i++) {
				if (fMassPtPionRejProj[i-1]) fMassPtPionRejProj[i-1]->Write();
        else printf("Projection for pt bin %d not found\n",i-1);
			}
		}
		// Saving QA FFTs
		for (Int_t i=0;i<fmaxBins-1;i++)
		{
			if (fFFTsumCorrSE_ProjFull[i]) {
				if(fObservable==0) fFFTsumCorrSE_ProjFull[i]->SetTitle(Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]));
				if(fObservable==1) fFFTsumCorrSE_ProjFull[i]->SetTitle(Form("%0.1f < z_{T}^{%s} < %0.1f",fArray_ZT_Bins[i],fTriggerName.Data(),fArray_ZT_Bins[i+1]));
//				if(fObservable==2) fFFTsumCorrSE_ProjFull[i]->SetTitle(Form("%0.1f < #xi^{%s} < %0.1f",fArray_XI_Bins[i],fTriggerName.Data(),fArray_XI_Bins[i+1]));
				if(fObservable==2) fFFTsumCorrSE_ProjFull[i]->SetTitle(Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]));
        fFFTsumCorrSE_ProjFull[i]->Write();
      }
    }
    // FIXME is something going wrong here?
    outputRootFile->Write();
	}
	outputRootFile->Close();
}

void PlotGHcorrelation2::DrawWIP(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size) {

	TLegend * leg = new TLegend(x,y,x+x_size,y+y_size);
//	leg->SetHeader("ALICE");
	leg->AddEntry(Histo,"ALICE","");
	leg->AddEntry(Histo,"Work in Progress","");
	leg->SetTextSize(0.06);
	leg->SetBorderSize(0);
	leg->SetFillColorAlpha(10,0);

	leg->Draw("SAME");
}


///
/// Funtion to set TH1 histograms to a similar style
///
//________________________________________________________________________
void PlotGHcorrelation2::SetTH1Histo(TH1 *Histo,TString Xtitle,TString Ytitle,Bool_t big)
{
	Histo->SetStats(0);
	Histo->SetTitle("");
	if(big==0)
	{
		Histo->GetYaxis()->SetTitleOffset(0.8); //1.4
    if (bNoYLabel) Histo->GetYaxis()->SetTitleOffset(0.71);
		Histo->GetXaxis()->SetTitleOffset(0.9); //1.4
		Histo->GetXaxis()->SetLabelSize(0.045);
		Histo->GetYaxis()->SetLabelSize(0.04);//0.05
		Histo->GetXaxis()->SetTitleSize(0.045);
		Histo->GetYaxis()->SetTitleSize(0.045);
    Histo->SetMarkerSize(2.5);
	}
	if(big==1)
	{
		Histo->GetYaxis()->SetTitleOffset(1.0);
    if (bNoYLabel) Histo->GetYaxis()->SetTitleOffset(0.71);
		Histo->GetYaxis()->SetLabelSize(0.07); //0.07
		Histo->GetYaxis()->SetTitleSize(0.08); 
		Histo->GetXaxis()->SetTitleOffset(0.36); //0.76 //0.82
		Histo->GetXaxis()->SetLabelSize(0.05); //0.07
		Histo->GetXaxis()->SetTitleSize(0.1);
    Histo->SetMarkerSize(1.0);
	  Histo->SetMarkerStyle(kFullSquare);
    Histo->SetLineColor(kBlack);
    Histo->SetMarkerColor(kBlack);
	}
	Histo->GetXaxis()->SetNdivisions(505);
	Histo->GetYaxis()->SetNdivisions(505);
	//..make nice font
	Histo->GetXaxis()->SetLabelFont(42);
	Histo->GetYaxis()->SetLabelFont(42);
	Histo->GetXaxis()->SetTitleFont(42);
	Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitle!="")Histo->GetXaxis()->SetTitle(Xtitle);
	if(Ytitle!="") {
    if (bEnableArbUnits) Histo->GetYaxis()->SetTitle(Form("%s (Arb. Units)",Ytitle.Data()));
    else Histo->GetYaxis()->SetTitle(Ytitle);
    if (bNoYLabel) {
      Histo->GetYaxis()->SetLabelColor(kWhite);
      Histo->GetYaxis()->SetLabelSize(0.1);
      //Histo->GetYaxis()->SetTitleOffset(0.71);
      //Histo->GetYaxis()->SetTitle("");
    }
    //else Histo->GetYaxis()->SetTitle(Ytitle);
  }

//	Histo->SetLineColor(kBlack);
//	Histo->SetMarkerColor(kBlack);
//	Histo->SetMarkerStyle(20);
	//Histo->SetMarkerSize(2.5);
}
///
/// Funtion to set TH2 histograms to a similar style
///
//________________________________________________________________________
void PlotGHcorrelation2::SetTH2Histo(TH2 *Histo,TString Xtitle,TString Ytitle,TString Ztitle,Int_t rows)
{
	Histo->SetStats(0);
	Histo->SetTitle("");

	//..ideally for 3x2 canvas
	if(rows==2) // This is the main one that gets used
	{
		Histo->GetXaxis()->SetTitleOffset(1.65); //1.7
		Histo->GetYaxis()->SetTitleOffset(1.65);
		Histo->GetZaxis()->SetTitleOffset(1.1);
    if (bNoYLabel) Histo->GetZaxis()->SetTitleOffset(0.65); // 0.55 //71
		Histo->GetXaxis()->SetLabelSize(0.03); //0.05
		Histo->GetYaxis()->SetLabelSize(0.03); //0.05
		//Histo->GetZaxis()->SetLabelSize(0.05);
		Histo->GetXaxis()->SetTitleSize(0.045);
		Histo->GetYaxis()->SetTitleSize(0.045);
		Histo->GetZaxis()->SetTitleSize(0.045);
	}
	//..ideally for 3x3 canvas
	if(rows==3)
	{
		Histo->GetXaxis()->SetTitleOffset(1.7);
		Histo->GetYaxis()->SetTitleOffset(1.7);
		Histo->GetZaxis()->SetTitleOffset(1.1);
		if (bNoYLabel) Histo->GetZaxis()->SetTitleOffset(0.65); //71
		Histo->GetXaxis()->SetLabelOffset(0.01);
		Histo->GetYaxis()->SetLabelOffset(0.01);
		//Histo->GetZaxis()->SetLabelOffset(0.01);
		Histo->GetXaxis()->SetLabelSize(0.06);
		Histo->GetYaxis()->SetLabelSize(0.06);
		//Histo->GetZaxis()->SetLabelSize(0.06);
		Histo->GetXaxis()->SetTitleSize(0.055);
		Histo->GetYaxis()->SetTitleSize(0.055);
		Histo->GetZaxis()->SetTitleSize(0.055);
	}
	//..For 2D eta phi
	if(rows==4)
	{
   // Histo->SetTitleSize(0.07);

		Histo->GetYaxis()->SetTitleOffset(0.7);
		Histo->GetYaxis()->SetLabelSize(0.07);
		Histo->GetYaxis()->SetTitleSize(0.06);
//		Histo->GetYaxis()->SetTitleSize(0.1);
		Histo->GetYaxis()->SetLabelOffset(0.01);

		Histo->GetXaxis()->SetTitleOffset(0.82);
		Histo->GetXaxis()->SetLabelSize(0.07);
		Histo->GetXaxis()->SetTitleSize(0.06);
//		Histo->GetXaxis()->SetTitleSize(0.1);
		Histo->GetXaxis()->SetLabelOffset(0.01);

		Histo->GetZaxis()->SetTitleOffset(1.1);
		Histo->GetZaxis()->SetTitleSize(0.055);
    if (bNoYLabel) Histo->GetZaxis()->SetTitleOffset(0.71);
	}

	Histo->GetXaxis()->CenterTitle();
	Histo->GetYaxis()->CenterTitle();
	Histo->GetXaxis()->SetNdivisions(505);
	Histo->GetYaxis()->SetNdivisions(505);
	//make nice font
	Histo->GetXaxis()->SetLabelFont(42);
	Histo->GetYaxis()->SetLabelFont(42);
	Histo->GetXaxis()->SetTitleFont(42);
	Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitle!="")Histo->GetXaxis()->SetTitle(Xtitle);
	if(Ytitle!="")Histo->GetYaxis()->SetTitle(Ytitle);
  if(Ztitle!=""){
    if (bEnableArbUnits) Histo->GetZaxis()->SetTitle(Form("%s (Arb. Units)",Ztitle.Data()));
    else Histo->GetZaxis()->SetTitle(Ztitle);
    if (bNoYLabel) {
      Histo->GetZaxis()->SetLabelColor(kWhite); // very official way of doing this
      Histo->GetZaxis()->SetLabelSize(0.1);
     // Histo->GetYaxis()->SetTitleOffset(0.71);
      //Histo->GetZaxis()->SetTitle("");
    }
//    else Histo->GetZaxis()->SetTitle(Ztitle);
  }
	//Histo->SetLineColorAlpha(kBlue+2,0.095);
	Histo->SetLineColorAlpha(kBlack,0.095);
	Histo->SetLineWidth(0.0000);
	//	Histo->GetYaxis()->SetRangeUser(-1.3,1.3);
}
///
/// Zoom the y axis into a reasonable range. minimum-border% and maximum+border%
///
//________________________________________________________________________
void PlotGHcorrelation2::ZoomYRange(TH1D *Histo,Double_t border,Double_t Range1,Double_t Range2)
{
	Double_t lowRange=Histo->GetXaxis()->GetFirst();
	Double_t highRange=Histo->GetXaxis()->GetLast();

	Double_t min;
	Double_t max;
	Double_t maxNew;
	Double_t minNew;
	if(Range1!=-1 && Range2!=-1)Histo->GetXaxis()->SetRangeUser(Range1,Range2);
	min=Histo->GetMinimum(0);
	max=Histo->GetBinContent(Histo->GetMaximumBin());
	Double_t range=max-min;
	if(range>0)
	{
		maxNew=min+range*(1.0+2*border);
		minNew=max-range*(1.0+border);
	}
	else
	{
		maxNew=max+(-1)*range*(1.0+2*border);
		minNew=min-(-1)*range*(1.0+border);
	}
	Histo->GetYaxis()->SetRangeUser(minNew,maxNew);
	if(Range1!=-1 && Range2!=-1)Histo->GetXaxis()->SetRange(lowRange,highRange);
	//cout<<"zoom into range: "<<min<<" - "<<max<<", "<<flush;
	//cout<<minNew<<" - "<<maxNew<<endl;
}
///
/// Zoom the z axis into a reasonable range for a Dphi-Deta plot with fluctuations
/// at large Deta
//________________________________________________________________________
void PlotGHcorrelation2::DrawEtaPhi2D(TH2 *Histo, TString sZTitle)
{
	double min = 0, max = 0;
	//  double initialLimitLow = HistoGetYaxis()->GetXmin();
	//  double initialLimitHigh = HistoGetYaxis()->GetXmax();
	//SetTH2Histo(Histo,Form("#Delta#varphi^{%s-h}",fTriggerName.Data()),Form("#Delta#eta^{%s-h}",fTriggerName.Data()),Form("#frac{1}{N^{%s_{can.}}} d^{2}N^{%s_{can.}-h}/d#Delta#eta#Delta#varphi",fTriggerName.Data(),fTriggerName.Data()),2); //rows=4
	//SetTH2Histo(Histo,Form("#Delta#varphi^{%s-h}",fTriggerName.Data()),Form("#Delta#eta^{%s-h}",fTriggerName.Data()),"#frac{1}{N_{trig}} d^{2}N^{assoc}/d#Delta#eta#Delta#varphi",2); //rows=4
	//SetTH2Histo(Histo,Form("#Delta#varphi^{%s-assoc}",fTriggerName.Data()),Form("#Delta#eta^{%s-assoc}",fTriggerName.Data()),"#frac{1}{N_{trig}} d^{2}N^{assoc}/d#Delta#eta#Delta#varphi",2); //rows=4

  TString sZaxisTitle="#frac{1}{N_{trig}} d^{2}N^{assoc}/d#Delta#eta d#Delta#varphi";
  if (sZTitle != "") {
    sZaxisTitle=sZTitle;
  }
	SetTH2Histo(Histo,"#Delta#varphi","#Delta#eta",sZaxisTitle.Data(),2); //rows=4
	Histo->GetYaxis()->SetRangeUser(-fDetaLimit,fDetaLimit);
//	Histo->GetYaxis()->SetRange(1,Histo->GetYaxis()->GetNbins());
	Histo->GetYaxis()->SetRangeUser(-fMaxDeltaEtaPlotRange,fMaxDeltaEtaPlotRange);
	min = Histo->GetBinContent(Histo->GetMinimumBin());
	max = Histo->GetBinContent(Histo->GetMaximumBin());

  Histo->SetLineWidth(1);

	//Histo->GetZaxis()->SetRangeUser(min,max);
	Histo->DrawCopy(fPlotOptions.Data());
	Histo->GetZaxis()->SetRangeUser(min,max);
}

///
///
/// Set a nice color sceme for plotting 2D ms
//________________________________________________________________________
void PlotGHcorrelation2::SetPlotStyle()
{
	const Int_t NRGBs = 5;
	const Int_t NCont = 99;//max possible?

	//Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 }; //log??
	//Double_t stops[NRGBs] = { 0.00, 0.18, 0.34, 0.61, 1.00 }; //log ELI test??
	//Double_t stops[NRGBs] = { 0.00, 0.25, 0.5, 0.75, 1.00 };  //linear
	Double_t stops[NRGBs] = { 0.00, 0.16, 0.39, 0.66, 1.00 };  //squeezing high and
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}
///
/// Plot a latex text into a histogram to add additonal infomation
///
//________________________________________________________________________
TLatex* PlotGHcorrelation2::PlotTopLegend(const char* label,Float_t x,Float_t y,Float_t size,Int_t color,Float_t angle)
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
///
/// 2 Gauss functions plus a polynomial
///
//________________________________________________________________________
Double_t PlotGHcorrelation2::PolyTwoGaussFitFunc(Double_t* x_val, Double_t* par)
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

	y = par0*(TMath::Gaus(x,CommomMean,par2,0))+par3*(TMath::Gaus(x,CommomMean,par5,0));
	y+=par6; // constant background
	//y+=(par6+par7*x+par8*x*x+par9*x*x*x+par10*x*x*x*x); // poly background
	return y;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TF1* PlotGHcorrelation2::allFitFuncVn(TString name,Int_t VnTerms,Double_t allPhiMin,Double_t allPhiMax)
{
	TF1* f1 = 0;
	if(VnTerms == 10)  f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 16);
	if(VnTerms == 7)   f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 11);
	if(VnTerms == 5)   f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 8);
	if(VnTerms == 4)   f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 7);
	if(VnTerms == 3)   f1 = new TF1 (name,JoelsVnFunctionValue, allPhiMin,allPhiMax, 5);

	return f1;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Double_t PlotGHcorrelation2::JoelsVnFunctionValue(const double * x, const double * p)
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
	//double phi = x[0]*TMath::Pi()/180.0;  //transform from deg to rad
	double phi = x[0];  //transform from deg to rad


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
	Part1  = 2.0*v1*v3*TMath::Cos(phi);//..let v1 be smaller than v3
	result = B*(1 + Part1 + Part2 + Part3 + Part4 + Part5 + Part6 + Part7 + Part8 + Part9 + Part10);

	return result;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void PlotGHcorrelation2::SetJoelsParLimits(TF1 *func, TH1 *histoToFit,Double_t par_V10[])
{
	par_V10[0] = (histoToFit->GetMinimum()+histoToFit->GetMaximum());
	func->SetParameters(&par_V10[0]);
	func->SetParNames("B","v1", "v2jet", "v2assoc", "v3", "v4jet", "v4assoc", "v5", "v6jet", "v6assoc", "v7");
	func->SetParName(11, "v8jet");
	func->SetParName(12, "v8assoc");
	func->SetParName(13, "v9");
	func->SetParName(14, "v10jet");
	func->SetParName(15, "v10assoc");

	for(Int_t g = 0; g < 16; g++)
	{
		func->ReleaseParameter(g);
		func->SetParameter(g,1e-4);
		if(g==0)func->SetParameter(g,1e+06); //!!! change when proper normalized
		func->SetParError(g,0.0);
	}

	// - completely arbirtary - set whats right for you
	func->SetParLimits(0,1e-6,histoToFit->GetMaximum());
	func->SetParLimits(1,1e-6,1e-5); //v1 should be smaller than v3;
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
//-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
Double_t PlotGHcorrelation2::FlowFunction(Double_t* x_val, Double_t* par)
{
	Double_t x, y, par0, par1, par2, par3;
	par0  = par[0]; //Background level
	par1  = par[1]; //v1
	par2  = par[2]; //v2
	par3  = par[3]; //0-1 (v3 smaller than v2)

	//x = x_val[0]*TMath::Pi()/180.0;  //transform from deg to rad
	x = x_val[0];  //transform from deg to rad

	y = par0*(1+2*par1*par2*par3*TMath::Cos(x)+2*par2*TMath::Cos(2*x)+2*par2*par3*TMath::Cos(3*x));
	return y;
}
///
/// A vertical line that can be plotted
///
//________________________________________________________________________
void PlotGHcorrelation2::PlotVerLineRange(Double_t x_val, Double_t yLow, Double_t yHigh, Int_t Line_Col)
{
	TLine* Zero_line = new TLine();
	Zero_line -> SetX1(x_val);
	Zero_line -> SetX2(x_val);
	Zero_line -> SetY1(yLow);
	Zero_line -> SetY2(yHigh);
	//cout << "x_val = " << x_val << ", Bin = " << Histo->FindBin(x_val) << ", Y2 = " << Histo->GetBinContent(Histo->FindBin(x_val)) << endl;
	Zero_line -> SetLineWidth(2);
	Zero_line -> SetLineStyle(2);
	Zero_line -> SetLineColor(Line_Col);
	Zero_line -> Draw();
	//delete Zero_line;
}
///
/// A vertical line that can be plotted
///
//________________________________________________________________________
void PlotGHcorrelation2::PlotVerLine3(Double_t x_val,TH1* Histo, Double_t y_fac, Int_t Line_Col)
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
///_______________________________________________________________________
///_______________________________________________________________________
///
/// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Here the plotting and histo manipulation functions start- - - - - - -
/// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
///________________________________________________________________________
///________________________________________________________________________

///
/// Draw Summary of SE plots
///
//________________________________________________________________________
void PlotGHcorrelation2::DrawSEplots()
{
	//..this function can be used to
	//..plot all SE bins together
	//..summed over all vertex bins to
	//..see the raw signal for each pT,zT or Xi bin

}
/// Analyzes the Pion Accepted, Rejected Histograms
///
//________________________________________________________________________
void PlotGHcorrelation2::AnalyzePionAccRej()
{

  // FIXME use the multicentrality one
	if(!fMassPtPionAcc || !fMassPtPionRej) {
		fprintf(stderr,"Missing Pion Accepted or Rejected Histograms; skipping.\n");
		return;
	}
	TH2D * fMassPtPionAll = (TH2D*) fMassPtPionRej->Clone("fMassPtPionAll");
	fMassPtPionAll->Add(fMassPtPionAcc);

	float massRangeMin = 0.05;
	float massRangeMax = 0.6;
	float ptRangeMin   = 4.;
	float ptRangeMax   = 22.;

//	fMassPtPionAcc->GetXaxis()->SetRangeUser(massRangeMin,massRangeMax);
	fMassPtPionAcc->GetYaxis()->SetRangeUser(ptRangeMin,ptRangeMax);
//	fMassPtPionRej->GetXaxis()->SetRangeUser(massRangeMin,massRangeMax);
	fMassPtPionRej->GetYaxis()->SetRangeUser(ptRangeMin,ptRangeMax);

  TCanvas * fCanvAccRej = new TCanvas("Canv_Acc_Rej","Accepted vs Rejected Pion Candidates",1400,700);
  if (fProducePi0AnalysisPlots) {

    fCanvAccRej->Divide(2,1);
    fCanvAccRej->cd(1);
    fMassPtPionAcc->Draw("COLZ");
    fCanvAccRej->cd(2);
    fMassPtPionRej->Draw("COLZ");

    fCanvAccRej->Print(Form("%s/%s.png",fOutputDir.Data(),fCanvAccRej->GetName()));
    fCanvAccRej->Print(Form("%s/CFiles/%s.C",fOutputDir.Data(),fCanvAccRej->GetName()));

    // Drawing Combined plot with lines
    fCanvAccRej->Clear();
  }
/*
	fMassPtPionAll->Rebin2D(5,2);
	fMassPtPionAll->GetXaxis()->SetRangeUser(0.,0.6);
	fMassPtPionAll->GetYaxis()->SetRangeUser(4.0,22.0);
	fMassPtPionAll->Draw(fPlotOptions.Data());
	// Copied and pasted
	Double_t fPi0MassFitParsValue[5] = {10.49,0.13852,-1.17e-4,2.861e-3,0};
	Double_t fPi0SigmaFitParsValue[5] = {8.34,9.90e-3,-1.09e-4,6.86e-4,0};

	float n_Sigma = 2;

	double d_m = fPi0MassFitParsValue[0];
	double e_m = fPi0MassFitParsValue[1];
	double m1_m = fPi0MassFitParsValue[2];
	double m2_m = fPi0MassFitParsValue[3];

	double d_s = fPi0SigmaFitParsValue[0];
	double e_s = fPi0SigmaFitParsValue[1];
	double m1_s = fPi0SigmaFitParsValue[2];
	double m2_s = fPi0SigmaFitParsValue[3];

	double x1_l,x2_l,x3_l,x4_l; // low mass points
	double x1_h,x2_h,x3_h,x4_h; // high mass points
	double y1,y2,y3,y4; // pt points

	y1 = 5.;
	y2 = min(d_m,d_s);
	y3 = max(d_m,d_s);
	y4 = ptRangeMax;

	double mass_1,mass_2,mass_3,mass_4;
	double sigma_1,sigma_2,sigma_3,sigma_4;

	mass_1 = (m1_m * y1 + e_m - m1_m * d_m);
	sigma_1 = (m1_s * y1 + e_s - m1_s * d_s);

	if (d_m < d_s) { // mass kink before sigma kink
		y2 = d_m;
		y3 = d_s;

		mass_2  = e_m;
		sigma_2 = m1_s * y2 + e_s - m1_s * d_s;
		mass_3  = m2_m * y3 + e_m - m2_m * d_m;
		sigma_3 = e_s;

	} else {
		y2 = d_s;
		y3 = d_m;

		mass_2  = m1_m * y2 + e_m - m1_m * d_m;
		sigma_2 = e_s;
		mass_3  = e_m;
		sigma_3 = m2_s * y3 + e_s - m2_s * d_s;

	}

	mass_4 = m2_m * y4 + e_m - m2_m * d_m;
	sigma_4 = m2_s * y4 + e_s - m2_s * d_s;

	x1_l = mass_1 - n_Sigma * sigma_1;
	x2_l = mass_2 - n_Sigma * sigma_2;
	x3_l = mass_3 - n_Sigma * sigma_3;
	x4_l = mass_4 - n_Sigma * sigma_4;

	x1_h = mass_1 + n_Sigma * sigma_1;
	x2_h = mass_2 + n_Sigma * sigma_2;
	x3_h = mass_3 + n_Sigma * sigma_3;
	x4_h = mass_4 + n_Sigma * sigma_4;

	TLine * linesLow[3];
	TLine * linesHigh[3];

	linesLow[0] = new TLine(x1_l,y1,x2_l,y2);
	linesLow[1] = new TLine(x2_l,y2,x3_l,y3);
	linesLow[2] = new TLine(x3_l,y3,x4_l,y4);

	linesHigh[0] = new TLine(x1_h,y1,x2_h,y2);
	linesHigh[1] = new TLine(x2_h,y2,x3_h,y3);
	linesHigh[2] = new TLine(x3_h,y3,x4_h,y4);

	for (int i = 0; i < 3; i++) {
		linesLow[i]->SetLineWidth(3);
		linesHigh[i]->SetLineWidth(3);

		linesLow[i]->Draw("SAME");
		linesHigh[i]->Draw("SAME");
	}
*/
//	fCanvAccRej->Print(Form("%s/Pion_Acc_Rej_Map.png",fOutputDir.Data()));
//	fCanvAccRej->Print(Form("%s/CFiles/Pion_Acc_Rej_Map.C",fOutputDir.Data()));

//	TH1D * fMassPtPionRejProj[kGammaNBINS];
//	TH1D * fMassPtPionAccProj[kGammaNBINS];
//	TH1D * fMassPtPionRejProj[kGammaNBINS];

	fCanvAccRej->Clear();
	//fCanvAccRej->Divide(3,3,0.001,0.001);
	fCanvAccRej->Divide(3,3,0,0);

	int nRebinAccRej=10;

  //vector<TPaveText *> fIntegralBoxes = {};
  TPaveText * fIntegralBoxes[kGammaNBINS];

	for (int i = 0; i < kGammaNBINS-1; i++) {
		fCanvAccRej->cd(i+1);
		double ptLow = fArray_G_Bins[i];
		double ptHigh = fArray_G_Bins[i+1];


		//FIXME add centrality info here.
    // fMassPtCentPionAcc
	/*	int ptBinLow = fMassPtPionAcc->GetYaxis()->FindBin(ptLow);
		int ptBinHigh = fMassPtPionAcc->GetYaxis()->FindBin(ptHigh);
    // FIXME does this need to be done with centrality in mind?
		fMassPtPionAccProj[i] = (TH1D * ) fMassPtPionAcc->ProjectionX(Form("fMassPtPionAccProj_%d",i),ptBinLow,ptBinHigh-1);
		fMassPtPionRejProj[i] = (TH1D * ) fMassPtPionRej->ProjectionX(Form("fMassPtPionRejProj_%d",i),ptBinLow,ptBinHigh-1); */
    // fCent
    int MinCentBin = fCent + 1;
    int MaxCentBin = fCent + 1;
    if (fCent < 0) {
      MinCentBin = 1;
      MaxCentBin = kCentBINS;
    }
		int ptBinLow = fMassPtCentPionAcc->GetYaxis()->FindBin(ptLow);
		int ptBinHigh = fMassPtCentPionAcc->GetYaxis()->FindBin(ptHigh);
    fMassPtCentPionAcc->GetZaxis()->SetRange(MinCentBin,MaxCentBin);
    fMassPtCentPionRej->GetZaxis()->SetRange(MinCentBin,MaxCentBin);
    printf("  Debug: MassPt Cent Range set to bins %d - %d (in my bins %d %d).\n",MinCentBin,MaxCentBin,MinCentBin-1,MaxCentBin-1);

    // FIXME double check the logic here

    // FIXME does this need to be done with centrality in mind?
//		fMassPtPionAccProj[i] = (TH1D * ) fMassPtPionAcc->ProjectionX(Form("fMassPtPionAccProj_%d",i),ptBinLow,ptBinHigh-1);
//		fMassPtPionRejProj[i] = (TH1D * ) fMassPtPionRej->ProjectionX(Form("fMassPtPionRejProj_%d",i),ptBinLow,ptBinHigh-1);
    TAxis * fMassPtAcc_PtAxis = fMassPtCentPionAcc->GetYaxis();
    double fMassPt_PtBinLowEdge = fMassPtAcc_PtAxis->GetBinLowEdge(ptBinLow);
    double fMassPt_PtBinUpEdge  = fMassPtAcc_PtAxis->GetBinUpEdge(ptBinHigh - 1);
    printf("  Debug: MassPt Set to bin range %d -%d (Pt Range %f - %f)\n",ptBinLow,ptBinHigh-1,fMassPt_PtBinLowEdge,fMassPt_PtBinUpEdge);

    fMassPtCentPionAcc->GetYaxis()->SetRange(ptBinLow,ptBinHigh-1);
    fMassPtCentPionRej->GetYaxis()->SetRange(ptBinLow,ptBinHigh-1);

		fMassPtPionAccProj[i] = (TH1D * ) fMassPtCentPionAcc->Project3D("x");
    fMassPtPionAccProj[i]->SetName(Form("fMassPtPionAccProj_%d",i));
		fMassPtPionRejProj[i] = (TH1D * ) fMassPtCentPionRej->Project3D("x");
    fMassPtPionRejProj[i]->SetName(Form("fMassPtPionRejProj_%d",i));

    double fIntegral = fMassPtPionAccProj[i]->Integral();
    printf("  Debug:  Found MassPt Acc Proj Integral = %f\n",fIntegral);
    TPaveText * fIntegralBox = new TPaveText(0.5,0.33,0.9,0.58,"NDC");
    fIntegralBox->AddText(Form("Trigger Count = %.1f",fIntegral));

    fIntegralBox->SetName(Form("PaveText_%d",i));
    fIntegralBoxes[i]=fIntegralBox;

		fMassPtPionAccProj[i]->Rebin(nRebinAccRej); //could rebin th2 instead
		fMassPtPionRejProj[i]->Rebin(nRebinAccRej);

		fMassPtPionAccProj[i]->SetLineColor(46);
		fMassPtPionRejProj[i]->SetLineColor(kBlack);
		fMassPtPionAccProj[i]->SetTitle(Form("%.0f #leq #it{p}_{T} < %.0f GeV/c",ptLow,ptHigh));
		fMassPtPionRejProj[i]->SetTitle(Form("%.0f #leq #it{p}_{T} < %.0f GeV/c",ptLow,ptHigh));

		fMassPtPionAccProj[i]->GetXaxis()->SetRangeUser(.05,.5);

    bool hasEntries= 0 < fMassPtPionAccProj[i]->GetEntries();

    if (hasEntries) {
			if (fGammaOrPi0 > 1) { // The sidebands should be drawn after rejects, since rejects include the pi0 peak
				fMassPtPionRejProj[i]->Draw();
				fMassPtPionAccProj[i]->Draw("SAME");
			} else {
				fMassPtPionAccProj[i]->Draw();
				fMassPtPionRejProj[i]->Draw("SAME");
				fMassPtPionAccProj[i]->Draw("SAME");
			}
    }


		PlotTopLegend(fMassPtPionAccProj[i]->GetTitle(),0.5,0.77,.1,kBlack);
		//PlotTopLegend(fMassPtPionAccProj[i]->GetTitle(),0.3,0.77,.15,kBlack);

    fIntegralBoxes[i]->Draw();
//		if (i == kGammaNBINS - 1) {
//			TLegend * legAccRej = new TLegend(.2,.1,.95,.90);
//			legAccRej->AddEntry(fMassPtPionAccProj[i],"Accepted","lp");
//			legAccRej->AddEntry(fMassPtPionRejProj[i],"Rejected","lp");
//			legAccRej->Draw("SAME");
//		}

	}
	fCanvAccRej->cd(kGammaNBINS);
	TLegend * legAccRej = new TLegend(.20,.25,.80,.75);
	legAccRej->AddEntry(fMassPtPionAccProj[0],"Accepted","lp");
	legAccRej->AddEntry(fMassPtPionRejProj[0],"Rejected","lp");
	legAccRej->Draw("SAME");

	fCanvAccRej->Print(Form("%s/Pion_Acc_Rej_Mass.png",fOutputDir.Data()));
	fCanvAccRej->Print(Form("%s/CFiles/Pion_Acc_Rej_Mass.C",fOutputDir.Data()));

	// Version With different centralities:
//	fCanvAccRej->Clear();

	int kCentColor[kCentBINS] = {40,41,42,45};
	TH1D * fMassPtCentPionAccProj[kGammaNBINS][kCentBINS];
	TH1D * fMassPtCentPionRejProj[kGammaNBINS][kCentBINS];
	for (int i = 0; i < kGammaNBINS-1; i++) {
		fCanvAccRej->cd(i+1);
		gPad->Clear();
		double ptLow = fArray_G_Bins[i];
		double ptHigh = fArray_G_Bins[i+1];
		int ptBinLow = fMassPtCentPionAcc->GetYaxis()->FindBin(ptLow);
		int ptBinHigh = fMassPtCentPionAcc->GetYaxis()->FindBin(ptHigh);

		for (int j = 0; j < kCentBINS; j++) {
			fMassPtCentPionAccProj[i][j] = (TH1D *) fMassPtCentPionAcc->ProjectionX(Form("fMassPtCentPionAccProj_%d_%d",j,i),ptBinLow,ptBinHigh-1,j+1,j+1);
			fMassPtCentPionAccProj[i][j]->SetTitle(Form("%.0f #leq #it{p}_{T} < %.0f GeV/c",ptLow,ptHigh));
			fMassPtCentPionRejProj[i][j] = (TH1D *) fMassPtCentPionRej->ProjectionX(Form("fMassPtCentPionRejProj_%d_%d",j,i),ptBinLow,ptBinHigh-1,j+1,j+1);
			// Drawing all
			fMassPtCentPionAccProj[i][j]->Sumw2();
			fMassPtCentPionRejProj[i][j]->Sumw2();

			fMassPtCentPionAccProj[i][j]->Rebin(nRebinAccRej);
			fMassPtCentPionRejProj[i][j]->Rebin(nRebinAccRej);

			fMassPtCentPionAccProj[i][j]->Add(fMassPtCentPionRejProj[i][j]);

			fMassPtCentPionAccProj[i][j]->SetLineColor(kCentColor[j]);
			fMassPtCentPionAccProj[i][j]->SetMarkerColor(kCentColor[j]);
			fMassPtCentPionAccProj[i][j]->SetMarkerStyle(kFullSquare);
		}
		// Scaling for Cent bin sizes
		fMassPtCentPionAccProj[i][0]->Scale(2.);

		fMassPtCentPionAccProj[i][0]->Draw("HIST");
		for (int j = 1; j < kCentBINS; j++) {
			fMassPtCentPionAccProj[i][j]->Draw("SAME HIST");
		} 
		PlotTopLegend(fMassPtCentPionAccProj[i][0]->GetTitle(),0.5,0.77,.1,kBlack);
	}
	fCanvAccRej->cd(kGammaNBINS);
	gPad->Clear();
	legAccRej->Clear();
	for (int j = 0; j < kCentBINS; j++) {
		if (j == 0) 
		legAccRej->AddEntry(fMassPtCentPionAccProj[0][j],Form("Cent. %.0f%% - %.0f%% (#times2)",fArray_cent_Bins[j],fArray_cent_Bins[j+1]),"lp");	
		else 
		legAccRej->AddEntry(fMassPtCentPionAccProj[0][j],Form("Cent. %.0f%% - %.0f%%",fArray_cent_Bins[j],fArray_cent_Bins[j+1]),"lp");	
	}
	legAccRej->Draw();
	fCanvAccRej->Print(Form("%s/Pion_Cent_Mass.png",fOutputDir.Data()));
	fCanvAccRej->Print(Form("%s/CFiles/Pion_Cent_Mass.C",fOutputDir.Data()));


}

/// In the case of poor statistic - for the moment - merge different ME histograms
/// for different bins
/// !!!!!! This needs more work
///
//________________________________________________________________________
void PlotGHcorrelation2::MergeMEplots()
{
	//..There are different way how to
	//..summ of ME histograms

	//..Way 1: add all z-vertex bins to improve statistic in the ME histogram
	//..Way 2: add all specific z-vertex bins to improve statistic in the ME histogram
	//..Way 3: add all ME histograms in many fObservable bins to improve statistic in the ME histogram
	for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
	{
		for (Int_t i = 0; i < fmaxBins; i++)
		{
			//..Add all ME for different z-vertexes together (Way1)
			fDetaDphi_ME_alt1[i]->Add(fDetaDphi_ME[i][j]);
			fDetaDphi_SE_alt1[i]->Add(fDetaDphi_SE[i][j]);

			//..Add all ME for different fObservables together (Way3)
			//..take a specific z-vertex bin
			fDetaDphi_ME_alt3[j]->Add(fDetaDphi_ME[i][j]);
			//.. Don't add SE for version 3 because we want to keep the individual fObservable bins
		}
	}
	//..very specific merging, add 4 bin on the outside and 2 on the inside(around 0) zVtx (Way2)
	for (Int_t i = 0; i < fmaxBins; i++)
	{
		Int_t bin=0;
		for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
		{
			if(fArray_zVtx_Bins[j]<fArray_zVtx_BinsAlt[bin+1] && fArray_zVtx_Bins[j]!=fArray_zVtx_BinsAlt[bin+1])
			{
				fDetaDphi_ME_alt2[i][bin]->Add(fDetaDphi_ME[i][j]);
				fDetaDphi_SE_alt2[i][bin]->Add(fDetaDphi_SE[i][j]);
			}
			else
			{
				bin++;
				fDetaDphi_ME_alt2[i][bin]->Add(fDetaDphi_ME[i][j]);
				fDetaDphi_SE_alt2[i][bin]->Add(fDetaDphi_SE[i][j]);
			}
		}
	}

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  Plot
	//..  The new merged SE distributions
	//. . . . . . . . . . . . . . . . . . . . . . . .

	TLegend *legME1;
	TLegend *legME2;
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fRaw_Plots_2D_alt2[i]->Divide(4,3,0.001,0.001);
		fRaw_Plots_2D_alt3[i]->Divide(4,3,0.001,0.001);

		for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
		//for(Int_t j=0;j<kNvertBins;j++)
		{
			if(j==0)legME1 = new TLegend(0.2,0.72,0.4,0.9);
			else    legME1 = new TLegend(0.2,0.85,0.4,0.9);
			legME1->AddEntry(fDetaDphi_ME[i][j],Form("%.0f < z_{vtx} < %.0f",fArray_zVtx_Bins[j],fArray_zVtx_Bins[j+1]),"");
			if(j==0 && fObservable==0)legME1->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.0f<#it{p}_{T}<%0.0f",fArray_G_Bins[i],fArray_G_Bins[i+1]),"");
			if(j==0 && fObservable==1)legME1->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//			if(j==0 && fObservable==2)legME1->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
			if(j==0 && fObservable==2)legME1->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#it{p}_{T}^{assoc}<%0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
			legME1->SetTextColor(kBlack);
			legME1->SetTextSize(0.075);
			legME1->SetBorderSize(0);
			legME1->SetFillColorAlpha(10, 0);

			if(j==0)legME2 = new TLegend(0.2,0.72,0.4,0.9);
			else    legME2 = new TLegend(0.2,0.85,0.4,0.9);
			if(j<8)legME2->AddEntry(fDetaDphi_ME[i][j],Form("%.0f < z_{vtx} < %.0f",fArray_zVtx_BinsAlt[j],fArray_zVtx_BinsAlt[j+1]),"");
			if(j==0 && fObservable==0)legME2->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.0f<#it{p}_{T}<%0.0f",fArray_G_Bins[i],fArray_G_Bins[i+1]),"");
			if(j==0 && fObservable==1)legME2->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//			if(j==0 && fObservable==2)legME2->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
			if(j==0 && fObservable==2)legME2->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#it{p}_{T}^{assoc}<%0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
			legME2->SetTextColor(kBlack);
			legME2->SetTextSize(0.075);
			legME2->SetBorderSize(0);
			legME2->SetFillColorAlpha(10, 0);

			//..plot the scaled 2D ME distribution
			//fMEPlots2DGamma[i]->cd(j+1);
			//DrawEtaPhi2D(fDetaDphi_ME[i][j]);
			//legME->Draw("same");

			if(j<8)
			{
				//..plot the 2D SE distribution
				fRaw_Plots_2D_alt2[i]->cd(j+1);
				DrawEtaPhi2D(fDetaDphi_SE_alt2[i][j]);
				legME2->Draw("same");
			}
			//..plot the 2D SE distribution
			fRaw_Plots_2D_alt3[i]->cd(j+1);
			DrawEtaPhi2D(fDetaDphi_SE_alt3[i][j]);
			legME1->Draw("same");
		}
		fRaw_Plots_2D_alt2[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fRaw_Plots_2D_alt2[i]->GetName()));
		fRaw_Plots_2D_alt2[i]->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fRaw_Plots_2D_alt2[i]->GetName()));
		fRaw_Plots_2D_alt3[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fRaw_Plots_2D_alt3[i]->GetName()));
		fRaw_Plots_2D_alt3[i]->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fRaw_Plots_2D_alt3[i]->GetName()));
	}
	//..Then you have to Normalize the Corrected SE plots. This is done in PlotCorrected2DHistograms()
}
/// **
/// Determine the normalization for the same event histograms
/// see how many triggers were in the collected statistic
/// to determina a normalization factor
/// **
//________________________________________________________________________
void PlotGHcorrelation2::NormalizeSEsignal()
{

}
///**
///Scale the ME background such that the plateau around 0,0 is at a height of 1.
///**
//________________________________________________________________________
void PlotGHcorrelation2::NormalizeMEsignal()
{
	TCanvas* DummyCan;
	// Place a switch here for the merged z-vertex bins?
	for(Int_t i=0;i<fmaxBins;i++)
	{
		DummyCan = new TCanvas();
		fMEPlots1DGamma[i]->Divide(4,3,0.001,0.001);
		DummyCan->Divide(4,3,0.001,0.001);

		for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
		//for(Int_t j=0;j<kNvertBins;j++)
		{
			//..project the 2D distribution in a narrow eta window (2nd and 3rd argument)
			//..and scaled the plateau to one
			ScaleMEbackground(fDetaDphi_ME[i][j],-fMEDEtaRangeForNorm,fMEDEtaRangeForNorm,fMEPlots1DGamma[i],j);

			if(fPlotMoreMEstrategies==1)
			{
				//..Scale also the ME background of the other 3 versions
				if(j==1)ScaleMEbackground(fDetaDphi_ME_alt1[i],-fMEDEtaRangeForNorm,fMEDEtaRangeForNorm,DummyCan,0);
				if(j<kNvertBins_alt)ScaleMEbackground(fDetaDphi_ME_alt2[i][j],-fMEDEtaRangeForNorm,fMEDEtaRangeForNorm,DummyCan,0);
				if(i==0)ScaleMEbackground(fDetaDphi_ME_alt3[j],-fMEDEtaRangeForNorm,fMEDEtaRangeForNorm,DummyCan,0);
			}
		}
		//..since it is normalized to 1 in a single z-Vtx bin
		//..adding kNvertBins of them demands dividing/scaling by kNvertBins.
		fMEPlots1DGamma[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fMEPlots1DGamma[i]->GetName()));
		fMEPlots1DGamma[i]->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fMEPlots1DGamma[i]->GetName()));
	}
}
///**
///Project and fit the correlation function to
///determine the exact scaling factor in order to scale the plateau around 0,0 to 1.
///**
//________________________________________________________________________
void PlotGHcorrelation2::ScaleMEbackground(TH2D* Histo,Double_t lowRange,Double_t highRange,TCanvas* Can,Int_t CanvasPad)
{
	//..project to the delta phi-axis. But only around y=0 (in range lowRange-highRange (high and low range are included [thats why -1 is subtracted])).
	TString ProjectionName;
	ProjectionName= Histo->GetName();
	ProjectionName+="_projX_range";

  double epsilon = 0.0001;

  int iLowRangeBin = Histo->GetYaxis()->FindBin(lowRange+epsilon);
  int iHighRangeBin = Histo->GetYaxis()->FindBin(highRange-epsilon);

	//TH1D *PprojX=Histo->ProjectionX((const char*)ProjectionName,iLowRangeBin,iHighRangeBin-1);
	TH1D *PprojX=Histo->ProjectionX((const char*)ProjectionName,iLowRangeBin,iHighRangeBin);
	//TH1D *PprojX=Histo->ProjectionX((const char*)ProjectionName,Histo->GetYaxis()->FindBin(lowRange),Histo->GetYaxis()->FindBin(highRange)-1);

//  Int_t nRebinMEForNorm = 5;
//  nRebinMEForNorm = 2; //FIXME test
  PprojX->Rebin(nRebinMEForNorm);

	SetTH1Histo(PprojX,Form("#Delta#varphi^{%s-h}",fTriggerName.Data()),Form("dN^{%s-h}",fTriggerName.Data()),1);
	ZoomYRange(PprojX,0.5);
	Can->cd(CanvasPad+1);
  PprojX->GetYaxis()->SetLabelSize(0.04);
  PprojX->GetYaxis()->SetTitleSize(0.04);
	PprojX->DrawCopy("E");

	//..fit the projection with a flat line to
	//..determine the plateau
	if(PprojX->GetEntries()==0.0)
	{
		printf("Empty Histogram, skipping\n");
		return;
	}
	Double_t fLineFitMin;// = (180. - 50.)*DTR;
	Double_t fLineFitMax;// = (180. + 50.)*DTR;


  Int_t nProjBinsX=PprojX->GetNbinsX();
	//..Determine max five bins in histogram
	Int_t startBin=0;
  //Int_t nMENormSideBins = 2; // 2 + 1 + 2 is our window size
  //nMENormSideBins = 3; // FIXME test

	Double_t sum;
	Double_t sumMax=0;
  // recall that index: 0 1 2    nBins nBins+1   (total of nBins+2)
  //               bin: u 1 2    nBins o
  
  // FIXME should this be i = nMENormSideBins + 1, since 0 is underflow
	//for(Int_t i=1+nMENormSideBins;i<=(nProjBinsX-nMENormSideBins);i++)
	for(Int_t i=1;i<=nProjBinsX;i++)
	{
		sum=0;
    //printf("Debug: ");
    for (Int_t j = i - nMENormSideBins; j <= i + nMENormSideBins; j++) {
      // want to map 0 to nPRojBinsX - 1, -1 to nProjBinsX - 2, nProjBinsX to nProjBinsX, nProjBinsX+1 to 1, nProjBinsX+2 to 2
      int localIndex = j % nProjBinsX; // ?. Want nProjBinsX to map to nProjBinsX
      if (localIndex <= 0) localIndex += nProjBinsX;
      //printf("%d, ",localIndex);
      // if ( i < nMENormSideBins or something)
      // Or apply modulus of nProjBinsX
      //sum+=PprojX->GetBinContent(j);
      if (localIndex >= 0 && localIndex <= nProjBinsX+1) {
        sum+=PprojX->GetBinContent(localIndex);
      }
    }
    //printf("\n");
		//sum+=PprojX->GetBinContent(i-2);
		//sum+=PprojX->GetBinContent(i-1);
		//sum+=PprojX->GetBinContent(i);
		//sum+=PprojX->GetBinContent(i+1);
		//sum+=PprojX->GetBinContent(i+2);
		if(sum>sumMax)
		{
			sumMax=sum;
			startBin=i;
		}
	}
  fMEIndexOfNormRegion = startBin;
	Double_t width= 0.5*PprojX->GetBinWidth(startBin);
	fLineFitMin   = PprojX->GetBinCenter(startBin-nMENormSideBins)-width;
	fLineFitMax   = PprojX->GetBinCenter(startBin+nMENormSideBins)+width;

  // Special cases when startBin - nMENormSideBins <= 0
  //               when startBin + nMENormSideBins >= nProjBinX

  TF1* LinFit = 0;
  TH1D * PprojX_Clone = 0;

  PprojX_Clone = (TH1D *) PprojX->Clone("PprojX_Clone");
  PprojX_Clone->SetLineColor(kGray);
  PprojX_Clone->SetMarkerColor(kGray);

//  printf("  startBin = %d, region: [%d,%d].  nProjBinsX = %d\n",startBin,startBin - nMENormSideBins,startBin + nMENormSideBins,nProjBinsX);
  if (startBin - nMENormSideBins <= 0 || startBin + nMENormSideBins > nProjBinsX) {
//    printf("   Found overflow,underflow of region\n");
//    PprojX_Clone->Draw("SAME");
 //   PprojX->Draw("SAME");


    LinFit = new TF1("pol0","pol0",fLineFitMin,fLineFitMax);
    LinFit->SetLineColor(kCyan); // FIXME for debugging 

    // Zero out the region for PprojX, but clone first, then draw clone  
    int zeroRegionMin = 1;
    int zeroRegionMax = nProjBinsX;
    if (startBin - nMENormSideBins <= 0) {
//      printf("  underflow ...\n");
      zeroRegionMin = startBin+nMENormSideBins+1;
      zeroRegionMax = startBin-nMENormSideBins+nProjBinsX-1;
      //zeroRegion(PprojX,startBin+nMENormSideBins+1,startBin-nMENormSideBins+nProjBinsX-1+1);
      LinFit->SetLineStyle(4); //FIXME for debugging
    }
    else {
//      printf("   overflow ...\n");
      zeroRegionMin = startBin+nMENormSideBins-nProjBinsX+1;
      zeroRegionMax = startBin-nMENormSideBins-1;
     // zeroRegion(PprojX,startBin+nMENormSideBins-nProjBinsX,startBin-nMENormSideBins-1);
    }
//    printf(" Zeroing region [%d,%d]\n",zeroRegionMin,zeroRegionMax);
    zeroRegion(PprojX,zeroRegionMin,zeroRegionMax);
    //PprojX_Clone->Draw("SAME");
    fLineFitMin = -TMath::Pi() / 2;
    fLineFitMax = 3*TMath::Pi() / 2;
    PprojX->Fit("pol0","Q","",fLineFitMin,fLineFitMax);//Q = quiet mode, no printout

  } else {


    int zeroRegionMin = 1;
    int zeroRegionMax = startBin-nMENormSideBins-1;
    if (zeroRegionMin <= zeroRegionMax) zeroRegion(PprojX,zeroRegionMin,zeroRegionMax);
    zeroRegionMin = startBin+nMENormSideBins+1;
    zeroRegionMax = nProjBinsX;
    if (zeroRegionMin <= zeroRegionMax) zeroRegion(PprojX,zeroRegionMin,zeroRegionMax);


    LinFit = new TF1("pol0","pol0",fLineFitMin,fLineFitMax);
    /*if(!Can)
    {
      cout<<"2 no canvas provided"<<endl;
      PprojX->Fit("pol0","NQ","",fLineFitMin,fLineFitMax);//Q = quiet mode, no printout
    }*/
    PprojX->Fit("pol0","Q","",fLineFitMin,fLineFitMax);//Q = quiet mode, no printout
  }
  PprojX_Clone->Draw("SAME");
  PprojX->Draw("SAME");

	//..determine/etabin yield (count bins over which it was integrated (first bin and last bin included))
	// FIXME should this have + 1? if high range == low range, then 1 bin is used (not 0)
	//Int_t nBins= Histo->GetYaxis()->FindBin(highRange)-Histo->GetYaxis()->FindBin(lowRange);
	//Int_t nBins= 1+Histo->GetYaxis()->FindBin(highRange)-Histo->GetYaxis()->FindBin(lowRange);
	Int_t nBins= 1 + iHighRangeBin - iLowRangeBin ; 

  //printf(" debugME nBins = %d, iLowRangeBin = %d, iHighRangeBin = %d\n",nBins,iLowRangeBin,iHighRangeBin);
  printf(" debugME LowRangeBin = %.3f HighRangeBin-1 = %.3f, HighRange = %.3f \n",Histo->GetYaxis()->GetBinLowEdge(iLowRangeBin),Histo->GetYaxis()->GetBinUpEdge(iHighRangeBin-1),Histo->GetYaxis()->GetBinUpEdge(iHighRangeBin));

//	cout<<"high range: "<<highRange<<", bin: "<<Histo->GetYaxis()->FindBin(highRange)<<"| low range: "<<lowRange<<", bin: "<<Histo->GetYaxis()->FindBin(lowRange)<<endl;
	////cout<<"high range: "<<highRange<<", bin: "<<Histo->GetYaxis()->FindBin(highRange)-1<<"| low range: "<<lowRange<<", bin: "<<Histo->GetYaxis()->FindBin(lowRange)<<endl;
//`	cout<<"Nbins integrated: "<<nBins<<endl;
	Double_t Scalef = LinFit->GetParameter(0)/(1.0*nBins*nRebinMEForNorm);

  fMEValueAtZero = PprojX->GetBinContent(PprojX->FindBin(0.));
  fMEValueAtPi   = PprojX->GetBinContent(PprojX->FindBin(TMath::Pi()));
  fMEMean1DValue = PprojX->Integral() / PprojX->GetNbinsX();
  fMEValueAt1DMax= PprojX->GetBinContent(PprojX->GetMaximumBin());
  fMEScaleValue = Scalef;
  
  // FIXME temporary test of (0,0) = 1
  //Scalef = fMEValueAtZero / (1.0 * nBins * nRebinMEForNorm);

	if (Scalef == 0)
	{
		fprintf(stderr,"Error: Mixed Event Scalef is 0 !!!\n");
		return;
	}
	Histo->Scale(1/Scalef);

  // Check the success of scaling
  fMEValueAtOrigin = Histo->GetBinContent(Histo->FindBin(0,0));
  // FIXME calculate Mean value over a small range in eta. -0.1,0.1 ?
  //fMEMean2DValue = Histo->Integral() / (Histo->GetNbinsX() * Hist->GetNbinsY());
  fMEValueAt2DMax   = Histo->GetBinContent(Histo->GetMaximumBin());
  printf("2D ME Hist value at 0,0 = %f, Value at Maximum = %f\n",fMEValueAtOrigin,fMEValueAt2DMax);

  // Find a way to keep track of the average values, scale values, origin values, max values
  // Value at pi? both 0 and pi should be where TPC and EMCAL modules align.
  // TTree?

	TLegend *legME = new TLegend(0.25,0.70,0.4,0.9);
	//legME->SetHeader("Fit to plateu in #Delta#varphi");
	legME->AddEntry(Histo,"#Delta#varphi projected in","pe");
	legME->AddEntry(Histo,Form("%0.2f<#Delta#eta^{%s-h}<%0.2f",lowRange,fTriggerName.Data(),highRange),"");
	legME->AddEntry(Histo,Form("%.1f < z_{vtx} < %.1f cm",fArray_zVtx_Bins[CanvasPad],fArray_zVtx_Bins[CanvasPad+1]),"");
	legME->AddEntry(Histo,Form("fit value: %0.2f/%i",LinFit->GetParameter(0),nBins),"");
	legME->SetTextColor(kCyan+3);
	legME->SetTextSize(0.05);
	legME->SetBorderSize(0);
	legME->SetFillColorAlpha(10, 0);
	legME->Draw("same");


  fMEScaleTree->Fill();
}
///**
/// Plot the 2D correlation functions for SE and ME
///**
//________________________________________________________________________
void PlotGHcorrelation2::Plot2DHistograms()
{
	TLegend *legME;
	for(Int_t i=0;i<fmaxBins;i++)
	{
    if (fPlotVtzBins) {
      printf("Creating 2D plots with separate z-vertex bins for obs bin %d ...\n",i);
      fMEPlots2DGamma[i]->Divide(4,3,0.001,0.001);
      fRaw_Plots_2D[i]->Divide(4,3,0.001,0.001);

      for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
      //for(Int_t j=0;j<kNvertBins;j++)
      {
        //..plot the scaled 2D ME distribution
        fMEPlots2DGamma[i]->cd(j+1);
        DrawEtaPhi2D(fDetaDphi_ME[i][j]);

        if(j==0)legME = new TLegend(0.2,0.72,0.4,0.9);
        else    legME = new TLegend(0.2,0.85,0.4,0.9);
        legME->AddEntry(fDetaDphi_ME[i][j],Form("%.0f < z_{vtx} < %.0f",fArray_zVtx_Bins[j],fArray_zVtx_Bins[j+1]),"");
        if(j==0 && fObservable==0)legME->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.0f<#it{p}_{T}<%0.0f GeV/c",fArray_G_Bins[i],fArray_G_Bins[i+1]),"");
        if(j==0 && fObservable==1)legME->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
  //			if(j==0 && fObservable==2)legME->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
        if(j==0 && fObservable==2)legME->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#it{p}_{T}^{assoc}<%0.1f GeV/c",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
        legME->SetTextColor(kBlack);
        legME->SetTextSize(0.075);
        legME->SetBorderSize(0);
        legME->SetFillColorAlpha(10, 0);
        legME->Draw("same");
        //..plot the 2D SE distribution
        fRaw_Plots_2D[i]->cd(j+1);
        DrawEtaPhi2D(fDetaDphi_SE[i][j]);
        legME->Draw("same");
      }
    } else {
      printf("Creating 2D plots without separate z-vertex bins for obs bin %d ...\n",i);
      legME = new TLegend(0.2,0.72,0.4,0.86);
      // add WIP title
      fMEPlots2DGamma[i]->cd();
      //fDetaDphi_ME_alt1[i]->GetZaxis()->SetTitle("d^2N^{assoc}/d#Delta#eta d#Delta#varphi (a.u.)");
      DrawEtaPhi2D(fDetaDphi_ME_alt1[i],"d^{2}N^{assoc}/d#Delta#eta d#Delta#varphi (a.u.)");
      legME->AddEntry(fDetaDphi_ME_alt1[i],"ALICE Work in Progress","");
      if(fObservable==0)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.0f<#it{p}_{T}<%0.0f GeV/c",fArray_G_Bins[i],fArray_G_Bins[i+1]),"");
      if(fObservable==1)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//			if(j==0 && fObservable==2)legME->AddEntry(fDetaDphi_ME[i][j],Form("All pads %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
      if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.1f<#it{p}_{T}^{assoc}<%0.1f GeV/c",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
      legME->SetTextColor(kBlack);
      legME->SetTextSize(0.075);
      legME->SetBorderSize(0);
      legME->SetFillColorAlpha(10, 0);
      legME->Draw("same");

      fRaw_Plots_2D[i]->cd();
      //fDetaDphi_SE_alt1[i]->GetZaxis()->SetTitle("d^2N^{assoc}/d#Delta#eta d#Delta#varphi (a.u.)");
      DrawEtaPhi2D(fDetaDphi_SE_alt1[i],"d^{2}N^{assoc}/d#Delta#eta d#Delta#varphi (a.u.)");
      legME->Draw("same");
    }
		fMEPlots2DGamma[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fMEPlots2DGamma[i]->GetName()));
		fMEPlots2DGamma[i]->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fMEPlots2DGamma[i]->GetName()));
		fRaw_Plots_2D[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fRaw_Plots_2D[i]->GetName()));
		fRaw_Plots_2D[i]->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fRaw_Plots_2D[i]->GetName()));
	}
  // Plot the example plotsa
  int iExampleObsBin = 3;
  int iExampleZVtx   = 5;
  TCanvas * cFancyExampleCanvas = new TCanvas("Example2DPlot","Example2DPlot",1000,600);
  cFancyExampleCanvas->Divide(2,1);
  cFancyExampleCanvas->cd(1);
  fDetaDphi_SE[iExampleObsBin][iExampleZVtx]->Draw("LEGO2");
  cFancyExampleCanvas->cd(2);
  fDetaDphi_ME[iExampleObsBin][iExampleZVtx]->Draw("LEGO2");
  printf("Will try to save a 2d plot of name %s\n",cFancyExampleCanvas->GetName());
  printf("    For %.0f < z_{vtx} < %.0f, obs bin %d\n",fArray_zVtx_Bins[iExampleZVtx],fArray_zVtx_Bins[iExampleZVtx+1],iExampleObsBin);
  cFancyExampleCanvas->Print(TString::Format("%s/%s.png",fOutputDir.Data(),cFancyExampleCanvas->GetName()));
//  cFancyExampleCanvas->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),cFancyExampleCanvas->GetName());

}
///**
/// Divide SE by ME and then plot the corrected 2D correlation functions (SE/ME)
///**
//________________________________________________________________________
void PlotGHcorrelation2::PlotCorrected2DHistograms()
{
	TLegend* legSE;
	TLegend* legSE_vtx;
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  First Step
	//..  Divide SE by ME for individual z-Vertexes and plot them
	//. . . . . . . . . . . . . . . . . . . . . . . .
	fPlots_2D_Corr_alt1->Divide(3,3,0.001,0.001);
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_2D_Corr[i]->Divide(4,3,0.001,0.001);
		//..correct the SE for each z-Vtx bin separatley
		for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
		//for(Int_t j=0;j<kNvertBins;j++)
		{
			fPlots_2D_Corr[i]->cd(j+1);
			fDetaDphi_SE[i][j]->Divide(fDetaDphi_ME[i][j]);
			DrawEtaPhi2D(fDetaDphi_SE[i][j]);

			if(j==0)legSE_vtx = new TLegend(0.2,0.72,0.4,0.9);
			else    legSE_vtx = new TLegend(0.2,0.85,0.4,0.9);
			legSE_vtx->AddEntry(fDetaDphi_SE[i][j],Form("%.0f < z_{vtx} < %.0f",fArray_zVtx_Bins[j],fArray_zVtx_Bins[j+1]),"");
			if(j==0 && fObservable==0)legSE_vtx->AddEntry(fDetaDphi_SE[i][j],Form("All pads %0.0f<#it{p}_{T}<%0.0f",fArray_G_Bins[i],fArray_G_Bins[i+1]),"");
			if(j==0 && fObservable==1)legSE_vtx->AddEntry(fDetaDphi_SE[i][j],Form("All pads %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//			if(j==0 && fObservable==2)legSE_vtx->AddEntry(fDetaDphi_SE[i][j],Form("All pads %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
			if(j==0 && fObservable==2)legSE_vtx->AddEntry(fDetaDphi_SE[i][j],Form("All pads %0.1f<#it{p}_{T}^{assoc}<%0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
			legSE_vtx->SetTextColor(kBlack);
			legSE_vtx->SetTextSize(0.075);
			legSE_vtx->SetBorderSize(0);
			legSE_vtx->SetFillColorAlpha(10, 0);
			legSE_vtx->Draw("same");
		}
		fPlots_2D_Corr[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_Corr[i]->GetName()));
		fPlots_2D_Corr[i]->Print(TString::Format("%s/%s.eps",fOutputDir.Data(),fPlots_2D_Corr[i]->GetName()));
		fPlots_2D_Corr[i]->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_2D_Corr[i]->GetName()));

		//..Do the same proceedure for different ME strategies
		if(fPlotMoreMEstrategies==1)
		{
			fPlots_2D_Corr_alt2[i]->Divide(3,3,0.001,0.001);
			fPlots_2D_Corr_alt1->cd(i+1);
			fDetaDphi_SE_alt1[i]->Divide(fDetaDphi_ME_alt1[i]);
			DrawEtaPhi2D(fDetaDphi_SE_alt1[i]);
			for(Int_t j=0;j<kNvertBins_alt;j++)
			{
				fPlots_2D_Corr_alt2[i]->cd(j+1);
				fDetaDphi_SE_alt2[i][j]->Divide(fDetaDphi_ME_alt2[i][j]);
				DrawEtaPhi2D(fDetaDphi_SE_alt2[i][j]);
			}
			fPlots_2D_Corr_alt3[i]->Divide(3,3,0.001,0.001);
			for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
			{
				fPlots_2D_Corr_alt3[i]->cd(j+1);
				fDetaDphi_SE_alt3[i][j]->Divide(fDetaDphi_ME_alt3[j]);
				DrawEtaPhi2D(fDetaDphi_SE_alt3[i][j]);
			}
			fPlots_2D_Corr_alt2[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_Corr_alt2[i]->GetName()));
			fPlots_2D_Corr_alt3[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_Corr_alt3[i]->GetName()));
		}
	}
	if(fPlotMoreMEstrategies==1)fPlots_2D_Corr_alt1->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_Corr_alt1->GetName()));

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  Second Step
	//..  Add all Corrected z-Vtx bins together to a single histogram for each fObservable
	//. . . . . . . . . . . . . . . . . . . . . . . .
	for(Int_t i=0;i<fmaxBins;i++)
	{
		for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
		{
			fsumCorrSE[i]->Add(fDetaDphi_SE[i][j]);
		}
		//..Do the same proceedure for different ME strategies
		if(fPlotMoreMEstrategies==1)
		{
			//..this is basically just passing over to another histogram because fDetaDphi_SE_alt1 is already summed up
			fsumCorrSE_alt1[i]->Add(fDetaDphi_SE_alt1[i]);
			for(Int_t j=0;j<kNvertBins_alt;j++)
			{
				fsumCorrSE_alt2[i]->Add(fDetaDphi_SE_alt2[i][j]);
			}
			for(Int_t j=fminZvtx;j<fmaxZvtx;j++)
			{
				fsumCorrSE_alt3[i]->Add(fDetaDphi_SE_alt3[i][j]);
			}
		}
	}
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  Third Step
	//..  Normalize the corrected correlation functions by the number of triggers
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..Normalize g-h pairs to No. of triggers
	//..possibly move this part into the NormalizeSEsignal() function
	cout<<"o Normalize SE"<<endl;
	fCanNormCheck->Divide(3,3,0.001,0.001);
	for(Int_t i=0;i<fmaxBins;i++)
	{
		// Normalize by N_{Triggers}
		double fScale = 1;
		if(fTrigger_SE[i])
		{
      // FIXME should this be fTriggerPt?
			//fScale = fTrigger_SE[i]->Integral();
      // FIXME this might only be correct for 
      int iIntegralPtBin = 1;
      if (fObservable==0) {
        iIntegralPtBin = i+1;
        fScale = fTriggerPt->Integral(iIntegralPtBin,iIntegralPtBin);
       // fScale = fTriggerPt->Integral(fPtMinBin,fPtMaxBin);
        printf(" DEBUG: fMassPtPionAccProj_%d gives %f.\n",iIntegralPtBin-1,fMassPtPionAccProj[iIntegralPtBin-1]->Integral());
      }
      else { // fObservable = 1 or 2
//        iIntegralPtBin = fPtMinBin
        printf(" DEBUG:  iIntegralPtBin = %d\n",iIntegralPtBin);
       // fScale = fTriggerPt->Integral(fPtMinBin,fPtMaxBin);
        // Since the source THnSparse is already limited to the proper range before it was projected, I can just take the
        // whole integral of fTriggerPt
        fScale = fTriggerPt->Integral();
        //fScale = fTriggerPt->Integral(iIntegralPtBin,iIntegralPtBin);
        // FIXME Make this loop over all the pt bins included
        printf(" DEBUG: fMassPtPionAccProj_%d gives %f.\n",fPtMinBin-1,fMassPtPionAccProj[fPtMinBin-1]->Integral());
      }

      // FIXME restrict to specific pt range
      printf(" DEBUG: fTriggerSE[%d] gives %f triggers, fTriggerPt gives %f triggers\n",i,fTrigger_SE[i]->Integral(),fTriggerPt->Integral(iIntegralPtBin,iIntegralPtBin));
			printf("  o Normalizing bin %d with N_{Trig} = %f\n",i,fScale);

			fCanNormCheck->cd(i+1);
			SetTH1Histo(fTrigger_SE[i],"z Vertex position (cm)","counts",0);
			fTrigger_SE[i]->DrawCopy("hist");

			if (fScale != 0.)
			{

				fsumCorrSE[i]->Scale(1./fScale);
				//..Do the same proceedure for different ME strategies
				if(fPlotMoreMEstrategies==1)
				{
					fsumCorrSE_alt1[i]->Scale(1./fScale);
					fsumCorrSE_alt2[i]->Scale(1./fScale);
					fsumCorrSE_alt3[i]->Scale(1./fScale);
				}
			}
			else cout<<"ERROR: Problem with trigger scaling. Scale is 0!"<<endl;
		}
		else
		{
			fprintf(stderr,"Missing Trigger Histogram for bin %d.\n",i);
		}
	}
  TCanvas * fLocalPlot = new TCanvas("Corr_2D_Indiv","Corr_2D_Indiv",1750,1400);
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  Fourth Step
	//..  Plot the corrected, summed, and normalized correlation functions for each fObservable bin
	//. . . . . . . . . . . . . . . . . . . . . . . .
	fPlots_2D_CorrSum->Divide(3,3,0.001,0.001);
	fPlots_2D_CorrSum_alt1->Divide(3,3,0.001,0.001);
	fPlots_2D_CorrSum_alt2->Divide(3,3,0.001,0.001);
	fPlots_2D_CorrSum_alt3->Divide(3,3,0.001,0.001);
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_2D_CorrSum->cd(i+1);
//		TVirtualPad * fLocalPad = fPlots_2D_CorrSum->cd(i+1);
		fsumCorrSE[i]->GetYaxis()->SetTitleOffset(0.5);

		// Scaling Z-axis to AwaySide max
		fsumCorrSE[i]->GetXaxis()->SetRangeUser(TMath::Pi()/2,3.*TMath::Pi()/2);
		Double_t awaySideMax = fsumCorrSE[i]->GetBinContent(fsumCorrSE[i]->GetMaximumBin());

		fsumCorrSE[i]->GetZaxis()->SetRangeUser(fsumCorrSE[i]->GetZaxis()->GetXmin(),awaySideMax*1.4);
		fsumCorrSE[i]->GetYaxis()->SetRangeUser(-fMaxDeltaEtaPlotRange,fMaxDeltaEtaPlotRange);
		fsumCorrSE[i]->GetXaxis()->SetRangeUser(-TMath::Pi()/2,3.*TMath::Pi()/2);
		DrawEtaPhi2D(fsumCorrSE[i]);

		//if(j==0)legSE = new TLegend(0.1,0.75,0.4,0.9);
		legSE = new TLegend(0.41,0.67,0.81,0.79); //0.55, 0.71, 0.87, 0.84
		if(fObservable==0)legSE->AddEntry(fsumCorrSE[i],Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
	  if(fObservable!=0)legSE->AddEntry((TObject*)0,Form("%0.1f < #it{p}_{T}^{%s} < %0.1f GeV/#it{c}",fArray_G_Bins[fPtMinBin-1],fTriggerName.Data(),fArray_G_Bins[fPtMaxBin]),"");
		if(fObservable==1)legSE->AddEntry(fsumCorrSE[i],Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//		if(fObservable==2)legSE->AddEntry(fsumCorrSE[i],Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
		if(fObservable==2)legSE->AddEntry(fsumCorrSE[i],Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
		legSE->SetTextColor(kBlack);
		//legSE->SetTextSize(0.08);
		legSE->SetTextSize(0.04);
		legSE->SetBorderSize(0);
//		legSE->SetFillColorAlpha(10, 0);
    legSE->SetFillStyle(0);
		legSE->Draw("same");
    // Draw individual plots
    printf("Trying to save individual plots ... \n");
 //   gPad->Print(TString::Format("%s/Test_Plot%d.png",fOutputDir.Data(),i));
  //  gPad->Print(TString::Format("%s/%s_Plot%d.png",fOutputDir.Data(),fPlots_2D_CorrSum->GetName(),i));
 //   gPad->Print(TString::Format("%s/%s_Plot%d.C",fOutputDir.Data(),fPlots_2D_CorrSum->GetName(),i));

    // Settings for Indivdual DPhiDEta 2D Plot

    fLocalPlot->cd();
		DrawEtaPhi2D(fsumCorrSE[i]);
    DrawAlicePerf(fsumCorrSE[i],0.2,0.78,0.12,0.12); //0.22,0.8,0.12,0.12
    legSE->Draw("SAME");
    fLocalPlot->Print(TString::Format("%s/%s_Plot%d.pdf",fOutputDir.Data(),fLocalPlot->GetName(),i));
    fLocalPlot->Print(TString::Format("%s/%s_Plot%d.png",fOutputDir.Data(),fLocalPlot->GetName(),i));
    fLocalPlot->Print(TString::Format("%s/%s_Plot%d.eps",fOutputDir.Data(),fLocalPlot->GetName(),i));
    fLocalPlot->Print(TString::Format("%s/%s_Plot%d.C",fOutputDir.Data(),fLocalPlot->GetName(),i));

		//..Do the same proceedure for different ME strategies
		if(fPlotMoreMEstrategies==1)
		{
			fPlots_2D_CorrSum_alt1->cd(i+1);
			DrawEtaPhi2D(fsumCorrSE_alt1[i]);
			legSE->Draw("same");
			fPlots_2D_CorrSum_alt2->cd(i+1);
			DrawEtaPhi2D(fsumCorrSE_alt2[i]);
			legSE->Draw("same");
			fPlots_2D_CorrSum_alt3->cd(i+1);
			DrawEtaPhi2D(fsumCorrSE_alt3[i]);
			legSE->Draw("same");
		}
	}
	fCanNormCheck->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fCanNormCheck->GetName()));
	fPlots_2D_CorrSum->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fPlots_2D_CorrSum->GetName()));
//	fPlots_2D_CorrSum->Print(TString::Format("%s/%s.eps",fOutputDir.Data(),fPlots_2D_CorrSum->GetName()));
	fPlots_2D_CorrSum->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_CorrSum->GetName()));
	fPlots_2D_CorrSum->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_2D_CorrSum->GetName()));
	fPlots_2D_CorrSum_alt1->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_CorrSum_alt1->GetName()));
	fPlots_2D_CorrSum_alt1->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_2D_CorrSum_alt1->GetName()));
	fPlots_2D_CorrSum_alt2->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_CorrSum_alt2->GetName()));
	fPlots_2D_CorrSum_alt2->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_2D_CorrSum_alt2->GetName()));
	fPlots_2D_CorrSum_alt3->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_CorrSum_alt3->GetName()));
	fPlots_2D_CorrSum_alt3->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_2D_CorrSum_alt3->GetName()));
}

///**
/// compare the different versions how to calculate the ME
/// background. A. differential in z-vertex, B. all z-vertices together
///**
//________________________________________________________________________
void PlotGHcorrelation2::PlotCompare2DHistoCorrection()
{
	TString ProjectionName;
	TLegend* legME= nullptr;
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  MIXED EVENT ALTERNATIVE -1-
	//..  add all Z-vtx bins for SE and ME and then divide them. (old simple method)
	//. . . . . . . . . . . . . . . . . . . . . . . .
	fPlots_2D_ME_alt1->Divide(3,3,0.001,0.001);
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_2D_ME_alt1->cd(i+1);
		DrawEtaPhi2D(fDetaDphi_ME_alt1[i]);

		legME = new TLegend(0.2,0.72,0.4,0.9);
		if(fObservable==0)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
		if(fObservable==1)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//		if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
		if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt1[i],Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
		legME->AddEntry(fDetaDphi_ME_alt1[i],Form("-10 < z-Vtx < 10"),"");
		legME->SetTextColor(kBlack);
		legME->SetTextSize(0.08);
		legME->SetBorderSize(0);
		//legME->SetFillColorAlpha(10, 0);
    legME->SetFillStyle(0);
		legME->Draw("same");
	}
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  MIXED EVENT ALTERNATIVE -2-
	//..  add all Z-vtx bins for SE and ME and then divide them. (old simple method)
	//. . . . . . . . . . . . . . . . . . . . . . . .
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_2D_ME_alt2[i]->Divide(3,3,0.001,0.001);
		for(Int_t j=0;j<kNvertBins_alt;j++)
		{
			fPlots_2D_ME_alt2[i]->cd(j+1);
			DrawEtaPhi2D(fDetaDphi_ME_alt2[i][j]);

			legME = new TLegend(0.2,0.72,0.4,0.9);
			//legME = new TLegend(0.25,0.83,0.4,0.88);
			if(fObservable==0)legME->AddEntry(fDetaDphi_ME_alt2[i][j],Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
			if(fObservable==1)legME->AddEntry(fDetaDphi_ME_alt2[i][j],Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
			//if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt2[i][j],Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
			if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt2[i][j],Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
			legME->AddEntry(fDetaDphi_ME_alt2[i][j],Form("%.0f < z-Vtx < %.0f",fArray_zVtx_BinsAlt[j],fArray_zVtx_BinsAlt[j+1]),"");
			legME->SetTextColor(kBlack);
			legME->SetTextSize(0.08);
			legME->SetBorderSize(0);
			//legME->SetFillColorAlpha(10, 0);
      legME->SetFillStyle(0);
			legME->Draw("same");
		}
	}
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  MIXED EVENT ALTERNATIVE -3-
	//..  add all Z-vtx bins over different fObservable bins, in total there will be
	//..  the same ME background for all bins but with 20 z-Vtx bins.
	//. . . . . . . . . . . . . . . . . . . . . . . .
	fPlots_2D_ME_alt3->Divide(3,3,0.001,0.001);
	for(Int_t j=0;j<kNvertBins;j++)
	{
		fPlots_2D_ME_alt3->cd(j+1);
		DrawEtaPhi2D(fDetaDphi_ME_alt3[j]);

		legME = new TLegend(0.2,0.72,0.4,0.9);
		//legME = new TLegend(0.25,0.83,0.4,0.88);
		if(fObservable==0)legME->AddEntry(fDetaDphi_ME_alt3[j],Form("for all #it{p}_{T}^{%s} bins",fTriggerName.Data()),"");
		if(fObservable==1)legME->AddEntry(fDetaDphi_ME_alt3[j],Form("for all z_{T} bins"),"");
//		if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt3[j],Form("for all #xi bins"),"");
		if(fObservable==2)legME->AddEntry(fDetaDphi_ME_alt3[j],Form("for all #it{p}_{T}^{assoc} bins"),"");
		legME->AddEntry(fDetaDphi_ME_alt3[j],Form("%.0f < z-Vtx < %.0f",fArray_zVtx_Bins[j],fArray_zVtx_Bins[j+1]),"");
		legME->SetTextColor(kBlack);
		legME->SetTextSize(0.08);
		legME->SetBorderSize(0);
		//legME->SetFillColorAlpha(10, 0);
    legME->SetFillStyle(0);
		legME->Draw("same");
	}


	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  Big Comparison
	//..  Now take all three alternative versions and
	//..  additionally the original corrected SE and compare their projections
	//. . . . . . . . . . . . . . . . . . . . . . . .

	//..compare the errors+yield of the different histograms
	fPlots_1D_MEcompare      ->Divide(3,3,0.001,0.001);
	fPlots_1D_MEcompareEta   ->Divide(3,3,0.001,0.001);
	fPlots_1D_MEcompareRatio ->Divide(3,3,0.001,0.001);
	fPlots_1D_MEcompareErrors->Divide(3,3,0.001,0.001);
	Double_t etaMin=-1.3;
	Double_t etaMax=1.3;
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_1D_MEcompare->cd(i+1);
		TH1D *PprojXSigOrig  =fsumCorrSE[i]     ->ProjectionX(Form("%s_ProjOrig",fsumCorrSE[i]->GetName()),fsumCorrSE[i]->GetYaxis()->FindBin(etaMin),fsumCorrSE[i]->GetYaxis()->FindBin(etaMax)); //..do not include the edges
		TH1D *PprojXSigAlt1  =fsumCorrSE_alt1[i]->ProjectionX(Form("%s_Proj_Alt1",fsumCorrSE_alt1[i]->GetName()),fsumCorrSE_alt1[i]->GetYaxis()->FindBin(etaMin),fsumCorrSE_alt1[i]->GetYaxis()->FindBin(etaMax)); //..do not include the edges
		TH1D *PprojXSigAlt2  =fsumCorrSE_alt2[i]->ProjectionX(Form("%s_Proj_Alt2",fsumCorrSE_alt2[i]->GetName()),fsumCorrSE_alt2[i]->GetYaxis()->FindBin(etaMin),fsumCorrSE_alt1[i]->GetYaxis()->FindBin(etaMax)); //..do not include the edges
		TH1D *PprojXSigAlt3  =fsumCorrSE_alt3[i]->ProjectionX(Form("%s_Proj_Alt3",fsumCorrSE_alt3[i]->GetName()),fsumCorrSE_alt3[i]->GetYaxis()->FindBin(etaMin),fsumCorrSE_alt1[i]->GetYaxis()->FindBin(etaMax)); //..do not include the edges
		PprojXSigOrig ->GetYaxis()->SetTitle("dN.../d");
		PprojXSigOrig ->SetMarkerStyle(kFullSquare);
	//	PprojXSigOrig ->SetMarkerSize(0.8);
		PprojXSigOrig ->SetMarkerStyle(kOpenSquare);
		PprojXSigOrig ->SetMarkerColor(kYellow-3);
		PprojXSigOrig ->SetLineColor(kYellow-3);
		PprojXSigAlt1 ->SetMarkerStyle(kFullSquare);
	//	PprojXSigAlt1 ->SetMarkerSize(0.8);
		PprojXSigAlt1 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt1 ->SetMarkerColor(kAzure-7);
		PprojXSigAlt1 ->SetLineColor(kAzure-7);
		PprojXSigAlt2 ->SetMarkerStyle(kFullSquare);
	//	PprojXSigAlt2 ->SetMarkerSize(0.8);
		PprojXSigAlt2 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt2 ->SetMarkerColor(kRed-7);
		PprojXSigAlt2 ->SetLineColor(kRed-7);
		PprojXSigAlt3 ->SetMarkerStyle(kFullSquare);
	//	PprojXSigAlt3 ->SetMarkerSize(0.8);
		PprojXSigAlt3 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt3 ->SetMarkerColor(kGreen-7);
		PprojXSigAlt3 ->SetLineColor(kGreen-7);
		PprojXSigOrig ->DrawCopy(" E");
		PprojXSigAlt1 ->DrawCopy("same E");
		PprojXSigAlt2 ->DrawCopy("same E");
		PprojXSigAlt3 ->DrawCopy("same E");
/*		PprojXSigOrig ->DrawNormalized(" E",1);
		PprojXSigAlt1 ->DrawNormalized("same E",1);
		PprojXSigAlt2 ->DrawNormalized("same E",1);
		PprojXSigAlt3 ->DrawNormalized("same E",1);
*/
		legME = new TLegend(0.55,0.72,0.7,0.9);
		if(fObservable==0)legME->AddEntry(PprojXSigOrig,Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/c",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
		if(fObservable==1)legME->AddEntry(PprojXSigOrig,Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//		if(fObservable==2)legME->AddEntry(PprojXSigOrig,Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
		if(fObservable==2)legME->AddEntry(PprojXSigOrig,Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
		legME->AddEntry(PprojXSigOrig,Form("%0.1f < #Delta#eta < %0.1f",etaMin,etaMax),"");
		legME->SetTextColor(kBlack);
		legME->SetTextSize(0.08);
		legME->SetBorderSize(0);
		//legME->SetFillColorAlpha(10, 0);
    legME->SetFillStyle(0);
		legME->Draw("same");

		if(i==(fmaxBins-1))
		{
			legME = new TLegend(0.2,0.5,0.4,0.9);
			legME->AddEntry(PprojXSigOrig,"Original corr. 10 zVtx bins","pl");
			legME->AddEntry(PprojXSigAlt1->Clone("PprojXSigAlt1"),"Alternative corr. 1 zVtx bin","pl");
			legME->AddEntry(PprojXSigAlt2,"Alternative corr. 8 zVtx bins","pl");
			legME->AddEntry(PprojXSigAlt3,"Alternative corr. 10 zVtx bins","pl");
			legME->AddEntry(PprojXSigAlt3,"summed over fObs","");
			legME->SetTextColor(kBlack);
			legME->SetTextSize(0.08);
			legME->SetBorderSize(0);
			//legME->SetFillColorAlpha(10, 0);
      legME->SetFillStyle(0);
		}

		//. . . . . . . . . . . . . . . . . . . . . . .
		//..compare the size of the errors
		fPlots_1D_MEcompareErrors->cd(i+1);
		TH1D* ErrorCompOrig = (TH1D*)PprojXSigOrig->Clone("ErrorCompOrig");
		TH1D* ErrorCompAlt1 = (TH1D*)PprojXSigAlt1->Clone("ErrorCompAlt1");
		TH1D* ErrorCompAlt2 = (TH1D*)PprojXSigAlt2->Clone("ErrorCompAlt2");
		TH1D* ErrorCompAlt3 = (TH1D*)PprojXSigAlt3->Clone("ErrorCompAlt3");
		ErrorCompOrig->Reset();
		ErrorCompAlt1->Reset();
		ErrorCompAlt2->Reset();
		ErrorCompAlt3->Reset();
		//FIXME check everthing for underflow overflow bins
		//..loop over the histogram and set the bin content as the bin error
		for(Int_t i=1; i<ErrorCompOrig->GetXaxis()->GetNbins()+1;i++)
		{
			ErrorCompOrig ->SetBinContent(i,PprojXSigOrig->GetBinError(i));
			ErrorCompAlt1 ->SetBinContent(i,PprojXSigAlt1->GetBinError(i));
			ErrorCompAlt2 ->SetBinContent(i,PprojXSigAlt2->GetBinError(i));
			ErrorCompAlt3 ->SetBinContent(i,PprojXSigAlt3->GetBinError(i));
		}
		ErrorCompAlt1->Divide(ErrorCompOrig);
		ErrorCompAlt2->Divide(ErrorCompOrig);
		ErrorCompAlt3->Divide(ErrorCompOrig);
		Double_t minRange,maxRange;
		maxRange= ErrorCompAlt2->GetMaximum();
		minRange= ErrorCompAlt3->GetMinimum();
		//cout<<"min and max found: "<<minRange<<" - "<<maxRange<<endl;

		//SetTH1Histo(ErrorCompAlt1,Form("#Delta#varphi^{%s-h}",fTriggerName.Data()),"Error of Var./Orig. Error",1);
		SetTH1Histo(ErrorCompAlt1,"#Delta#varphi","Error of Var./Orig. Error",1);
		ErrorCompAlt1->SetLineColor(kAzure-7);
		ErrorCompAlt1 ->GetYaxis()->SetRangeUser(minRange*0.9,maxRange*1.05);
		ErrorCompAlt1 ->DrawCopy(" hist");
		ErrorCompAlt2 ->DrawCopy("same hist");
		ErrorCompAlt3 ->DrawCopy("same hist");

		//. . . . . . . . . . . . . . . . . . . . . . .
		//..Plot the ratio to compare the shape
		fPlots_1D_MEcompareRatio->cd(i+1);
//		PprojXSigAlt1 ->Add(PprojXSigOrig,-1);
		PprojXSigAlt1 ->Divide(PprojXSigOrig);
		PprojXSigAlt1 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt1 ->SetMarkerColor(kAzure-7);
		PprojXSigAlt1 ->SetLineColor(kAzure-7);
//		PprojXSigAlt2 ->Add(PprojXSigOrig,-1);
		PprojXSigAlt2 ->Divide(PprojXSigOrig);
//		PprojXSigAlt3 ->Add(PprojXSigOrig,-1);
		PprojXSigAlt3 ->Divide(PprojXSigOrig);
		PprojXSigAlt1 ->GetYaxis()->SetRangeUser(0.98,1.01);
		if(i>1)PprojXSigAlt1 ->GetYaxis()->SetRangeUser(0.9,1.05);
		if(i>2)PprojXSigAlt1 ->GetYaxis()->SetRangeUser(0.7,1.1);
		if(i>3)PprojXSigAlt1 ->GetYaxis()->SetRangeUser(0.4,1.5);
		SetTH1Histo(PprojXSigAlt1,"#Delta#varphi","#Delta#varphi Variat./#Delta #varphi Orig.",1);
		//SetTH1Histo(PprojXSigAlt1,Form("#Delta#varphi^{%s-h}",fTriggerName.Data()),"#Delta#varphi Variat./#Delta #varphi Orig.",1);
//		SetTH1Histo(PprojXSigAlt1,Form("#Delta#varphi^{%s-h}",fTriggerName.Data()),"#Delta#varphi (Variat.-Orig.)/Orig.",1);
		PprojXSigAlt1 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt1 ->SetMarkerColor(kAzure-7);
		PprojXSigAlt1 ->SetLineColor(kAzure-7);
		PprojXSigAlt1 ->DrawCopy("E");
		PprojXSigAlt2 ->DrawCopy("same E");
		PprojXSigAlt3 ->DrawCopy("same E");
	}

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..  Big Comparison
	//..  Eta part
	//. . . . . . . . . . . . . . . . . . . . . . . .

	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_1D_MEcompareEta->cd(i+1);
		TH1D *PprojXSigOrig  =fsumCorrSE[i]     ->ProjectionY(Form("%s_ProjOrigEta",fsumCorrSE[i]->GetName()));
		TH1D *PprojXSigAlt1  =fsumCorrSE_alt1[i]->ProjectionY(Form("%s_ProjEta_Alt1",fsumCorrSE_alt1[i]->GetName()));
		TH1D *PprojXSigAlt2  =fsumCorrSE_alt2[i]->ProjectionY(Form("%s_ProjEta_Alt2",fsumCorrSE_alt2[i]->GetName()));
		TH1D *PprojXSigAlt3  =fsumCorrSE_alt3[i]->ProjectionY(Form("%s_ProjEta_Alt3",fsumCorrSE_alt3[i]->GetName()));
		PprojXSigOrig ->GetYaxis()->SetTitle("dN.../d");
		PprojXSigOrig ->SetMarkerStyle(kFullSquare);
//		PprojXSigOrig ->SetMarkerSize(0.8);
//		PprojXSigOrig ->SetMarkerStyle(kOpenSquare);
//		PprojXSigOrig ->SetMarkerColor(kYellow-3);
//		PprojXSigOrig ->SetLineColor(kYellow-3);
		PprojXSigOrig ->SetMarkerColor(kBlack);
		PprojXSigOrig ->SetLineColor(kBlack);
//		PprojXSigAlt1 ->SetMarkerStyle(kFullSquare);
//		PprojXSigAlt1 ->SetMarkerSize(0.8);
		PprojXSigAlt1 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt1 ->SetMarkerColor(kAzure-7);
		PprojXSigAlt1 ->SetLineColor(kAzure-7);
//		PprojXSigAlt2 ->SetMarkerStyle(kFullSquare);
//		PprojXSigAlt2 ->SetMarkerSize(0.8);
		PprojXSigAlt2 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt2 ->SetMarkerColor(kRed-7);
		PprojXSigAlt2 ->SetLineColor(kRed-7);
//		PprojXSigAlt3 ->SetMarkerStyle(kFullSquare);
//		PprojXSigAlt3 ->SetMarkerSize(0.8);
		PprojXSigAlt3 ->SetMarkerStyle(kOpenSquare);
		PprojXSigAlt3 ->SetMarkerColor(kGreen-7);
		PprojXSigAlt3 ->SetLineColor(kGreen-7);

		//..for performance
		PprojXSigOrig->GetXaxis()->SetRangeUser(-0.8,0.8);
		/*
		Double_t lowRange =PprojXSigOrig->GetXaxis()->GetFirst();
		Double_t highRange=PprojXSigOrig->GetXaxis()->GetLast();
		Double_t min = PprojXSigOrig->GetMinimum();
		PprojXSigOrig->GetXaxis()->SetRangeUser(-0.8,0.8);
		Double_t range1= PprojXSigOrig->GetMaximum()-min;
		range1= PprojXSigOrig->GetMaximum()-PprojXSigOrig->GetMinimum();
		PprojXSigOrig->GetYaxis()->SetRangeUser(PprojXSigOrig->GetMinimum()-range1*0.4,PprojXSigOrig->GetMinimum()+range1*1.4);
		PprojXSigOrig->GetXaxis()->SetRange(lowRange,highRange);
*/

		PprojXSigOrig ->DrawCopy(" E");
		PprojXSigAlt1 ->DrawCopy("same E");
		PprojXSigAlt2 ->DrawCopy("same E");
		PprojXSigAlt3 ->DrawCopy("same E");
/*		PprojXSigOrig ->DrawNormalized(" E",1);
		PprojXSigAlt1 ->DrawNormalized("same E",1);
		PprojXSigAlt2 ->DrawNormalized("same E",1);
		PprojXSigAlt3 ->DrawNormalized("same E",1);
*/
		legME = new TLegend(0.55,0.72,0.7,0.9);
		if(fObservable==0)legME->AddEntry(PprojXSigOrig,Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
		if(fObservable==1)legME->AddEntry(PprojXSigOrig,Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
//		if(fObservable==2)legME->AddEntry(PprojXSigOrig,Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
		if(fObservable==2)legME->AddEntry(PprojXSigOrig,Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
		//legME->AddEntry(PprojXSigOrig,Form("%0.1f < #Delta#eta < %0.1f",etaMin,etaMax),"");
		legME->SetTextColor(kBlack);
		legME->SetTextSize(0.08);
		legME->SetBorderSize(0);
		//legME->SetFillColorAlpha(10, 0);
    legME->SetFillStyle(0);
		legME->Draw("same");

		if(i==(fmaxBins-1))
		{
			legME = new TLegend(0.2,0.45,0.4,0.9);
			legME->AddEntry(PprojXSigOrig,"Original corr. 10 zVtx bins","pl");
			legME->AddEntry(PprojXSigAlt1->Clone("PprojXSigAlt1"),"Alternative corr. 1 zVtx bin","pl");
			legME->AddEntry(PprojXSigAlt2,"Alternative corr. 8 zVtx bins","pl");
			legME->AddEntry(PprojXSigAlt3,"Alternative corr. 10 zVtx bins","pl");
			legME->AddEntry(PprojXSigAlt3,"summed over fObs","");
			legME->SetTextColor(kBlack);
			legME->SetTextSize(0.08);
			legME->SetBorderSize(0);
			//legME->SetFillColorAlpha(10, 0);
      legME->SetFillStyle(0);
		}
	}
	//..plot the legends into the canvas
	fPlots_1D_MEcompare->cd(8);
	legME->Draw("");
	fPlots_1D_MEcompareRatio->cd(8);
	legME->Draw("");
	fPlots_1D_MEcompareErrors->cd(8);
	legME->Draw("");
	fPlots_1D_MEcompareEta->cd(8);
	legME->Draw("");

	//..Save the compare Canvases
	fPlots_2D_ME_alt1->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_ME_alt1->GetName()));
	fPlots_2D_ME_alt1->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_2D_ME_alt1->GetName()));
	for(Int_t i=0;i<fmaxBins;i++)
	{
		fPlots_2D_ME_alt2[i]->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_ME_alt2[i]->GetName()));
	}
	fPlots_2D_ME_alt3->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_2D_ME_alt3->GetName()));
	fPlots_1D_MEcompare->Print(TString::Format("%s/%s.pdf",fOutputDir.Data(),fPlots_1D_MEcompare->GetName()));
	fPlots_1D_MEcompare->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_1D_MEcompare->GetName()));
	fPlots_1D_MEcompare->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_1D_MEcompare->GetName()));
	fPlots_1D_MEcompareEta->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_1D_MEcompareEta->GetName()));
	fPlots_1D_MEcompareEta->Print(TString::Format("%s/CFiles/%s.C",fOutputDir.Data(),fPlots_1D_MEcompareEta->GetName()));
	fPlots_1D_MEcompareRatio->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_1D_MEcompareRatio->GetName()));
	fPlots_1D_MEcompareErrors->Print(TString::Format("%s/%s.png",fOutputDir.Data(),fPlots_1D_MEcompareErrors->GetName()));
}
///**
/// Project the corrected histogram on the eta and phi axis to determine the NS width
///**
//________________________________________________________________________
void PlotGHcorrelation2::DetermineEtaWidths(TH2D* corrHistoSE[])
{
	//cout<<"inside of DetermineWidths()"<<endl;
	TF1* GaussFunc  = new TF1("GaussFunc",&PlotGHcorrelation2::PolyTwoGaussFitFunc,-100*DTR,300*DTR,11);
	TF1* GaussFunc1 = new TF1("GaussFunc1",&PlotGHcorrelation2::PolyTwoGaussFitFunc,-100*DTR,300*DTR,11);
	TF1* GaussFunc2 = new TF1("GaussFunc2",&PlotGHcorrelation2::PolyTwoGaussFitFunc,-100*DTR,300*DTR,11);
	TLegend* legEta;
	TLegend* legEta2;

	Bool_t plotVersion; //..two versions of plotting (0) without background subtracted, (1) with background subtracted
	plotVersion=0; //(0) without background subtracted
	//plotVersion=1;  // (1) with background subtracted
    Double_t BkGRange1=80*DTR;
    Double_t BkGRange2=4.6;
    Double_t SigRange1=-80*DTR;
    Double_t SigRange2=80*DTR;



  TCanvas * fLocalPlotDEta = new TCanvas("Corr1D_DEta_ProjFull_Indiv","Corr1D_DEta_ProjFull_Indiv",1750,1400);

	fCanCorr1D->Divide(3,3,0.001,0.001);
	fCanCorr1D_Sub->Divide(3,3,0.001,0.001);
	//for(Int_t i=0;i<fmaxBins;i++)
	for(Int_t i=0;i<fmaxBins;i++) //for the poster
	{
		fCanCorr1D->cd(i+1);
//		if(i>0)fCanCorr1D->cd(i);
//		if(i>2)fCanCorr1D->cd(i+1);
		//..project to the y-axis. (But only the near side.->Changed that)
//		corrHistoSE[i]->GetYaxis()->SetRangeUser(-1.6,1.6);
		corrHistoSE[i]->GetYaxis()->SetRangeUser(-fMaxDeltaEtaPlotRange,fMaxDeltaEtaPlotRange);
		TH1D *projY=corrHistoSE[i]->ProjectionY((const char*)Form("%s_projY",corrHistoSE[i]->GetName()),corrHistoSE[i]->GetXaxis()->FindBin(SigRange1),corrHistoSE[i]->GetXaxis()->FindBin(SigRange2));

    if (i <= kRebinDEtaThreshold) {
      projY->Rebin(kRebinDEta);
      projY->Scale(1./kRebinDEta);
    }

		//TH1D *projY=corrHistoSE[i]->ProjectionY((const char*)Form("%s_projY",corrHistoSE[i]->GetName()));
		if(projY->GetEntries()==0)continue;
		//SetTH1Histo(projY,Form("#Delta#eta^{%s-h}",fTriggerName.Data()),Form("dN^{%s-h}/N^{%s}",fTriggerName.Data(),fTriggerName.Data()),1);
//		SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),Form("dN^{%s_{can.}-h}/N^{%s_{can.}} #bullet factor",fTriggerName.Data(),fTriggerName.Data()),1);
//		SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),Form("factor #bullet dN^{%s_{can.}-h}/N^{%s_{can.}}",fTriggerName.Data(),fTriggerName.Data()),1);
//		SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),Form("dN^{%s_{can.}-h}/N^{%s_{can.}} (arb.units)",fTriggerName.Data(),fTriggerName.Data()),1);
//		SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),Form("1/N^{%s_{can.}} dN^{%s_{can.}-h}/d#Delta#eta^{%s_{can.}-h} (arb.units)",fTriggerName.Data(),fTriggerName.Data(),fTriggerName.Data()),1);
//		SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),Form("1/N^{%s_{can.}} dN^{%s_{can.}-h}/d#Delta#eta (arb.units)",fTriggerName.Data(),fTriggerName.Data()),1);
//		SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),Form("#frac{1}{N^{%s_{can.}}} dN^{%s_{can.}-h}/d#Delta#eta",fTriggerName.Data(),fTriggerName.Data()),1);
		//SetTH1Histo(projY,Form("#Delta#eta^{%s_{can.}-h}",fTriggerName.Data()),"#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#eta}",1);
		SetTH1Histo(projY,"#Delta#eta","#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#eta}",1);
		ZoomYRange(projY,0.3);
		//ZoomYRange(projY,0.3,-1.2,1.2);

		//..for performance
		// FIXME
//		projY->GetXaxis()->SetRangeUser(-0.8,0.8);
/*
		Double_t lowRange =projY->GetXaxis()->GetFirst();
		Double_t highRange=projY->GetXaxis()->GetLast();
		Double_t min = projY->GetMinimum();
		projY->GetXaxis()->SetRangeUser(-0.8,0.8);
		Double_t range1= projY->GetMaximum()-min;
		projY->Scale(1.0/(projY->GetMaximum()+range1*0.1));
		range1= projY->GetMaximum()-projY->GetMinimum();
		projY->GetYaxis()->SetRangeUser(projY->GetMinimum()-range1*0.1,1);
		projY->GetXaxis()->SetRange(lowRange,highRange);
*/

		projY->GetXaxis()->SetLabelSize(0.06);
		projY->GetYaxis()->SetLabelSize(0.06);
		projY->GetXaxis()->SetTitleSize(0.06);
		projY->GetYaxis()->SetTitleSize(0.04); //0.06
		projY->GetYaxis()->SetTitleOffset(1.1); //1.4
    if (bNoYLabel) projY->GetYaxis()->SetTitleOffset(0.71);
		projY->GetXaxis()->SetTitleOffset(1.1);
		projY->GetXaxis()->CenterTitle(kFALSE);
		projY->DrawCopy("E");

		fDeta_Proj[i] = (TH1D*) projY->Clone(Form("fDeta_Proj_%d",i));
		fDeta_Proj[i]->SetDirectory(0);

		//..try to model the sometimes squewed background underneath the delta eta peak
		TH1D *projY_BKG=corrHistoSE[i]->ProjectionY(Form("%s_ProjBackground",corrHistoSE[i]->GetName()),corrHistoSE[i]->GetXaxis()->FindBin(BkGRange1),corrHistoSE[i]->GetXaxis()->FindBin(BkGRange2));
    SetTH1Histo(projY_BKG,"","",1);
		projY_BKG->SetLineColor(kCyan-2);
		projY_BKG->SetMarkerColor(kCyan-2);
//		projY_BKG->SetMarkerSize(0.7);
		projY_BKG->SetMarkerStyle(kFullSquare);
    // somewhere, add projY_BKG to fDeta_AwaySide array

    if (i <= kRebinDEtaThreshold) {
      projY_BKG->Rebin(kRebinDEta);
      projY_BKG->Scale(1./kRebinDEta);
    }

		//..Normalize histograms to each other (excluding the NS eta peak)
    // M: is this necessary? a constant difference could just be fit?
		Double_t intA = projY->Integral(projY->FindBin(-1.2),projY->FindBin(-0.6));
		Double_t intB = projY_BKG->Integral(projY_BKG->FindBin(-1.2),projY_BKG->FindBin(-0.6));
		intA += projY->Integral(projY->FindBin(0.6),projY->FindBin(1.2));
		intB += projY_BKG->Integral(projY_BKG->FindBin(0.6),projY_BKG->FindBin(1.2));
    printf(" dEta   Scaling the away-side dEta projection by %f/%f=%f\n",intA,intB,intA/intB);

		projY_BKG->Scale(intA/intB);
		projY_BKG->DrawCopy("same E");

    // FIXME
    //  save projY_BKG fDeta_AwaySide
    // want to undo the rescale when fitting? or just clone before scale

	//	DrawAlicePerf(projY,0.22,0.81,0.12,0.12);

//		legEta = new TLegend(0.25,0.71,0.4,0.92);
		//legEta = new TLegend(0.12,0.44,0.35,0.71);  //0.24,0.6,0.4,0.78
		legEta = new TLegend(0.65,0.64,0.85,0.78);  //0.24,0.6,0.4,0.78
		if(fObservable==0)legEta->AddEntry(fsumCorrSE[i],Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
	  if(fObservable!=0)legEta->AddEntry((TObject*)0,Form("%0.1f < #it{p}_{T}^{%s} < %0.1f GeV/#it{c}",fArray_G_Bins[fPtMinBin-1],fTriggerName.Data(),fArray_G_Bins[fPtMaxBin]),"");
		if(fObservable==1)legEta->AddEntry(fsumCorrSE[i],Form("%0.2f < z_{T} < %0.2f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
		//if(fObservable==2)legEta->AddEntry(fsumCorrSE[i],Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
		if(fObservable==2)legEta->AddEntry(fsumCorrSE[i],Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");
		legEta->AddEntry(projY,Form("#Delta#varphi [%0.1f , %0.1f]",SigRange1,SigRange2),"pe");
		legEta->AddEntry(projY_BKG,Form("#Delta#varphi [%0.1f , %0.1f] (scaled)",BkGRange1,BkGRange2),"pe");
		legEta->SetTextColor(kBlack);
		legEta->SetTextSize(0.041); // 0.035
		legEta->SetBorderSize(0);
		//legEta->SetFillColorAlpha(10, 0);
    legEta->SetFillStyle(0);
		legEta->Draw("same");
		//DrawAlicePerf(projY,0.22,0.81,0.12,0.12);
		DrawAlicePerf(projY,0.62,0.81,0.12,0.12);

    // Now draw individually
    fLocalPlotDEta->cd();

    SetTH1Histo(projY,"","",0); // Settings for indiv plot
    SetTH1Histo(projY_BKG,"","",0); // Settings for indiv plot
    //projY->SetMarkerSize(2.5);
    //projY_BKG->SetMarkerSize(2.5);

		projY->DrawCopy("E");
		projY_BKG->DrawCopy("same E");

    gPad->SetTickx();
    gPad->SetTicky();

    // Settings for Individual DEta plot
		DrawAlicePerf(projY,0.13,0.73,0.12,0.12);
    legEta->Draw("SAME");

    fLocalPlotDEta->Print(TString::Format("%s/%s_Plot%d.pdf",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));
    fLocalPlotDEta->Print(TString::Format("%s/%s_Plot%d.png",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));
//    fLocalPlotDEta->Print(TString::Format("%s/%s_Plot%d.eps",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));
    fLocalPlotDEta->Print(TString::Format("%s/CFiles/%s_Plot%d.C",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));


		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		//..Subtract the modelled background
    // Drawing in the multiplot
		fCanCorr1D_Sub->cd(i+1);

		projY->Add(projY_BKG,-1);
		//ZoomYRange(projY,0.3,-1.2,1.2);
		ZoomYRange(projY,0.3);

    // FIXME
    SetTH1Histo(projY,"","",1); // Settings for large plot
    projY->SetMarkerSize(0.8);
		projY->DrawCopy("E");

		//..fit, draw and save widths
		FitGaussAndDraw(i+1,projY,GaussFunc,GaussFunc1,GaussFunc2,0);

    bool fDebugEtaFit = 1;

		fDeta_ProjSub[i] = (TH1D*) projY->Clone(Form("fDeta_NearSideProjSub_%d",i));
		fDeta_ProjSub[i]->SetDirectory(0);

    fDeta_AwaySide[i] = (TH1D*) projY_BKG->Clone(Form("fDeta_AwaySideProj_%d",i));
    fDeta_AwaySide[i]->SetDirectory(0);
    fDeta_AwaySide[i]->Scale(intB/intA); // Undoing nearside matching scale

    if (fDebugEtaFit) legEta2 = new TLegend(0.25,0.45,0.4,0.92); //..Bkg subtracted
		else legEta2 = new TLegend(0.25,0.67,0.4,0.92); //..Bkg subtracted

		legEta2->AddEntry(projY,Form("Bkg. sub. proj. in #Delta#varphi [#pi/2-3/2#pi]"),"pe");
		if(fObservable==0)legEta2->AddEntry(fsumCorrSE[i],Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[i],fTriggerName.Data(),fArray_G_Bins[i+1]),"");
		if(fObservable==1)legEta2->AddEntry(fsumCorrSE[i],Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]),"");
		//if(fObservable==2)legEta2->AddEntry(fsumCorrSE[i],Form("%0.1f < #xi < %0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]),"");
		if(fObservable==2)legEta2->AddEntry(fsumCorrSE[i],Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[i],fArray_HPT_Bins[i+1]),"");

    legEta2->AddEntry(GaussFunc1,"Sum of two Gauss","l");

	  //legEta2->AddEntry(GaussFunc1,Form("#mu_{1}: %0.3f",GaussFunc->GetParameter(1)),"l");
		//legEta2->AddEntry(GaussFunc,Form("#mu_{1}: %0.3f",GaussFunc->GetParameter(1)),"l");
    if (fDebugEtaFit) legEta2->AddEntry(GaussFunc2,Form("Y_{1}: %0.3e #pm %.3e",GaussFunc->GetParameter(0),GaussFunc->GetParError(0)),"l");
    if (fDebugEtaFit) legEta2->AddEntry(GaussFunc1,Form("Y_{2}/Y_{1}: %0.3e #pm %.3e",GaussFunc->GetParameter(3),GaussFunc->GetParError(3)),"l");


		legEta2->AddEntry(GaussFunc1,Form("#sigma_{1}: %0.3f #pm %.3f",GaussFunc->GetParameter(2),GaussFunc->GetParError(2)),"l");
		legEta2->AddEntry(GaussFunc1,Form("#sigma_{2}/#sigma_{1}: %0.3f #pm %.3f",GaussFunc->GetParameter(5),GaussFunc->GetParError(5)),"l");
		legEta2->AddEntry(GaussFunc,Form("#chi^{2}/NDF: %0.3f",GaussFunc->GetChisquare()/GaussFunc->GetNDF()),"");

		//legEta2->AddEntry(GaussFunc,Form("#sigma_{1}: %0.3f",GaussFunc->GetParameter(2)),"l");
		legEta2->SetTextColor(kBlack);
		legEta2->SetTextSize(0.05);
		legEta2->SetBorderSize(0);
		//legEta2->SetFillColorAlpha(10, 0);
    legEta2->SetFillStyle(0);
		legEta2->Draw("same");
//		DrawWIP(fsumCorrSE[i],0.55,0.5,0.35,0.1);


    // Now Draw individually
    fLocalPlotDEta->cd();
    SetTH1Histo(projY,"#Delta#eta","#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#eta}",0); // Settings for indiv plot
		projY->DrawCopy("E");
		//FitGaussAndDraw(i+1,projY,GaussFunc,GaussFunc1,GaussFunc2,0);
    GaussFunc->Draw("SAME");
    GaussFunc1->Draw("SAME");
    GaussFunc2->Draw("SAME");
		legEta2->Draw("same");

    fLocalPlotDEta->Print(TString::Format("%s/%s_FitPlot%d.pdf",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));
    fLocalPlotDEta->Print(TString::Format("%s/%s_FitPlot%d.png",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));
    fLocalPlotDEta->Print(TString::Format("%s/CFiles/%s_FitPlot%d.C",fOutputDir.Data(),fLocalPlotDEta->GetName(),i));

	}
}
//________________________________________________________________________
void PlotGHcorrelation2::DeterminePhiWidths()
{
	//DeterminePhiWidths(TH2D* corrHistoSE[])
	//..!!!Do not do that here - will not yield resonable results
	//..!!! peak in delta phi is jet+flow ->flow widens the peak
	/*fCanCorr1D->cd(i*2+2);
	TH1D *projX=corrHistoSE[i]->ProjectionX();
	SetTH1Histo(projX,"","dN^{#gamma-h}/dN^{#gamma}",1);
	ZoomYRange(projX);
	projX->DrawCopy("E");
	FitGaussAndDraw(i+1,projX,GaussFunc,GaussFunc1,GaussFunc2,1);
	 */
}
///***
/// Fit a Gauss function to the near side peak
/// to determine its width
///***
//________________________________________________________________________
void PlotGHcorrelation2::FitGaussAndDraw(Int_t bin,TH1D* corrProjHistoSE,TF1* Func1, TF1* Func2,TF1* Func3,Bool_t EtaPhi)
{
  // Idea: set Widthmin = bin width

  double fWidthMin=0.01;

  double fWidthBinScale = 0.25;
  fWidthMin = fWidthBinScale * corrProjHistoSE->GetBinWidth(corrProjHistoSE->FindBin(0.0));

  double fWidthMax=0.5;
  printf("Debug FitGausAndDraw Setting width min,max to %f,%f\n",fWidthMin,fWidthMax);

	for(Int_t g = 0; g < 11; g++)
	{
		Func1->ReleaseParameter(g);
		Func1->SetParameter(g,0);
		Func1->SetParError(g,0.0);
	}
	Func1->SetParameter(1,0); //..start value for mean1
	Func1->SetParameter(4,0); //..start value for mean2

  // Naming
  Func1->SetParName(0,"amp_1");
  Func1->SetParName(1,"mu_1");
  Func1->SetParName(2,"sigma_1");
  Func1->SetParName(3,"amp_1");
  Func1->SetParName(4,"mu_1");
  Func1->SetParName(5,"sigma_2/sigma_1");
  Func1->SetParName(6,"a");

	if(EtaPhi==0)//..Eta case
	{
    Bool_t bFixDEtaPeak=1;

		//- - - - - - - - - - - - - - - -
		//..Flat background
		//..for eta -> estimate the background level on eta=-0.7&+0.7
		Double_t backgroundLevel=corrProjHistoSE->Integral(corrProjHistoSE->FindBin(-0.7),corrProjHistoSE->FindBin(-0.5));
		backgroundLevel+=corrProjHistoSE->Integral(corrProjHistoSE->FindBin(0.5),corrProjHistoSE->FindBin(0.7));
		//..Divide by number of bins
		Double_t bins=corrProjHistoSE->FindBin(-0.5)-corrProjHistoSE->FindBin(-0.7);
		bins+=corrProjHistoSE->FindBin(0.7)-corrProjHistoSE->FindBin(0.5);
		backgroundLevel*=1/bins;
		//cout<<"Integral: "<<backgroundLevel*bins<<", n bins= "<<bins<<", average bgk: "<<backgroundLevel<<endl;
		Func1->SetParameter(6,backgroundLevel);
		Func1->SetParLimits(6,backgroundLevel*0.9,backgroundLevel*1.1);  //..allow a variation of +-10%

		//- - - - - - - - - - - - - - - -
		//..big, narrow gaussian
		//..etimate amplitude by 0 heigth and background level
		//Double_t amplEst=corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(0));
		Double_t amplEst=TMath::Abs(corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(0)));
		amplEst+=TMath::Abs(corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(0)-1));
		amplEst+=TMath::Abs(corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(0)+1));
    // FIXME at lowest statistics, amplEst can be negative!
	//	amplEst-=3*backgroundLevel;
    amplEst = amplEst / 3.0;
    double widthEst = 0.45;

    printf("Debug eta width fit: using amp est %f\n",amplEst);
		//amplEst-=backgroundLevel;
		Func1->SetParameter(0,amplEst);    //..amplitude
		Func1->SetParameter(2,widthEst);       //..width
		Func1->SetParLimits(0,amplEst*0.9,amplEst*1.1);
    if (bFixDEtaPeak) { 
      Func1->FixParameter(1,0.0);
    } else {
      Func1->SetParLimits(1,-0.1,0.1);  //..mean limits
    }
		Func1->SetParLimits(2,fWidthMin,fWidthMax);  //..width limits

		//- - - - - - - - - - - - - - - -
		//..small, wide gaussian
		Func1->SetParameter(3,0.05);       //..amplitude (ratio to amplitude of main peak)
		Func1->SetParameter(5,1.1);        //..width
		Func1->SetParLimits(3,0.05,0.5);   //..amplitude limits 5-50% of the main peak (ampl2 = param0*param3)
    if (bFixDEtaPeak) { 
      Func1->FixParameter(4,0.0);   //..fixing mean parameter at 0
    } else {
      Func1->SetParLimits(4,-0.1,0.1);   //..mean limits
    }

    if (bin == 0) { // Disabling lowest bin, where we always fail
      Func1->FixParameter(2,0.5);


    }



		Func1->SetParLimits(5,1.05,3.0);   //..width limits 105%-300% of the big one (width2 = param2*param5)

		// Fixing the background to zero for the eta case
		Func1->FixParameter(6,0);

		//- - - - - - - - - - - - - - - -
		//..plot range for eta projection
		Func1->SetRange(-1.0,1.0);
		Func2->SetRange(-1.0,1.0);
		Func3->SetRange(-1.0,1.0);
	}

	if(EtaPhi==1)//..phi case
	{
		//- - - - - - - - - - - - - - - -
		//..Flat background
		//..for phi -> estimate the background level on phi=-85&+85
		Double_t backgroundLevel=corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(-85));
		backgroundLevel+=corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(+85));
		backgroundLevel*=0.5;
		//cout<<"Background level="<<backgroundLevel<<endl;
		Func1->SetParameter(6,backgroundLevel);
		Func1->SetParLimits(6,backgroundLevel*0.9,backgroundLevel*1.1);  //..allow a variation of +-10%

		//- - - - - - - - - - - - - - - -
		//..big, narrow gaussian
		//..etimate amplitude by 0 heigth and background level
		Double_t amplEst=corrProjHistoSE->GetBinContent(corrProjHistoSE->FindBin(0));
		amplEst-=backgroundLevel;
		Func1->SetParameter(0,amplEst);   //..amplitude
		Func1->SetParameter(2,0.15);      //..width
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
  //
	//Func1->FixParameter(3,0); //..Use only 1 gaussian for fitting

	Func1->SetLineColor(kDEtaFitColor);
	Func1->SetLineStyle(kDEtaFitStyle); //2



	TString Name= Func1->GetName();
	if(EtaPhi==0)corrProjHistoSE->Fit(Name,"","",-0.7,0.7); //Eta case Q = quiet mode, no printout
	if(EtaPhi==1)corrProjHistoSE->Fit(Name,"Q","",-70,70);   //Phi case Q = quiet mode, no printout
	if(EtaPhi==1)cout<<"Error - change to rad"<<endl;
	//..width that is used to define an eta range not contaminated by the near side peak
	//..you can define the width in multiple ways.
	//..THINK ABOUT THAT!
	//Width   =Func1->GetParameter(5)*Func1->GetParameter(2);//bigger width in delta phi
  // FIXME parameter 2 is usually the smaller width
	Double_t Width   =Func1->GetParameter(2); //bigger width in delta phi
	Double_t ErrWidth=Func1->GetParError(2);
	if(EtaPhi==0)//..eta case
	{
		fEtaWidth->SetBinContent(bin,Width);
		fEtaWidth->SetBinError(bin,ErrWidth);
	}
	if(EtaPhi==1)//..phi case
	{
		fPhiWidth->SetBinContent(bin,Width);
		fPhiWidth->SetBinError(bin,ErrWidth);
	}
	//  cout<<"Width gauss2: "<<Func1->GetParameter(5)<<" times of width 1"<<endl;

	for(Int_t g = 0; g < 11; g++)
	{
		Func2->SetParameter(g,Func1->GetParameter(g));
		Func3->SetParameter(g,Func1->GetParameter(g));
	}
	//..small, wide gaussian
	//..due to the fact that param 0 and param 3 are proportional
	//..we do a little hack here. Setting param0 to 0 is necessary
	//..to see only the small wide gaussian. If we set param0 to 0
	//..however, param3 will become 0 by defualt. We can however
	//..set param0 to a negligibly small value x and multiply param3
	//..by the inverse of x. (normally param3 is in the range 0-1, but we omit this for this specific case)
	Double_t Shrinkage=0.00001;
	//Func2
	Func2->SetParameter(0,Shrinkage);
	Func2->SetParameter(3,1.0*Func1->GetParameter(3)*Func1->GetParameter(0)/Shrinkage);
	Func2->SetLineColor(kDEtaFitWideColor); //kPink-9);
  Func2->SetLineStyle(kDEtaFitWideStyle);

	if(Func1->GetParameter(3)!=0)Func2 ->DrawCopy("same"); //..only when the small-broad gaussian is not set to 0

	//..big, narrow gaussian (green)
	Func3->SetParameter(3,0);
	//Func3->SetLineColor(kCyan-2);
	Func3->SetLineColor(kDEtaFitThinColor);
	Func3->SetLineStyle(kDEtaFitThinStyle);
	Func3 ->DrawCopy("same");

	//..flat background
  // This was why the parameters appeared to be 0 later.
  // Do I need the flat background drawn? This would then be Func4
	//Func1->SetParameter(0,0);
	//Func1->SetParameter(3,0);
	//Func1 ->DrawCopy("same");

	//..Draw legend
	/*Obsolete as of Feb28
	TString TopLegendText;
	Double_t x_Pos;
	if(EtaPhi==0)x_Pos=0.2;//Eta case
	if(EtaPhi==1)x_Pos=0.6;//Phi case
	TopLegendText=Form("#mu_{1}: %0.3f",Func1->GetParameter(1));
	PlotTopLegend((const char*)TopLegendText,x_Pos,0.85,0.07,kGreen-2);
	TopLegendText=Form("#sigma_{1}: %0.3f",Func1->GetParameter(2));
	PlotTopLegend((const char*)TopLegendText,x_Pos,0.78,0.07,kGreen-2);

	if(Func1->GetParameter(3)!=0)
	{
		TopLegendText=Form("#mu_{2}: %0.3f",Func1->GetParameter(4));
		PlotTopLegend((const char*)TopLegendText,x_Pos,0.71,0.07,kPink-9);
		TopLegendText=Form("#sigma_{2}: %0.2f",Func1->GetParameter(5)*Func1->GetParameter(2));
		PlotTopLegend((const char*)TopLegendText,x_Pos,0.64,0.07,kPink-9);
	}*/
}
///**
/// Draw the determined widths of the NS peaks in the
/// different bins of the analysis
///**
//________________________________________________________________________
void PlotGHcorrelation2::DrawWidths()
{
	//..Draw a canvas with the eta and phi widths
//	fCanvWidth->Divide(2);
	//fCanvWidth->cd(1);
	fCanvWidth->cd();
	if(fObservable==0)SetTH1Histo(fEtaWidth,Form("E_{%s}",fTriggerName.Data()),"#sigma of #Delta#eta");
	//	if(fObservable==0)SetTH1Histo(fEtaWidth,"E_{#gamma}","#sigma of #Delta#eta");
	if(fObservable==1)SetTH1Histo(fEtaWidth,"z_{T}","#sigma of #Delta#eta");
	//if(fObservable==2)SetTH1Histo(fEtaWidth,"#xi","#sigma of #Delta#eta");
	if(fObservable==2)SetTH1Histo(fEtaWidth,"#it{p}_{T}^{assoc} (GeV/#it{c})","#sigma of #Delta#eta");
	//if(fObservable==0)fEtaWidth->GetXaxis()->SetRangeUser(0,1);
	if(fObservable==1)fEtaWidth->GetXaxis()->SetRangeUser(0,1);
	if(fObservable==2)fEtaWidth->GetXaxis()->SetRangeUser(0,12);
	fEtaWidth->SetMarkerColor(17);
	fEtaWidth->SetMarkerStyle(kFullSquare);
	fEtaWidth->SetMarkerSize(1.5);
	fEtaWidth->DrawCopy("E");
	fEtaWidth->SetLineColor(kRed-1);
	fEtaWidth->SetMarkerColor(kRed-1); //kTeal-3
	fEtaWidth->DrawCopy("same E");
  gPad->SetLogy();
	/*fCanvWidth->cd(2);
	if(fObservable==0)SetTH1Histo(fPhiWidth,"E_{#gamma}","#sigma of #Delta#phi");
	if(fObservable==1)SetTH1Histo(fPhiWidth,"z_{T}","#sigma of #Delta#phi");
	if(fObservable==2)SetTH1Histo(fPhiWidth,"#xi","#sigma of #Delta#phi");
	//if(fObservable==0)fPhiWidth->GetXaxis()->SetRangeUser(0,1);
	if(fObservable==1)fPhiWidth->GetXaxis()->SetRangeUser(0,1);
	if(fObservable==2)fPhiWidth->GetXaxis()->SetRangeUser(0,2);
	fPhiWidth->SetMarkerColor(kTeal-3);
	fPhiWidth->SetMarkerStyle(20);
	fPhiWidth->SetMarkerSize(1.1);
	fPhiWidth->DrawCopy("E");
	fPhiWidth->SetLineColor(kTeal-3);
	fPhiWidth->SetMarkerColor(kSpring+8);
	fPhiWidth->SetMarkerSize(0.8);
	fPhiWidth->DrawCopy("same E");*/
}
///**
/// Calls the project and fit function for every bin
///**
//________________________________________________________________________
void PlotGHcorrelation2::FitEtaSides(TH2D* corrHistoSE[],Double_t sigma1,TCanvas* can1,TCanvas* can2,TCanvas* can3)
{
	Double_t width=0;

	can1->Divide(4,4,0.001,0.0012);
	can2->Divide(3,3,0.001,0.0012);
	can3->Divide(3,3,0.001,0.0012);
	for(Int_t i=0;i<fmaxBins;i++) // why was this fmaxbins - 1? just for zt?
	{

    // FIXME add adjustable scale.? Or should that be done in FitEtaSide?

		width=fEtaWidth->GetBinContent(i+1);
		FitEtaSide(corrHistoSE[i],width,sigma1,can1,can2,can3,i);
    // FIXME save the width to a TTree to pass down.
    // array branch with a branch nObsBins.
	}
}
///**
/// After you have determined the projection ranges
/// for the near side fit you can update the eta distribution
/// with the determined ranges
///**
//________________________________________________________________________
void PlotGHcorrelation2::UpdateEtaCanvas(Int_t bin,Double_t edgeLow,Double_t edgeSigL,Double_t edgeSigR,Double_t edgeHigh)
{
	TVirtualPad *pad=fCanCorr1D_Sub->cd(bin+1);
	//..plot 3 sigma range (enough to split off near side peak)
	Double_t lowYrange=0;
	Double_t highYrange=0;

	TString name;
	name=fsumCorrSE[bin]->GetName();
	name+="_projY"; //..!!THis is hard coded stuff - not very elegant!!

	TH1* test=(TH1*)pad->GetPrimitive(name);
	//cout<<"found histo with name: "<<test->GetName()<<endl;

	//..after adding dividing etc. min and max value are unfortunatley
	//..recalculated. This is important for
	/*Double_t currentMin = projY->GetMinimumStored();
	 Double_t currentMax = projY->GetMaximumStored();
	 if((currentMin == -1111) || (currentMax == -1111))
	 {
		 projY->TH1::GetMinimumAndMaximum(currentMin, currentMax);
		 projY->SetMinimum(currentMin);
		 projY->SetMaximum(currentMax);
	 }
	 */
	//..works for unsubtracted
	lowYrange =test->GetMinimumStored();
	highYrange=test->GetMaximumStored();
	//..
	//highYrange=test->GetBinContent(test->FindBin(0))*2;
	//lowYrange =-test->GetBinContent(test->FindBin(0));
	//test->GetYaxis()->SetRangeUser(lowYrange,highYrange);
	//cout<<"high range: "<<highYrange<<", low range: "<<lowYrange<<endl;

	PlotVerLineRange(edgeSigL,lowYrange,highYrange,17);
	PlotVerLineRange(edgeSigR,lowYrange,highYrange,17);

	PlotVerLineRange(edgeLow,lowYrange,highYrange,17);
	PlotVerLineRange(edgeHigh,lowYrange,highYrange,17);
}
///**
/// This function finds the last good bin in eta direction that
/// adds statistic without blowing up the errors.
/// We look for the error/value (relative error) to not grow when
/// adding an extra eta bin.
///**
//________________________________________________________________________
void PlotGHcorrelation2::FindLastGoodBin(TH2D* corrHistoSE,Int_t startSearch,Int_t direction,Int_t *lastGoodBin)
{
	if(direction!=-1 && direction!=1) cout<<"Major error in FindLastGoodBin. direction needs to be +1/-1"<<endl;
	//..Project to the delta eta axis
	//..unzoom whatever was zoomed before
	corrHistoSE->GetYaxis()->UnZoom();
	//..sum over all DeltaPhi bins and project onto the eta axis
	TH1D *PprojEtaAxis = corrHistoSE->ProjectionY("Test projection");
	Double_t cont,error;
	Double_t contSum=0;
	Double_t contSumOld=0;
	Double_t errorSum=0;
	Double_t errorSumOld=0;

	Double_t relErrOld=100000000000;
	Double_t relErrAddingNew=0;

	//..Find the lowest/highest good bin
	//..look at the projection onto the eta axis to
	//..see whether there is a bin which has very large accumulated errors
	for(Int_t binLoop=0;binLoop<25;binLoop++)
	{
		cont =PprojEtaAxis->GetBinContent(startSearch+(direction*binLoop));
		if(cont==0)continue;
		error=PprojEtaAxis->GetBinError(startSearch+(direction*binLoop));
		contSum+=cont;
		errorSum=sqrt(pow(errorSumOld,2)+pow(error,2));
		//cout<<"------loop: "<<binLoop<<", poision"<<PprojEtaAxis->GetBinCenter(startSearch+(direction*binLoop))<<"-------"<<endl;
		//cout<<"content: "<<cont<<", err: "<<error<<", rel. err: "<<error/cont*1000<<"*10-3"<<endl;
		//cout<<"Sum content: "<<contSum<<", rel. err: "<<errorSum/contSum*1000<<"*10-3"<<endl;

		//..calculate how the relative error would be
		//..adding this additional bin to the integral
		if(contSum!=0)relErrAddingNew=errorSum/contSum;
		else cout<<"Error in FindLastGoodBin :: Divinding by 0"<<endl;

		//..If the relative error gets smaller
		//..whenn adding the new bin continue with
		//..finding the lowest good bin
		if(relErrAddingNew<relErrOld)
		{
			errorSumOld=errorSum;
			contSumOld =contSum;
			relErrOld  =errorSumOld/contSumOld;
		}
		else
		{
			//cout<<"----This bin will not be added!------"<<endl;
			//..if the relative error of the integated projection
			//..would get bigger adding the new eta bin
			//..Dont do it! The previous bin was the last bin over
			//..which one should integrate in the projeciton
			*lastGoodBin=startSearch+(direction*(binLoop-1));
			break;
		}
	}
}
///
///
///
//________________________________________________________________________
void PlotGHcorrelation2::FitEtaSide(TH2D* corrHistoSE,Double_t width,Double_t Sigma2,TCanvas* Can,TCanvas* Can2,TCanvas* Can3,Int_t CanvasPad)
{
	TString TopLegendText;
	TString ProjectionName;
	//TF1* BackgroundFunction = new TF1("FlowFunction",FlowFunction,-100,300,4);
	TF1* BackgroundFunction = new TF1("FlowFunction",FlowFunction,-100*DTR,300*DTR,4);
	Int_t vN=3;
	//TF1 *allFit             = allFitFuncVn("JoelsFlowFunction",vN,-100,300);
	TF1 *allFit             = allFitFuncVn("JoelsFlowFunction",vN,-100*DTR,300*DTR);

  TCanvas * fLocalPlotBoth = new TCanvas("Corr1D_ProjBoth_Indiv","Corr1D_ProjBoth_Indiv",1750,1400);
  TCanvas * fLocalPlotFarDeltaEta = new TCanvas("Corr1D_FarDeltaEta_Indiv","Corr1D_ProjFarDeltaEta_Indiv",1750,1400);

  // Could use the fact that CanvasPad corresponds to the observable bin to use bin dependent fixed widths
  Double_t SignalEtaWidth = width*Sigma2;

  double epsilon = 0.000001; // To avoid asymmetry from binning

  // Base on iFixedDEtaCutIndex 
  // 0 - > sigma?
  // 1 - > 0.8, 1.35
  // 2 - > 0.7, 1.35
  // 3 - > 0.8, 1.2
  // 4 - > 0.7, 1.2
  // 5 - > 0.8, 1.45
  // 6 - > 0.7, 1.45
  // 7 - > 0.9, 1.45
  // 8 - > 0.9, 1.5
  switch (iFixedDEtaCutIndex) {
    case 8:
      fMaxDEtaSignalRange = 0.9 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.5  - epsilon;
      break;
    case 7:
      fMaxDEtaSignalRange = 0.9 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.45  - epsilon;
      break;
    case 6:
      fMaxDEtaSignalRange = 0.7 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.45  - epsilon;
      break;
    case 5:
      fMaxDEtaSignalRange = 0.8 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.45 - epsilon;
      break;
    case 4:
      fMaxDEtaSignalRange = 0.7 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.2  - epsilon;
      break;
    case 3:
      fMaxDEtaSignalRange = 0.8 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.2  - epsilon;
      break;
    case 2:
      fMaxDEtaSignalRange = 0.7 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.35  - epsilon;
      break;
    case 1:
      fMaxDEtaSignalRange = 0.8 + epsilon;
      fMinDEtaSignalRange = fMaxDEtaSignalRange;
      fMaxDeltaEtaRange = 1.35  - epsilon;
      break;
    default:
    case 0:
      break;
  }




  if (SignalEtaWidth > fMaxDEtaSignalRange) SignalEtaWidth = fMaxDEtaSignalRange;
  if (SignalEtaWidth < fMinDEtaSignalRange) SignalEtaWidth = fMinDEtaSignalRange;


  if (fFixedDEtaCut > 0) {
    SignalEtaWidth = fFixedDEtaCut;
    //SignalEtaWidth = fFixedDEtaCut+0.00001; // magic
  }
	//Double_t SignalEtaRange = 2.*width*Sigma2; 
  Double_t SignalEtaRange = 2*SignalEtaWidth;
  // fMaxDEtaSignalRange

  // FIXME make sigma2 depend on Z? or set a maximum on SignalEtaRange?
  // Rewrite following code to strictly use SignalEtaRange/2. and apply a minimum

	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//.. .. select the eta range in which you want to project your signal.. ..
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	Int_t SignalCenter     = corrHistoSE->GetYaxis()->FindBin(0.0);
	Int_t SignalEdgeLow    = corrHistoSE->GetYaxis()->FindBin(0-SignalEtaWidth);
	Int_t SignalEdgeHigh   = corrHistoSE->GetYaxis()->FindBin(0+SignalEtaWidth);
	//..Test the eta sideband window
	Double_t variation=0.3;
	Int_t SignalEdgeLowT1  = corrHistoSE->GetYaxis()->FindBin(0-SignalEtaWidth*(1+variation));//..change sigma range to check how much signal is inside this eta band
	Int_t SignalEdgeHighT1 = corrHistoSE->GetYaxis()->FindBin(0+SignalEtaWidth*(1+variation));//..change sigma range to check how much signal is inside this eta band
	Int_t SignalEdgeLowT2  = corrHistoSE->GetYaxis()->FindBin(0-SignalEtaWidth*(1-variation));//..change sigma range to check how much signal is inside this eta band
	Int_t SignalEdgeHighT2 = corrHistoSE->GetYaxis()->FindBin(0+SignalEtaWidth*(1-variation));//..change sigma range to check how much signal is inside this eta band
	Int_t LowestBin=0;
	Int_t HighestBin=0;


	// FIXME

//	Float_t fMaxDeltaEtaRange = 1.2-0.001; // epsilon = 0.001 to avoid bin effects
	if (fUseFindLastGoodBin) {
		//..This function finds the last good bin in the delta eta axis for
		//..a good projection without overblown errors.
		FindLastGoodBin(corrHistoSE,SignalEdgeLow,-1,&LowestBin);
		FindLastGoodBin(corrHistoSE,SignalEdgeHigh,+1,&HighestBin);
	} else {
		LowestBin = corrHistoSE->GetYaxis()->FindBin(-fMaxDeltaEtaRange);
		HighestBin = corrHistoSE->GetYaxis()->FindBin(fMaxDeltaEtaRange);
	}

	//..check that the mean+-sigma is not larger or smaller than the histogram range
	if(SignalEdgeLow<LowestBin || SignalEdgeHigh>HighestBin)
	{
		cout<<"Error: Problem detected!"<<endl;
		cout<<"In Histo: "<<corrHistoSE->GetName()<<endl;
		cout<<"Signal range is reaching outside the histogram boundaries - please correct"<<endl;
		cout<<"bins lowes"<<LowestBin<<", edge "<<SignalEdgeLow<<", center"<<SignalCenter <<", uppedge"<< SignalEdgeHigh<<", highestbin "<< HighestBin<<endl;
		SignalEdgeLow=LowestBin;
		SignalEdgeHigh=HighestBin;
	}
//	Double_t fitRangeL1=corrHistoSE->GetYaxis()->GetBinCenter(LowestBin)-0.5*corrHistoSE->GetYaxis()->GetBinWidth(LowestBin);
//	Double_t fitRangeL2=corrHistoSE->GetYaxis()->GetBinCenter(SignalEdgeLow)+0.5*corrHistoSE->GetYaxis()->GetBinWidth(SignalEdgeLow);
//	Double_t fitRangeR1=corrHistoSE->GetYaxis()->GetBinCenter(SignalEdgeHigh)-0.5*corrHistoSE->GetYaxis()->GetBinWidth(SignalEdgeHigh);
//	Double_t fitRangeR2=corrHistoSE->GetYaxis()->GetBinCenter(HighestBin)+0.5*corrHistoSE->GetYaxis()->GetBinWidth(HighestBin);
  // Alternate Version
  
	Double_t fitRangeL1=corrHistoSE->GetYaxis()->GetBinLowEdge(LowestBin);
	Double_t fitRangeL2=corrHistoSE->GetYaxis()->GetBinUpEdge(SignalEdgeLow);
	Double_t fitRangeR1=corrHistoSE->GetYaxis()->GetBinLowEdge(SignalEdgeHigh);
	Double_t fitRangeR2=corrHistoSE->GetYaxis()->GetBinUpEdge(HighestBin);
  
  
  cout<<"In Histo: "<<corrHistoSE->GetName()<<endl;
	cout<<"low range:"<<fitRangeL1<<"-"<<fitRangeL2<<", high range: "<<fitRangeR1<<"-"<<fitRangeR2<<endl;
	//..draw the projection ranges into the eta canvas
	UpdateEtaCanvas(CanvasPad,fitRangeL1,fitRangeL2,fitRangeR1,fitRangeR2);

	Double_t SideBandEtaRange = fitRangeL2 - fitRangeL1 + fitRangeR2 - fitRangeR1;
	Double_t RealSignalEtaRange = fitRangeR1 - fitRangeL2;  // Accounting for binning effects

	cout<<"Signal   range ="<<SignalEtaRange<<"."<<endl;
	cout<<"RealSig  range ="<<RealSignalEtaRange<<"."<<endl;
	cout<<"Sideband range ="<<SideBandEtaRange<<"."<<endl;

	//-------------------------------------------------------------------------------------------
	//.. determine a scale factor to fit NS to away side. This is necessary since the flow .. ..
	//.. function is only fit to the side band but needs to be scaled to the whole eta range
	//-------------------------------------------------------------------------------------------
  /*
	Int_t nBinsSB = (SignalEdgeLow-LowestBin)+(HighestBin-SignalEdgeHigh)+2; //..add 2 because the projection INCLUDES lowes and highest bin
	Int_t nBinsNS =  SignalEdgeHigh-SignalEdgeLow-1;//..subtract 1 because projection is done until the edges not with them
	fscaleFactorSBtoNS[CanvasPad]=1.0*nBinsSB/nBinsNS;  //IS this used? -MO
	//cout<<"^-^"<<"bins SB: "<<nBinsSB<<", bins NS: "<<nBinsNS<<", scaling: "<<fscaleFactorSBtoNS[CanvasPad]<<endl;
*/

	//-------------------------------------------------------------------------------------------
	//.. Project the 2D histogram to the Delta Phi axis ..
	//-------------------------------------------------------------------------------------------
	ProjectionName=Form("%s_Pproj%1fSig",corrHistoSE->GetName(),Sigma2*10);
	TH1D *PprojXSig  =corrHistoSE->ProjectionX(ProjectionName,SignalEdgeLow+1,SignalEdgeHigh-1); //..do not include the edges
  // FIXME is this where some delta eta norm is going bad??
	ProjectionName=Form("%s_Pproj%1fSigSide1",corrHistoSE->GetName(),Sigma2*10);
	TH1D *PprojXSide1=corrHistoSE->ProjectionX(ProjectionName,LowestBin,SignalEdgeLow);      //..include the edges
	ProjectionName=Form("%s_Pproj%1fSigSide2",corrHistoSE->GetName(),Sigma2*10);
	TH1D *PprojXSide2=corrHistoSE->ProjectionX(ProjectionName,SignalEdgeHigh,HighestBin);    //..include the edges
	//.. Test part increase by 10%
	ProjectionName=Form("%s_Pproj%1fSigSide1_Test1",corrHistoSE->GetName(),Sigma2*10*(1+variation));
	TH1D *PprojXSide1_T1=corrHistoSE->ProjectionX(ProjectionName,LowestBin,SignalEdgeLowT1);   //..include the edges
	ProjectionName=Form("%s_Pproj%1fSigSide2_Test1",corrHistoSE->GetName(),Sigma2*10*(1+variation));
	TH1D *PprojXSide2_T1=corrHistoSE->ProjectionX(ProjectionName,SignalEdgeHighT1,HighestBin); //..include the edges
	//.. Test part reduce by 10%
	ProjectionName=Form("%s_Pproj%1fSigSide1_Test2",corrHistoSE->GetName(),Sigma2*10*(1-variation));
	TH1D *PprojXSide1_T2=corrHistoSE->ProjectionX(ProjectionName,LowestBin,SignalEdgeLowT2);   //..include the edges
	ProjectionName=Form("%s_Pproj%1fSigSide2_Test2",corrHistoSE->GetName(),Sigma2*10*(1-variation));
	TH1D *PprojXSide2_T2=corrHistoSE->ProjectionX(ProjectionName,SignalEdgeHighT2,HighestBin); //..include the edges
    //.. Do a full projection for an intermediate result
	ProjectionName=Form("%s_FullProjection",corrHistoSE->GetName());
	//TH1D *PprojXFull=corrHistoSE->ProjectionX(ProjectionName,corrHistoSE->GetYaxis()->FindBin(-maxDeltaEtaRange),corrHistoSE->GetYaxis()->FindBin(maxDeltaEtaRange));   //..include the edges
	TH1D *PprojXFull=corrHistoSE->ProjectionX(ProjectionName,LowestBin,HighestBin);   //..include the edges

		
  // The step of normalizing by the eta range
  printf("   Scaling full dEta range by 1 / %f\n",(fitRangeR2-fitRangeL1));
	PprojXFull->Scale(1./(fitRangeR2-fitRangeL1));  // FIXME is this the correct eta range?
//	PprojXFull->Scale(0.5/maxDeltaEtaRange); // 1 / (2 * maxDeltaEtaRange)  // this would not account for bins
  printf("      Scaling again for bin width 1 / %f\n",PprojXFull->GetXaxis()->GetBinWidth(3));
	PprojXFull->Scale(1./PprojXFull->GetXaxis()->GetBinWidth(3)); // Normalize to dDeltaPhi bin
	PprojXFull->SetName(Form("dPhi_ObsBin%d_Full",CanvasPad));

	fsumCorrSE_ProjFull[CanvasPad] = PprojXFull;


	//-------------------------------------------------------------------------------------------
	//.. Draw NS signal region ..
	//-------------------------------------------------------------------------------------------
	Can->cd(CanvasPad*2+1);
	SetTH1Histo(PprojXSig,"",Form("d^{2}N^{%s-h}/N^{%s}d#Delta#eta d#Delta#phi",fTriggerName.Data(),fTriggerName.Data()),1);
	//FIXME not sure if this is the right place for this
	//PprojXSig->Scale(1./SignalEtaRange);
  printf("   Scaling signal eta range by 1 / %f\n",RealSignalEtaRange);
	PprojXSig->Scale(1./RealSignalEtaRange);
  printf("      Scaling again for bin width 1 / %f\n",PprojXSig->GetXaxis()->GetBinWidth(3));
	PprojXSig->Scale(1./PprojXSig->GetXaxis()->GetBinWidth(3));

	PprojXSig->SetName(Form("dPhi_ObsBin%d_NearEta",CanvasPad));
	fsumCorrSE_NearEta[CanvasPad] = PprojXSig;

  PprojXSig->SetMarkerStyle(kProjNearEtaStyle);
  PprojXSig->SetMarkerColor(kProjNearEtaColor);
  PprojXSig->SetLineColor(kProjNearEtaColor);

  // FIXME trying to figure out what's going on with the eta range
	//ZoomYRange(PprojXSig,0.4);
	ZoomYRange(PprojXSig,0.4,-1.5,1.5);
	PprojXSig->DrawCopy("E");
	TLegend* leg1;
	leg1 = new TLegend(0.25,0.68,0.4,0.92); //..Bkg subtracted

  TLegend* legBoth = new TLegend(0.6,0.6,0.95,0.92);

	leg1->AddEntry(PprojXSig,"Projection in range:","pe");
	leg1->AddEntry(PprojXSig,Form("%0.2f < #Delta#eta #leq %0.2f",fitRangeL2,fitRangeR1),"");



	if(fObservable==0) {
    leg1->AddEntry(PprojXSig,Form("      %0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[CanvasPad],fTriggerName.Data(),fArray_G_Bins[CanvasPad+1]),"");
    legBoth->AddEntry(PprojXSig,Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[CanvasPad],fTriggerName.Data(),fArray_G_Bins[CanvasPad+1]),"");
  }
	if(fObservable==1) {
    leg1->AddEntry(PprojXSig,Form("       %0.1f < z_{T} < %0.1f",fArray_ZT_Bins[CanvasPad],fArray_ZT_Bins[CanvasPad+1]),"");
    legBoth->AddEntry(PprojXSig,Form("%0.1f < z_{T} < %0.1f",fArray_ZT_Bins[CanvasPad],fArray_ZT_Bins[CanvasPad+1]),"lp");
  }
	//if(fObservable==2)leg1->AddEntry(PprojXSig,Form("       %0.1f < #xi < %0.1f",fArray_XI_Bins[CanvasPad],fArray_XI_Bins[CanvasPad+1]),"");
	if(fObservable==2) {
    leg1->AddEntry(PprojXSig,Form("       %0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[CanvasPad],fArray_HPT_Bins[CanvasPad+1]),"");
    legBoth->AddEntry(PprojXSig,Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[CanvasPad],fArray_HPT_Bins[CanvasPad+1]),"");
  }
	legBoth->AddEntry(PprojXSig,"Projection in range:","");
	legBoth->AddEntry(PprojXSig,Form("%0.2f < #Delta#eta #leq %0.2f",fitRangeL2,fitRangeR1),"pl");
	legBoth->SetBorderSize(0);
	legBoth->SetFillColorAlpha(10, 0);
  legBoth->SetFillStyle(0);


	leg1->SetTextColor(kBlack);
	leg1->SetTextSize(0.07);
	leg1->SetBorderSize(0);
	leg1->SetFillColorAlpha(10, 0);
  leg1->SetFillStyle(0);
	leg1->Draw("same");

	// Temporary Label (hopefully!)
//	DrawWIP(PprojXSide1,0.45,0.5,0.31,0.13);

  fLocalPlotBoth->cd();
	SetTH1Histo(PprojXSig,"",Form("d^{2}N^{%s-h}/N^{%s}d#Delta#eta d#Delta#phi",fTriggerName.Data(),fTriggerName.Data()),0);
	PprojXSig->DrawCopy("E");


	//-------------------------------------------------------------------------------------------
	//.. Draw the eta SideBand ..
	//-------------------------------------------------------------------------------------------
	Can->cd(CanvasPad*2+2);
	PprojXSide1  ->Add(PprojXSide2);
	SetTH1Histo(PprojXSide1,"",Form("d^{2}N^{%s-h}/N^{%s}d#Delta#eta d#Delta#phi",fTriggerName.Data(),fTriggerName.Data()),1);
	// FIXME not sure if this is the right place to do this
//	PprojXSide1->Scale(1./SideBandEtaRange);
  printf("   Scaling far eta range by 1 / %f\n",SideBandEtaRange);
	PprojXSide1->Scale(1./SideBandEtaRange); 
  printf("      Scaling again for bin width 1 / %f\n",PprojXSide1->GetXaxis()->GetBinWidth(3));
	PprojXSide1->Scale(1./PprojXSide1->GetXaxis()->GetBinWidth(3));

	PprojXSide1->SetName(Form("dPhi_ObsBin%d_FarEta",CanvasPad));
	fsumCorrSE_FarEta[CanvasPad] = PprojXSide1;

	//ZoomYRange(PprojXSide1,0.4);
	ZoomYRange(PprojXSide1,0.4,-1.5,1.5);
	PprojXSide1->SetMarkerStyle(kProjFarEtaStyle);
	PprojXSide1->SetMarkerColor(kProjFarEtaColor);
	PprojXSide1->SetLineColor(kProjFarEtaColor);
	PprojXSide1->DrawCopy("E");

	TLegend* leg2;
	leg2 = new TLegend(0.25,0.68,0.4,0.92); //..Bkg subtracted
	leg2->AddEntry(PprojXSide1,"Projection in range:","pe");
  // FIXME for some reason this legend doesn't seem to get marker style/color for PprojXSide1
  // Maybe hey get changed elsewhere?
	leg2->AddEntry(PprojXSide1,Form("%0.2f < #Delta#eta #leq %0.2f",fitRangeL1,fitRangeL2),"lp");
	leg2->AddEntry(PprojXSide1,Form("+ %0.2f < #Delta#eta #leq %0.2f",fitRangeR1,fitRangeR2),"");
	leg2->SetTextColor(kBlack);
	leg2->SetTextSize(0.07);
	leg2->SetBorderSize(0);
	leg2->SetFillColorAlpha(10, 0);
  leg2->SetFillStyle(0);
	leg2->Draw("same");

	legBoth->AddEntry(PprojXSide1,Form("%0.2f < #Delta#eta #leq %0.2f",fitRangeL1,fitRangeL2),"lp");
	legBoth->AddEntry(PprojXSide1,Form("+ %0.2f < #Delta#eta #leq %0.2f",fitRangeR1,fitRangeR2),"");
  fLocalPlotBoth->cd();
	SetTH1Histo(PprojXSide1,"",Form("d^{2}N^{%s-h}/N^{%s}d#Delta#eta d#Delta#phi",fTriggerName.Data(),fTriggerName.Data()),0);
	PprojXSide1->DrawCopy("SAME");
	legBoth->Draw("same");
  DrawAlicePerf(0,0.06,0.8,0.4,0.1);
  fLocalPlotBoth->Print(Form("%s/%s_%d.pdf",fOutputDir.Data(),fLocalPlotBoth->GetName(),CanvasPad));
  fLocalPlotBoth->Print(Form("%s/%s_%d.png",fOutputDir.Data(),fLocalPlotBoth->GetName(),CanvasPad));

  fLocalPlotFarDeltaEta->cd();
	PprojXSide1->DrawCopy("");
  fLocalPlotFarDeltaEta->Print(Form("%s/%s_%d.pdf",fOutputDir.Data(),fLocalPlotFarDeltaEta->GetName(),CanvasPad));
  fLocalPlotBoth->Print(Form("%s/%s_%d.png",fOutputDir.Data(),fLocalPlotBoth->GetName(),CanvasPad));

	//-------------------------------------------------------------------------------------------
	//..Ratio of 2 sigma ranges to get a feel what a good range could be
	//-------------------------------------------------------------------------------------------
	Can2->cd(CanvasPad+1);
	PprojXSide1_T1->Add(PprojXSide2_T1);
	PprojXSide1_T2->Add(PprojXSide2_T2);
	//..Scale the sideband projections to the away-side jet hump to compare the near side region
	//	Double_t Integral1=PprojXSide1  ->Integral(PprojXSide1->FindBin(160),PprojXSide1->FindBin(200));
	//	Double_t Integral2=PprojXSide1_T->Integral(PprojXSide1_T->FindBin(160),PprojXSide1_T->FindBin(200));
	Double_t Integral1=PprojXSide1  ->Integral(PprojXSide1->FindBin(160*DTR),PprojXSide1->FindBin(200*DTR));
	Double_t IntegralT1=PprojXSide1_T1->Integral(PprojXSide1_T1->FindBin(160*DTR),PprojXSide1_T1->FindBin(200*DTR));
	Double_t IntegralT2=PprojXSide1_T2->Integral(PprojXSide1_T2->FindBin(160*DTR),PprojXSide1_T2->FindBin(200*DTR));
	SetTH1Histo(PprojXSide1,"",Form("dN^{%s-h}/N^{%s}",fTriggerName.Data(),fTriggerName.Data()),1);
	ZoomYRange(PprojXSide1,0.4);
	PprojXSide1->SetLineColor(kBlack);
	PprojXSide1->SetMarkerColor(kBlack);
	PprojXSide1->DrawCopy("E");
	PprojXSide1_T1->SetLineColor(kAzure+1);
	PprojXSide1_T1->SetMarkerColor(kAzure+1);
	PprojXSide1_T1->SetMarkerStyle(kOpenSquare);
	PprojXSide1_T1->Scale(Integral1/IntegralT1);
	PprojXSide1_T1->DrawCopy("same E");
	PprojXSide1_T2->SetLineColor(kOrange+1);
	PprojXSide1_T2->SetMarkerColor(kOrange+1);
	PprojXSide1_T2->SetMarkerStyle(kOpenSquare);
	PprojXSide1_T2->Scale(Integral1/IntegralT2);
	PprojXSide1_T2->DrawCopy("same E");

	TLegend* leg3;
	leg3 = new TLegend(0.25,0.68,0.4,0.92); //..Bkg subtracted
	leg3->AddEntry(PprojXSide1,Form("Proj. %0.1f#sigma from NS peak",Sigma2),"pe");
	leg3->AddEntry(PprojXSide1_T1,Form("Proj. %0.1f#sigma from NS peak",Sigma2*(1+variation)),"pe");
	leg3->AddEntry(PprojXSide1_T2,Form("Proj. %0.1f#sigma from NS peak",Sigma2*(1-variation)),"pe");
	if(fObservable==0)leg3->AddEntry(PprojXSig,Form("      %0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[CanvasPad],fTriggerName.Data(),fArray_G_Bins[CanvasPad+1]),"");
	if(fObservable==1)leg3->AddEntry(PprojXSig,Form("       %0.1f < z_{T} < %0.1f",fArray_ZT_Bins[CanvasPad],fArray_ZT_Bins[CanvasPad+1]),"");
	if(fObservable==2)leg3->AddEntry(PprojXSig,Form("       %0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[CanvasPad],fArray_HPT_Bins[CanvasPad+1]),"");
//	if(fObservable==0)leg3->AddEntry(PprojXSig,Form("      %0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[CanvasPad],fTriggerName.Data(),fArray_G_Bins[CanvasPad+1]),"");
//	if(fObservable==1)leg3->AddEntry(PprojXSig,Form("       %0.1f < z_{T} < %0.1f",fArray_ZT_Bins[CanvasPad],fArray_ZT_Bins[CanvasPad+1]),"");
//	if(fObservable==2)leg3->AddEntry(PprojXSig,Form("       %0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[CanvasPad],fArray_HPT_Bins[CanvasPad+1]),"");
	leg3->SetTextColor(kBlack);
	leg3->SetTextSize(0.05);
	leg3->SetBorderSize(0);
	//leg3->SetFillColorAlpha(10, 0);
  leg3->SetFillStyle(0);
	leg3->Draw("same");

	//..Ratio of 2 sigma ranges to get a feel what a good range could be
	Can3->cd(CanvasPad+1);
	SetTH1Histo(PprojXSide1_T1,"",Form("Ratio (%0.1f#sigma#pm%0.1f%%)/#bf{%0.1f#sigma}",Sigma2,variation*10,Sigma2),1);
	//PprojXSide1_T->Add(PprojXSide1,-1);
	PprojXSide1_T1->Divide(PprojXSide1);
	PprojXSide1_T1->SetLineColor(kAzure+1);
	PprojXSide1_T1->SetMarkerColor(kAzure+1);
	PprojXSide1_T1->SetMarkerStyle(kOpenSquare);
	PprojXSide1_T1->GetYaxis()->SetTitleSize(0.08);
  PprojXSide1_T1->GetYaxis()->SetTitleOffset(1.3); //FIXME
  if (bNoYLabel) PprojXSide1_T1->GetYaxis()->SetTitleOffset(0.71);
	PprojXSide1_T1->DrawCopy("E");
	PprojXSide1_T2->Divide(PprojXSide1);
	PprojXSide1_T2->DrawCopy("same E");

  TCanvas * fLocalPlot = new TCanvas("Corr1D_ProjFull_Indiv","Corr1D_ProjFull_Indiv",1750,1400);


	//-------------------------------------------------------------------------------------------
	//.. Draw as an intermediate result the full projection ..
	//-------------------------------------------------------------------------------------------
	fCanProjFull->cd(CanvasPad+1);
//	if(CanvasPad>0)fCanProjFull->cd(CanvasPad);
//	if(CanvasPad>2)fCanProjFull->cd(CanvasPad+1);
	//SetTH1Histo(PprojXFull,Form("#Delta#varphi^{%s_{can.}-h}",fTriggerName.Data()),Form("factor #bullet dN^{%s_{can.}-h}/N^{%s_{can.}}",fTriggerName.Data(),fTriggerName.Data()),1);
//	SetTH1Histo(PprojXFull,Form("#Delta#varphi^{%s_{can.}-h}",fTriggerName.Data()),Form("1/N^{%s_{can.}} dN^{%s_{can.}-h}/d#Delta#varphi^{%s_{can.}-h} (arb.units)",fTriggerName.Data(),fTriggerName.Data(),fTriggerName.Data()),1);
//	SetTH1Histo(PprojXFull,Form("#Delta#varphi^{%s_{can.}-h}",fTriggerName.Data()),Form("1/N^{%s_{can.}} dN^{%s_{can.}-h}/d#Delta#varphi (arb.units)",fTriggerName.Data(),fTriggerName.Data()),1);
	//SetTH1Histo(PprojXFull,Form("#Delta#varphi^{%s_{can.}-h}",fTriggerName.Data()),Form("#frac{1}{N^{%s_{can.}}} d^{2}N^{%s_{can.}-h}/d#Delta#eta#Delta#varphi",fTriggerName.Data(),fTriggerName.Data()),1);

	SetTH1Histo(PprojXFull,"#Delta#varphi","#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}",1);
	//SetTH1Histo(PprojXFull,Form("#Delta#varphi^{trig-assoc}",fTriggerName.Data()),"#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}",1);
	//SetTH1Histo(PprojXFull,Form("#Delta#varphi^{%s_{can.}-h}",fTriggerName.Data()),"#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}",1);

	//PprojXFull->SetLineColor(kAzure+1);
	//PprojXFull->SetMarkerColor(kAzure+1);
	Double_t range= PprojXFull->GetMaximum()-PprojXFull->GetMinimum();
//	PprojXFull->Scale(1.0/(PprojXFull->GetMaximum()+range*0.4));
  // FIXME does some scaling need to be done here for the used eta range??
	range= PprojXFull->GetMaximum()-PprojXFull->GetMinimum();
	//PprojXFull->GetYaxis()->SetRangeUser(PprojXFull->GetMinimum()-range*0.1,PprojXFull->GetMaximum()+range*0.4);
	PprojXFull->GetYaxis()->SetRangeUser(PprojXFull->GetMinimum()-range*0.1,PprojXFull->GetMaximum()+range*0.1);

	PprojXFull->SetMarkerStyle(kProjFullStyle);
  PprojXFull->SetMarkerColor(kProjFullColor);
  PprojXFull->SetLineColor(kProjFullColor);

	PprojXFull->GetXaxis()->SetLabelSize(0.045);
	PprojXFull->GetYaxis()->SetLabelSize(0.06);
	PprojXFull->GetXaxis()->SetTitleSize(0.06);
	PprojXFull->GetYaxis()->SetTitleSize(0.06);
	PprojXFull->GetYaxis()->SetTitleOffset(1.4);
  if (bNoYLabel) PprojXFull->GetYaxis()->SetTitleOffset(0.71);
	PprojXFull->GetXaxis()->SetTitleOffset(1.1);
	PprojXFull->GetXaxis()->CenterTitle(kFALSE);
	PprojXFull->DrawCopy("E");

  gPad->SetTickx();
  gPad->SetTicky();

	DrawAlicePerf(PprojXFull,0.22,0.8,0.12,0.12);
	TLegend* leg4=nullptr;
	if(fObservable==0)leg4 = new TLegend(0.48,0.78,0.83,0.84); //..Bkg subtracted
	if(fObservable!=0)leg4 = new TLegend(0.48,0.55,0.83,0.71); //..Bkg subtracted
	if(fObservable==1)leg4->AddEntry(PprojXFull,"0.15 < #it{p}_{T}^{assoc} < 30 GeV/#it{c}","");
	if(fObservable!=0)leg4->AddEntry(PprojXFull,Form("%0.1f < #it{p}_{T}^{%s} < %0.1f GeV/#it{c}",fArray_G_Bins[fPtMinBin-1],fTriggerName.Data(),fArray_G_Bins[fPtMaxBin]),"");
	//leg4->AddEntry(PprojXFull,Form("Proj. over #Delta#eta [-1.2,1.2]"),"");
	if(fObservable==0)leg4->AddEntry(PprojXFull,Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fArray_G_Bins[CanvasPad],fTriggerName.Data(),fArray_G_Bins[CanvasPad+1]),"pe");
	if(fObservable==1)leg4->AddEntry(PprojXFull,Form("%0.2f < #it{z}_{T} < %0.2f",fArray_ZT_Bins[CanvasPad],fArray_ZT_Bins[CanvasPad+1]),"pe");
	//if(fObservable==2)leg4->AddEntry(PprojXFull,Form("%0.2f < #it{#xi} < %0.2f",fArray_XI_Bins[CanvasPad],fArray_XI_Bins[CanvasPad+1]),"pe");
	if(fObservable==2)leg4->AddEntry(PprojXFull,Form("%0.1f < #it{p}_{T}^{assoc} < %0.1f GeV/#it{c}",fArray_HPT_Bins[CanvasPad],fArray_HPT_Bins[CanvasPad+1]),"pe");
	//if(fObservable!=0)leg4->AddEntry(PprojXFull,"0.15 < #it{p}_{T}^{assoc} < 30 GeV/#it{c}","");
	leg4->SetTextColor(kBlack);
	leg4->SetTextSize(0.045);
	leg4->SetBorderSize(0);
	//leg4->SetFillColorAlpha(10, 0);
  leg4->SetFillStyle(0);
	leg4->Draw("same");


  // Do the same for an individual plot
  fLocalPlot->cd();

  // Settings for individual DPhi plot
	SetTH1Histo(PprojXFull,"","",0); // Settings for indiv plot
  PprojXFull->SetMarkerSize(2.5);
	PprojXFull->DrawCopy("E");
	DrawAlicePerf(PprojXFull,0.48,0.73,0.12,0.12); //0.22,0.8,0.12,0.12
  leg4->Draw("SAME");

  fLocalPlotBoth->cd();

  fLocalPlot->Print(TString::Format("%s/%s_Plot%d.pdf",fOutputDir.Data(),fLocalPlot->GetName(),CanvasPad));
  fLocalPlot->Print(TString::Format("%s/%s_Plot%d.png",fOutputDir.Data(),fLocalPlot->GetName(),CanvasPad));
//  fLocalPlot->Print(TString::Format("%s/%s_Plot%d.eps",fOutputDir.Data(),fLocalPlot->GetName(),CanvasPad));
  fLocalPlot->Print(TString::Format("%s/%s_Plot%d.C",fOutputDir.Data(),fLocalPlot->GetName(),CanvasPad));







	//-------------------------------------------------------------------------------------------
	//.. fit with the flow function ..
	//-------------------------------------------------------------------------------------------
	/*	//.. Fit Range
	Double_t limits=0;
	limits= 90;                      //..pi
	//limits=1.25*180.0/TMath::Pi(); //..joels paper range 1 (71.6 deg too small)
	limits=1.57*180.0/TMath::Pi();   //..joels paper range 2 (90 deg)

	// arbirtary - set for you, I'm just copying and generalizing stuff here
	//double par_V10[15] = {450, 0.8e-1, 5e-2, 1e-4, 1e-2, 1e-2,1e-4, 1e-3, 1e-3, 1e-5, 1e-3, 1e-3, 1e-5, 1e-3, 1e-3};
	double par_V10[15] = {450, 0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0};
	SetJoelsParLimits(allFit,PprojXSide1,par_V10);
	BackgroundFunction->SetParNames("B", "comb v1", "comb v2", "comb v3");

	for(Int_t g = 0; g < 4; g++)
	{
		BackgroundFunction->ReleaseParameter(g);
		BackgroundFunction->SetParameter(g,10e-3);
		if(g==0)BackgroundFunction->SetParameter(g,2.5*10e6);//!!! this should change once everything is properly nrmalized
		BackgroundFunction->SetParError(g,0.0);
	}
	//..flat line
 	//BackgroundFunction->FixParameter(1,0.0);
	//BackgroundFunction->FixParameter(2,0.0);
	//BackgroundFunction->FixParameter(3,0.0);

	BackgroundFunction->SetParLimits(1,0,1);     //..v1 is usually very small and is positive or negative depending on rapidity
	BackgroundFunction->SetParLimits(2,0.0,10);  //..do not allow negative values
	BackgroundFunction->SetParLimits(3,0.0,1);   //..do not allow negative values, v3<v2 (v3 strength fraction of v2 strength)

    //..Do the actual fit
	TString funcName    = BackgroundFunction->GetName();
	TString JoelfuncName= allFit->GetName();
	//cout<<" . . . . . . . . . . . . . . "<<endl;
	//cout<<"Fit with function 1: "<<funcName<<endl;
	TFitResultPtr r  = PprojXSide1->Fit(funcName,"","",-limits,limits);//Q = quiet mode, no printout
	//..get status of the fit and if not CONVERGED than print out detailed info
	if(gMinuit->fCstatu.Contains("CONVERGED")==0)
	{
		cout<<">>> Error: Fourier fit at bin "<<CanvasPad<<" did not converge"<<endl;
		gMinuit->mnprin(2,0);
	}
	//cout<<" . . . . . . . . . . . . . . "<<endl;
	//cout<<"Fit with function 2: "<<funcName<<endl;
	//..Do a first fit to get the rough order of magnitude for
	//..the starting values set quiet mode and do not plot
	PprojXSide1->Fit(JoelfuncName,"0Q","",-limits,limits);//Q = quiet mode, no printout
	//if(gMinuit->fCstatu.Contains("CONVERGED")==0)
	//{
	//	cout<<">>> Error: Fourier fit at bin "<<CanvasPad<<" did not converge"<<endl;
	//	gMinuit->mnprin(2,0);
	//}
	for(Int_t g = 0; g < 5; g++)
	{
		Int_t oOfMag=round(log10(allFit->GetParameter(g)));
		//..set the order of magnitude as the new starting value
		allFit->SetParameter(g,pow(10,oOfMag));
	}
	//..fit again and hope that it converges
	TFitResultPtr r2 = PprojXSide1->Fit(JoelfuncName,"SQ","",-limits,limits);//Q = quiet mode, no printout
	if(gMinuit->fCstatu.Contains("CONVERGED")==0)
	{
		cout<<">>> Error: Fourier fit at bin "<<CanvasPad<<" did not converge"<<endl;
		gMinuit->mnprin(2,0);
	}

    //-------------------------------------------------------------------------------------------
	//.. Draw the fits ..
	//-------------------------------------------------------------------------------------------
	allFit->SetLineColor(kAzure-7);//kRed-9
	allFit->SetLineWidth(3);
	allFit->SetRange(-100,300);
	allFit->DrawCopy("same");
	BackgroundFunction->SetLineStyle(9);
	BackgroundFunction->SetLineColor(kYellow-3);
	BackgroundFunction->SetLineWidth(3);
	BackgroundFunction->SetRange(-100,300);
	BackgroundFunction->DrawCopy("same");

//	DrawSupComponentsJoel(allFit,3);
//	DrawSupComponentsEli(BackgroundFunction,3);
	 */

	//..Add the scaled fit also to the NS
	/*	Can->cd(CanvasPad*2+1);
	BackgroundFunction->SetParameter(0,BackgroundFunction->GetParameter(0)/fscaleFactorSBtoNS[CanvasPad]);
	BackgroundFunction->DrawCopy("same");
	allFit->SetParameter(0,allFit->GetParameter(0)/fscaleFactorSBtoNS[CanvasPad]);
	allFit->DrawCopy("same");

	//-------------------------------------------------------------------------------------------
	//.. Return the background subtracted correlation spectrum ..
	//-------------------------------------------------------------------------------------------
	ProjectionName = PprojXSide1->GetName();
	ProjectionName+="_backgroundSubtracted";
	BackgroundFunction->SetLineColor(kWhite);
	TH1D* BackgroundSubtraction = (TH1D*)PprojXSide1->Clone(ProjectionName);
	BackgroundSubtraction->Add(BackgroundFunction,-1);
//	BackgroundSubtraction->GetFunction(funcName)->Delete();//..to not draw the fit function at any point later
	 */
	//  return BackgroundSubtraction;
}
///
/// Plot different components of v_n fit
//  for Joel's fit function
//________________________________________________________________________
void PlotGHcorrelation2::DrawSupComponentsJoel(TF1 *func, Int_t harmonic)
{
	//..lookup of harmonic vs. parameters to loop over
	Int_t loop=0;
	if(harmonic==3)loop=5;
	if(harmonic==4)loop=7;
	if(harmonic==5)loop=8;
	if(harmonic==7)loop=11;
	if(harmonic==6 || harmonic>7)cout<<"Error in PlotGHcorrelation2::DrawSupComponentsJoel"<<endl;

	TF1 *allFitSingle[11];

	for(Int_t i=0;i<11;i++)
	{
		//allFitSingle[i]= allFitFuncVn(Form("JoelsFlowFunction%d",i),harmonic,-100,300);
		allFitSingle[i]= allFitFuncVn(Form("JoelsFlowFunction%d",i),harmonic,-100*DTR,300*DTR);
		allFitSingle[i]->SetParName(i,func->GetParName(i));
	}

	//.. present single components to the total spectrum?
	for(Int_t j=0;j<loop;j++)
	{
		//cout<<"Orig. Parameter value: p("<<j<<") ="<<func->GetParameter(j)<<endl;
	}

	for(Int_t component=0;component<loop;component++)
	{
		for(Int_t j=0;j<loop;j++)
		{
			allFitSingle[component]->SetParameter(j,0);
		}
		allFitSingle[component]->SetParameter(component,func->GetParameter(component));
		if(component==1)
		{
			allFitSingle[component]->SetParameter(1,func->GetParameter(2));
			//..due to the fact that param 1 and param 4 are proportional
			//..we do a little hack here. Setting param4 to 0 is neseccary
			//..to see only the small wiede gaussian. If we set param4 to 0
			//..however, param1 will become 0 by defualt. We can however
			//..set param4 to a negligibly small value x and multiply param1
			//..by the inverse of x. (normally param3 is in the range 0-1, but we omit this for this specific case)
			Double_t Shrinkage=0.00000000000001;
			//..
			allFitSingle[component]->SetParameter(4,Shrinkage);//..that is v3
			allFitSingle[component]->SetParameter(1,1.0*func->GetParameter(1)*func->GetParameter(4)/Shrinkage);
		}
		if(component==2 || component==3)
		{
			allFitSingle[component]->SetParameter(2,func->GetParameter(2));
			allFitSingle[component]->SetParameter(3,func->GetParameter(3));
		}
		if(component==5 || component==6)
		{
			allFitSingle[component]->SetParameter(5,func->GetParameter(5));
			allFitSingle[component]->SetParameter(6,func->GetParameter(6));
		}
		allFitSingle[component]->SetParameter(0,func->GetParameter(0));
		if(component<6)allFitSingle[component]->SetLineColor(fColorSceme[component]);//..not defined for >5
		allFitSingle[component]->SetLineStyle(1);
		for(Int_t j=0;j<loop;j++)
		{
			//cout<<"Parameter value: p("<<j<<") ="<<allFitSingle[component]->GetParameter(j)<<endl;
		}
		if(component!=3 || component!=6)allFitSingle[component]->Draw("same");
		//cout<<"Drawn: "<<allFitSingle[component]->GetName()<<"(Round: "<<component<<")"<<endl;
		//cout<<"func value at 0= "<<allFitSingle[component]->Eval(0)<<", value at pi"<<allFitSingle[component]->Eval(90)<<endl;
	}

	TLegend *legSubC = new TLegend(0.80,0.65,0.85,0.90);
	//legSubC->SetHeader("ALICE performance (27.1.2017)");
	for(Int_t component=0;component<loop;component++)
	{
		if(component!=3 || component!=6)legSubC->AddEntry(allFitSingle[component],allFitSingle[component]->GetParName(component),"pl");
	}
	legSubC->SetTextSize(0.05);
	legSubC->SetBorderSize(0);
	//legSubC->SetFillColorAlpha(10, 0);
  legSubC->SetFillStyle(0);
	legSubC->Draw("same");

}
///
/// Plot different components of v_n fit
//  for Eliane's fit function
//________________________________________________________________________
void PlotGHcorrelation2::DrawSupComponentsEli(TF1 *func, Int_t harmonic)
{
	TF1 *allFitSingle[11];
	for(Int_t i=0;i<11;i++)
	{
		//allFitSingle[i]= new TF1(Form("ElianesFlowFunction%d",i),FlowFunction,-100,300,4);
		allFitSingle[i]= new TF1(Form("ElianesFlowFunction%d",i),FlowFunction,-100*DTR,300*DTR,4);
		allFitSingle[i]->SetParName(i,func->GetParName(i));
	}

	for(Int_t component=0;component<harmonic+1;component++)
	{
		for(Int_t j=0;j<(harmonic+1);j++)
		{
			allFitSingle[component]->SetParameter(j,0);
		}
		allFitSingle[component]->SetParameter(component,func->GetParameter(component));
		allFitSingle[component]->SetParameter(0,func->GetParameter(0));
		if(component==1)
		{
			//..due to the fact that param 1 and param 4 are proportional
			//..we do a little hack here. Setting param4 to 0 is neseccary
			//..to see only the small wiede gaussian. If we set param4 to 0
			//..however, param1 will become 0 by defualt. We can however
			//..set param4 to a negligibly small value x and multiply param1
			//..by the inverse of x. (normally param3 is in the range 0-1, but we omit this for this specific case)
			Double_t Shrinkage=0.00000000000001;
			//..
			allFitSingle[component]->SetParameter(2,Shrinkage);//..that is v2
			allFitSingle[component]->SetParameter(3,Shrinkage);//..that is v3
			allFitSingle[component]->SetParameter(1,1.0*func->GetParameter(1)*func->GetParameter(2)*func->GetParameter(3)/(Shrinkage*Shrinkage));
		}
		if(component==3)
		{
			//..due to the fact that param 1 and param 4 are proportional
			//..we do a little hack here. Setting param4 to 0 is neseccary
			//..to see only the small wiede gaussian. If we set param4 to 0
			//..however, param1 will become 0 by defualt. We can however
			//..set param4 to a negligibly small value x and multiply param1
			//..by the inverse of x. (normally param3 is in the range 0-1, but we omit this for this specific case)
			Double_t Shrinkage=0.00000000000001;
			//..
			allFitSingle[component]->SetParameter(2,Shrinkage);//..that is v2
			allFitSingle[component]->SetParameter(3,1.0*func->GetParameter(2)*func->GetParameter(3)/(Shrinkage));
		}
		if(component<6)allFitSingle[component]->SetLineColor(fColorSceme[component]);//..not defined for >5
		allFitSingle[component]->SetLineStyle(9);
		for(Int_t j=0;j<(harmonic+1);j++)
		{
			//cout<<"Parameter value: p("<<j<<") ="<<allFitSingle[component]->GetParameter(j)<<endl;
		}
		allFitSingle[component]->Draw("same");
		//cout<<"Drawn: "<<allFitSingle[component]->GetName()<<"(Round: "<<component<<")"<<endl;
		//cout<<"func value at 0= "<<allFitSingle[component]->Eval(0)<<", value at pi"<<allFitSingle[component]->Eval(90)<<endl;
	}
}
///**
///Draw a legend box with standard info
///**
//________________________________________________________________________
void PlotGHcorrelation2::DrawAliceInfoBox(TObject* Histo)
{
	//	\newcommand{\PbPb}{\ensuremath{\mbox{Pb--Pb}} }
	//	\newcommand{\sqrtSnn}{\ensuremath{\sqrt{s_{\mathrm{NN}}}}}

	TLegend *leg = new TLegend(0.60,0.60,0.9,0.85);
	leg->SetHeader("ALICE performance (27.1.2017)");
	leg->AddEntry(Histo,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pe");
	leg->AddEntry(Histo,"#it{E}_{T}^{ cluster}> 10 GeV/#it{c}","");
	leg->AddEntry(Histo,"#it{p}_{T}^{ h}      > 0.15 GeV/#it{c}","");
	leg->SetTextSize(0.05);
	leg->SetBorderSize(0);
	//leg->SetFillColorAlpha(10, 0);
  leg->SetFillStyle(0);

	leg->Draw("same");
}
///
/// Draw an ALICE performance legend entry
//
//________________________________________________________________________
void PlotGHcorrelation2::DrawAlicePerf(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size)
{
  const char *kMonthList[12] = {"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec."};
  const char *kCentList[5] = {"0-80%","0-10%","10-30%","30-50%","50-80%"}; // index=fCent+1
  TLegend * leg  = new TLegend(x,y,x+x_size,y+y_size);
  TDatime * time = new TDatime();
  const char * month = kMonthList[time->GetMonth()-1];

  //leg->SetHeader(Form("ALICE Performance - %d %s %d",time->GetDay(),month,time->GetYear()));
  //if (fPerformance) leg->AddEntry(Histo,Form("ALICE Performance %d %s %d",time->GetDay()-1,month,time->GetYear()),"");  
  //else leg->AddEntry(Histo,Form("Work in Progress %d %s %d",time->GetDay()-1,month,time->GetYear()),""); 
  if (fPerformance) leg->AddEntry(Histo,"ALICE Performance","");  
  else leg->AddEntry(Histo,Form("Work in Progress %d %s %d",time->GetDay(),month,time->GetYear()),""); 
 // leg->AddEntry(Histo,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 0-90%","");
  leg->AddEntry(Histo,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, %s",kCentList[fCent+1]),"");
  
  // Info about event plane


  leg->SetTextSize(0.041); // 0.045
  leg->SetBorderSize(0);
  //leg->SetFillColorAlpha(10,0);
  leg->SetFillStyle(0);

  leg->Draw("SAME");
}


// Copied Flow Analysis from Phase 1
void PlotGHcorrelation2::ProduceDeltaPsiPlots() {
//  if (hPtEPAnglePionAcc) hPtEPAnglePionAcc->Rebin2D(nRebinDeltaPsi,1);
  Int_t kUsedPi0TriggerPtBins = 5; // How many Pt bins do we actually use

  std::vector<double> Pi0PtBins = {5,7,9,11,14,17};

  if (hPtEPAnglePionAcc==0) {
    fprintf(stderr,"Missing hPtEPAnglePionAcc\n");
    return;
  }

  printf("about to make some event plane projections\n");


  // In MC, so far, we only have hPtEPAnglePionAcc

  if (fPtEPAnglePionAccCent == 0) {
    printf("Missing fPtEPAnglePionAcc\n");
  } else {
    // Set the Centrality range
    fPtEPAnglePionAccCent->GetZaxis()->SetRange(fCent+1,fCent+1);
  }

  if (fPtEP3AnglePionAccCent) {
    fPtEP3AnglePionAccCent->GetZaxis()->SetRange(fCent+1,fCent+1);
  } else {
    printf("Missing hPtEP3AnglePionAcc\n");
  }
  if (fPtEP4AnglePionAccCent) {
    fPtEP4AnglePionAccCent->GetZaxis()->SetRange(fCent+1,fCent+1);
  } else {
    printf("Missing hPtEP4AnglePionAcc\n");
  }

  // For each of the above, project onto DeltaPsi 6 pt bins
  for (int i = 0; i < kUsedPi0TriggerPtBins; i++) {
    double fMinPt = Pi0PtBins[i];
    double fMaxPt = Pi0PtBins[i+1];

    // Note that the yaxis has bins of 500 MeV.
    int iMinBin = hPtEPAnglePionAcc->GetYaxis()->FindBin(fMinPt);
    //int iMaxBin = hPtEPAnglePionAcc->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound
    int iMaxBin = hPtEPAnglePionAcc->GetYaxis()->FindBin(fMaxPt-0.001); // Want the bin with fMaxPt as an upper bound
    
    printf("Projecting Event Plane Histograms in bin from %.1f to %.1f\n",hPtEPAnglePionAcc->GetYaxis()->GetBinLowEdge(iMinBin),hPtEPAnglePionAcc->GetYaxis()->GetBinUpEdge(iMaxBin));

    TString sFormat = "%s_Proj_%d";
    TString sPtRange = Form("%.0f #leq #it{p}_{T} < %.0f GeV/#it{c}",fMinPt,fMaxPt);

    // Accepted Pi0s vs event plane (all centralities in the wagon)

    TH1F * hLocalPtEPAnglePionAcc_Proj = (TH1F *) hPtEPAnglePionAcc->ProjectionX(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i),iMinBin,iMaxBin);
    hLocalPtEPAnglePionAcc_Proj->Sumw2();
    hLocalPtEPAnglePionAcc_Proj->SetTitle(Form("#pi_{0}^{Cand} #Delta#Psi_{EP} (%s)",sPtRange.Data()));
    hLocalPtEPAnglePionAcc_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
    hPtEPAnglePionAcc_Proj.push_back(hLocalPtEPAnglePionAcc_Proj);

    if (fPtEPAnglePionAccCent == 0) continue;


    int iMinBinCent = fPtEPAnglePionAccCent->GetYaxis()->FindBin(fMinPt);
    int iMaxBinCent = fPtEPAnglePionAccCent->GetYaxis()->FindBin(fMaxPt-0.001);

    printf("Projecting the Trigger Vs EP (centrality selection) range %d - %d\n",iMinBinCent,iMaxBinCent);

    // Just the centrality selected by the task
    TH1F * fLocalPtEPAnglePionAccCent_Proj = (TH1F *) fPtEPAnglePionAccCent->ProjectionX(Form(sFormat.Data(),fPtEPAnglePionAccCent->GetName(),i),iMinBinCent, iMaxBinCent, fCent+1, fCent+1);
    fLocalPtEPAnglePionAccCent_Proj->Sumw2();
    fLocalPtEPAnglePionAccCent_Proj->SetTitle(Form("#pi_{0}^{Cand} #Delta#Psi_{EP} (%s)",sPtRange.Data()));
    fLocalPtEPAnglePionAccCent_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
    fPtEPAnglePionAccCent_Proj.push_back(fLocalPtEPAnglePionAccCent_Proj);

    if (fPtEP3AnglePionAccCent) {
      TH1F * fLocalPtEP3AnglePionAccCent_Proj = (TH1F *) fPtEP3AnglePionAccCent->ProjectionX(Form(sFormat.Data(),fPtEP3AnglePionAccCent->GetName(),i),iMinBinCent, iMaxBinCent, fCent+1, fCent+1);
      fLocalPtEP3AnglePionAccCent_Proj->Sumw2();
      fLocalPtEP3AnglePionAccCent_Proj->SetTitle(Form("#pi_{0}^{Cand} #Delta#Psi_{EP,3} (%s)",sPtRange.Data()));
      fLocalPtEP3AnglePionAccCent_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
      fPtEP3AnglePionAccCent_Proj.push_back(fLocalPtEP3AnglePionAccCent_Proj);
    }
    if (fPtEP3AnglePionAccCent) {
      TH1F * fLocalPtEP4AnglePionAccCent_Proj = (TH1F *) fPtEP4AnglePionAccCent->ProjectionX(Form(sFormat.Data(),fPtEP4AnglePionAccCent->GetName(),i),iMinBinCent, iMaxBinCent, fCent+1, fCent+1);
      fLocalPtEP4AnglePionAccCent_Proj->Sumw2();
      fLocalPtEP4AnglePionAccCent_Proj->SetTitle(Form("#pi_{0}^{Cand} #Delta#Psi_{EP,4} (%s)",sPtRange.Data()));
      fLocalPtEP4AnglePionAccCent_Proj->GetYaxis()->SetTitle("N_{#pi_{0}^{Cand}}");
      fPtEP4AnglePionAccCent_Proj.push_back(fLocalPtEP4AnglePionAccCent_Proj);
    }
  }

  printf("Finished making projections of event plane trigger histograms\n");

}

void PlotGHcorrelation2::MeasureVn() {


}



///
/// More advanced QA of corrected DPhi correlations: look for residual TPC sector effect
//
//________________________________________________________________________
void PlotGHcorrelation2::DPhiQA() {
	cout<<"o Doing DeltaPhi QA with FFTs: "<<endl;
//  TH1D * fsumCorrSE_ProjFull[10];  

  for (Int_t i=0;i<fmaxBins;i++)
  {
    printf("step %d\n",i);
    TH1D * hFFTTemp = 0;
    //fFFTsumCorrSE_ProjFull[i] = (TH1D *) fsumCorrSE_ProjFull[i]->FFT(fFFTsumCorrSE_ProjFull[i],"MAG");
    hFFTTemp = (TH1D *) fsumCorrSE_ProjFull[i]->FFT(hFFTTemp,"MAG");
    hFFTTemp->SetName(Form("QA_FFT_ProjFull%d",i));
    fFFTsumCorrSE_ProjFull[i] = hFFTTemp;
  }
}



/*
///
/// Plot the 2D correlation functions
//
//________________________________________________________________________
void PlotGHcorrelation2::PlotProjectedCorr(Int_t version,TH2D* corrHistoSE[],TH2D* corrHistoSEClone[])
{
	Double_t SigmaEtaSigmaPhi[4];
	TCanvas *corrPlotsWidth[version] = new TCanvas("Corr1D_GPlots","Corrected_G 1D Plots 1) Full",...,1400,1400);
	TCanvas *CanvProj = new TCanvas("Corr1D_GPlots_SB","Corrected_G 1D Plots 2) SB",1400,1400);
	TCanvas *CanvProj2 = new TCanvas("Corr1D_GPlots_SB2","Corrected_G 1D Plots 2) SB",1400,1400);

	Int_t maxCanvasPad=0;
	TString histoName=rawHistoSE[0]->GetName();
	//kan ich auch mit version machen
	if(histoName.Contains("DEtaDPhiG")) {maxCanvasPad = kGammaNBINS;}
	if(histoName.Contains("DEtaDPhiZT")){maxCanvasPad = kZtNBINS;	}
	if(histoName.Contains("DEtaDPhiXI")){maxCanvasPad = kXiNBINS;	}

	corrPlotsWidth[version] ->Divide(4,4,0.001,0.001);
	corrPlotsNSfit[version] ->Divide(4,4);
	corrPlotsNSfit2[version]->Divide(4,4);

	for(Int_t i=0;i<maxCanvasPad;i++)
	{
		//..perform a 2D gaussian fit to determine the near side width in eta and phi
		//..this is useful to determine the range of eta not influenced by the jet signal
		DetermineWidths(corrHistoSE[i],SigmaEtaSigmaPhi,Canv,i);
		//..plot large delta eta for checking out the background (determine flow coefficients for g_dir + g_decay mixture)
		//..get the width from the previous function

		FitEtaSides(corrHistoSEClone[i],SigmaEtaSigmaPhi,2.0,corrPlotsNSfit2[version],i);

		FitEtaSides(corrHistoSE[i],SigmaEtaSigmaPhi,3.5,corrPlotsNSfit[version],i);

		if(version==0)
		{
			fEtaWidth_gamma->SetBinContent(i+1,SigmaEtaSigmaPhi[0]);
			fEtaWidth_gamma->SetBinError(i+1,SigmaEtaSigmaPhi[1]);
			fPhiWidth_gamma->SetBinContent(i+1,SigmaEtaSigmaPhi[2]);
			fPhiWidth_gamma->SetBinError(i+1,SigmaEtaSigmaPhi[3]);
		}
		if(version==1)
		{
			fEtaWidth_Zt->SetBinContent(i+1,SigmaEtaSigmaPhi[0]);
			fEtaWidth_Zt->SetBinError(i+1,SigmaEtaSigmaPhi[1]);
			fPhiWidth_Zt->SetBinContent(i+1,SigmaEtaSigmaPhi[2]);
			fPhiWidth_Zt->SetBinError(i+1,SigmaEtaSigmaPhi[3]);
		}
		if(version==2)
		{
			fEtaWidth_Xi->SetBinContent(i+1,SigmaEtaSigmaPhi[0]);
			fEtaWidth_Xi->SetBinError(i+1,SigmaEtaSigmaPhi[1]);
			fPhiWidth_Xi->SetBinContent(i+1,SigmaEtaSigmaPhi[2]);
			fPhiWidth_Xi->SetBinError(i+1,SigmaEtaSigmaPhi[3]);
		}
	}
    //..Draw a canvas with the eta and phi widths
	TCanvas *Canv_Width = new TCanvas("Canv_Width","Width of NS peak vs E_g",1400,700);
	Canv_Width->Divide(2);
	Canv_Width->cd(1);
	SetTH1Histo(EtaWidth_gamma,"E_{#gamma}","#sigma of #Delta#eta");
	EtaWidth_gamma->SetMarkerColor(17);
	EtaWidth_gamma->SetMarkerStyle(20);
	EtaWidth_gamma->SetMarkerSize(1.1);
	EtaWidth_gamma->DrawCopy("E");
	EtaWidth_gamma->SetLineColor(colorSceme[0]);
	EtaWidth_gamma->SetMarkerColor(colorSceme[0]);
	EtaWidth_gamma->SetMarkerSize(0.8);
	EtaWidth_gamma->DrawCopy("same E");

	Canv_Width->cd(2);
	SetTH1Histo(PhiWidth_gamma,"E_{#gamma}","#sigma of #Delta#phi");
	PhiWidth_gamma->SetMarkerColor(17);
	PhiWidth_gamma->SetMarkerStyle(20);
	PhiWidth_gamma->SetMarkerSize(1.1);
	PhiWidth_gamma->DrawCopy("E");
	PhiWidth_gamma->SetLineColor(colorSceme[0]);
	PhiWidth_gamma->SetMarkerColor(colorSceme[0]);
	PhiWidth_gamma->SetMarkerSize(0.8);
	PhiWidth_gamma->DrawCopy("same E");


	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 	cout<<"Draw Projections and fit near side with for ZT"<<endl;
	TCanvas *Canv2 = new TCanvas("Corr1D_ZTPlots","Corrected_ZT 1D Plots 1)Full",1400,1400);
	Canv2->Divide(4,4,0.001,0.001);
	TCanvas *Canv2Proj = new TCanvas("Corr1D_ZTPlots_SB","Corrected_ZT 1D Plots 2) SB",1400,1400);
	Canv2Proj->Divide(4,4);
	for(Int_t i=0;i<7;i++)
	{
		DetermineWidths(dEdP_ZT[i],SigmaEtaSigmaPhi,Canv2,i);
		dEdP_ZT_BgSub[i] = FitEtaSides(dEdP_ZT[i],SigmaEtaSigmaPhi,3.5,Canv2Proj,i);

		EtaWidth_Zt->SetBinContent(i+1,SigmaEtaSigmaPhi[0]);
		EtaWidth_Zt->SetBinError(i+1,SigmaEtaSigmaPhi[1]);
		PhiWidth_Zt->SetBinContent(i+1,SigmaEtaSigmaPhi[2]);
		PhiWidth_Zt->SetBinError(i+1,SigmaEtaSigmaPhi[3]);
	}

	//..Draw a canvas with the eta and phi widths
	TCanvas *Canv_Width2 = new TCanvas("Canv_Width2","Width of NS peak vs Zt",1400,700);
	Canv_Width2->Divide(2);
	Canv_Width2->cd(1);
	SetTH1Histo(EtaWidth_Zt,"z_{T}","#sigma of #Delta#eta");
	EtaWidth_Zt->GetXaxis()->SetRangeUser(0,1);
	EtaWidth_Zt->SetMarkerColor(17);
	EtaWidth_Zt->SetMarkerStyle(20);
	EtaWidth_Zt->SetMarkerSize(1.1);
	EtaWidth_Zt->DrawCopy("E");
	EtaWidth_Zt->SetLineColor(colorSceme[0]);
	EtaWidth_Zt->SetMarkerColor(colorSceme[0]);
	EtaWidth_Zt->SetMarkerSize(0.8);
	EtaWidth_Zt->DrawCopy("same E");

	Canv_Width2->cd(2);
	SetTH1Histo(PhiWidth_Zt,"z_{T}","#sigma of #Delta#phi");
	PhiWidth_Zt->GetXaxis()->SetRangeUser(0,1);
	PhiWidth_Zt->SetMarkerColor(17);
	PhiWidth_Zt->SetMarkerStyle(20);
	PhiWidth_Zt->SetMarkerSize(1.1);
	PhiWidth_Zt->DrawCopy("E");
	PhiWidth_Zt->SetLineColor(colorSceme[0]);
	PhiWidth_Zt->SetMarkerColor(colorSceme[0]);
	PhiWidth_Zt->SetMarkerSize(0.8);
	PhiWidth_Zt->DrawCopy("same E");

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	cout<<"Draw Projections and fit near side with for Xi"<<endl;
	TCanvas *Canv3 = new TCanvas("Corr1D_XIPlots","Corrected_XI 1D Plots 1) Full",1400,1400);
	Canv3->Divide(4,4,0.001,0.001);
	TCanvas *Canv3Proj = new TCanvas("Corr1D_XIPlots_SB","Corrected_XI 1D Plots 2) SB",1400,1400);
	Canv3Proj->Divide(4,4);
	for(Int_t i=0;i<8;i++)
	{
		DetermineWidths(dEdP_XI[i],SigmaEtaSigmaPhi,Canv3,i);
		dEdP_XI_BgSub[i] = FitEtaSides(dEdP_XI[i],SigmaEtaSigmaPhi,3.5,Canv3Proj,i);

		EtaWidth_Xi->SetBinContent(i+1,SigmaEtaSigmaPhi[0]);
		EtaWidth_Xi->SetBinError(i+1,SigmaEtaSigmaPhi[1]);
		PhiWidth_Xi->SetBinContent(i+1,SigmaEtaSigmaPhi[2]);
		PhiWidth_Xi->SetBinError(i+1,SigmaEtaSigmaPhi[3]);
	}

	//..Draw a canvas with the eta and phi widths
	TCanvas *Canv_Width3 = new TCanvas("Canv_Width3","Width of NS peak vs Xi",1400,700);
	Canv_Width3->Divide(2);
	Canv_Width3->cd(1);
	SetTH1Histo(EtaWidth_Xi,"#xi","#sigma of #Delta#eta");
	EtaWidth_Xi->GetXaxis()->SetRangeUser(0,2);
	EtaWidth_Xi->SetMarkerColor(17);
	EtaWidth_Xi->SetMarkerStyle(20);
	EtaWidth_Xi->SetMarkerSize(1.1);
	EtaWidth_Xi->DrawCopy("E");
	EtaWidth_Xi->SetLineColor(colorSceme[0]);
	EtaWidth_Xi->SetMarkerColor(colorSceme[0]);
	EtaWidth_Xi->SetMarkerSize(0.8);
	EtaWidth_Xi->DrawCopy("same E");

	Canv_Width3->cd(2);
	SetTH1Histo(PhiWidth_Xi,"#xi","#sigma of #Delta#phi");
	PhiWidth_Xi->GetXaxis()->SetRangeUser(0,2);
	PhiWidth_Xi->SetMarkerColor(17);
	PhiWidth_Xi->SetMarkerStyle(20);
	PhiWidth_Xi->SetMarkerSize(1.1);
	PhiWidth_Xi->DrawCopy("E");
	PhiWidth_Xi->SetLineColor(colorSceme[0]);
	PhiWidth_Xi->SetMarkerColor(colorSceme[0]);
	PhiWidth_Xi->SetMarkerSize(0.8);
	PhiWidth_Xi->DrawCopy("same E");

}
 */
/*  |
   |
   |
  \  /
   \/
 to be updated

 */

