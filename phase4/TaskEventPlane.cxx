#if defined(__CINT__)
#define _SYS_TYPES_H_
#endif
#define DTR TMath::DegToRad()
// --- ROOT system ---
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TMultiGraph.h>
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

// C,C++ Stuff
#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>

#include "TaskEventPlane.h"
#include "TaskEventPlaneMathTools.cxx"
#include "TaskEventPlaneGraphicsTools.cxx"

using namespace std;

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;
/// \cond CLASSIMP
ClassImp(TaskEventPlane);

static const Bool_t bUnifiedNorm = true; // True if each input is normalized to the total triggers, not just the triggers within that EP bin



void SetGraphColorStyle(TGraph * graph, Int_t iColor, Int_t iMarkerStyle, Int_t iMarkerSize = -1){
  graph->SetLineColor(iColor);
  graph->SetMarkerColor(iColor);
  graph->SetMarkerStyle(iMarkerStyle);
  if (iMarkerSize != -1) graph->SetMarkerSize(iMarkerSize);
}

void SetTH1ColorStyle(TH1 * hist, Int_t iColor, Int_t iMarkerStyle){
  hist->SetLineColor(iColor);
  hist->SetMarkerColor(iColor);
  hist->SetMarkerStyle(iMarkerStyle);
}


///
/// Default constructor
///
//________________________________________________________________________

TaskEventPlane::TaskEventPlane():TObject(),fIsMCGenMode(0),fOutputFile(0),iRPFMode(0),fObservable(-1),fObservableName(),nObsBins(0),fObsBins(),
fAllTriggerPt(0),fEPBinTriggerPt(0),
fFullDPhiProjAll(0),fFullDPhiProj(0),
fNearEtaDPhiProjAll(0),fNearEtaDPhiProj(0),
fFarEtaDPhiProjAll(0),fFarEtaDPhiProj(0),
fFullDPhiProjAll_Sub(0),fFullDPhiProj_Sub(0),
fNearEtaDPhiProjAll_Sub(0),fNearEtaDPhiProj_Sub(0),
fFarEtaDPhiProjAll_Sub(0),fFarEtaDPhiProj_Sub(0),
fFullDPhiProj_Rescale(0),
fNearEtaDPhiProj_Rescale(0),
fFarEtaDPhiProj_Rescale(0)
{

  gROOT->SetBatch(kTRUE);
	
	fDebugLevel = 0;

	fInputFileAll = 0;
	for (Int_t i = 0; i < 3; i++) fInputFileEvt[i] = 0;
	fPlotOptions = "LEGO2";
	fOutputDir = "output";
	fSavePlots = 1;

	fPlaneLabels.push_back("In Plane");
	fPlaneLabels.push_back("Mid Plane");
	fPlaneLabels.push_back("Out of Plane");
	fPlaneLabels.push_back("All Angles");

	SetStyle();
}

void TaskEventPlane::SetStyle() {
//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	//gStyle->SetPalette(53);  //standard is 1
	gStyle->SetCanvasColor(10);
	//  TGaxis::SetMaxDigits(4);  //..ELI I don't remember why I wanted 4 changed to 2
	gStyle->SetPadTopMargin(0.07);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	//  gStyle->SetPadRightMargin(0.045);
	gStyle->SetPadRightMargin(0.08);
	gStyle->SetPadLeftMargin(0.21);
	gStyle->SetFrameFillColor(10);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(0.07,"X");
//	gStyle->SetTitleSize(5.0,"Y");
	TGaxis::SetMaxDigits(3); //2
	gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas

  //wierd fun settings
  /*
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameFillColor(kGray);
  gStyle->SetAxisColor(0,"XYZ");
  gStyle->SetCanvasColor(kGray);
  gStyle->SetLabelColor(0,"XYZ");
  gStyle->SetStatColor(0);
  gStyle->SetStatTextColor(0);
  gStyle->SetLegendFillColor(kGray);
  gStyle->SetTitleColor(kWhite,"XYZ");
  gStyle->SetTitleColor(kWhite,"T");
  gStyle->SetTitleFillColor(kGray);
  gStyle->SetTitleTextColor(kWhite);
  gStyle->SetTextColor(kWhite);
*/

}

/**
 * Sets the In-plane, mid-plane, out-of-plane titles 
 */
void TaskEventPlane::SetEPTitles(vector<vector<TH1D *>> HistArray) {
  for (UInt_t i = 0; i < HistArray.size(); i++) {
    vector<TH1D *> localArray = HistArray[i];
    for (UInt_t j = 0; j < localArray.size(); j++) {
      localArray[j]->SetTitle(fEPBinTitles[j]);
      localArray[j]->SetTitleSize(0.06);
    } 
  }
}

void TaskEventPlane::LoadHistograms() {
	cout<<"Loading Histograms"<<endl;

  cout<<"Loading ALICE Published Flow Measurements"<<endl;
  TFile * fFlowMeasurementsFile = new TFile(sFlowGraphPath.Data(),"OPEN");
  if (!fFlowMeasurementsFile) {

    fprintf(stderr,"Error: could not find ALICE Published flow measurements files.\n");
  } else {
    fV2Graph2D = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v2Graph2D");
    fV3Graph2D = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v3Graph2D");
    fV4Graph2D = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v4Graph2D");
    fV2Graph2DErrUp = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v2Graph2DErrUp");
    fV3Graph2DErrUp = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v3Graph2DErrUp");
    fV4Graph2DErrUp = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v4Graph2DErrUp");
    fV2Graph2DErrDown = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v2Graph2DErrDown");
    fV3Graph2DErrDown = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v3Graph2DErrDown");
    fV4Graph2DErrDown = (TGraph2DErrors *) fFlowMeasurementsFile->Get("v4Graph2DErrDown");
  }
  if (!fV2Graph2D || !fV3Graph2D || !fV4Graph2D) {
    fprintf(stderr,"Error: Missing one or more 2D vN graphs\n");
  }
  if (!fV2Graph2DErrUp || !fV3Graph2DErrUp || !fV4Graph2DErrUp) {
    fprintf(stderr,"Error: Missing one or more 2D vN Error Up range graphs\n");
  }
  if (!fV2Graph2DErrDown || !fV3Graph2DErrDown || !fV4Graph2DErrDown) {
    fprintf(stderr,"Error: Missing one or more 2D vN Error Down range graphs\n");
  }



	// Loading our observable settings:
	TH1D * VariableInfo = (TH1D *) fInputFileAll->Get("VariableInfo");
	if (VariableInfo) {
		fObservable = VariableInfo->GetBinContent(1);
	} else {
		cout<<"No Variable Input TH1D found in Input File All."<<endl;
		cout<<"Using default values for Observable Info and name."<<endl;
		//fObservable=1; // the Z_t
		fObservable=2; // the pTa
	}
	if (fObservable==0)      fObservableName = "p_{T} (GeV/#it{c}";
	else if (fObservable==1) fObservableName = "z_{T}";
	else if (fObservable==2) fObservableName = "p_{T}^{a} (GeV/#it{c})";

  hHistTrackPsiEPPtCent = (TH3F *) fInputFileAll->Get("fHistTrackPsiEPPtCent");
  if (!hHistTrackPsiEPPtCent) {
    fprintf(stderr,"Missing Th3F fHistTrackPsiEPPtCent\n");
  }
  hHistTrackPsiEP3PtCent = (TH3F *) fInputFileAll->Get("fHistTrackPsiEP3PtCent");
  if (!hHistTrackPsiEP3PtCent) {
    fprintf(stderr,"Missing Th3F fHistTrackPsiEP3PtCent\n");
  }
  hHistTrackPsiEP4PtCent = (TH3F *) fInputFileAll->Get("fHistTrackPsiEP4PtCent");
  if (!hHistTrackPsiEP4PtCent) {
    fprintf(stderr,"Missing Th3F fHistTrackPsiEP4PtCent\n");
  }

  // Load the phase 1 v_n information (soon to be obsolute

  if (fIsMCGenMode) {
    gTrigger_Bv = (TGraphErrors *) fInputFileAll->Get("Trigger_Bv");
    gTrigger_V2 = (TGraphErrors *) fInputFileAll->Get("Trigger_V2");
    gTrigger_V4 = (TGraphErrors *) fInputFileAll->Get("Trigger_V4");
    gTrigger_V6 = (TGraphErrors *) fInputFileAll->Get("Trigger_V6");
    gTrack_Bv = (TGraphErrors *) fInputFileAll->Get("Track_Bv");
    gTrack_V2 = (TGraphErrors *) fInputFileAll->Get("Track_V2");
    gTrack_V4 = (TGraphErrors *) fInputFileAll->Get("Track_V4");
    gTrack_V6 = (TGraphErrors *) fInputFileAll->Get("Track_V6");
  } else { // Get them from data file
    // The data file is either presubtracted, or eventually it will be sideband subtracted
    gTrigger_V2 = (TGraphErrors *) fInputFileAll->Get("TriggerFlowPostSub_V2");
    gTrigger_V4 = (TGraphErrors *) fInputFileAll->Get("TriggerFlowPostSub_V4");
    gTrigger_V6 = (TGraphErrors *) fInputFileAll->Get("TriggerFlowPostSub_V6");
    

    // Better yet, recalculate these, using the track vs ep histograms

    // Loading the phase3 flow histograms


    // Track Vn will be calculated here in FitFlow


  }






  if (!gTrigger_V2 || !gTrigger_V4) {
    fprintf(stderr,"Missing Trigger Flow fits from phase 1\n");
  }
  if (!gTrack_V2 || !gTrack_V4) {
    fprintf(stderr,"Missing Track Flow fits from phase 1\n");
  }

  // Load our trigger PtBins;

  if (!fIsMCGenMode) {
    fAllTriggerPt = (TH1D *) fInputFileAll->Get("fTriggerPt");
  } else {
    fAllTriggerPt = (TH1D *) fInputFileAll->Get("leadingJetPtEP_0");
    if (!fAllTriggerPt) {
      fprintf(stderr,"Could not find leadingJetPtEP_0\n");
      return;
    }
    // Restrict to proper pt range FIXME. Probably only needed for MC at this point
    // FIXME check that I don't need to do this in data
  }

  int iMinTriggerPtBin = 1;
  int iMaxTriggerPtBin = 1;
  //PtBins[6]
  if (fIsMCGenMode) {
    iMinTriggerPtBin = fAllTriggerPt->FindFixBin(0.001 + PtBins[iPtBin-1]);
    iMaxTriggerPtBin = fAllTriggerPt->FindFixBin(-0.001 + PtBins[iPtBin]);
    printf("For MC Gen, restricting trigger pt range to bins %d %d, range %f %f\n",iMinTriggerPtBin,iMaxTriggerPtBin,fAllTriggerPt->GetXaxis()->GetBinLowEdge(iMinTriggerPtBin),fAllTriggerPt->GetXaxis()->GetBinUpEdge(iMaxTriggerPtBin));

    fAllTriggerPt->GetXaxis()->SetRange(iMinTriggerPtBin,iMaxTriggerPtBin);
  }


  for (Int_t i = 0; i < kNEPBins; i++) {
    //fEPBinTriggerPt.push_back((TH1D *) fInputFileEvt[i]->Get("fTriggerPt"));
    if (!fIsMCGenMode) {
      fEPBinTriggerPt.push_back((TH1D *) fInputFileEvt[i]->Get("fTriggerPtWithinEPBin"));
    } else {
      fEPBinTriggerPt.push_back((TH1D *) fInputFileAll->Get(Form("leadingJetPtEP_%d",i+1)));
      // Restrict to proper pt range FIXME
      fEPBinTriggerPt[i]->GetXaxis()->SetRange(iMinTriggerPtBin,iMaxTriggerPtBin);
    }
  }
  printf("HOLA\n");

  // Load the track pt histograms
  fTrackPtProjectionSE = (TH1D *) fInputFileAll->Get("TrackPtProjectionSE");
  fTrackPtProjectionME = (TH1D *) fInputFileAll->Get("TrackPtProjectionME");
  fTrackPtFromTrackPsi = (TH1D *) fInputFileAll->Get("TrackPtFromTrackPsi");

	// This loop will also count the number of bins in our observable
	for (Int_t i = 0; i < 13; i++) { // max num should be >= all possibilities
		TH1D * fLocal = 0;
    TString fLocalName = Form("%s_FullDPhi_ObsBin%d",fDPhiHistName.Data(),i);
//		TString fLocalName = Form("dPhi_ObsBin%d_Full",i);
    if (fIsMCGenMode) fLocalName = Form("Proj_PtBin%d_EP-1_FullDPhi_ObsBin%d",iPtBin,i);
		fLocal = (TH1D *) fInputFileAll->Get(fLocalName);
		if (!fLocal ) {
			if (fDebugLevel) printf("  Histo %s not found.  End of search.",fLocalName.Data());
			break;
		}
		fLocal->SetName(Form("%s_All",fLocalName.Data()));
		fFullDPhiProjAll.push_back(fLocal);

		nObsBins++;
	}
  // Building the vector for the Python results
  fPyBkgParGraphs={};
  fPyRPSParGraphs={};
  for (Int_t i = 0; i < nObsBins; i++) {
    fPyBkgParGraphs.push_back(0);
    fPyRPSParGraphs.push_back(0);
  }

  // For MC Gen
  // Example projection names
  // pi0Hadron_EP_3_pi0Pt_14_17_particlePt_250_420_FarEtaProj
  // pi0Hadron_EP_3_pi0Pt_14_17_particlePt_250_420_Projection

	printf("  Observable Bins: %d\n",nObsBins);
	for (Int_t i = 0; i < nObsBins; i++) { 
		vector<TH1D *> fLocalVector = {};
		for (Int_t j = 0; j < kNEPBins; j++) {
			TH1D * fLocal = 0;
			TString fLocalName = Form("%s_FullDPhi_ObsBin%d",fDPhiHistName.Data(),i);
//			TString fLocalName = Form("dPhi_ObsBin%d_Full",i);
      if (fIsMCGenMode) {
        fLocalName = Form("Proj_PtBin%d_EP%d_FullDPhi_ObsBin%d",iPtBin,j,i);
        fLocal = (TH1D *) fInputFileAll->Get(fLocalName);
      } else {
        fLocal = (TH1D *) fInputFileEvt[j]->Get(fLocalName);
      }

			if (!fLocal ) {
				printf("  Histo %s not found.\n",fLocalName.Data());
				break;
			}
			fLocal->SetName(Form("%s_Evt%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fFullDPhiProj.push_back(fLocalVector);
	}
	
	// Now loading Near Eta Projections
	//   All Evt Plane
	for (Int_t i = 0; i < nObsBins; i++) {
		TH1D * fLocal = 0;
    TString fLocalName = Form("%s_NearEtaDPhi_ObsBin%d",fDPhiHistName.Data(),i);
    if (fIsMCGenMode) fLocalName = Form("Proj_PtBin%d_EP-1_NearEtaDPhi_ObsBin%d",iPtBin,i);
//		TString fLocalName = Form("dPhi_ObsBin%d_NearEta",i);
		fLocal = (TH1D *) fInputFileAll->Get(fLocalName);
		if (!fLocal) {
			printf("Could not find histogram %s in all Evt Plane\n",fLocalName.Data());
		}
		fLocal->SetName(Form("%s_All",fLocalName.Data()));
		fNearEtaDPhiProjAll.push_back(fLocal);
	}
	//   Indiv Evt Plane
	for (Int_t i = 0; i < nObsBins; i++) { 
		vector<TH1D *> fLocalVector = {};
		for (Int_t j = 0; j < kNEPBins; j++) {
			TH1D * fLocal = 0;
      TString fLocalName = Form("%s_NearEtaDPhi_ObsBin%d",fDPhiHistName.Data(),i);
//			TString fLocalName = Form("dPhi_ObsBin%d_NearEta",i);
      if (fIsMCGenMode) {
        //fLocalName = Form("Proj_PtBin%d_EP%d_FullDPhi_ObsBin%d",iPtBin,j,i); // Had this by mistake
        fLocalName = Form("Proj_PtBin%d_EP%d_NearEtaDPhi_ObsBin%d",iPtBin,j,i);
        fLocal = (TH1D *) fInputFileAll->Get(fLocalName);
      } else {
        fLocal = (TH1D *) fInputFileEvt[j]->Get(fLocalName);
      }
			if (!fLocal ) {
				printf("  Histo %s not found in file %s.\n",fLocalName.Data(),fInputFileEvt[j]->GetName());
				break;
			}
			fLocal->SetName(Form("%s_Evt%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fNearEtaDPhiProj.push_back(fLocalVector);
	}

	// Now loading Far Eta Projections
	//   All Evt Plane
	for (Int_t i = 0; i < nObsBins; i++) {
		TH1D * fLocal = 0;
    TString fLocalName = Form("%s_FarEtaDPhi_ObsBin%d",fDPhiHistName.Data(),i);
    if (fIsMCGenMode) fLocalName = Form("Proj_PtBin%d_EP-1_FarEtaDPhi_ObsBin%d",iPtBin,i);
//		TString fLocalName = Form("dPhi_ObsBin%d_FarEta",i);
		fLocal = (TH1D *) fInputFileAll->Get(fLocalName);
		if (!fLocal) {
			printf("Could not find histogram %s in all Evt Plane\n",fLocalName.Data());
		}
		fLocal->SetName(Form("%s_All",fLocalName.Data()));
		fFarEtaDPhiProjAll.push_back(fLocal);
	}
	//   Indiv Evt Plane
	for (Int_t i = 0; i < nObsBins; i++) { 
		vector<TH1D *> fLocalVector = {};
		for (Int_t j = 0; j < kNEPBins; j++) {
			TH1D * fLocal = 0;
      TString fLocalName = Form("%s_FarEtaDPhi_ObsBin%d",fDPhiHistName.Data(),i);
//			TString fLocalName = Form("dPhi_ObsBin%d_FarEta",i);
      if (fIsMCGenMode) {
        fLocalName = Form("Proj_PtBin%d_EP%d_FarEtaDPhi_ObsBin%d",iPtBin,j,i);
        fLocal = (TH1D *) fInputFileAll->Get(fLocalName);
      } else {
        fLocal = (TH1D *) fInputFileEvt[j]->Get(fLocalName);
      }
			if (!fLocal ) {
				printf("  Histo %s not found in file %s.\n",fLocalName.Data(),fInputFileEvt[j]->GetName());
				break;
			}
			fLocal->SetName(Form("%s_Evt%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fFarEtaDPhiProj.push_back(fLocalVector);
	}

  SetEPTitles(fFullDPhiProj);
  SetEPTitles(fNearEtaDPhiProj);
  SetEPTitles(fFarEtaDPhiProj);

}

void TaskEventPlane::LoadPublishedFlow() {
	cout<<"Loading Published Flow"<<endl;

  // Build gAliTrack_V2,3,4
  // Build from fV2Graph2D
  // x axis is pt
  // y axis is cent

  // need to know effective centrality
  // I have iCent Bin = 0,1,2,3a
  // These are just the average centrality.
  // This does not account for the bias of triggers 
  // towards more central events from NBinary, NPart
  double fEffCentArray[4] = {5.,20.,40.,65.};

  // To do: vary the centrality used in the interpolation
  //
  // With the plusminus graphs and sampling two additional centralities
  // that gives 3*3 = 9 values.

  // How to actually calculate the variance? TProfile

  // If I had the charged particle spectrum in pt vs centrality , I 
  // would use that here to weight the vN function in centrality.

  // I do have the charged particle spectra, but not in fine centrality bins.


  
  double fCentBinWidth = fCentArray[iCentBin+1] - fCentArray[iCentBin];

  double fRangeFactor = 0.25;

  double fEffectiveCent = fEffCentArray[iCentBin];
  double fEffCentLow = fEffCentArray[iCentBin] - fRangeFactor * fCentBinWidth;
  double fEffCentHigh = fEffCentArray[iCentBin] + fRangeFactor * fCentBinWidth;
  





  // What pt bins to use? I can just use my default, since i'm already
  // doing some interpolating no matter what.
  // Get those pt bins from that histogram I worked so hard to propagate
  TAxis * fAnalysisPtAxis = fTrackPtProjectionSE->GetXaxis();



  int nBins = fAnalysisPtAxis->GetNbins();

  gAliTrack_V2 = new TGraphErrors(nBins);
  gAliTrack_V2->SetName("AliTrack_V2");
  gAliTrack_V2->SetMarkerStyle(kFullSquare);
  gAliTrack_V2->SetMarkerColor(46);
  gAliTrack_V2->SetLineColor(46);
  gAliTrack_V2->SetTitle("v_{2} (charged tracks)");
  gAliTrack_V2->GetXaxis()->SetTitle("p_{T}^{a} (GeV/#it{c})");
  gAliTrack_V2->GetYaxis()->SetTitle("v_{2}");


  gAliTrack_V3 = new TGraphErrors(nBins);
  gAliTrack_V3->SetName("AliTrack_V3");
  gAliTrack_V3->SetMarkerStyle(kFullSquare);
  gAliTrack_V3->SetMarkerColor(9);
  gAliTrack_V3->SetLineColor(9);
  gAliTrack_V3->SetTitle("v_{3} (charged tracks)");
  gAliTrack_V3->GetXaxis()->SetTitle("p_{T}^{a} (GeV/#it{c})");
  gAliTrack_V3->GetYaxis()->SetTitle("v_{3}");

  gAliTrack_V4 = new TGraphErrors(nBins);
  gAliTrack_V4->SetName("AliTrack_V4");
  gAliTrack_V4->SetMarkerStyle(kFullSquare);
  gAliTrack_V4->SetMarkerColor(8);
  gAliTrack_V4->SetLineColor(8);
  gAliTrack_V4->SetTitle("v_{4} (charged tracks)");
  gAliTrack_V4->GetXaxis()->SetTitle("p_{T}^{a} (GeV/#it{c})");
  gAliTrack_V4->GetYaxis()->SetTitle("v_{4}");


  for (int i = 0; i < nBins; i++) {
    // some of these may be unnecessary
    double fPtChMin = fAnalysisPtAxis->GetBinLowEdge(i+1);
    double fPtChCenter = fAnalysisPtAxis->GetBinCenter(i+1);
    double fPtChErr = 0.5 * fAnalysisPtAxis->GetBinWidth(i+1);
    double fPtChMax = fAnalysisPtAxis->GetBinUpEdge(i+1);
    printf("  PtCh bin %d is [%f,%f)\n",i,fPtChMin,fPtChMax);   

    double fInterpolateV2 = fV2Graph2D->Interpolate(fPtChCenter,fEffectiveCent);
    //double fInterpolateV2Err = 0; // FIXME
    double fInterpolateV2Err = 0.5 * TMath::Abs(fV2Graph2DErrUp->Interpolate(fPtChCenter,fEffectiveCent) - fV2Graph2DErrDown->Interpolate(fPtChCenter,fEffectiveCent));
    // would be better to do use all three values, 

    printf("     2D interpolation finds V2 = %f \\pm %f\n",fInterpolateV2,fInterpolateV2Err);

    gAliTrack_V2->SetPoint(i,fPtChCenter,fInterpolateV2);
    gAliTrack_V2->SetPointError(i,fPtChErr,fInterpolateV2Err);

    double fInterpolateV3 = fV3Graph2D->Interpolate(fPtChCenter,fEffectiveCent);
    //double fInterpolateV3Err = 0; // FIXME
    double fInterpolateV3Err = 0.5 * TMath::Abs(fV3Graph2DErrUp->Interpolate(fPtChCenter,fEffectiveCent) - fV3Graph2DErrDown->Interpolate(fPtChCenter,fEffectiveCent));
    gAliTrack_V3->SetPoint(i,fPtChCenter,fInterpolateV3);
    gAliTrack_V3->SetPointError(i,fPtChErr,fInterpolateV3Err);
    
    double fInterpolateV4 = fV4Graph2D->Interpolate(fPtChCenter,fEffectiveCent);
    //double fInterpolateV4Err = 0; // FIXME
    double fInterpolateV4Err = 0.5 * TMath::Abs(fV4Graph2DErrUp->Interpolate(fPtChCenter,fEffectiveCent) - fV4Graph2DErrDown->Interpolate(fPtChCenter,fEffectiveCent));
    gAliTrack_V4->SetPoint(i,fPtChCenter,fInterpolateV4);
    gAliTrack_V4->SetPointError(i,fPtChErr,fInterpolateV4Err);



  }
}

void TaskEventPlane::FitFlow() {
	cout<<"Fitting flow histograms"<<endl;

  if (fIsMCGenMode) {
    printf("Since this is MCGen mode, we will trust the fitting done in the event generator analysis.\n");
    if (!gTrigger_Bv || !gTrigger_V2 || !gTrigger_V4) {
      fprintf(stderr, "Missing a flow graph from MC analysis\n");

    }
    return;
  }
  // FIXME are these loaded yet?
  Double_t fEPRes_R2 = fEPRes[1];
  Double_t fEPRes_R3 = fEPRes[2];
  Double_t fEPRes_R4 = fEPRes[3];

  Double_t fEPRes_R6 = fEPRes[5];

  Double_t fEP3Res_R2 = fEP3Res[1];
  Double_t fEP3Res_R3 = fEP3Res[2];
  Double_t fEP3Res_R4 = fEP3Res[3];

  Double_t fEP3Res_R6 = fEP3Res[5];

  Double_t fEP4Res_R2 = fEP4Res[1];
  Double_t fEP4Res_R3 = fEP4Res[2];
  Double_t fEP4Res_R4 = fEP4Res[3];

  Double_t fEP4Res_R6 = fEP4Res[5];

  // 
  if (!hHistTrackPsiEPPtCent) {
    fprintf(stderr,"Missing Th3F fHistTrackPsiEPPtCent\n");
  } else {
    printf("Successfully found the Track Th3\n");
  }


  // This requires the pion vs event plane histogram
  // Creating the track Vn graphs

  gTrigger_Bv = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_Bv->SetName("Trigger_Bv");
  gTrigger_Bv->SetTitle("B Value from V_{n} fit");
  gTrigger_Bv->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_Bv->GetYaxis()->SetTitle("B");
  gTrigger_Bv_Presub = (TGraphErrors *) gTrigger_Bv->Clone("Trigger_Bv_Presub");

  gTrigger_V2 = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_V2->SetName("Trigger_V2");
  gTrigger_V2->SetTitle("Calculated #tilde{v}_{2}^{trigger} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrigger_V2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_V2->GetYaxis()->SetTitle("#tilde{v}_{2}");
  gTrigger_V2_Presub = (TGraphErrors *) gTrigger_V2->Clone("Trigger_V2_Presub");

  gTrigger_V4 = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_V4->SetName("Trigger_V4");
  gTrigger_V4->SetTitle("Calculated #tilde{v}_{4}^{trigger} (Event Plane method)");
  gTrigger_V4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_V4->GetYaxis()->SetTitle("#tilde{v}_{4}");
  gTrigger_V4_Presub = (TGraphErrors *) gTrigger_V4->Clone("Trigger_V4_Presub");

  gTrigger_V6 = new TGraphErrors(kUsedPi0TriggerPtBins);
  gTrigger_V6->SetName("Trigger_V6");
  gTrigger_V6->SetTitle("Calculated #tilde{v}_{6}^{trigger} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrigger_V6->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrigger_V6->GetYaxis()->SetTitle("#tilde{v}_{6}");
  gTrigger_V6_Presub = (TGraphErrors *) gTrigger_V6->Clone("Trigger_V6_Presub");



  // Values for the 2nd order event plane

  gTrack_Bv = new TGraphErrors(kNTrackPtBins);
  gTrack_Bv->SetName("Track_Bv");
  gTrack_Bv->SetTitle("B Value from V_{n} fit");
  gTrack_Bv->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_Bv->GetYaxis()->SetTitle("B");

  gTrack_V2 = new TGraphErrors(kNTrackPtBins);
  gTrack_V2->SetName("Track_V2");
  gTrack_V2->SetTitle("Calculated #tilde{v}_{2}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrack_V2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V2->GetYaxis()->SetTitle("#tilde{v}_{2}");
  gTrack_V2->SetMarkerStyle(kOpenSquare);

  //gTrack_V3 = new TGraphErrors(kNTrackPtBins);
  //gTrack_V3->SetName("Track_V3");
  //gTrack_V3->SetTitle("Calculated #tilde{v}_{3}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  //gTrack_V3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  //gTrack_V3->GetYaxis()->SetTitle("#tilde{v}_{3}");

  gTrack_V4 = new TGraphErrors(kNTrackPtBins);
  gTrack_V4->SetName("Track_V4");
  gTrack_V4->SetTitle("Calculated #tilde{v}_{4}^{Track} (Event Plane method)");
  gTrack_V4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V4->GetYaxis()->SetTitle("#tilde{v}_{4}");

  gTrack_V6 = new TGraphErrors(kNTrackPtBins);
  gTrack_V6->SetName("Track_V6");
  gTrack_V6->SetTitle("Calculated #tilde{v}_{6}^{Track} (Event Plane method)"); // n in the title for the purpose of the drawn graph
  gTrack_V6->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V6->GetYaxis()->SetTitle("#tilde{v}_{6}");

  // Values for the 3rd order event plane
  gTrack_Bv_EP3 = new TGraphErrors(kNTrackPtBins);
  gTrack_Bv_EP3->SetName("Track_Bv_EP3");
  gTrack_Bv_EP3->SetTitle("B Value from V_{3} fit");
  gTrack_Bv_EP3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_Bv_EP3->GetYaxis()->SetTitle("B");

  gTrack_V3_EP3 = new TGraphErrors(kNTrackPtBins);
  gTrack_V3_EP3->SetName("Track_V3_EP3");
  gTrack_V3_EP3->SetTitle("Calculated #tilde{v}_{3}^{Track} (Event Plane method)");
  gTrack_V3_EP3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V3_EP3->GetYaxis()->SetTitle("#tilde{v}_{3}");
  gTrack_V3_EP3->SetMarkerStyle(kOpenSquare);
  //  will need EPR{3,3}

  // Values for the 4th order event plane
  //   
  gTrack_Bv_EP4 = new TGraphErrors(kNTrackPtBins);
  gTrack_Bv_EP4->SetName("Track_Bv_EP4");
  gTrack_Bv_EP4->SetTitle("B Value from V_{4} fit");
  gTrack_Bv_EP4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_Bv_EP4->GetYaxis()->SetTitle("B");

  gTrack_V4_EP4 = new TGraphErrors(kNTrackPtBins);
  gTrack_V4_EP4->SetName("Track_V4_EP4");
  gTrack_V4_EP4->SetTitle("Calculated #tilde{v}_{4}^{Track} (Event Plane method)");
  gTrack_V4_EP4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gTrack_V4_EP4->GetYaxis()->SetTitle("#tilde{v}_{4}");
  //  will need EPR{4,4}



  // Copied from phase1
  hHistTrackPsiEPPtCent->GetZaxis()->SetRange(iCentBin+1,iCentBin+1);
  TH2F * hHistTrackPsiEPPt = (TH2F *) hHistTrackPsiEPPtCent->Project3D("yx");
  hHistTrackPsiEPPt->SetTitle("Track #it{p}_{T} vs #Delta#Psi_{EP}");

  for (int i = 0; i < kNTrackPtBins; i++) {
    printf("Doing flow pt things for track pt bin %d\n",i);
    double fMinPt = fTrackPtBins[i];
    double fMaxPt = fTrackPtBins[i+1];

    int iMinBin = hHistTrackPsiEPPtCent->GetYaxis()->FindBin(fMinPt);
    int iMaxBin = hHistTrackPsiEPPtCent->GetYaxis()->FindBin(fMaxPt) - 1;

    printf("Projecting Tracks Event Plane Histograms in bins from %d to %d (%.2f,%.2f)\n",iMinBin,iMaxBin,fMinPt,fMaxPt);
    //printf("Projecting Tracks Event Plane Histograms in bin from %.1f to %.1f\n",hPtEPAnglePionAcc->GetYaxis()->GetBinLowEdge(iMinBin),hPtEPAnglePionAcc->GetYaxis()->GetBinUpEdge(iMaxBin));

    TString sFormat = "%s_Proj_%d";
    TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

//(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i),iMinBin,iMaxBin);
    hHistTrackPsiEPPtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
    TH1F * hLocalPtEPAngleTrack_Proj = (TH1F *) hHistTrackPsiEPPtCent->Project3D("xe");
    hLocalPtEPAngleTrack_Proj->SetName(Form(sFormat.Data(),hLocalPtEPAngleTrack_Proj->GetName(),i));
//    hLocalPtEPAngleTrack_Proj->SetName(Form(sFormat.Data(),hPtEPAnglePionAcc->GetName(),i));
    //hLocalPtEPAnglePionAcc_Proj->Sumw2();
    hLocalPtEPAngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP} (%s)",sPtRange.Data()));
    hLocalPtEPAngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
    hPtEPAngleTrack_Proj.push_back(hLocalPtEPAngleTrack_Proj);
  }

  if (hHistTrackPsiEP3PtCent) {
    printf("Found the TH3 for tracks relative to 3rd order event plane\n");

    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      int iMinBin = hHistTrackPsiEP3PtCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hHistTrackPsiEP3PtCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

      TString sFormat = "%s_EP3Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

      hHistTrackPsiEP3PtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
      TH1F * hLocalPtEP3AngleTrack_Proj = (TH1F *) hHistTrackPsiEP3PtCent->Project3D("xe");
      hLocalPtEP3AngleTrack_Proj->SetName(Form(sFormat.Data(),hHistTrackPsiEP3PtCent->GetName(),i));
      //hLocalPtEP3AnglePionAcc_Proj->Sumw2();
      hLocalPtEP3AngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP,3} (%s)",sPtRange.Data()));
      hLocalPtEP3AngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
      hPtEP3AngleTrack_Proj.push_back(hLocalPtEP3AngleTrack_Proj);

    }
  } else {
    printf("Did not find the TH3 for tracks relative to 3rd order event plane\n");
  }


  if (hHistTrackPsiEP4PtCent) {
    printf("Found the TH3 for tracks relative to 4th order event plane (which is not reconstructed well in V0.\n");

    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      int iMinBin = hHistTrackPsiEP4PtCent->GetYaxis()->FindBin(fMinPt);
      int iMaxBin = hHistTrackPsiEP4PtCent->GetYaxis()->FindBin(fMaxPt) - 1; // Want the bin with fMaxPt as an upper bound

      TString sFormat = "%s_EP4Proj_%d";
      TString sPtRange = Form("%.2f #leq #it{p}_{T} < %.2f GeV/#it{c}",fMinPt,fMaxPt);

      hHistTrackPsiEP4PtCent->GetYaxis()->SetRange(iMinBin,iMaxBin);
      TH1F * hLocalPtEP4AngleTrack_Proj = (TH1F *) hHistTrackPsiEP4PtCent->Project3D("xe");
      hLocalPtEP4AngleTrack_Proj->SetName(Form(sFormat.Data(),hHistTrackPsiEP4PtCent->GetName(),i));
      //hLocalPtEP3AnglePionAcc_Proj->Sumw2();
      hLocalPtEP4AngleTrack_Proj->SetTitle(Form("Track #Delta#Psi_{EP,4} (%s)",sPtRange.Data()));
      hLocalPtEP4AngleTrack_Proj->GetYaxis()->SetTitle("N_{Tracks}");
      hPtEP4AngleTrack_Proj.push_back(hLocalPtEP4AngleTrack_Proj);
    }
  }



  TCanvas * cVn = new TCanvas("cVn","cVn");


  printf("Starting the fitting of tracks vs the event plane\n");
  TH1F * hTrackEP = 0;
  for (int i = 0; i < kNTrackPtBins; i++) {
    printf("Fitting track pt bin %d\n",i);
    double fMinPt = fTrackPtBins[i];
    double fMaxPt = fTrackPtBins[i+1];

    printf("Hello\n");
    // FIXME
    hTrackEP = hPtEPAngleTrack_Proj[i];
    if (!hTrackEP) {
      printf("Missing the track vs EP histogram\n");
      return;
    }

    TF1 * fitTrackEP = new TF1(Form("Track_VnFit_%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x)+2*[3]*TMath::Cos(6*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
    fitTrackEP->SetParameter(0,hTrackEP->Integral("width") / (TMath::Pi() / 2));
    fitTrackEP->SetParameter(1,0.01);
    fitTrackEP->SetParameter(2,0.001);
    fitTrackEP->SetParameter(3,0.);

    fitTrackEP->SetParLimits(1,0.0,0.5);
    fitTrackEP->SetParLimits(2,0.0,0.25);
    fitTrackEP->SetParLimits(3,0.0,0.25);

    fitTrackEP->SetParName(0,"B");
    fitTrackEP->SetParName(1,"v_2");
    fitTrackEP->SetParName(2,"v_4");
    fitTrackEP->SetParName(3,"v_6");

    hTrackEP->Fit(fitTrackEP);
    fitTrackEP->SetLineColor(kCyan);

    hTrackEP->Draw();
    fitTrackEP->Draw("SAME");

    printf("Using EPR2 = %f, EPR4 = %f\n",fEPRes_R2,fEPRes_R4);

    printf("Found v2(track) = %f \\pm %f\n",fitTrackEP->GetParameter(1)/fEPRes_R2,fitTrackEP->GetParError(1)/fEPRes_R2);
    printf("Found v4(track) = %f \\pm %f\n",fitTrackEP->GetParameter(3)/fEPRes_R4,fitTrackEP->GetParError(3)/fEPRes_R4);


    gTrack_Bv->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(0));
    gTrack_V2->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(1)/fEPRes_R2);
    gTrack_V4->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(2)/fEPRes_R4);
    gTrack_V6->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(3)/fEPRes_R6);

    gTrack_Bv->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(0));
    gTrack_V2->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(1)/fEPRes_R2);
    gTrack_V4->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(2)/fEPRes_R4);
    gTrack_V6->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(3)/fEPRes_R6);

    cVn->Print(Form("%s/EPStudy_Track_Pt_%.2f_%.2f.pdf",fOutputDir.Data(),fMinPt,fMaxPt));
    cVn->Print(Form("%s/CFiles/EPStudy_Track_Pt_%.2f_%.2f.C",fOutputDir.Data(),fMinPt,fMaxPt));

  }

  // FIXME calculate v3(ep3) v4(ep4)
  if (hHistTrackPsiEP3PtCent) {
    printf("Calculating V3\n");

    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      hTrackEP = hPtEP3AngleTrack_Proj[i];
      if (!hTrackEP) {
        printf("Missing the v3 track vs EP3 histograms\n");
        return;
      }
      printf("Doing the v3 thing with histogram %s (%s)\n",hTrackEP->GetName(),hTrackEP->GetTitle());
      TF1 * fitTrackEP = new TF1(Form("Track_V3Fit_%d",i),"[0]*(1+2*[1]*TMath::Cos(3*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
      fitTrackEP->SetParameter(0,hTrackEP->Integral("width") / (TMath::Pi() / 2));
      fitTrackEP->SetParameter(1,0.01);

      fitTrackEP->SetParName(0,"B");
      fitTrackEP->SetParName(1,"v_3");

      hTrackEP->Fit(fitTrackEP);
      fitTrackEP->SetLineColor(kCyan);

      hTrackEP->Draw();
      fitTrackEP->Draw("SAME");
//      gTrack_V3_EP3->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(1)/fEPRes_R3);
//      gTrack_V3_EP3->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(1)/fEPRes_R3);
      gTrack_V3_EP3->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(1)/fEP3Res_R3);
      gTrack_V3_EP3->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(1)/fEP3Res_R3);

      cVn->Print(Form("%s/EPStudy_Track_Pt_%.2f_%.2f_v3.pdf",fOutputDir.Data(),fMinPt,fMaxPt));
      cVn->Print(Form("%s/CFiles/EPStudy_Track_Pt_%.2f_%.2f_v3.C",fOutputDir.Data(),fMinPt,fMaxPt));
    }
  }


  if (hHistTrackPsiEP4PtCent) {
    printf("Calculating V4 (EP4), which won't work well.\n");

    for (int i = 0; i < kNTrackPtBins; i++) {
      double fMinPt = fTrackPtBins[i];
      double fMaxPt = fTrackPtBins[i+1];

      hTrackEP = hPtEP4AngleTrack_Proj[i];
      if (!hTrackEP) return;
      printf("Doing the v4 thing with histogram %s (%s)\n",hTrackEP->GetName(),hTrackEP->GetTitle());
      TF1 * fitTrackEP = new TF1(Form("Track_V4Fit_%d",i),"[0]*(1+2*[1]*TMath::Cos(4*x))",hTrackEP->GetXaxis()->GetXmin(),hTrackEP->GetXaxis()->GetXmax());
      fitTrackEP->SetParameter(0,hTrackEP->Integral("width") / (TMath::Pi() / 2));
      fitTrackEP->SetParameter(1,0.01);

      fitTrackEP->SetParName(0,"B");
      fitTrackEP->SetParName(1,"v_{4,EP4}");

      hTrackEP->Fit(fitTrackEP);
      fitTrackEP->SetLineColor(kCyan);

      hTrackEP->Draw();
      fitTrackEP->Draw("SAME");
      gTrack_V4_EP4->SetPoint(i,(fMinPt+fMaxPt)/2.,fitTrackEP->GetParameter(1)/fEPRes_R4);
      gTrack_V4_EP4->SetPointError(i,(fMaxPt-fMinPt)/2.,fitTrackEP->GetParError(1)/fEPRes_R4);

      cVn->Print(Form("%s/EPStudy_Track_Pt_%.2f_%.2f_v4.pdf",fOutputDir.Data(),fMinPt,fMaxPt));
      cVn->Print(Form("%s/CFiles/EPStudy_Track_Pt_%.2f_%.2f_v4.C",fOutputDir.Data(),fMinPt,fMaxPt));
    }
  }

  // Fitting the Vn Values

  // Parameter Fit Functions
  // v2 = V2FP_0 * TMath::Landau(pt,V2FP_1,V2FP_2,false);
  // range? Could check the last point of a graph. but not a Py graph
  // magic number for now
  // Example params:
  //          V2FP_0 = 5.39182e-01;//   8.36928e-03
  //         V2FP_1 = 3.70690e+00;//   1.94541e-01
  //        V2FP_2 = 1.54252e+00;//   9.17354e-02
  //        V4FP_0 = 3.89071e-01;//   2.13836e-02
  //        V4FP_1 = 3.43959e+00;//   1.90255e-01
  //        V4FP_2 = 1.23961e+00;//   8.31609e-02

  

  TF1 * fV2A_Fit = new TF1("V2A_Fit","[0]*(TMath::Landau(x,[1],[2],0))",0.5,18.);
  PrepLandauFit(fV2A_Fit);
  fV2A_Fit->SetLineColor(kRed-1);
  fV2A_Fit->SetLineStyle(2);
  TF1 * fV3A_Fit = new TF1("V3A_Fit","[0]*(TMath::Landau(x,[1],[2],0))",0.5,18.);
  PrepLandauFit(fV3A_Fit);
  fV3A_Fit->SetLineColor(kRed-1);
  fV3A_Fit->SetLineStyle(2);
  TF1 * fV4A2_Fit = new TF1("V4A2_Fit","[0]*(TMath::Landau(x,[1],[2],0))",0.5,18.);
  PrepLandauFit(fV4A2_Fit);
  fV4A2_Fit->SetLineColor(kRed-1);
  fV4A2_Fit->SetLineStyle(2);
  TF1 * fV6A2_Fit = new TF1("V6A2_Fit","[0]*(TMath::Landau(x,[1],[2],0))",0.5,18.);
  PrepLandauFit(fV6A2_Fit);
  fV6A2_Fit->SetLineColor(kRed-1);
  fV6A2_Fit->SetLineStyle(2);

  cout<<"Fitting V2A_Fit"<<endl;
  gTrack_V2->Fit(fV2A_Fit,"R");
  PrintLandauFit(fV2A_Fit,2);

  cout<<"Fitting V3A_Fit"<<endl;
  gTrack_V3_EP3->Fit(fV3A_Fit,"R");
  PrintLandauFit(fV3A_Fit,3);
  
  cout<<"Fitting V4A2_Fit"<<endl;
  gTrack_V4->Fit(fV4A2_Fit,"R");
  PrintLandauFit(fV4A2_Fit,4);

  cout<<"Fitting V6A2_Fit"<<endl;
  gTrack_V6->Fit(fV6A2_Fit,"R");
  PrintLandauFit(fV6A2_Fit,6);



  cVn->Clear();
  // Draw and Compare flow from different sources
  TLegend * legCompareSource = new TLegend(0.7,0.6,0.9,0.85);
  legCompareSource->AddEntry(gTrack_V2,"This analysis {EP2}","lp");
  // Could add in additonal notes, like the deltaEta cut
  // or how the interpolation is done.
  legCompareSource->AddEntry(gAliTrack_V2,"ALICE Published {2}","lp");
  TMultiGraph * mgV2 = new TMultiGraph();
  mgV2->Add(gTrack_V2);
  mgV2->Add(gAliTrack_V2);
  mgV2->GetXaxis()->SetTitle(gAliTrack_V2->GetXaxis()->GetTitle());
  mgV2->GetYaxis()->SetTitle(gAliTrack_V2->GetYaxis()->GetTitle());
  mgV2->Draw("ALP");

 // gTrack_V2->Draw("ALP");
 // gAliTrack_V2->Draw("LP SAME");
  legCompareSource->Draw("SAME");
  cVn->Print(Form("%s/FlowCmp_v2_Sources.pdf",fOutputDir.Data()));
  cVn->Print(Form("%s/FlowCmp_v2_Sources.png",fOutputDir.Data()));

  cVn->Clear();
  legCompareSource->Clear();
  TMultiGraph * mgV3 = new TMultiGraph();
  // FIXME missing y analysis of v3?
  legCompareSource->AddEntry(gTrack_V3_EP3,"This analysis","lp");
  legCompareSource->AddEntry(gAliTrack_V3,"ALICE Published {2}","lp");
//  gTrack_V3_EP3->Draw("ALP");
//  gAliTrack_V3->Draw("ALP");
//  gAliTrack_V3->Draw("LP SAME");
  mgV3->Add(gTrack_V3_EP3);
  mgV3->Add(gAliTrack_V3);
  mgV3->GetXaxis()->SetTitle(gAliTrack_V3->GetXaxis()->GetTitle());
  mgV3->GetYaxis()->SetTitle(gAliTrack_V3->GetYaxis()->GetTitle());
  mgV3->Draw("ALP");

  legCompareSource->Draw("SAME");
  cVn->Print(Form("%s/FlowCmp_v3_Sources.pdf",fOutputDir.Data()));
  cVn->Print(Form("%s/FlowCmp_v3_Sources.png",fOutputDir.Data()));

  cVn->Clear();
  legCompareSource->Clear();
  legCompareSource->AddEntry(gTrack_V4,"This analysis {EP,2}","lp");
  legCompareSource->AddEntry(gAliTrack_V4,"ALICE Published {2}","lp");

  TMultiGraph * mgV4 = new TMultiGraph();
  mgV4->Add(gTrack_V4);
  mgV4->Add(gAliTrack_V4);
  mgV4->GetXaxis()->SetTitle(gAliTrack_V4->GetXaxis()->GetTitle());
  mgV4->GetYaxis()->SetTitle(gAliTrack_V4->GetYaxis()->GetTitle());
  mgV4->Draw("ALP");

 // gTrack_V4->Draw("ALP");
 // gAliTrack_V4->Draw("LP SAME");
  legCompareSource->Draw("SAME");
  cVn->Print(Form("%s/FlowCmp_v4_Sources.pdf",fOutputDir.Data()));
  cVn->Print(Form("%s/FlowCmp_v4_Sources.png",fOutputDir.Data()));

}



void TaskEventPlane::InitArrays() {
	cout<<"Initializing Arrays ..."<<endl;
	
  // Setting Event Plane Resolutions
  //Int_t iEPRSet;   Int_t iCentBin; 
  //Double_t kEPRes[6] = {0.0,0.885,0.605,0.245,0.0,0.1};

//  Double_t fEPRes[4][6] = {0}; // [Cent][N]
//  fEPRes = {0}; // [Cent][N]
  //             {R_1,   R_2,   R_3,   R_4,  R_5,  R_6}
  // fEPRes[0] = {  0,  0.73,  0.62, 0.275,  0.0,  0.0};  
  // fEPRes[1] = {  0, 0.885, 0.605, 0.245,  0.0,  0.0};
  // fEPRes[2] = {  0,  0.85,  0.49,  0.21,  0.0,  0.0};
  // fEPRes[3] = {  0,  0.58,  0.22,  0.08,  0.0,  0.0};

  // FIXME add error bars

  // EPR Set EP2, from QnVector MB (Train unknown, prior to T59)
//  Double_t fEPRes_Set_0[4][6]  = {
 //                   {  0.765960, 0.619163,  0.509267, 0.348666, 0.318429, 0.187868},  
 //                   {  0.858157, 0.822691, 0.692985, 0.580624,  0.502229,  0.375755},
 //                   {  0.832549,  0.771133,  0.639423,  0.507014,  0.439729,  0.305388},
 //                   {  0.704550,  0.445893,  0.380824,  0.196809,  0.211605,  0.084895}};
  // EPR Set EP2, from QnVector MB (T59)
  Double_t fEPRes_Set_0[4][6]  = {
  //{      R_1,      R_2,      R_3,      R_4,      R_5,      R_6}
    { 0.765625, 0.619035, 0.509148, 0.348548, 0.318585, 0.187958},
    { 0.858052, 0.822638, 0.693049, 0.580605, 0.502134, 0.375592},
    { 0.832519, 0.771184, 0.639522, 0.507076, 0.439909, 0.305442},
    { 0.704443, 0.445964, 0.380898, 0.196861, 0.211537, 0.084192}
  };
  Double_t fEPRes_Set_0_Err[4][6]  = {
    { 0.000434, 0.000548, 0.000690, 0.000979, 0.001097, 0.001653},
    { 0.000262, 0.000233, 0.000343, 0.000435, 0.000517, 0.000681},
    { 0.000275, 0.000272, 0.000380, 0.000495, 0.000581, 0.000791},
    { 0.000271, 0.000424, 0.000516, 0.000847, 0.000876, 0.001653}
  };

/*  Double_t fEPRes_Set_1[4][6]  = {
                    {  0,  0.6192508430757114,  0., 0.34878092755772117,  0.0,  0.18777865138044672},  
                    {  0, 1., 0, 0.,  0.0,  0.0},   // Don't have Cent1, Cent3
                    {  0,  0.7703651242647157,  0.,  0.5046126852106662,  0.0,  0.3020062445564112},
                    {  0,  1.,  0.,  0.,  0.0,  0.0}  };
*/

  // EPR Set EP2, from QnVector EGA (T59)
  Double_t fEPRes_Set_1[4][6]  = {
                    {  0.767842, 0.622406, 0.525535, 0.390981, 0.352418, 0.236852},
                    {  0.781089, 0.651664, 0.549970, 0.423973, 0.376066, 0.259773},
                    {  0.779045, 0.646402, 0.544922, 0.414330, 0.369320, 0.252558},
                    {  0.762544, 0.610570, 0.515899, 0.376056, 0.343955, 0.226055}};
  Double_t fEPRes_Set_1_Err[4][6]  = {
    { 0.000986, 0.001280, 0.001574, 0.002121, 0.002375, 0.003367},
    { 0.000894, 0.001117, 0.001381, 0.001821, 0.002049, 0.002843},
    { 0.001326, 0.001662, 0.002064, 0.002738, 0.003081, 0.004337},
    { 0.001927, 0.002536, 0.003097, 0.004162, 0.004653, 0.006587} 
  };

  // Read off from some plot
  Double_t fEPRes_Set_2[4][6]  = {
                    {  0,  0.73,  0.62, 0.275,  0.0,  0.0},  
                    {  0, 0.885, 0.605, 0.245,  0.0,  0.0},
                    {  0,  0.85,  0.49,  0.21,  0.0,  0.0},
                    {  0,  0.58,  0.22,  0.08,  0.0,  0.0}  };

  // MCGen (EPR = 1)
  Double_t fEPRes_Set_3[4][6]  = {
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0},  
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0}  };
  Double_t fEPRes_Set_3_Err[4][6]  = {
                    {  0, 0.0, 0.0, 0.0, 0.0, 0.0},  
                    {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {  0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {  0, 0.0, 0.0, 0.0, 0.0, 0.0}  };

  // EPR Set EP3 from QnVector MB (T59)
  Double_t fEP3Res_Set_0[4][6]  = { // FIXME
                    {  0.845080, 0.525658, 0.393520, 0.391907, 0.259402, 0.139690},
                    {  0.842722, 0.513905, 0.372819, 0.375413, 0.246513, 0.125048},  
                    {  0.833544, 0.461511, 0.265899, 0.299722, 0.200492, 0.065516},
                    {  0.825777, 0.407792, 0.087980, 0.204440, 0.162711, 0.014950}};

  // EPR Set EP3 from QnVector EGA (T59)
  Double_t fEP3Res_Set_1[4][6]  = { // FIXME
                    {  0.835607, 0.469750, 0.275954, 0.307249, 0.209391, 0.087959},
                    {  0.834803, 0.466532, 0.271198, 0.304266, 0.209051, 0.092139},
                    {  0.834278, 0.461766, 0.257766, 0.294599, 0.200117, 0.075660},
                    {  0.832765, 0.452335, 0.234578, 0.280758, 0.201374, 0.092235}};


  // EPR Set EP4 from QnVector MB (T59)
  Double_t fEP4Res_Set_0[4][6]  = {
                    {  0.926312, 0.712967, 0.384993, 0.141703, 0.308486, 0.569209},
                    {  0.925901, 0.711299, 0.380526, 0.133776, 0.308111, 0.574361},
                    {  0.924596, 0.706943, 0.369987, 0.095097, 0.306386, 0.600262},
                    {  0.923124, 0.703032, 0.365919, 0.029369, 0.330021, 0.652974}};

  // EPR Set EP4 from QnVector EGA (T59)
  Double_t fEP4Res_Set_1[4][6]  = {
                    {  0.924481, 0.707366, 0.374302, 0.101455, 0.312070, 0.598036},
                    {  0.924401, 0.706881, 0.372533, 0.097955, 0.312920, 0.605455},
                    {  0.924679, 0.708054, 0.375002, 0.096597, 0.312343, 0.609336},
                    {  0.924669, 0.708298, 0.375728, 0.090440, 0.311984, 0.603048}};
 
  // MC (1.0)
  Double_t fEP3Res_Set_3[4][6]  = {
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0 }};
    Double_t fEP3Res_Set_3_Err[4][6]  = {
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
  Double_t fEP4Res_Set_3[4][6]  = {
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0.0, 1.0, 1.0, 1.0, 1.0, 1.0 }};
    Double_t fEP4Res_Set_3_Err[4][6]  = {
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};



  switch (iEPRSet) {
    default:
    case 0: // From QnVector MB
      printf("Using EPR Set 0 (MB)\n");
      memcpy(fEPRes,fEPRes_Set_0[iCentBin], sizeof(fEPRes));
      memcpy(fEPRes_Err,fEPRes_Set_0_Err[iCentBin], sizeof(fEPRes_Err));
      memcpy(fEP3Res,fEP3Res_Set_0[iCentBin], sizeof(fEP3Res));
      memcpy(fEP4Res,fEP4Res_Set_0[iCentBin], sizeof(fEP4Res));


     /* fEPRes[0] = {  0,  0.73,  0.62, 0.275,  0.0,  0.0};  
      fEPRes[1] = {  0, 0.885, 0.605, 0.245,  0.0,  0.0};
      fEPRes[2] = {  0,  0.85,  0.49,  0.21,  0.0,  0.0};
      fEPRes[3] = {  0,  0.58,  0.22,  0.08,  0.0,  0.0}; */
      break;
    case 1: // From QnVector EGA 
      printf("Using EPR Set 1 (EGA)\n");
      memcpy(fEPRes,fEPRes_Set_1[iCentBin], sizeof(fEPRes));
      memcpy(fEPRes_Err,fEPRes_Set_1_Err[iCentBin], sizeof(fEPRes_Err));
      memcpy(fEP3Res,fEP3Res_Set_1[iCentBin], sizeof(fEP3Res));
      memcpy(fEP4Res,fEP4Res_Set_1[iCentBin], sizeof(fEP4Res));

      //old// From Raymond's analysis. Cent0M with QnVector correction
      break;
    case 2: // read off of that one graph. TPC R_n values
      memcpy(fEPRes,fEPRes_Set_2[iCentBin], sizeof(fEPRes));
      break;
    case 3: // Full resolution. Ideal for MC
      printf("Using EPR Set 3 (MCGen, All EPRs set to 1)\n");
      memcpy(fEPRes,fEPRes_Set_3[iCentBin], sizeof(fEPRes));
      memcpy(fEPRes_Err,fEPRes_Set_3_Err[iCentBin], sizeof(fEPRes_Err));
      memcpy(fEP3Res,fEPRes_Set_3[iCentBin], sizeof(fEP3Res));
      memcpy(fEP4Res,fEPRes_Set_3[iCentBin], sizeof(fEP4Res));
 
  }

  for (Int_t i = 0; i < RPF_Functor::kTotalNumberOfRn; i++) {
    printf("Loading resolution R_%d = %f\n",i+1,fEPRes[i]);
  }

  // FIXME Make TGraphErrors of the RPF here, draw.
  
  fEP2RGraph = new TGraphErrors(6);
  fEP3RGraph = new TGraphErrors(6);
  for (int i = 0; i < 6; i++) {
    printf("Setting point %d\n",i);
    fEP2RGraph->SetPoint(i,i+1,fEPRes[i]);
    fEP2RGraph->SetPointError(i,0,fEPRes_Err[i]);

    fEP3RGraph->SetPoint(i,i+1,fEP3Res[i]);
    fEP3RGraph->SetPointError(i,0,0); // Haven't inlcuded EP3 errors
  }
  fEP2RGraph->SetName("EP2RGraph");
  fEP2RGraph->SetTitle("Event Plane Resolutions (w.r.t. 2nd order EP)");
  fEP2RGraph->GetXaxis()->SetTitle("n (order)");
  fEP2RGraph->GetYaxis()->SetTitle("R_{n,2}^{V0}");
  fEP2RGraph->SetMarkerStyle(kFullSquare);
  fEP2RGraph->SetMarkerColor(kAzure);
  fEP2RGraph->SetLineColor(kAzure);

  TCanvas * cCanvasEPR = new TCanvas("CanvasEPR","CanvasEPR");
  fEP2RGraph->Draw("ALP");
  cCanvasEPR->Print(Form("%s/EventPlane2_Res.pdf",fOutputDir.Data()));
  cCanvasEPR->Print(Form("%s/EventPlane2_Res.png",fOutputDir.Data()));

  fEP3RGraph->SetName("EP3RGraph");
  fEP3RGraph->SetTitle("Event Plane Resolutions (w.r.t. 3rd order EP)");
  fEP3RGraph->GetXaxis()->SetTitle("n (order)");
  fEP3RGraph->GetYaxis()->SetTitle("R_{n,3}^{V0}");
  fEP3RGraph->SetMarkerStyle(kFullSquare);
  fEP3RGraph->SetMarkerColor(kOrange+1);
  fEP3RGraph->SetLineColor(kOrange+1);

  fEP3RGraph->Draw("ALP");
  cCanvasEPR->Print(Form("%s/EventPlane3_Res.pdf",fOutputDir.Data()));
  cCanvasEPR->Print(Form("%s/EventPlane3_Res.png",fOutputDir.Data()));






	Double_t fZtStep = 1.0/(7 - 1.0);
	Double_t fXiStep = 2.5/(8 - 1.0);
	
	Double_t * ObsArray = 0;

  Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,20,23,30,60};
  Double_t array_ZT_BinsValue[kZtNBINS+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
  Double_t array_XI_BinsValue[kXiNBINS+1]   ={-100,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};
  Double_t array_HPT_BinsValue[kNoHPtBins+1]={0.15,0.4,0.8,1.45,2.5,4.2,6.95,11.4,18.6};


	if (fObservable == 0) ObsArray = array_G_BinsValue;
	else if (fObservable == 1) ObsArray = array_ZT_BinsValue;
	//else if (fObservable == 2) ObsArray = array_XI_BinsValue;
	else if (fObservable == 2) ObsArray = array_HPT_BinsValue;

	for (Int_t i = 0; i <= nObsBins; i++) {
		fObsBins.push_back(ObsArray[i]);
	}
}


/**
 *	Draw the in,mid,out, and all histograms in row
 */
void TaskEventPlane::DrawRawOmniPlots() {
	cout<<"Drawing Omni Plots"<<endl;
	
  // Note: Rescaling the All event histograms here

 // Bool_t bUnifiedNorm = true; // True if each input is normalized to the total triggers, not just the triggers within that EP bin

  // FIXME why am I scaling these down by 1/3?

	// Full Projections
	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fDPhiSet = {}; 
		for (Int_t j = 0; j < kNEPBins; j++) {
			fDPhiSet.push_back(fFullDPhiProj[i][j]);
		}
    if (bUnifiedNorm) fFullDPhiProjAll[i]->Scale(1./kNEPBins); // Don't do this if normalizing per triggers within the same EPbin
		fDPhiSet.push_back(fFullDPhiProjAll[i]);
		DrawOmniPlots_Type(fDPhiSet,Form("FullProj_ObsBin%d",i));
	}

	// Near Eta Projections
	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fDPhiSet = {}; 
		for (Int_t j = 0; j < kNEPBins; j++) {
			fDPhiSet.push_back(fNearEtaDPhiProj[i][j]);
		}
    if (bUnifiedNorm) fNearEtaDPhiProjAll[i]->Scale(1./kNEPBins);
		fDPhiSet.push_back(fNearEtaDPhiProjAll[i]);
		DrawOmniPlots_Type(fDPhiSet,Form("NearEta_ObsBin%d",i));
	}

	// Far Eta Projections
	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fDPhiSet = {}; 
		for (Int_t j = 0; j < kNEPBins; j++) {
			fDPhiSet.push_back(fFarEtaDPhiProj[i][j]);
		}
    if (bUnifiedNorm) fFarEtaDPhiProjAll[i]->Scale(1./kNEPBins);
		fDPhiSet.push_back(fFarEtaDPhiProjAll[i]);
		DrawOmniPlots_Type(fDPhiSet,Form("FarEta_ObsBin%d",i));
	}
	cout<<"Finished Drawing Omni Plots"<<endl;
}


/**
 *	Draw the in,mid,out, and all histograms in row
 * draw these for the histograms that have been rescaled after subtraction
 */
void TaskEventPlane::DrawRescaleOmniPlots() {
	cout<<"Drawing Omni Plots for Rescaled histograms"<<endl;

  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    for (Int_t i = 0; i < nObsBins; i++) {
      printf("Rescale  Obs Bin %d \n",i);
      vector<TH1D *> fDPhiSet = {}; 
      for (Int_t j = 0; j < kNEPBins; j++) {
        printf("  [%d,%d] Debug:\n",i,j);
        printf("    fFullDPhiProj_Rescale array has length %d\n",(int) fFullDPhiProj_Rescale[iV].size());
        if (!fFullDPhiProj_Rescale[iV][i][j]) {
          printf("Missing the rescaled projection somehow!!\n");
        }
        printf("  [%d,%d] Including histogram of name %s\n",i,j,fFullDPhiProj_Rescale[iV][i][j]->GetName());
        fDPhiSet.push_back(fFullDPhiProj_Rescale[iV][i][j]);
      }
      fDPhiSet.push_back(fFullDPhiProjAll[i]);
      DrawOmniPlots_Type(fDPhiSet,Form("Rescale_RPFMethod%d_FullProj_ObsBin%d",iV,i));
    }


    for (Int_t i = 0; i < nObsBins; i++) {
      printf("Rescale  Obs Bin %d \n",i);
      vector<TH1D *> fDPhiSet = {}; 
      for (Int_t j = 0; j < kNEPBins; j++) {
        printf("  [%d,%d] Debug:\n",i,j);
        printf("    fNearEtaDPhiProj_Rescale array has length %d\n",(int) fNearEtaDPhiProj_Rescale[iV].size());
        if (!fNearEtaDPhiProj_Rescale[iV][i][j]) {
          printf("Missing the rescaled projection somehow!!\n");
        }
        printf("  [%d,%d] Including histogram of name %s\n",i,j,fNearEtaDPhiProj_Rescale[iV][i][j]->GetName());
        fDPhiSet.push_back(fNearEtaDPhiProj_Rescale[iV][i][j]);
      }
      fDPhiSet.push_back(fNearEtaDPhiProjAll[i]);
      DrawOmniPlots_Type(fDPhiSet,Form("Rescale_RPFMethod%d_NearEtaProj_ObsBin%d",iV,i));
    }

  }
  // Not much reason to draw the far eta projections

  cout<<"Done drawing rescaled omni plots"<<endl;
}

void TaskEventPlane::DrawOmniPlots_Type(vector<TH1D *> fHists, TString fLabel, vector<TF1 *> fFits) {

	//Bool_t bIncludeFit = (((void *) fFits )!= NULL);
	Bool_t bIncludeFit = (fFits.size() != 0);

	TString canvasName = "cRawOmni";
	if (bIncludeFit) canvasName = "cFitOmni";

	TCanvas * cRawOmni = new TCanvas(canvasName.Data(),canvasName.Data(),900,250);
	
	cRawOmni->Divide(4,1,0,0);
  // Might need to be done later, in each subpad
  cRawOmni->SetGridx(kEnableGridX);
  cRawOmni->SetGridy(kEnableGridY);

	Double_t fCommonMin=0.1;
	Double_t fCommonMax=0.9;
	FindCommonMinMax(fHists,&fCommonMin,&fCommonMax); // not applying to sum, for now

	Int_t nHists = fHists.size();
	for (Int_t j = 0; j < nHists; j++) {
		cRawOmni->cd(j+1);
		TLegend * leg = new TLegend(0.55,0.82,0.76,0.93);
 //   if ( j == 3) {
//      leg->SetHeader("Title");
 //   }
		leg->AddEntry(fHists[j],fPlaneLabels[j].c_str(),"");
		fHists[j]->GetYaxis()->SetRangeUser(fCommonMin,fCommonMax);
    fHists[j]->UseCurrentStyle();
		fHists[j]->Draw();
		if (bIncludeFit) {
			fFits[j]->Draw("SAME");
		}	

		leg->SetTextSize(0.06);
		leg->SetBorderSize(0);
////		leg->SetFillColorAlpha(10,0);

		leg->Draw("SAME");

		if (j == 0) {
			gPad->SetLeftMargin(0.14);
			fHists[j]->GetYaxis()->SetTitleOffset(1.0);
		}
		gPad->SetBottomMargin(0.18);
		if (j == kNEPBins) gPad->SetRightMargin(0.13);
	}

	TString name = "OmniPlotRaw";
	if (bIncludeFit) name = "FitOmniPlot";

	cRawOmni->Print(Form("%s/%s_%s.pdf",fOutputDir.Data(),name.Data(),fLabel.Data()));
	cRawOmni->Print(Form("%s/%s_%s.png",fOutputDir.Data(),name.Data(),fLabel.Data()));
	cRawOmni->Print(Form("%s/CFiles/%s_%s.C",fOutputDir.Data(),name.Data(),fLabel.Data()));
	cRawOmni->Clear();
	delete cRawOmni;
}

/**
 * Do a prelimary calculation of observables
 */
void TaskEventPlane::PrelimCalculation() {
  cout<<"Conducting preliminary calculations of observables"<<endl;

  fPrelimASYieldsInc_Array = {};
  fPrelimASYieldsOutOverIn_Array = {};

  fPrelimASYieldsEP_Array = {};

  fPrelimNSYieldsInc_Array = {};
  fPrelimNSYieldsOutOverIn_Array = {};

  fPrelimNSYieldsEP_Array = {};

  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    PrelimCalculation_Step(iV);
  }

}

void TaskEventPlane::PrelimCalculation_Step(Int_t iV) {

  TGraphErrors * fPrelimASYieldsInc  = new TGraphErrors(nObsBins);
  fPrelimASYieldsInc->SetName(Form("PrelimASYieldsInc_RPFMethod%d",iV));
  fPrelimASYieldsInc->SetTitle("Preliminary Yield (Away-Side)");
  fPrelimASYieldsInc->GetXaxis()->SetTitle(fObservableName.Data());
  fPrelimASYieldsInc->SetMarkerStyle(kFullSquare);
  fPrelimASYieldsInc_Array.push_back(fPrelimASYieldsInc);

  TGraphErrors * fPrelimASYieldsOutOverIn  = new TGraphErrors(nObsBins);
  fPrelimASYieldsOutOverIn->SetName(Form("PrelimASYieldOutOverIn_RPFMethod%d",iV));
  fPrelimASYieldsOutOverIn->SetTitle("Preliminary OutOverIn Ratio (Away-Side)");
  fPrelimASYieldsOutOverIn->GetXaxis()->SetTitle(fObservableName.Data());
  fPrelimASYieldsOutOverIn->SetMarkerStyle(kFullSquare);
  fPrelimASYieldsOutOverIn_Array.push_back(fPrelimASYieldsOutOverIn);

  vector<TGraphErrors *> fPrelimASYieldsEP = {};
  for (Int_t k = 0; k < kNEPBins; k++) {
    TGraphErrors * fPrelimASYields = new TGraphErrors(nObsBins);
    fPrelimASYields->SetName(Form("PrelimASYields_RPFMethod%d_EP%d",iV,k));
    fPrelimASYields->SetTitle(Form("Preliminary ASYield (%s)",fEPBinTitles[k].Data()));
    fPrelimASYields->GetXaxis()->SetTitle(fObservableName.Data());

    fPrelimASYields->SetLineColor(kEPColorList[k+1]);
    fPrelimASYields->SetMarkerColor(kEPColorList[k+1]);
    fPrelimASYields->SetMarkerStyle(kFullSquare);

    fPrelimASYieldsEP.push_back(fPrelimASYields);
  }
  fPrelimASYieldsEP_Array.push_back(fPrelimASYieldsEP);

  // Prepping Near Side
  TGraphErrors * fPrelimNSYieldsInc  = new TGraphErrors(nObsBins);
  fPrelimNSYieldsInc->SetName(Form("PrelimNSYieldsInc_RPFMethod%d",iV));
  fPrelimNSYieldsInc->SetTitle("Preliminary Yield (Near-Side)");
  fPrelimNSYieldsInc->GetXaxis()->SetTitle(fObservableName.Data());
  fPrelimNSYieldsInc->SetMarkerStyle(kFullSquare);
  fPrelimNSYieldsInc_Array.push_back(fPrelimNSYieldsInc);

  TGraphErrors * fPrelimNSYieldsOutOverIn  = new TGraphErrors(nObsBins);
  fPrelimNSYieldsOutOverIn->SetName(Form("PrelimNSYieldOutOverIn_RPFMethod%d",iV));
  fPrelimNSYieldsOutOverIn->SetTitle("Preliminary OutOverIn Ratio (Away-Side)");
  fPrelimNSYieldsOutOverIn->GetXaxis()->SetTitle(fObservableName.Data());
  fPrelimNSYieldsOutOverIn->SetMarkerStyle(kFullSquare);
  fPrelimNSYieldsOutOverIn_Array.push_back(fPrelimNSYieldsOutOverIn);

  vector<TGraphErrors *> fPrelimNSYieldsEP = {};
  for (Int_t k = 0; k < kNEPBins; k++) {
    TGraphErrors * fPrelimNSYields = new TGraphErrors(nObsBins);
    fPrelimNSYields->SetName(Form("PrelimNSYields_RPFMethod%d_EP%d",iV,k));
    fPrelimNSYields->SetTitle(Form("Preliminary NSYield (%s)",fEPBinTitles[k].Data()));
    fPrelimNSYields->GetXaxis()->SetTitle(fObservableName.Data());

    fPrelimNSYields->SetLineColor(kEPColorList[k+1]);
    fPrelimNSYields->SetMarkerColor(kEPColorList[k+1]);
    fPrelimNSYields->SetMarkerStyle(kFullSquare);

    fPrelimNSYieldsEP.push_back(fPrelimNSYields);
  }
  fPrelimNSYieldsEP_Array.push_back(fPrelimNSYieldsEP);

  double fPrelimIntegrationRange = TMath::Pi() / 3;
  double fIntMin = TMath::Pi() - fPrelimIntegrationRange;
  double fIntMax = TMath::Pi() + fPrelimIntegrationRange;
  // Away Side

  // Starting with the NearEta region
	for (Int_t i = 0; i < nObsBins; i++) {
    // General bins
    double fXValue = i + 0.5;
    double fXError = 0.5;
    if (fObservable == 1) {
      fXValue = (i+0.5)/6.;
      fXError = 0.5 /6.;
    }
    TH1D * fInclusiveHist = fNearEtaDPhiProjAll_Sub[iV][i];

//    TH1D * fInPlaneHist = fNearEtaDPhiProj_Rescale[i][0];
//    TH1D * fMidPlaneHist = fNearEtaDPhiProj_Rescale[i][1];
//    TH1D * fOutPlaneHist = fNearEtaDPhiProj_Rescale[i][2];

    double fYield_Inc_err = 0;
    double fYield_Inc = fInclusiveHist->IntegralAndError(fInclusiveHist->FindBin(fIntMin),fInclusiveHist->FindBin(fIntMax),fYield_Inc_err,"width");
    fPrelimASYieldsInc->SetPoint(i,fXValue,fYield_Inc);
    fPrelimASYieldsInc->SetPointError(i,fXError,fYield_Inc_err);
    double fYieldEPArray[3] = {0,0,0};
    double fYieldEPArray_err[3] = {0,0,0};

    for (Int_t j = 0; j < kNEPBins; j++) {
      TH1D * fHist = fNearEtaDPhiProj_Rescale[iV][i][j];
      double fYield_EP_err = 0;
      double fYield_EP = fHist->IntegralAndError(fHist->FindBin(fIntMin),fHist->FindBin(fIntMax),fYield_EP_err,"width");
      fPrelimASYieldsEP[j]->SetPoint(i,fXValue,fYield_EP);
      fPrelimASYieldsEP[j]->SetPointError(i,fXError,fYield_EP_err);
      fYieldEPArray[j] = fYield_EP;
      fYieldEPArray_err[j] = fYield_EP_err;
    }
    double fOutOverIn = 0;
    double fOutOverIn_err = 0;
    if (fYieldEPArray[0] != 0) {
      fOutOverIn = fYieldEPArray[2] / fYieldEPArray[0];
      if (fOutOverIn != 0) fOutOverIn_err = fOutOverIn * TMath::Sqrt(TMath::Power(fYieldEPArray_err[2]/fYieldEPArray[2],2) + TMath::Power(fYieldEPArray_err[0]/fYieldEPArray[0],2));
    }
    fPrelimASYieldsOutOverIn->SetPoint(i,fXValue,fOutOverIn);
    fPrelimASYieldsOutOverIn->SetPointError(i,fXError,fOutOverIn_err);
  }


  fPrelimIntegrationRange = TMath::Pi() / 3;
  fIntMin = 0. - fPrelimIntegrationRange;
  fIntMax = 0. + fPrelimIntegrationRange;

  // Near Side
	for (Int_t i = 0; i < nObsBins; i++) {
    // General bins
    double fXValue = i + 0.5;
    double fXError = 0.5;
    if (fObservable == 1) {
      fXValue = (i+0.5)/6.;
      fXError = 0.5 /6.;
    }
    TH1D * fInclusiveHist = fNearEtaDPhiProjAll_Sub[iV][i];

//    TH1D * fInPlaneHist = fNearEtaDPhiProj_Rescale[i][0];
//    TH1D * fMidPlaneHist = fNearEtaDPhiProj_Rescale[i][1];
//    TH1D * fOutPlaneHist = fNearEtaDPhiProj_Rescale[i][2];

    double fYield_Inc_err = 0;
    double fYield_Inc = fInclusiveHist->IntegralAndError(fInclusiveHist->FindBin(fIntMin),fInclusiveHist->FindBin(fIntMax),fYield_Inc_err,"width");
    fPrelimNSYieldsInc->SetPoint(i,fXValue,fYield_Inc);
    fPrelimNSYieldsInc->SetPointError(i,fXError,fYield_Inc_err);
    double fYieldEPArray[3] = {0,0,0};
    double fYieldEPArray_err[3] = {0,0,0};

    for (Int_t j = 0; j < kNEPBins; j++) {
      TH1D * fHist = fNearEtaDPhiProj_Rescale[iV][i][j];
      double fYield_EP_err = 0;
      double fYield_EP = fHist->IntegralAndError(fHist->FindBin(fIntMin),fHist->FindBin(fIntMax),fYield_EP_err,"width");
      fPrelimNSYieldsEP[j]->SetPoint(i,fXValue,fYield_EP);
      fPrelimNSYieldsEP[j]->SetPointError(i,fXError,fYield_EP_err);
      fYieldEPArray[j] = fYield_EP;
      fYieldEPArray_err[j] = fYield_EP_err;
    }
    double fOutOverIn = 0;
    double fOutOverIn_err = 0;
    if (fYieldEPArray[0] != 0) {
      fOutOverIn = fYieldEPArray[2] / fYieldEPArray[0];
      if (fOutOverIn != 0) fOutOverIn_err = fOutOverIn * TMath::Sqrt(TMath::Power(fYieldEPArray_err[2]/fYieldEPArray[2],2) + TMath::Power(fYieldEPArray_err[0]/fYieldEPArray[0],2));
    }
    fPrelimNSYieldsOutOverIn->SetPoint(i,fXValue,fOutOverIn);
    fPrelimNSYieldsOutOverIn->SetPointError(i,fXError,fOutOverIn_err);
  }



}


/**
  * Produce RPF Fit objects from the inputted parameters from the Python methods
  */ 
void TaskEventPlane::FormatPythonRPFs() {
  cout<<"Formatting Python Results.."<<endl;

  Double_t fNumTriggersInclusive  = fAllTriggerPt->Integral();
  Double_t Min = -0.5 * TMath::Pi();
  Double_t Max = 1.5 * TMath::Pi();

  Double_t fExtraScale = 1; // Rescale just for Monte Carlo Generators
  //if (fIsMCGenMode) fExtraScale = 1./fMCRescaleFactor;
  if (fIsMCGenMode) fExtraScale = 0.5/fMCRescaleFactor;
  // Does the 0.5 factor need to be applied in data as well?

  for (int k = 1; k < nRPFMethods; k++) {
    //Int_t nPar = 1 + fGlobalFit->GetNpar();
 //   Int_t nPar = fPyBkgParGraphs.size();
    // FIXME
    //Int_t nPar = 6;
    Int_t nPar = 7; // updated to include v1
    printf("Format debug: it appears nPar = %d\n",nPar);

    vector<TGraphErrors *> fParamGraphArray;

    if (k == 1) {
      fParamGraphArray = fPyBkgParGraphs;
    } else if (k == 2) {
      fParamGraphArray = fPyRPSParGraphs;
    }
    vector<TF1 *> fFitArray_Local = {};

    for (Int_t i = 0; i < nObsBins; i++) {

      RPF_Functor *fFitFunctor = new RPF_Functor();
      // Set Event Plane Resolutions
      // FIXME
      for (Int_t l = 0; l < RPF_Functor::kTotalNumberOfRn; l++) {
        fFitFunctor->SetEPRes(l,fEPRes[l]);
      }
      TString lName = Form("Obs_%d_PythonVersion%d_Fit",i,k);
      TF1 * lFit = new TF1(lName.Data(),fFitFunctor,Min,Max,nPar);


//      for (Int_t p = 0; p < nPar; p++) { 
      lFit->SetParName(0,"EventPlanePar");
      for (Int_t p = 0; p < nPar-1; p++) { 
        //lFit->SetParName(p,Form("Par%d",p)); // FIXME
        lFit->SetParName(p+1,Form("Par%d",p)); // FIXME

        double fParValue_Local = fParamGraphArray[p]->GetY()[i];
        double fParError_Local = fParamGraphArray[p]->GetEY()[i];
        printf("From Python RPFversion %d, got parameter %d = %e \\pm %e\n",k,p,fParValue_Local,fParError_Local);
        if (p == 0) { // B value is normalized differently in python versions
          // FIXME double check the 3 or pi
          //fParValue_Local = fExtraScale * fParValue_Local / (3.*fNumTriggersInclusive);
         // fParError_Local = fExtraScale * fParError_Local / (3.*fNumTriggersInclusive);
         // fParValue_Local = fExtraScale * fParValue_Local / (2.*fNumTriggersInclusive);
         // fParError_Local = fExtraScale * fParError_Local / (2.*fNumTriggersInclusive);
          fParValue_Local = fExtraScale * fParValue_Local / (fNumTriggersInclusive);
          fParError_Local = fExtraScale * fParError_Local / (fNumTriggersInclusive);
          // updating the TGraphErrors:
          fParamGraphArray[p]->SetPoint(i,fParamGraphArray[p]->GetX()[i],fParValue_Local);
          fParamGraphArray[p]->SetPointError(i,fParamGraphArray[p]->GetEX()[i],fParError_Local);
          printf("  B value renormalized to %e \\pm %e\n",fParValue_Local,fParError_Local);
        }
        //lFit->SetParameter(p,fParValue_Local);
        //lFit->SetParError(p,fParError_Local);
        lFit->SetParameter(p+1,fParValue_Local);
        lFit->SetParError(p+1,fParError_Local);
      }
      fFitArray_Local.push_back(lFit);
    }
    fRPFFits.push_back(fFitArray_Local);
  }

  cout<<"  ... Done."<<endl;
}

/**
  * Plot Comparisons of parameters from different methods and flow measurements
  */
void TaskEventPlane::CompareParameters() {
  cout<<"Comparing Parameters"<<endl;


  // Produce TGraphErrors for the parameters we can get from flow fits
  TGraphErrors * fGraphFlowV2T = new TGraphErrors(1);
  fGraphFlowV2T->SetName("GraphFlowV2T");
  TGraphErrors * fGraphFlowV4T = new TGraphErrors(1);
  fGraphFlowV4T->SetName("GraphFlowV4T");
  TGraphErrors * fGraphFlowV6T = new TGraphErrors(1);
  fGraphFlowV6T->SetName("GraphFlowV6T");

  double fObsMin = 0.;
  double fObsMax = 1.;

  if (fObservable==2) {
    fObsMin = 0.150;
    fObsMax = 18.6;
  }

  // Testing Estimating Vn of associates for z=p_t^h/p_t^\pi observable
  
  /*
   Have iPtBin 
   fAllTriggerPt 
  */

  TCanvas * cTest = new TCanvas("cTest","cTest");
  
  fAllTriggerPt->Draw();

  cTest->Print("Test.pdf");

  // FIXME next error somewhere after here  // FIXME next error somewhere after here

  if (gTrigger_V2) {
    fGraphFlowV2T->SetPoint(0,0.5*(fObsMax+fObsMin),gTrigger_V2->GetY()[iPtBin]);
    fGraphFlowV2T->SetPointError(0,0.5*(fObsMax-fObsMin),gTrigger_V2->GetEY()[iPtBin]);
  }
  if (gTrigger_V4) {
    fGraphFlowV4T->SetPoint(0,0.5*(fObsMax+fObsMin),gTrigger_V4->GetY()[iPtBin]);
    fGraphFlowV4T->SetPointError(0,0.5*(fObsMax-fObsMin),gTrigger_V4->GetEY()[iPtBin]);
  }
  if (gTrigger_V6) {
    fGraphFlowV6T->SetPoint(0,0.5*(fObsMax+fObsMin),gTrigger_V6->GetY()[iPtBin]);
    fGraphFlowV6T->SetPointError(0,0.5*(fObsMax-fObsMin),gTrigger_V6->GetEY()[iPtBin]);
  }




  SetGraphColorStyle(fGraphFlowV2T,kFlowFitColor,kFlowFitMarkerStyle);
  fGraphFlowV2T->SetFillColor(kFlowFitColor);
  fGraphFlowV2T->SetFillStyle(kFlowFillStyle);

  SetGraphColorStyle(fGraphFlowV4T,kFlowFitColor,kFlowFitMarkerStyle);
  fGraphFlowV4T->SetFillColor(kFlowFitColor);
  fGraphFlowV4T->SetFillStyle(kFlowFillStyle);

//  vector<TGraphErrors *> fParamGraphArray = {};

  TCanvas * cComparison = new TCanvas("Comparison","Comparison",900,600);
  cComparison->cd(1);

  TMultiGraph * mg1 = new TMultiGraph();
  TLegend * lCmp = new TLegend(0.55,0.75,0.95,0.95);
  lCmp->AddEntry(fChiSqGraph,"RPF1 Bkg-Only Fit","lp");
  lCmp->AddEntry(fPyBkgChiSqGraph,"RPF2 Bkg-Only Fit","lp");
  if (nRPFMethods > 2) lCmp->AddEntry(fPyRPSChiSqGraph,"RPF2 RP-Signal Fit","lp");

  SetGraphColorStyle(fChiSqGraph,kCMethodColor,kCMethodMarkerStyle);
  SetGraphColorStyle(fPyBkgChiSqGraph,kPyBkgColor,kPyBkgMarkerStyle);
  if (nRPFMethods > 2) SetGraphColorStyle(fPyRPSChiSqGraph,kPyRPDepColor,kPyRPDepMarkerStyle);

  mg1->Add(fChiSqGraph);
  mg1->Add(fPyBkgChiSqGraph);
  if (nRPFMethods > 2) mg1->Add(fPyRPSChiSqGraph);
  mg1->SetTitle("#chi/N_{dof}");

  if (fObservable == 1) mg1->GetXaxis()->SetTitle("z_{T}");
  if (fObservable == 2) mg1->GetXaxis()->SetTitle("p_{T}^{a}");

  mg1->GetYaxis()->SetTitle("#chi/N_{dof}");
  mg1->Draw("ALP");
  lCmp->Draw("SAME");

  cComparison->Print(Form("%s/Cmp_ChiSq.pdf",fOutputDir.Data()));
  cComparison->Print(Form("%s/CFiles/Cmp_ChiSq.C",fOutputDir.Data()));


  int nParams = 9; 
  int maxParamsFromPython = 6;

  int nRPF1Pars = (int) fParGraphs.size();
  int nRPF2Pars = (int) fPyBkgParGraphs.size();

  printf("  RPF1 ParGraphs: %d\n",(int)fParGraphs.size());
  printf("  RPF2  ParGraphs: %d\n",(int)fPyBkgParGraphs.size());

  // FIXME this whole part needs to be debugged;
  nParams = min(nRPF1Pars,nRPF2Pars);






  for (int i = 0; i < nParams; i++) {
 // for (int i = 0; i < (int) fPyBkgParGraphs.size(); i++) {
    printf("Drawing the comparison graphs for i = %d / %d\n",i,nParams);
//    printf("Drawing the comparison graphs for i = %d / %d\n",i,(int) fPyBkgParGraphs.size());
//    printf("    graph name: %s\n",fPyBkgParGraphs[i]->GetName());
//    mg1->Clear();
    cComparison->Clear();
    TMultiGraph * mg2 = new TMultiGraph();
    TLegend * lCmp2 = new TLegend(0.55,0.75,0.95,0.95);
    if (i < nRPF1Pars) {
      lCmp2->AddEntry(fParGraphs[i],"RPF1 Bkg-Only Fit","lp");
    }
    if (i < maxParamsFromPython) {
      lCmp2->AddEntry(fPyBkgParGraphs[i],"RPF2 Bkg-Only Fit","lp");
      if (nRPFMethods > 2) lCmp2->AddEntry(fPyRPSParGraphs[i],"RPF2 RP-Signal Fit","lp");
    }

    int kCMethodMarkerSize = 2;
    int kPyBkgMarkerSize = 2;
    int kPyRPDepMarkerStyle = 2;

    if (i < nRPF1Pars) {
      SetGraphColorStyle(fParGraphs[i],kCMethodColor,kCMethodMarkerStyle,kCMethodMarkerSize);
    }
    if (i < maxParamsFromPython) {
      SetGraphColorStyle(fPyBkgParGraphs[i],kPyBkgColor,kPyBkgMarkerStyle,kPyBkgMarkerSize);
      if (nRPFMethods > 2) SetGraphColorStyle(fPyRPSParGraphs[i],kPyRPDepColor,kPyRPDepMarkerStyle,kPyRPDepMarkerStyle);
    }
/*
    fParGraphs[i]->SetLineColor(kCMethodColor);
    fParGraphs[i]->SetMarkerColor(kCMethodColor);
    fParGraphs[i]->SetMarkerStyle(kCMethodMarkerStyle);
    fPyBkgParGraphs[i]->SetLineColor(kPyBkgColor);
    fPyBkgParGraphs[i]->SetMarkerColor(kPyBkgColor);
    fPyBkgParGraphs[i]->SetMarkerStyle(kPyBkgMarkerStyle);
    fPyRPSParGraphs[i]->SetLineColor(kPyRPDepColor);
    fPyRPSParGraphs[i]->SetMarkerColor(kPyRPDepColor);
    fPyRPSParGraphs[i]->SetMarkerStyle(kPyRPDepMarkerStyle);
*/
    TString sParName = fParNames[i];
    TString sParTitle = fParTitles[i];
    /*TString sParName = "ChiSquare";
    TString sParTitle = "#Chi^{2}/NDF";
    if (i > 0) {
      sParName  = fParNames[i-1];
      sParTitle = fParTitles[i-1];
    }*/
    printf("  Starting to add things to multigraph\n");

    mg2->SetTitle(fParGraphs[i]->GetTitle());
    mg2->GetXaxis()->SetTitle(fParGraphs[i]->GetXaxis()->GetTitle());
    mg2->GetYaxis()->SetTitle(fParGraphs[i]->GetYaxis()->GetTitle());
    // I think the first graph is chisqaure
    //if (i > 0) mg2->Add(fParGraphs[i-1]);
    mg2->Add(fParGraphs[i]);
    if (i < maxParamsFromPython) {
      mg2->Add(fPyBkgParGraphs[i]);
      if (nRPFMethods > 2) mg2->Add(fPyRPSParGraphs[i]);
    }
    mg2->Draw("ALP");

    printf("Drawing has started\n");

    double yMin = fmin(0.,mg2->GetYaxis()->GetXmin());
  //  printf("Debug: Setting mg graph min to %f\n",yMin);
    mg2->GetYaxis()->SetLimits(yMin,mg2->GetYaxis()->GetXmax());
    cComparison->Modified();

    // FIXME will have to update the indices when adding v1
    if (i==1) {
    //  fGraphFlowV2T->Draw("SAME 2");
      mg2->Add(fGraphFlowV2T);
      lCmp2->AddEntry(fGraphFlowV2T,"Flow Fit Parameter","lp");
    }
    if (i==2 && fObservable==2) {
      SetGraphColorStyle(gTrack_V2,kFlowFitColor,kFlowFitMarkerStyle);
      gTrack_V2->SetFillColor(kFlowFitColor);
      gTrack_V2->SetFillStyle(kFlowFitMarkerStyle);
      gTrack_V2->Draw("SAME P");
      lCmp2->AddEntry(gTrack_V2,"Flow Fit Parameter (This analysis)","flp");
      
      if(!fIsMCGenMode) {
        gAliTrack_V2->Draw("SAME P");
        lCmp2->AddEntry(gAliTrack_V2,"ALICE Flow Value (Interpolated)","flp");
      }
    }

    if (i==3) { //V3
      // Setting V3 plot limits to experimental limits (~0.005 should be the highest value)
      mg2->GetYaxis()->SetLimits(-0.0025,0.0075);

    }

    if (i==4) {
      //fGraphFlowV4T->Draw("SAME 2");
      mg2->Add(fGraphFlowV4T);
      lCmp2->AddEntry(fGraphFlowV4T,"Flow Fit Parameter","lp");
    }
    if (i==5 && fObservable==2) {
      SetGraphColorStyle(gTrack_V4,kFlowFitColor,kFlowFitMarkerStyle);
      gTrack_V4->SetFillColor(kFlowFitColor);
      gTrack_V4->SetFillStyle(kFlowFitMarkerStyle);
      gTrack_V4->Draw("SAME P");
      lCmp2->AddEntry(gTrack_V4,"Flow Fit Parameter","lp");

      if (!fIsMCGenMode) {
        gAliTrack_V4->Draw("SAME P");
        lCmp2->AddEntry(gAliTrack_V4,"ALICE Flow Value (Interpolated)","flp");
      }
    }
//    mg2->Draw("SAME ALP");
    if (i==7 && fObservable==2) {
      mg2->Add(fGraphFlowV6T);
      lCmp2->AddEntry(fGraphFlowV6T,"Flow Fit Parameter","lp");
    }
    if (i==8 && fObservable==2 && gTrack_V6) {
      SetGraphColorStyle(gTrack_V6,kFlowFitColor,kFlowFitMarkerStyle);
      gTrack_V6->SetFillColor(kFlowFitColor);
      gTrack_V6->SetFillStyle(kFlowFitMarkerStyle);
      mg2->Add(gTrack_V6);
  //    gTrack_V6->Draw("SAME P");
      lCmp2->AddEntry(gTrack_V6,"Flow Fit Parameter","lp");
    }


    lCmp2->Draw("SAME");

    cComparison->Print(Form("%s/Cmp_Par_%s.pdf",fOutputDir.Data(),sParName.Data()));
    cComparison->Print(Form("%s/CFiles/Cmp_Par_%s.C",fOutputDir.Data(),sParName.Data()));
  }
  cout<<"Done with comparison"<<endl;
}


/**
  * Produce the RPF functions for individual EP bins
  */
void TaskEventPlane::ProcessFitParams() {
	cout<<"Processing Fit Parameters"<<endl;
		
  // potential alternative for new loop is to add a configuration of which is the global fit to use
  // not computationally efficient for comparisons
  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    vector<vector<TF1*>> IndivFits_Version = {};
    vector<vector<TH1D *>> ResidualsIndiv_Version = {};
    for (Int_t i = 0; i < nObsBins; i++) {
      TF1 * fGlobalFit = fRPFFits[iV][i]; // Choose where to get the global fit
      Int_t nPar = 1 + fGlobalFit->GetNpar();
      Double_t Min = -0.5 * TMath::Pi();
      Double_t Max = 1.5 * TMath::Pi();
      vector<TF1 *> IndivFits = {};
      vector<TH1D *> IndivResiduals = {};
      for (Int_t j = 0; j <= kNEPBins; j++) {
        TString lName = Form("RPFMethod_%d_Obs_%d_EP_%d_Fit",iV,i,j);
        if (j == kNEPBins) lName = Form("RPFMethod_%d_Obs_%d_EP_All_Fit",iV,i);
        RPF_Functor_Single * fFitSingleFunctor = new RPF_Functor_Single();
        // Set the proper resolution
        for (Int_t l = 0; l < RPF_Functor::kTotalNumberOfRn; l++) {
          fFitSingleFunctor->SetEPRes(l,fEPRes[l]);
        }
        TF1 * lFit = new TF1(lName.Data(),fFitSingleFunctor,Min,Max,nPar);
        // why was this second one commented out?
        // for good reason
//        TF1 * lFit = new TF1(lName.Data(),RPFFunction_Single,Min,Max,nPar); // How does this compile when it is getting a class as an arguments???

        if (!fFarEtaDPhiProj[i][j]) printf("Could not Find fFarEtaDPhiProj[%d][%d]\n",i,j);
        TH1D * lResidual = 0;
        if (j < kNEPBins) { 
          lResidual = (TH1D *) fFarEtaDPhiProj[i][j]->Clone(Form("FarEtaDPhi_Res_RPFMethod%d_ObsBin%d_EP%d",iV,i,j));
        } else {
          lResidual = (TH1D *) fFarEtaDPhiProjAll[i]->Clone(Form("FarEtaDPhi_Res_RPFMethod%d_ObsBin%d_EPAll",iV,i));
        }

        // Loading the Parameters from the Fit
//        for (Int_t k = 0; k < nPar-1; k++) {
        for (Int_t k = 1; k < nPar; k++) { // k = 0 is just Event plane
          TString tParName = fGlobalFit->GetParName(k);
          double tParValue = fGlobalFit->GetParameter(k);
          double tParError = fGlobalFit->GetParError(k);
          //printf("  ProcessFitParams Debug: Setting par. %s to %f \\pm %f\n",tParName.Data(),tParValue,tParError);
          lFit->SetParName(k,tParName);
          lFit->SetParameter(k,tParValue);
          lFit->SetParError(k,tParError);

         /* lFit->SetParName(k,fGlobalFit->GetParName(k));
          lFit->SetParameter(k,fGlobalFit->GetParameter(k));
          lFit->SetParError(k,fGlobalFit->GetParError(k));*/
        }

       /* lFit->SetParName(nPar-1,"iEP");
        if (j == kNEPBins) lFit->SetParameter(nPar-1,-1);
        else lFit->SetParameter(nPar-1,j);*/
        lFit->SetParName(0,"iEP");
        if (j == kNEPBins) lFit->SetParameter(0,-1);
        else lFit->SetParameter(0,j);

        lFit->SetLineColor(kViolet-5);

        lResidual->Add(lFit,-1);
        lResidual->Divide(lFit);
        //lResidual->GetYaxis()->SetTitle("Data - Fit");
        lResidual->GetYaxis()->SetTitle("(Data - Fit)/Fit");

        IndivFits.push_back(lFit);
        IndivResiduals.push_back(lResidual);
      }
      IndivFits_Version.push_back(IndivFits);
      ResidualsIndiv_Version.push_back(IndivResiduals);
//      fRPFFits_Indiv.push_back(IndivFits);
//      fRPF_Residuals_Indiv.push_back(IndivResiduals);
    }
    fRPFFits_Indiv.push_back(IndivFits_Version);
    fRPF_Residuals_Indiv.push_back(ResidualsIndiv_Version);
  }
}

void TaskEventPlane::PlotFitParams() {
	cout<<"Drawing RPF Fits"<<endl;

	TCanvas * cFitParams = new TCanvas("cFitParams","cFitParams",900,600);
//	Int_t nParams = fRPFFits[0]->GetNpar();
	Int_t nParams = fRPFFits[0][0]->GetNpar() - 1; // skipping the event plane one
  if (iOverallMode > 0) { // Forcing the number of parameters to avoid issues later
    nParams = 6;
  }

	// Would it be better to use a histogram??
	// currently using obs bin numbers as x-axis
	// FIXME add in way to get real x-values for each observable

	//for (Int_t i = 0; i < nParams; i++) {
	for (Int_t i = 0; i < nParams; i++) {
		//TString sParName = fRPFFits[0][0]->GetParName(i);
		TString sParName = fRPFFits[0][0]->GetParName(i+1);
		TGraphErrors * fParamGraph = new TGraphErrors(nObsBins);
    fParamGraph->SetName(Form("ParGraph%d",i));
    fParamGraph->SetTitle(sParName.Data());
		fParamGraph->SetMarkerStyle(kFullSquare);
		fParamGraph->GetXaxis()->SetTitle(fObservableName.Data());
		fParamGraph->GetYaxis()->SetTitle(sParName.Data());

		for (Int_t j = 0; j < nObsBins; j++) {
			Double_t fMeanObs = (fObsBins[j+1] + fObsBins[j] )/ 2.;
			Double_t fWidthObs = (fObsBins[j+1] - fObsBins[j]) / 2;
//			fParamGraph->SetPoint(j,fMeanObs,fRPFFits[0][j]->GetParameter(i));
//			fParamGraph->SetPointError(j,fWidthObs,fRPFFits[0][j]->GetParError(i));
      if (i < fRPFFits[0][j]->GetNpar()) {
        fParamGraph->SetPoint(j,fMeanObs,fRPFFits[0][j]->GetParameter(i+1));
        fParamGraph->SetPointError(j,fWidthObs,fRPFFits[0][j]->GetParError(i+1));
      }
		}
    fParamGraph->UseCurrentStyle();
		fParamGraph->Draw("AP");
		fParGraphs.push_back(fParamGraph);
		cFitParams->Print(Form("%s/Method0_FitParam_Param%d.pdf",fOutputDir.Data(),i));
		cFitParams->Print(Form("%s/CFiles/Method0_FitParam_Param%d.C",fOutputDir.Data(),i));
//		cFitParams->Print(Form("%s/FitParam_%s.pdf",fOutputDir.Data(),sParName.Data()));
//		cFitParams->Print(Form("%s/FitParam_%s.C",fOutputDir.Data(),sParName.Data()));
 //   fOutputFile->Add(fParamGraph);//FIXME
	}

	// Drawing ChiSq / NDF
	fChiSqGraph = new TGraphErrors(nObsBins);
  fChiSqGraph->SetName("ChiSqGraph");
  fChiSqGraph->SetTitle("#chi^{2}/N_{d.o.f.}");
	fChiSqGraph->GetXaxis()->SetTitle(fObservableName.Data());
	fChiSqGraph->GetYaxis()->SetTitle("#chi^{2}/N_{d.o.f.}");
	for (Int_t j = 0; j < nObsBins; j++) {
		Double_t lChiSquare = fRPFFits[0][j]->GetChisquare();
		Double_t lNDF       = fRPFFits[0][j]->GetNDF();
		Double_t lChiSqOverNDF = -1;
		if (lNDF > 0) lChiSqOverNDF = lChiSquare / lNDF;
		Double_t fMeanObs = (fObsBins[j+1] + fObsBins[j] )/ 2.;
		Double_t fWidthObs = (fObsBins[j+1] - fObsBins[j]) / 2;
		fChiSqGraph->SetPoint(j,fMeanObs,lChiSqOverNDF);
		fChiSqGraph->SetPointError(j,fWidthObs,0);
	}
	cFitParams->Clear();
  fChiSqGraph->UseCurrentStyle();
	fChiSqGraph->SetMarkerStyle(kFullSquare);
	fChiSqGraph->Draw("AP");
//  fArrayChiSqGraph.push_back(fChiSqGraph);
	cFitParams->Print(Form("%s/Fit_ChiSqOverNDF.pdf",fOutputDir.Data()));
	cFitParams->Print(Form("%s/CFiles/Fit_ChiSqOverNDF.C",fOutputDir.Data()));
}


void TaskEventPlane::DoRPFThing() {
	cout<<"Beginning RPF Procedure"<<endl;
	
	// For loop over obsbins
	// Using far eta windows  fFarEtaDPhiProj[0,1,2]
	// Need to Construct histogram 

  Double_t fV2T_Fixed = -1;

  // FIXME I could see this causing an issue
  fRPFFits = {{}}; // initialize one array within one array

	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fDPhiSet = fFarEtaDPhiProj[i];
		TString fLabel = Form("ObsBin%d",i);

    DoRPFThing_Step(fDPhiSet,fLabel,i,fV2T_Fixed); 

    if (bFixV2T && i == 0) { // Extract the V2T from the first bin
//      fV2T_Fixed = fRPFFits[0][i]->GetParameter(1);  
      fV2T_Fixed = fRPFFits[0][i]->GetParameter(3);   // V1 Fix
      printf("Found V_2,T = %f in first Zt bin\n",fV2T_Fixed);
    }
	}

	// Preprocess Fit Parameters 
	// PreprocessFitParams();

	// print plots of the fit parameters
	PlotFitParams();	

	// Produce RPF fits for each EP bin
//	ProcessFitParams();

	// Print out plots of the fit
//	DrawFitOmniPlots();

	cout<<"Finished RPF Procedure"<<endl;
	return;
}

void TaskEventPlane::DrawFitOmniPlots() {

  // FIXME need to add loop for methods
	cout<<"Drawing DPhi Projections with Fits (Omni Plots)"<<endl;

  for (Int_t iV = 0; iV < nRPFMethods; iV++) {

    // Far Eta Projections
    for (Int_t i = 0; i < nObsBins; i++) {
      vector<TH1D *> fDPhiSet = {}; 
      for (Int_t j = 0; j < kNEPBins; j++) {
        fDPhiSet.push_back(fFarEtaDPhiProj[i][j]);
      }
      fDPhiSet.push_back(fFarEtaDPhiProjAll[i]);
  //		DrawOmniPlots_Type(fDPhiSet,Form("FarEta_ObsBin%d",i),fRPFFits_Indiv[i]);
      DrawOmniPlots_Type(fDPhiSet,Form("FarEta_RPFMethod%d_ObsBin%d",iV,i),fRPFFits_Indiv[iV][i]);
    }
    
    // Full Projections
    for (Int_t i = 0; i < nObsBins; i++) {
      vector<TH1D *> fDPhiSet = {}; 
      for (Int_t j = 0; j < kNEPBins; j++) {
        fDPhiSet.push_back(fFullDPhiProj[i][j]);
      }
      fDPhiSet.push_back(fFullDPhiProjAll[i]);
  //		DrawOmniPlots_Type(fDPhiSet,Form("Full_ObsBin%d",i),fRPFFits_Indiv[i]);
      DrawOmniPlots_Type(fDPhiSet,Form("Full_RPFMethod%d_ObsBin%d",iV,i),fRPFFits_Indiv[iV][i]);
    }

    // Near Eta
    if (fDebugLevel > 0) {
      for (Int_t i = 0; i < nObsBins; i++) {
        vector<TH1D *> fDPhiSet = {}; 
        for (Int_t j = 0; j < kNEPBins; j++) {
          fDPhiSet.push_back(fNearEtaDPhiProj[i][j]);
        }
        fDPhiSet.push_back(fNearEtaDPhiProjAll[i]);
  //			DrawOmniPlots_Type(fDPhiSet,Form("NearEta_ObsBin%d",i),fRPFFits_Indiv[i]);
        DrawOmniPlots_Type(fDPhiSet,Form("NearEta_RPFMethod%d_ObsBin%d",iV,i),fRPFFits_Indiv[iV][i]);
      }
    }
  }
}

/**
 * Does the RPF fit for the input set of three (in plane, mid plane, out of plane) histograms
 */
void TaskEventPlane::DoRPFThing_Step(vector<TH1D *> fHists, TString fLabel, Int_t iObsBin, Double_t fV2T_Fixed = -1) {


  printf("DEBUGFlow DOING RPFTHING_STEP bin %d, ObservableType = %d\n",iObsBin,fObservable);
  printf("DEBUGFlow ModeTrigger = %d  ModeV3 = %d  ModeAssoc = %d\n",iFlowTermModeTrigger,bFixV3To0,iFlowTermModeAssoc);


	Int_t nHists = fHists.size();
	
	// Merge all event planes horizontally
	TH1D * fMergedHist = MergeEvtPlanesForRPF(fHists,fLabel);
	
  // Create RPF Functor with appropriate event plane resolution set
  RPF_Functor *fFitFunctor = new RPF_Functor();
  // Set Event Plane Resolutions
  for (Int_t i = 0; i < RPF_Functor::kTotalNumberOfRn; i++) {
    fFitFunctor->SetEPRes(i,fEPRes[i]);
  }

  Int_t indexForTriggerFlow = iPtBin;
  if (fObservable == 0) {
    printf("ERROR, I didn't code for this yet\n"); // FIXME
  }

  if (!gTrigger_V2) {
    fprintf(stderr, "I'm about to crash because of missing gTrigger_V2 object. Send help.\n");

  }


  // Set the absolute maxima

  // Always allowing negative v1, for now
  fFitFunctor->SetV1Range(-fV1_AbsMax,fV1_AbsMax);

  if (bAllowNegativeVn) {
    fFitFunctor->SetV2TRange(-fV2T_AbsMax,fV2T_AbsMax);
    fFitFunctor->SetV2ARange(-fV2A_AbsMax,fV2A_AbsMax);

    fFitFunctor->SetV3Range(-fV3_AbsMax,fV3_AbsMax);

    fFitFunctor->SetV4TRange(-fV4T_AbsMax,fV4T_AbsMax);
    fFitFunctor->SetV4ARange(-fV4A_AbsMax,fV4A_AbsMax);
  } else {
    fFitFunctor->SetV2TRange(0,fV2T_AbsMax);
    fFitFunctor->SetV2ARange(0,fV2A_AbsMax);

    fFitFunctor->SetV3Range(0,fV3_AbsMax);

    fFitFunctor->SetV4TRange(0,fV4T_AbsMax);
    fFitFunctor->SetV4ARange(0,fV4A_AbsMax);


  }



  double fV2T = gTrigger_V2->GetY()[indexForTriggerFlow];
  double fV2Te= gTrigger_V2->GetEY()[indexForTriggerFlow];

  double fV4T = gTrigger_V4->GetY()[indexForTriggerFlow];
  double fV4Te= gTrigger_V4->GetEY()[indexForTriggerFlow];

  double fV2A = 0;
  double fV2Ae = 0.2;

  double fV4A = 0;
  double fV4Ae = 0.2;

  double fV3TV3A = 0;
  double fV3TV3Ae = 0.03;

  printf("fObservable = %d\n",fObservable);

  printf("Beginning FlowFinder\n");

  if (fObservable == 1) { // Zt
    // Determine pTa from pT bin and Zt (preferably with average values

  } else if (fObservable == 2) { // pTA

    //int iPTABin = fTrackPtProjectionSE->GetXaxis()->FindFixBin
    double fPtAMin = -1;
    double fPtAMax = -1;
    if (fTrackPtProjectionSE) {
      fPtAMin = fTrackPtProjectionSE->GetXaxis()->GetBinLowEdge(iObsBin+1);
      fPtAMax = fTrackPtProjectionSE->GetXaxis()->GetBinUpEdge(iObsBin+1);
    } else {
      fprintf(stderr,"DoRPFThingStep:  MISSING Track ProjectionSE\n");
    }
    // FIXME simple using middle bin
    double fPtAValue = 0.5 * (fPtAMin + fPtAMax);
    printf("DEBUGFlow evaluating V2 at pT %f  (bin %d should be [%f , %f) )\n",fPtAValue,iObsBin,fPtAMin,fPtAMax);


    // iFlowFinderMode or iFlowSource
    // Here, I get the v2,v4 of charged particles from existing graphs
    // choice 0: my measurements
    // choice 1: ALICE published measurements
    // note: the ALICE measurements might not describe triggered data

    // gTrack_V2, gTrack_V4 have my decent flow measurments
    // gAliTrack_V2,V3,V4 have ALICE official measurements (with interpolation)
    // The v4 is not the same (possibly due to ep4 vs ep2 being only partially correlated).

    // could use iFlowFinderMode to try weighted evaulations

    if (iFlowSource == 0) { //local fits


      // Could replace this with call to TaskEventPlane::GetFlowVNAFromObsBin()

      fV2A = gTrack_V2->Eval(fPtAValue);
      // get error from slope or something
      double fV2A_min = gTrack_V2->Eval(fPtAMin);
      double fV2A_max = gTrack_V2->Eval(fPtAMax);
      fV2Ae = 0.5 * TMath::Abs(fV2A_max - fV2A_min);

      fV4A = gTrack_V4->Eval(fPtAValue);
      double fV4A_min = gTrack_V4->Eval(fPtAMin);
      double fV4A_max = gTrack_V4->Eval(fPtAMax);
      fV4Ae = 0.5 * TMath::Abs(fV4A_max - fV4A_min);

    } else if (iFlowSource == 1) {

      if (fIsMCGenMode) {
        fprintf(stderr,"Error: iFlowSource 1 is being attempted in MCGen mode. This should not happen yet\n");

      }

      fV2A = gAliTrack_V2->Eval(fPtAValue);
      // get error from slope or something
      double fV2A_min = gAliTrack_V2->Eval(fPtAMin);
      double fV2A_max = gAliTrack_V2->Eval(fPtAMax);
      fV2Ae = 0.5 * TMath::Abs(fV2A_max - fV2A_min);

      fV4A = gAliTrack_V4->Eval(fPtAValue);
      double fV4A_min = gAliTrack_V4->Eval(fPtAMin);
      double fV4A_max = gAliTrack_V4->Eval(fPtAMax);
      fV4Ae = 0.5 * TMath::Abs(fV4A_max - fV4A_min);

      // Default option for V3:

      double fV3A = gAliTrack_V3->Eval(fPtAValue);
      double fV3A_min = gAliTrack_V3->Eval(fPtAMin);
      double fV3A_max = gAliTrack_V3->Eval(fPtAMax);

      // FIXME
      fV3TV3A = fV3A;
      fV3TV3Ae = TMath::Abs(fV3A_max - fV3A_min);


    }




  }

  printf("DEBUGFlow Found Measured V2T = %f +- %f\n",fV2T,fV2Te);
  printf("DEBUGFlow Found Measured V4T = %f +- %f\n",fV4T,fV4Te);

  printf("DEBUGFlow Found Measured V2A = %f +- %f\n",fV2A,fV2Ae);
  printf("DEBUGFlow Found Measured V4A = %f +- %f\n",fV4A,fV4Ae);


  printf("DEBUGFlow Estimated V3 = %f +- %f\n",fV3TV3A,fV3TV3Ae);




  // Pass the FitFunctor information from the Vn graphs
    // FIXME need to get index for trigger pt, track pt

  fFitFunctor->SetInitV1(0);

  // Set the initial values based on flow findings.
  fFitFunctor->SetInitV2T(fV2T);
  fFitFunctor->SetInitV2A(fV2A);
  fFitFunctor->SetInitV3(0.02); // Room for improvement
  fFitFunctor->SetInitV4T(fV4T);
  fFitFunctor->SetInitV4A(fV4A);

  switch (iFlowTermModeTrigger) {
    case 2:
      fFitFunctor->SetV2TRange(fV2T - fV2Te, fV2T + fV2Te);
      fFitFunctor->SetV4TRange(fV4T - fV4Te, fV4T + fV4Te);
      break;
    case 1://
      // gTrigger_V2
      fFitFunctor->SetFixedV2T(fV2T);
      fFitFunctor->SetFixedV4T(fV4T);
      break;
    default:
    case 0:
      break;
  }

  if (bFixV3To0) {
    fFitFunctor->SetFixedV3(0.);
  }


  switch (iFlowTermModeAssoc) {
    case 2:  // Range
      fFitFunctor->SetV2ARange(fV2A - fV2Ae, fV2A + fV2Te);
      fFitFunctor->SetV4ARange(fV4A - fV4Ae, fV4A + fV4Te);
      break;
    case 1:  // Fixed Value
      fFitFunctor->SetFixedV2A(fV2A);
      fFitFunctor->SetFixedV4A(fV4A);
      printf("Fixing v2a to %f and v4a to %f\n",fV2A,fV4A);
      break; 
    default:
    case 0:
      break; 
  }

  switch (iFlowV5Mode) {

    default:
    case 0:
      fFitFunctor->SetFixedV5(0);
      break;
    case 1:
      break;
  }

  switch (iFlowV6TMode) {

    default:
    case 0:
      fFitFunctor->SetFixedV6T(0);
      break;
    case 1:
      break;
  }

  switch (iFlowV6AMode) {

    default:
    case 0:
      fFitFunctor->SetFixedV6A(0);
      break;
    case 1:
      break;
  }


  fFitFunctor->DebugPrint();
	// Fit with RPF
  switch (iOverallMode) {
    case 2:
      cout<<"Using FarEtaAve: Set B to average of far delta eta regions on the nearside from the three EP."<<endl;
      break;
    case 1:
      cout<<"Using ZYAM for background. To be Coded."<<endl;
      break;
    case 0:
    default:
      cout<<"Using Reaction Plane Fit."<<endl;
      break;
  }
	TF1 * fit = FitRPF(fMergedHist,fFitFunctor,fLabel,fV2T_Fixed, iOverallMode);

  fit->SetLineColor(kFitLineColor);
	TCanvas * MergeCanvas = new TCanvas("MergeCanvas","MergeCanvas",600,300);
	MergeCanvas->cd();

  fMergedHist->UseCurrentStyle();
	fMergedHist->Draw();
	fit->Draw("SAME");

	fRPFFits[0].push_back(fit);
//	fRPFFits.push_back(fit);

	MergeCanvas->Print(Form("%s/Merge_%s.pdf",fOutputDir.Data(),fLabel.Data()));
	MergeCanvas->Print(Form("%s/CFiles/Merge_%s.C",fOutputDir.Data(),fLabel.Data()));
	delete MergeCanvas;

	return;
}

void TaskEventPlane::SubtractBackground() {
	cout<<"Beginning background subtraction"<<endl;
  
  // FIXME add different RPFs 

  // FIXME Add plots
  // Compare separate and full subtraction

  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    vector<TH1D *> fFullDPhiProjAll_Sub_PerMethod = {};
    vector<vector<TH1D *>> fFullDPhiProj_Sub_PerMethod = {};
    vector<TH1D *> fNearEtaDPhiProjAll_Sub_PerMethod = {};
    vector<vector<TH1D *>> fNearEtaDPhiProj_Sub_PerMethod = {};
    vector<TH1D *> fFarEtaDPhiProjAll_Sub_PerMethod = {};
    vector<vector<TH1D *>> fFarEtaDPhiProj_Sub_PerMethod = {};


    for (Int_t i = 0; i < nObsBins; i++) {
      // Full Projections
      TH1D * fFullDPhiProjAll_Sub_Local = (TH1D *) fFullDPhiProjAll[i]->Clone(Form("dPhi_RPFMethod%d_Full_AllEP_RPFSub_ObsBin%d",iV,i));
      fFullDPhiProjAll_Sub_Local->SetDirectory(0);
      //fFullDPhiProjAll_Sub_Local->Add(fRPFFits_Indiv[i][kNEPBins],-1);
      fFullDPhiProjAll_Sub_Local->Add(fRPFFits_Indiv[iV][i][kNEPBins],-1);

      vector<TH1D *> fFullDPhiProj_Sub_Local = {};
      for (Int_t j = 0; j < kNEPBins; j++) {
        TH1D * fFullDPhiProjEP_Sub_Local = (TH1D *) fFullDPhiProj[i][j]->Clone(Form("dPhi_RPFMethod%d_Full_EP%d_RPFSub_ObsBin%d",iV,j,i));
        fFullDPhiProjEP_Sub_Local->SetDirectory(0);
  //      fFullDPhiProjEP_Sub_Local->Add(fRPFFits_Indiv[i][j],-1);
        fFullDPhiProjEP_Sub_Local->Add(fRPFFits_Indiv[iV][i][j],-1);
        fFullDPhiProj_Sub_Local.push_back(fFullDPhiProjEP_Sub_Local);
      }

      // NearEta Projections
      TH1D * fNearEtaDPhiProjAll_Sub_Local = (TH1D *) fNearEtaDPhiProjAll[i]->Clone(Form("dPhi_RPFMethod%d_NearEta_AllEP_RPFSub_ObsBin%d",iV,i));
      fNearEtaDPhiProjAll_Sub_Local->SetDirectory(0);
      //fNearEtaDPhiProjAll_Sub_Local->Add(fRPFFits_Indiv[i][kNEPBins],-1);
      fNearEtaDPhiProjAll_Sub_Local->Add(fRPFFits_Indiv[iV][i][kNEPBins],-1);

      vector<TH1D *> fNearEtaDPhiProj_Sub_Local = {};
      for (Int_t j = 0; j < kNEPBins; j++) {
        TH1D * fNearEtaDPhiProjEP_Sub_Local = (TH1D *) fNearEtaDPhiProj[i][j]->Clone(Form("dPhi_RPFMethod%d_NearEta_EP%d_RPFSub_ObsBin%d",iV,j,i));
        fNearEtaDPhiProjEP_Sub_Local->SetDirectory(0);
  //      fNearEtaDPhiProjEP_Sub_Local->Add(fRPFFits_Indiv[i][j],-1);
        fNearEtaDPhiProjEP_Sub_Local->Add(fRPFFits_Indiv[iV][i][j],-1);
        fNearEtaDPhiProj_Sub_Local.push_back(fNearEtaDPhiProjEP_Sub_Local);
      }

      // FarEta Projections
      TH1D * fFarEtaDPhiProjAll_Sub_Local = (TH1D *) fFarEtaDPhiProjAll[i]->Clone(Form("dPhi_RPFMethod%d_FarEta_AllEP_RPFSub_ObsBin%d",iV,i));
      fFarEtaDPhiProjAll_Sub_Local->SetDirectory(0);
  //    fFarEtaDPhiProjAll_Sub_Local->Add(fRPFFits_Indiv[i][kNEPBins],-1);
      fFarEtaDPhiProjAll_Sub_Local->Add(fRPFFits_Indiv[iV][i][kNEPBins],-1);

      vector<TH1D *> fFarEtaDPhiProj_Sub_Local = {};
      for (Int_t j = 0; j < kNEPBins; j++) {
        TH1D * fFarEtaDPhiProjEP_Sub_Local = (TH1D *) fFarEtaDPhiProj[i][j]->Clone(Form("dPhi_RPFMethod%d_FarEta_EP%d_RPFSub_ObsBin%d",iV,j,i));
        fFarEtaDPhiProjEP_Sub_Local->SetDirectory(0);
  //      fFarEtaDPhiProjEP_Sub_Local->Add(fRPFFits_Indiv[i][j],-1);
        fFarEtaDPhiProjEP_Sub_Local->Add(fRPFFits_Indiv[iV][i][j],-1);
        fFarEtaDPhiProj_Sub_Local.push_back(fFarEtaDPhiProjEP_Sub_Local);
      }

      /*fFullDPhiProjAll_Sub.push_back(fFullDPhiProjAll_Sub_Local);    
      fFullDPhiProj_Sub.push_back(fFullDPhiProj_Sub_Local);
      fNearEtaDPhiProjAll_Sub.push_back(fNearEtaDPhiProjAll_Sub_Local);    
      fNearEtaDPhiProj_Sub.push_back(fNearEtaDPhiProj_Sub_Local);
      fFarEtaDPhiProjAll_Sub.push_back(fFarEtaDPhiProjAll_Sub_Local);    
      fFarEtaDPhiProj_Sub.push_back(fFarEtaDPhiProj_Sub_Local);*/

      fFullDPhiProjAll_Sub_PerMethod.push_back(fFullDPhiProjAll_Sub_Local);    
      fFullDPhiProj_Sub_PerMethod.push_back(fFullDPhiProj_Sub_Local);
      fNearEtaDPhiProjAll_Sub_PerMethod.push_back(fNearEtaDPhiProjAll_Sub_Local);    
      fNearEtaDPhiProj_Sub_PerMethod.push_back(fNearEtaDPhiProj_Sub_Local);
      fFarEtaDPhiProjAll_Sub_PerMethod.push_back(fFarEtaDPhiProjAll_Sub_Local);    
      fFarEtaDPhiProj_Sub_PerMethod.push_back(fFarEtaDPhiProj_Sub_Local);
    }
    fFullDPhiProjAll_Sub.push_back(fFullDPhiProjAll_Sub_PerMethod);
    fFullDPhiProj_Sub.push_back(fFullDPhiProj_Sub_PerMethod);
    fNearEtaDPhiProjAll_Sub.push_back(fNearEtaDPhiProjAll_Sub_PerMethod);    
    fNearEtaDPhiProj_Sub.push_back(fNearEtaDPhiProj_Sub_PerMethod);
    fFarEtaDPhiProjAll_Sub.push_back(fFarEtaDPhiProjAll_Sub_PerMethod);    
    fFarEtaDPhiProj_Sub.push_back(fFarEtaDPhiProj_Sub_PerMethod);
  }
  printf("Subtraction Finished.\n");
}

void TaskEventPlane::DrawOmniSandwichPlots() {
  cout<<"Drawing the Big Sandwich Plots"<<endl;
  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    for (Int_t i = 0; i < nObsBins; i++) {
      DrawOmniSandwichPlots_Step(iV,i);
    }
  }

/*{
		vector<TH1D *> fDPhiSet = {}; 
		for (Int_t j = 0; j < kNEPBins; j++) {
			fDPhiSet.push_back(fFullDPhiProj[i][j]);
		}
    fDPhiSet.push_back();
  }*/
}

/** Draws the fancy sandwich plots
  * iV = RPF Method, iObsBin = Observable bin
  */
void TaskEventPlane::DrawOmniSandwichPlots_Step(Int_t iV, Int_t iObsBin) {
  cout<<"Drawing the Big Sandwich Plot for bin "<<iObsBin<<endl;

  bool bEnableComponentRow = false; // FIXME disabled for now

  TCanvas * cOmniSandwich = new TCanvas("cOmniSandwich","cOmniSandwich");
  cOmniSandwich->SetGridx(kEnableGridX);
  cOmniSandwich->SetGridy(kEnableGridY);
  if (bEnableComponentRow) cOmniSandwich->Divide(kNEPBins + 1,4,0,0);
  else cOmniSandwich->Divide(kNEPBins + 1,3,0,0);
  

  TLegend * legTop = new TLegend(0.4,0.45,1.0,0.9);
 
  TLegend * legFunctions = new TLegend(0.4,0.45,1.0,0.9);

  TF1 * fZeroFunction = new TF1("ZeroFunction","0*x",-TMath::Pi()/2,3*TMath::Pi()/2);
  fZeroFunction->SetLineColor(1);

  // Drawing the top plot: signal (or full) and bkg-dom region, fit
  vector<TH1D *> fTopHistList = {fNearEtaDPhiProjAll[iObsBin],fFarEtaDPhiProjAll[iObsBin]};
  for (Int_t j = 0; j < kNEPBins; j++) {
    fTopHistList.push_back(fNearEtaDPhiProj[iObsBin][j]);
    fTopHistList.push_back(fFarEtaDPhiProj[iObsBin][j]);
  }

  printf("HELLO\n");

  Double_t fCommonMin = 0.1;
  Double_t fCommonMax = 0.9;
  FindCommonMinMax(fTopHistList,&fCommonMin,&fCommonMax);
  for (Int_t j = 0; j <= kNEPBins; j++) {
    cOmniSandwich->cd(j+1);
    TH1D * histSignal = 0;
    TH1D * histBkg    = 0;
    TF1  * RPF_Fit    = fRPFFits_Indiv[iV][iObsBin][j]; // j = kNEPBins is All
//    TF1  * RPF_Fit    = fRPFFits_Indiv[iObsBin][j]; // j = kNEPBins is All
    printf("  Using fit %s\n (title %s)\n",RPF_Fit->GetName(),RPF_Fit->GetTitle());
    if (j >= kNEPBins) { 
      histSignal = fNearEtaDPhiProjAll[iObsBin];
      histBkg    = fFarEtaDPhiProjAll[iObsBin];


      histSignal->SetTitle("All EP Angles (#times 1/3)");

    } else {
      histSignal = fNearEtaDPhiProj[iObsBin][j];
      histBkg    = fFarEtaDPhiProj[iObsBin][j];
    }
    if (!histSignal) {printf("Missing histSignal\n"); return;}
    if (!histBkg) {printf("Missing histBkg\n"); return;}

    histSignal->GetYaxis()->SetRangeUser(fCommonMin,fCommonMax);
    histBkg->GetYaxis()->SetRangeUser(fCommonMin,fCommonMax);
    //histBkg->SetMarkerColor(kRed);
    histBkg->SetMarkerSize(kOmniMarkerSize);
    histSignal->UseCurrentStyle();
    histSignal->SetLineColor(kBlue);
    histSignal->SetMarkerColor(kBlue);
    histSignal->SetMarkerStyle(kFullCircle);
    histSignal->SetMarkerSize(kOmniMarkerSize);
    histSignal->Draw();
    RPF_Fit->Draw("SAME");
    histBkg->Draw("SAME");
    histSignal->Draw("SAME");
    // FIXME Draw the title of the histSignal histogram as an TPaveBox or something
    if (j == kNEPBins) {
      legTop->AddEntry(histSignal,"Signal+Bkg","lp");
      legTop->AddEntry(histBkg,"Bkg Dominated","lp");
      legTop->AddEntry(RPF_Fit,"RPF Bkg","l");
      legTop->Draw("SAME");
    }
  }
  printf("  Finished the top row\n");
  vector<TF1 *> fEvenArray= {};
  double fComponentRowMin = 0.1;
  double fComponentRowMax = 0.8;

  vector<TH1D*> fBkgHistArray = {};
  for (Int_t j = 0; j < kNEPBins; j++) {
    fBkgHistArray.push_back(fFarEtaDPhiProj[iObsBin][j]);
  }
  fBkgHistArray.push_back(fFarEtaDPhiProjAll[iObsBin]);
  FindCommonMinMax(fBkgHistArray,&fComponentRowMin,&fComponentRowMax);

  // Creating functors just for the component plots
  // Not sure if I need more than one functor
  RPF_Functor_Single *fFitFunctorV2V4 = new RPF_Functor_Single();
  RPF_Functor_Single *fFitFunctorV3 = new RPF_Functor_Single();

  if (bEnableComponentRow) {

    for (Int_t l = 0; l < RPF_Functor::kTotalNumberOfRn; l++) {
      fFitFunctorV2V4->SetEPRes(l,fEPRes[l]);
      fFitFunctorV3->SetEPRes(l,fEPRes[l]);
    }


    printf("===  FitFunctorV2V4 Debug  ==========================\n");
    for (Int_t l = 0; l < RPF_Functor::kTotalNumberOfRn; l++) {
      printf("    EPR[%d] = %f\n",l,fFitFunctorV2V4->GetEPRes(l));
    }
    printf("===  FitFunctorV3 Debug  ==========================\n");
    for (Int_t l = 0; l < RPF_Functor::kTotalNumberOfRn; l++) {
      printf("    EPR[%d] = %f\n",l,fFitFunctorV3->GetEPRes(l));
    }



    Double_t Min = -0.5 * TMath::Pi();
    Double_t Max = 1.5 * TMath::Pi();
    //Int_t nPar = 6;
    Int_t nPar = 7; // V1 Fix

    for (Int_t j = 0; j <= kNEPBins; j++) {
      cOmniSandwich->cd(j+1+1*(kNEPBins+1));
      TH1D * histSignal = 0;
      TH1D * histBkg    = 0;
      TF1  * RPF_Fit    = fRPFFits_Indiv[iV][iObsBin][j]; // j = kNEPBins is All

      if (j >= kNEPBins) { 
        histSignal = fNearEtaDPhiProjAll[iObsBin];
        histBkg    = fFarEtaDPhiProjAll[iObsBin];
      } else {
        histSignal = fNearEtaDPhiProj[iObsBin][j];
        histBkg    = fFarEtaDPhiProj[iObsBin][j];
      }

      //double fMinBkg = histBkg->GetBinContent(histBkg->GetMinimumBin());
      //double fMaxBkg = histBkg->GetBinContent(histBkg->GetMaximumBin());
      //histBkg->GetYaxis()->SetRangeUser(fMinBkg - 0.1*(fMaxBkg - fMinBkg), fMaxBkg + 0.1 * (fMaxBkg - fMinBkg));
      
      histBkg->GetYaxis()->SetRangeUser(fComponentRowMin,fComponentRowMax);
      histBkg->Draw();
      //histBkg->Draw("AXIS"); // Draw axes only?
      // What range to set for the yaxis?

      // How can I get the v2_eff, v3_eff, v4_eff?
      // Can split into the even part, and separate v3,v1,v5
      TString fName_RPF_Evens = Form("RPF_Even_EP%d",j);
      //TF1 * fRPF_Fit_Evens = (TF1 *) RPF_Fit->Clone(Form("RPF_Even_EP%d",j));
      TF1 * fRPF_Fit_Evens = new TF1(fName_RPF_Evens.Data(),fFitFunctorV2V4,Min,Max,nPar);



      //TF1 fRPF_Fit_Evens = *RPF_Fit;
//      fRPF_Fit_Evens->Copy(*RPF_Fit);
      //(TF1 *) RPF_Fit->Copy(Form("RPF_Even_EP%d",j));
  
  //    fRPF_Fit_Evens->SetName(Form("RPF_Even_EP%d",j));
      // Copy relevant parameters
      fRPF_Fit_Evens->SetParameter(0,RPF_Fit->GetParameter(0)); // EP
      fRPF_Fit_Evens->SetParameter(1,RPF_Fit->GetParameter(1)); // B
      fRPF_Fit_Evens->SetParameter(2,0); // V1
      fRPF_Fit_Evens->SetParameter(3,RPF_Fit->GetParameter(3)); // Vt2
      fRPF_Fit_Evens->SetParameter(4,RPF_Fit->GetParameter(4)); // Va2
      fRPF_Fit_Evens->SetParameter(5,0);                        // V3
      fRPF_Fit_Evens->SetParameter(6,RPF_Fit->GetParameter(6)); // Vt4
      fRPF_Fit_Evens->SetParameter(7,RPF_Fit->GetParameter(7)); // Va4
      fRPF_Fit_Evens->SetParameter(8,0);                        // V5
      fRPF_Fit_Evens->SetParameter(9,RPF_Fit->GetParameter(9)); // Vt6
      fRPF_Fit_Evens->SetParameter(10,RPF_Fit->GetParameter(10)); // Va6



      fRPF_Fit_Evens->SetLineColor(kGreen+1);

      TString fName_RPF_V3 = Form("RPF_V3_EP%d",j);
      //TF1 * fRPF_Fit_V3 = (TF1 *) RPF_Fit->Clone(Form("RPF_V3_EP%d",j));
      TF1 * fRPF_Fit_V3 = new TF1(fName_RPF_V3.Data(),fFitFunctorV3,Min,Max,nPar);

      fRPF_Fit_V3->SetName(Form("RPF_V3_EP%d",j));
      fRPF_Fit_V3->SetParameter(0,RPF_Fit->GetParameter(0));  // EP
      fRPF_Fit_V3->SetParameter(1,RPF_Fit->GetParameter(1));  // B
      fRPF_Fit_V3->SetParameter(2,0);  // V1
      fRPF_Fit_V3->SetParameter(3,0);  // Vt2
      fRPF_Fit_V3->SetParameter(4,0); // Va2
      fRPF_Fit_V3->SetParameter(5,RPF_Fit->GetParameter(5)); // V3
      fRPF_Fit_V3->SetParameter(6,0); // Vt4
      fRPF_Fit_V3->SetParameter(7,0); // Va4
      fRPF_Fit_V3->SetParameter(8,0); // V5
      fRPF_Fit_V3->SetParameter(9,0); // Vt6
      fRPF_Fit_V3->SetParameter(10,0); //VtA
      fRPF_Fit_V3->SetLineColor(kAzure);
      fRPF_Fit_V3->SetLineStyle(3);
      fRPF_Fit_V3->SetNpx(30);

      //fRPF_Fit_Evens->Draw("SAME");
      fRPF_Fit_V3->Draw("SAME");
      fRPF_Fit_Evens->Draw("SAME");
      fEvenArray.push_back(fRPF_Fit_Evens);
      if (j==kNEPBins) {
        // Draw a legend?
        legFunctions->AddEntry(histBkg,"Background Region","lp");
        legFunctions->AddEntry(fRPF_Fit_Evens,"v2,v4","l");
        legFunctions->AddEntry(fRPF_Fit_V3,"v3","l");
        legFunctions->Draw("SAME");
      }
    }
    printf("  Finished the component row\n");
  }

  // Drawing the (data - fit) / fit
  for (Int_t j = 0; j <= kNEPBins; j++) {
    if (bEnableComponentRow) cOmniSandwich->cd(j+1+2*(kNEPBins+1));
    else cOmniSandwich->cd(j+1+(kNEPBins+1));

    fRPF_Residuals_Indiv[iV][iObsBin][j]->SetLineColor(kBlack);
    fRPF_Residuals_Indiv[iV][iObsBin][j]->SetMarkerColor(kBlack);
    fRPF_Residuals_Indiv[iV][iObsBin][j]->SetMarkerStyle(kOpenCircle);
    fRPF_Residuals_Indiv[iV][iObsBin][j]->SetMarkerSize(kOmniMarkerSize);

    fRPF_Residuals_Indiv[iV][iObsBin][j]->Draw();
    fZeroFunction->Draw("SAME");
    fRPF_Residuals_Indiv[iV][iObsBin][j]->Draw("SAME");
    fRPF_Residuals_Indiv[iV][iObsBin][j]->GetYaxis()->SetRangeUser(-0.5,0.5);
  }
  printf("  Finished the 2nd row\n");

  // Drawing the Total - RPF_Background
  // Signal Region (Near Side)
  vector<TH1D *> fBottomHistList = {fNearEtaDPhiProjAll_Sub[iV][iObsBin]};
  for (Int_t j = 0; j < kNEPBins; j++) {
    fBottomHistList.push_back(fNearEtaDPhiProj_Sub[iV][iObsBin][j]);
  }

  FindCommonMinMax(fBottomHistList,&fCommonMin,&fCommonMax);
  for (Int_t j = 0; j <= kNEPBins; j++) {
    if (bEnableComponentRow) cOmniSandwich->cd(j+1+3*(kNEPBins+1));
    else cOmniSandwich->cd(j+1+2*(kNEPBins+1));

    TH1D * histTotalMinusBkg = 0;
    if (j >= kNEPBins) {
      histTotalMinusBkg = fNearEtaDPhiProjAll_Sub[iV][iObsBin];
    } else {
      histTotalMinusBkg = fNearEtaDPhiProj_Sub[iV][iObsBin][j];
    }
    histTotalMinusBkg->GetYaxis()->SetRangeUser(fCommonMin,fCommonMax);
    histTotalMinusBkg->UseCurrentStyle();
    histTotalMinusBkg->SetLineColor(kBlue-3);
    histTotalMinusBkg->SetMarkerColor(kBlue-3);
    histTotalMinusBkg->SetMarkerStyle(kFullCircle);
    histTotalMinusBkg->SetMarkerSize(kOmniMarkerSize);
    histTotalMinusBkg->Draw();
    fZeroFunction->Draw("SAME");
    histTotalMinusBkg->Draw("SAME");
  }

  // FIXME also draw Full (NearEta + FarEta)?

  printf("Finished drawing things for this sandwich plot bin\n");

  cOmniSandwich->Print(Form("%s/OmniSandwich_RPFMethod%d_ObsBin%d.pdf",fOutputDir.Data(),iV,iObsBin));
  cOmniSandwich->Print(Form("%s/CFiles/OmniSandwich_RPFMethod%d_ObsBin%d.C",fOutputDir.Data(),iV,iObsBin));
  delete cOmniSandwich;
}

void TaskEventPlane::Rescale() {
  cout<<"Rescaling the subtracted correlations by N_{t,inclusive}/N_{t,within EP bin}"<<endl;
 
  // Starting up them arrays 
  fFullDPhiProj_Rescale = {};
  fNearEtaDPhiProj_Rescale = {};
  fFarEtaDPhiProj_Rescale = {};

  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    fFullDPhiProj_Rescale.push_back({});
    fNearEtaDPhiProj_Rescale.push_back({});
    fFarEtaDPhiProj_Rescale.push_back({});

    for (Int_t iObsBin = 0; iObsBin < nObsBins; iObsBin++) {
      for (Int_t iRegion = 0; iRegion < 3; iRegion++) {
        RescaleRegion(iV,iObsBin,iRegion);     
      }
    }
  }

  // Are the presubtracted inclusive histograms in need to rescaling ?

  // Why is this rescaling needed?

  Double_t fExtraScale = 1;
//  if (fIsMCGenMode) fExtraScale = 1./fMCRescaleFactor;


  //Also rescale the inclusive ones that were scaled down by 1/3
  if (bUnifiedNorm) {
    printf("Rescaling out that 1/3 factor in the inclusive bin\n");
    for (Int_t iV = 0; iV < nRPFMethods; iV++) {
      for (Int_t i = 0; i < nObsBins; i++) {
        fFullDPhiProjAll_Sub[iV][i]->Scale(fExtraScale*1.0*kNEPBins);
        fNearEtaDPhiProjAll_Sub[iV][i]->Scale(fExtraScale*1.0*kNEPBins);
        fFarEtaDPhiProjAll_Sub[iV][i]->Scale(fExtraScale*1.0*kNEPBins);
      }
    }
  }
}

void TaskEventPlane::RescaleRegion(Int_t iV, Int_t iObsBin, Int_t iRegion) {
  // fFullEtaDPhiProj[iObsBin][iEPBin]
  // fNearEtaDPhiProj[iObsBin][iEPBin]
  // fFarEtaDPhiProj[iObsBin][iEPBin]

  vector<vector<TH1D *>> fRegionDPhiProj_Sub;
  switch (iRegion) {
    default:
    case 0:
      fRegionDPhiProj_Sub = fFullDPhiProj_Sub[iV];
    break;
    case 1:
      fRegionDPhiProj_Sub = fNearEtaDPhiProj_Sub[iV];
    break;
    case 2:
      fRegionDPhiProj_Sub = fFarEtaDPhiProj_Sub[iV];
  }
  vector<TH1D *> fRegionDPhiProj_Rescale = {};

  Double_t fExtraScale = 1;
//  if (fIsMCGenMode) fExtraScale = 1./fMCRescaleFactor;

  for (Int_t iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    TH1D * hLocalDPhiProj_Rescale = (TH1D *) fRegionDPhiProj_Sub[iObsBin][iEPBin]->Clone(Form("%s_Rescale",fRegionDPhiProj_Sub[iObsBin][iEPBin]->GetName()));
    // Calculate the ReScale
    if (!fAllTriggerPt || !fEPBinTriggerPt[iEPBin]) {
      fprintf(stderr,"Missing a trigger distribution!!!");
      return ;
    }   
    Double_t fNumTriggersInclusive  = fAllTriggerPt->Integral();
    Double_t fNumTriggersInEvtPlane = fEPBinTriggerPt[iEPBin]->Integral();
    printf("   Getting info from histograms %s (%s) and [iEPBin=%d] %s (%s)\n",fAllTriggerPt->GetName(),fAllTriggerPt->GetTitle(),iEPBin,fEPBinTriggerPt[iEPBin]->GetName(),fEPBinTriggerPt[iEPBin]->GetTitle());
    printf("Will rescale obs bin %d, region bin %d by %f / %f\n",iObsBin,iRegion,fExtraScale * fNumTriggersInclusive,fNumTriggersInEvtPlane);
    if (fNumTriggersInEvtPlane > 0) hLocalDPhiProj_Rescale->Scale(fExtraScale * fNumTriggersInclusive / fNumTriggersInEvtPlane);
    fRegionDPhiProj_Rescale.push_back(hLocalDPhiProj_Rescale);
  }

  // Save the array of rescaled
  switch (iRegion) {
    default:
    case 0:
      fFullDPhiProj_Rescale[iV].push_back(fRegionDPhiProj_Rescale);
    break;
    case 1:
      fNearEtaDPhiProj_Rescale[iV].push_back(fRegionDPhiProj_Rescale);
    break;
    case 2:
      fFarEtaDPhiProj_Rescale[iV].push_back(fRegionDPhiProj_Rescale);
  }

}

/**
  * Returns the flow value for the track v2
  * Currently only has valid code for pTA observable
  */
double TaskEventPlane::GetFlowVNAFromObsBin(int N, int iObsBin) {
  double fVNA = 0;
  if (fObservable == 2) {
    double fPtAMin = -1;
    double fPtAMax = -1;
    if (fTrackPtProjectionSE) {
      fPtAMin = fTrackPtProjectionSE->GetXaxis()->GetBinLowEdge(iObsBin+1);
      fPtAMax = fTrackPtProjectionSE->GetXaxis()->GetBinUpEdge(iObsBin+1);

      double fPtAValue = 0.5 * (fPtAMin + fPtAMax);

      switch (N) {
        case 2:
          fVNA = gTrack_V2->Eval(fPtAValue);
          break;
        case 4:
          fVNA = gTrack_V4->Eval(fPtAValue);
          break;
        default:
          printf("Invalid N\n");

      }
    } else {
      fprintf(stderr,"FlowVNAFromObsBin: MISSING Track ProjectionSE\n");
    }
  } else {
    return 0;
    // For Observable 1, can use ptTrigger * zT to get estimate p2A
  }
  return fVNA;
}

/**
  * Returns an error for the flow value for the track vN
  * This is in a separate function to play nice with python
  * This error is just based on the size of the Observable pt Bin
  * Currently only has valid code for pTA observable.
  * For the Obs1 (zT) this uncertainty is much more complicated
  */
double TaskEventPlane::GetFlowVNAeFromObsBin(int N, int iObsBin) {
  double fVNAe = 0;
  if (fObservable == 2) {
    double fPtAMin = -1;
    double fPtAMax = -1;
    if (fTrackPtProjectionSE) {
      fPtAMin = fTrackPtProjectionSE->GetXaxis()->GetBinLowEdge(iObsBin+1);
      fPtAMax = fTrackPtProjectionSE->GetXaxis()->GetBinUpEdge(iObsBin+1);

      double fPtAValue = 0.5 * (fPtAMin + fPtAMax);

      //fV2A = gTrack_V2->Eval(fPtAValue);
      // get error from slope or something
      double fVNA_min = 0;
      double fVNA_max = 0;
      switch (N) {
        case 2:
          fVNA_min = gTrack_V2->Eval(fPtAMin);
          fVNA_max = gTrack_V2->Eval(fPtAMax);
          fVNAe = 0.5 * TMath::Abs(fVNA_max - fVNA_min);
          break;
        case 4:
          fVNA_min = gTrack_V4->Eval(fPtAMin);
          fVNA_max = gTrack_V4->Eval(fPtAMax);
          fVNAe = 0.5 * TMath::Abs(fVNA_max - fVNA_min);
          break;
        default:
          printf("Invalid N\n");
      }
    } else {
      fprintf(stderr,"GetFlowVNAe: MISSING Track ProjectionSE\n");
    }
  } else {
    return 0;
    // For Observable 1, can use ptTrigger * zT to get estimate p2A
  }
  return fVNAe;
}




void TaskEventPlane::Run_Part1() {
	cout<<"Beginning Task Event Plane"<<endl;
	if (fDebugLevel) Debug(0);

	LoadHistograms();

  // Temp FIXME
//	if (fDebugLevel) Debug(2);

	InitArrays();

  if (!fIsMCGenMode) LoadPublishedFlow();

  FitFlow();

	if (fDebugLevel) Debug(1);


  // FIXME add in option to just use ZYAM ?

	if (fSavePlots) DrawRawOmniPlots();

	// RPF Fit
	DoRPFThing();


	return;
}

void TaskEventPlane::Run_Part2() {
  cout<<"Beginning part 2 of phase 4"<<endl;

  // Debug
  printf("PyBkg graphs has length %d and PyRPS has length %d\n",(int) fPyBkgParGraphs.size(),(int) fPyRPSParGraphs.size());

  int nPyBkgGraphs = 0;
  int nPyRPSGraphs = 0;

//  for (int i = 0; i < nObsBins; i++) {
  for (int i = 0; i < 6; i++) {
    TGraphErrors * fPyBkgPar = fPyBkgParGraphs[i];
    if (fPyBkgPar == 0) fprintf(stderr,"Missing  a graph\n");
    else {
      printf("  graph has name %s and title %s\n",fPyBkgPar->GetName(),fPyBkgPar->GetTitle());
      nPyBkgGraphs++;
    }

    TGraphErrors * fPyRPSPar = fPyRPSParGraphs[i];
    if (fPyRPSPar == 0) fprintf(stderr,"Missing a RPS graph\n");
    else {
      printf("    graph has name %s and title %s\n",fPyRPSPar->GetName(),fPyRPSPar->GetTitle());
      nPyRPSGraphs++;
    }
  }
  if (nPyBkgGraphs > 0) nRPFMethods++; 
  if (nPyRPSGraphs > 0) nRPFMethods++; 

  printf("Plotting %d RPF Methods...\n",nRPFMethods);

  // Produce RPF Fits with the parameters from python methods
  FormatPythonRPFs();

  // Debug only
  printf("Debugging RPFFits\n");
  for (int i = 0; i < (int) fRPFFits.size(); i++) {
    printf("  RPF Method %d\n",i);
    for (int j = 0; j < (int) fRPFFits[i].size(); j++) {
      TF1 * TheFit = fRPFFits[i][j];
      if (!TheFit) fprintf(stderr,"Missing a fit!\n");
      else {
        printf("    Fit name = %s\n",TheFit->GetName());
//        printf("        Par 1 (V2T) = %f\n",TheFit->GetParameter(1));
        printf("        Par 1 (V2T) = %f\n",TheFit->GetParameter(2));
      }
    }
  }

  // Compare parameters from different methods
  CompareParameters();

	// Produce RPF fits for each EP bin
	ProcessFitParams();

	// Print out plots of the fit
	DrawFitOmniPlots();

	if (fDebugLevel) Debug(2);

	// Subtract Background from All Angles
	SubtractBackground();

  // Draw the 4x3 plot comparing presubtraction, fit, and postsubtraction
  DrawOmniSandwichPlots();

  // Rescale by N_{triggers in EP bin) to get the final plots, draw new OmniPlots
  // For MCGen, this is where the MC Rescale done for the Python code is undone, currently
  Rescale();
  DrawRescaleOmniPlots();

	// Calculate Yields? Or do in a separate task?
  PrelimCalculation();   
  
  SaveOutput();


  return;
}

void TaskEventPlane::SaveOutput() {
  if (!fOutputFile) {
    fprintf(stderr,"Error: No output file set!\n");
    return;
  }

  if (fEP2RGraph) fOutputFile->Add(fEP2RGraph);
  if (fEP3RGraph) fOutputFile->Add(fEP3RGraph);

  // Saving the ALICE published info used
  if (gAliTrack_V2) fOutputFile->Add(gAliTrack_V2);
  if (gAliTrack_V3) fOutputFile->Add(gAliTrack_V3);
  if (gAliTrack_V4) fOutputFile->Add(gAliTrack_V4);




  // Parameters
  printf("Saving Initial Histograms...\n");

  // Save some initial histograms?
  for (Int_t i = 0; i < nObsBins; i++) {
    fOutputFile->Add(fFullDPhiProjAll[i]);
    fOutputFile->Add(fFarEtaDPhiProjAll[i]);
    fOutputFile->Add(fNearEtaDPhiProjAll[i]);
  }
  for (Int_t i = 0; i < nObsBins; i++) {
    for (Int_t j = 0; j < kNEPBins; j++) {
      fOutputFile->Add(fFullDPhiProj[i][j]);
      fOutputFile->Add(fFarEtaDPhiProj[i][j]);
      fOutputFile->Add(fNearEtaDPhiProj[i][j]);
    }
  }

  // Fit Parameter Graphs
//  for (const auto& fChiSqGraph : fArrayChiSqGraph) {
  fOutputFile->Add(fChiSqGraph);
//  }
  for (const auto& fParamGraph : fParGraphs) {
    fOutputFile->Add(fParamGraph);
  }  

  // If the Python RPF Graphs are loaded properly, save them
  fOutputFile->Add(fPyBkgChiSqGraph);
  for (const auto& fParamGraph : fPyBkgParGraphs) {
    if (fParamGraph != 0) fOutputFile->Add(fParamGraph);
  }  
  fOutputFile->Add(fPyRPSChiSqGraph);
  for (const auto& fParamGraph : fPyRPSParGraphs) {
    if (fParamGraph != 0) fOutputFile->Add(fParamGraph);
  }  


  // Subtracted Histograms
  for (Int_t iV = 0; iV < nRPFMethods; iV++) {
    //fOutputFile->Add(fFullDPhiProjAll_Sub[0]);
    for (Int_t i = 0; i < nObsBins; i++) {
      if (fFullDPhiProjAll_Sub[iV][i]) {
        printf("Adding Histogram %s\n",fFullDPhiProjAll_Sub[iV][i]->GetName());
        fOutputFile->Add(fFullDPhiProjAll_Sub[iV][i]);
      }
      if (fFullDPhiProj_Sub[iV][i][0]) {
        for (Int_t j = 0; j < kNEPBins; j++) {
          printf("Adding Histogram %s\n",fFullDPhiProj_Sub[iV][i][j]->GetName());
          fOutputFile->Add(fFullDPhiProj_Sub[iV][i][j]);
          printf("Adding Histogram %s\n",fFullDPhiProj_Rescale[iV][i][j]->GetName());
          fOutputFile->Add(fFullDPhiProj_Rescale[iV][i][j]);
        }
      }
    }
    for (Int_t i = 0; i < nObsBins; i++) {
      if (fNearEtaDPhiProjAll_Sub[iV][i]) {
        printf("Adding Histogram %s\n",fNearEtaDPhiProjAll_Sub[iV][i]->GetName());
        fOutputFile->Add(fNearEtaDPhiProjAll_Sub[iV][i]);
      }
      if (fNearEtaDPhiProj_Sub[iV][i][0]) {
        for (Int_t j = 0; j < kNEPBins; j++) {
          printf("Adding Histogram %s\n",fNearEtaDPhiProj_Sub[iV][i][j]->GetName());
          fOutputFile->Add(fNearEtaDPhiProj_Sub[iV][i][j]);
          printf("Adding Histogram %s\n",fNearEtaDPhiProj_Rescale[iV][i][j]->GetName());
          fOutputFile->Add(fNearEtaDPhiProj_Rescale[iV][i][j]);
        }
      }
    }
    for (Int_t i = 0; i < nObsBins; i++) {
      if (fFarEtaDPhiProjAll_Sub[iV][i]) {
        printf("Adding Histogram %s\n",fFarEtaDPhiProjAll_Sub[iV][i]->GetName());
        fOutputFile->Add(fFarEtaDPhiProjAll_Sub[iV][i]);
      }
      if (fFarEtaDPhiProj_Sub[iV][i][0]) {
        for (Int_t j = 0; j < kNEPBins; j++) {
          printf("Adding Histogram %s\n",fFarEtaDPhiProj_Sub[iV][i][j]->GetName());
          fOutputFile->Add(fFarEtaDPhiProj_Sub[iV][i][j]);
          printf("Adding Histogram %s\n",fFarEtaDPhiProj_Rescale[iV][i][j]->GetName());
          fOutputFile->Add(fFarEtaDPhiProj_Rescale[iV][i][j]);
        }
      }
    }


    fOutputFile->Add(fPrelimNSYieldsInc_Array[iV]);
    fOutputFile->Add(fPrelimASYieldsInc_Array[iV]);
    for (Int_t j = 0; j < kNEPBins; j++) {
      fOutputFile->Add(fPrelimNSYieldsEP_Array[iV][j]);
      fOutputFile->Add(fPrelimASYieldsEP_Array[iV][j]);
    }
    fOutputFile->Add(fPrelimNSYieldsOutOverIn_Array[iV]);
    fOutputFile->Add(fPrelimASYieldsOutOverIn_Array[iV]);
  }

  fOutputFile->Write();
}

void TaskEventPlane::Debug(Int_t stage = 0) {
	printf("Debug : Level %d, Stage %d\n",fDebugLevel,stage);
	TCanvas * cProjectionCmp = 0;
	switch (stage) {
		case 2:
			cProjectionCmp = new TCanvas ("ProjectionCmp","ProjectionCmp",900,250);
			for (Int_t i = 0; i < nObsBins; i++) {
				cProjectionCmp->Divide(kNEPBins+1,1);
				for (Int_t j = -1; j < kNEPBins; j++) {
					TH1D *hFull, *hNear, *hFar;					
					if (j < 0) { // all event planes
						cProjectionCmp->cd(kNEPBins+1);
						hFull = fFullDPhiProjAll[i];
						hNear = fNearEtaDPhiProjAll[i];
						hFar  = fFarEtaDPhiProjAll[i];
					} else {
						cProjectionCmp->cd(j+1);
						hFull = fFullDPhiProj[i][j];
						hNear = fNearEtaDPhiProj[i][j];
						hFar  = fFarEtaDPhiProj[i][j];
					}
					// Need eta ranges to add together properly
//					TH1D * hSum = (TH1D *) hNear->Clone(Form("TempSum_%d_%d",i,j));			
//					hSum->Add(hFar);
//					hSum->Scale(0.5);
//          hFull->UseCurrentStyle();

					hFull->SetLineColor(kBlack);
					hNear->SetLineColor(kAzure);
					hFar->SetLineColor(kRed+1);

					hFull->SetMarkerColor(kBlack);
					hNear->SetMarkerColor(kAzure);
					hFar->SetMarkerColor(kRed+1);

          hFull->SetMarkerStyle(kFullSquare);
          hNear->SetMarkerStyle(kFullCircle);
          hFar->SetMarkerStyle(kFullCircle);


	//				hSum->SetMarkerColor(kGreen+1);
	//				hSum->SetMarkerStyle(kOpenSquare);
          
					hFull->Draw();
					hNear->Draw("SAME");
					hFar->Draw("SAME");
	//				hSum->Draw("SAME");
				}

        // For the AllEP, do a comparison between the sum of all 3 and the inclusive case
        TH1D * hTempSum = (TH1D *) fFullDPhiProj[i][0]->Clone(Form("TempSum_ObsBin_%d",i));
        hTempSum->Add(fFullDPhiProj[i][1]);
        hTempSum->Add(fFullDPhiProj[i][2]);
        hTempSum->Scale(1./kNEPBins);

        hTempSum->SetLineColor(kGreen+1);
        hTempSum->SetMarkerColor(kGreen+1);
        hTempSum->SetMarkerStyle(kOpenSquare);

        cProjectionCmp->cd(kNEPBins+1);
        hTempSum->Draw("SAME");

				cProjectionCmp->Print(Form("%s/Debug_Proj_Cmp_%d.pdf",fOutputDir.Data(),i));
				cProjectionCmp->Print(Form("%s/CFiles/Debug_Proj_Cmp_%d.C",fOutputDir.Data(),i));
				cProjectionCmp->Clear();
			}
			return;
		case 1:
			printf("Observable: %d (%s), with bins:\n\t",fObservable,fObservableName.Data());
			for (Int_t i = 0; i < nObsBins; i++) {
				printf("%f,",fObsBins[i]);
			}
			printf("%f\n",fObsBins[nObsBins]);
			for (Int_t i = 0; i < nObsBins; i++) {
				TH1D * hLocal = fFullDPhiProjAll[i];
				printf("Histogram %s has %.1f entries\n",hLocal->GetName(),hLocal->GetEntries());
			}
			printf("2D Array has %d entries in the first dim\n", (int) fFullDPhiProj.size());
			for (Int_t i = 0; i < nObsBins; i++) {
				printf("Looking at full projections for obs %d\n",i);	
				vector<TH1D *> fLocalVector = fFullDPhiProj[i];
				for (Int_t j = 0; j < (int) fLocalVector.size(); j++) {
					TH1D * fLocal = fLocalVector[j];
					printf("Histogram %s has %.1f entries\n",fLocal->GetName(),fLocal->GetEntries());
				}
			}
			for (Int_t i = 0; i < nObsBins; i++) {
				TH1D * hLocal = fNearEtaDPhiProjAll[i];
				printf("Histogram %s has %.1f entries\n",hLocal->GetName(),hLocal->GetEntries());
			}
			printf("2D Array has %d entries in the first dim\n", (int) fFullDPhiProj.size());
			for (Int_t i = 0; i < nObsBins; i++) {
				printf("Looking at full projections for obs %d\n",i);	
				vector<TH1D *> fLocalVector = fNearEtaDPhiProj[i];
				for (Int_t j = 0; j < (int) fLocalVector.size(); j++) {
					TH1D * fLocal = fLocalVector[j];
					printf("Histogram %s has %.1f entries\n",fLocal->GetName(),fLocal->GetEntries());
				}
			}
			for (Int_t i = 0; i < nObsBins; i++) {
				TH1D * hLocal = fFarEtaDPhiProjAll[i];
				printf("Histogram %s has %.1f entries\n",hLocal->GetName(),hLocal->GetEntries());
			}
			printf("2D Array has %d entries in the first dim\n", (int) fFullDPhiProj.size());
			for (Int_t i = 0; i < nObsBins; i++) {
				printf("Looking at full projections for obs %d\n",i);	
				vector<TH1D *> fLocalVector = fFarEtaDPhiProj[i];
				for (Int_t j = 0; j < (int) fLocalVector.size(); j++) {
					TH1D * fLocal = fLocalVector[j];
					printf("Histogram %s has %.1f entries\n",fLocal->GetName(),fLocal->GetEntries());
				}
			}

			return;
		case 0:
		default:
      printf(" MC Mode: %d\n",fIsMCGenMode);
			printf("Input File All: %p\n",fInputFileAll);
			if (fInputFileAll) printf("          Name : %s\n",fInputFileAll->GetName());
      if (!fIsMCGenMode) {
        for (Int_t i = 0; i < 3; i++) {
          printf("Input File   %d: %p\n",i, fInputFileEvt[i]);
          printf("          Name : %s\n",fInputFileEvt[i]->GetName());
        }
      }
			return;
	}
}


