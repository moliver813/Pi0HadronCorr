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



void SetGraphColorStyle(TGraph * graph, Int_t iColor, Int_t iMarkerStyle){
  graph->SetLineColor(iColor);
  graph->SetMarkerColor(iColor);
  graph->SetMarkerStyle(iMarkerStyle);
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
	if (fObservable==0)      fObservableName = "p_{T}";
	else if (fObservable==1) fObservableName = "z_{T}";
	else if (fObservable==2) fObservableName = "p_{T}^{a}";

  // Load the phase 1 v_n information
  gTrigger_Bv = (TGraphErrors *) fInputFileAll->Get("Trigger_Bv");
  gTrigger_V2 = (TGraphErrors *) fInputFileAll->Get("Trigger_V2");
  gTrigger_V4 = (TGraphErrors *) fInputFileAll->Get("Trigger_V4");
  gTrigger_V6 = (TGraphErrors *) fInputFileAll->Get("Trigger_V6");
  gTrack_Bv = (TGraphErrors *) fInputFileAll->Get("Track_Bv");
  gTrack_V2 = (TGraphErrors *) fInputFileAll->Get("Track_V2");
  gTrack_V4 = (TGraphErrors *) fInputFileAll->Get("Track_V4");
  gTrack_V6 = (TGraphErrors *) fInputFileAll->Get("Track_V6");

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
        fLocalName = Form("Proj_PtBin%d_EP%d_FullDPhi_ObsBin%d",iPtBin,j,i);
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

  Double_t fEPRes_Set_0[4][6]  = {
                    {  0.765960, 0.619163,  0.509267, 0.348666, 0.318429, 0.187868},  
                    {  0.858157, 0.822691, 0.692985, 0.580624,  0.502229,  0.375755},
                    {  0.832549,  0.771133,  0.639423,  0.507014,  0.439729,  0.305388},
                    {  0.704550,  0.445893,  0.380824,  0.196809,  0.211605,  0.084895}};

  Double_t fEPRes_Set_1[4][6]  = {
                    {  0,  0.6192508430757114,  0., 0.34878092755772117,  0.0,  0.18777865138044672},  
                    {  0, 1., 0, 0.,  0.0,  0.0},   // Don't have Cent1, Cent3
                    {  0,  0.7703651242647157,  0.,  0.5046126852106662,  0.0,  0.3020062445564112},
                    {  0,  1.,  0.,  0.,  0.0,  0.0}  };

  Double_t fEPRes_Set_2[4][6]  = {
                    {  0,  0.73,  0.62, 0.275,  0.0,  0.0},  
                    {  0, 0.885, 0.605, 0.245,  0.0,  0.0},
                    {  0,  0.85,  0.49,  0.21,  0.0,  0.0},
                    {  0,  0.58,  0.22,  0.08,  0.0,  0.0}  };

  Double_t fEPRes_Set_3[4][6]  = {
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0},  
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0},
                    {  0, 1.0, 1.0, 1.0, 1.0, 1.0}  };

  switch (iEPRSet) {
    default:
    case 0: // From QnVector
      memcpy(fEPRes,fEPRes_Set_0[iCentBin], sizeof(fEPRes));
     /* fEPRes[0] = {  0,  0.73,  0.62, 0.275,  0.0,  0.0};  
      fEPRes[1] = {  0, 0.885, 0.605, 0.245,  0.0,  0.0};
      fEPRes[2] = {  0,  0.85,  0.49,  0.21,  0.0,  0.0};
      fEPRes[3] = {  0,  0.58,  0.22,  0.08,  0.0,  0.0}; */
      break;
    case 1: // From Raymond's analysis. Cent0M with QnVector correction
      memcpy(fEPRes,fEPRes_Set_1[iCentBin], sizeof(fEPRes));
      break;
    case 2: // read off of that one graph. TPC R_n values
      memcpy(fEPRes,fEPRes_Set_2[iCentBin], sizeof(fEPRes));
      break;
    case 3: // Full resolution. Ideal for MC
      memcpy(fEPRes,fEPRes_Set_3[iCentBin], sizeof(fEPRes));
  }

  for (Int_t i = 0; i < RPF_Functor::kTotalNumberOfRn; i++) {
    printf("Loading resolution R_%d = %f\n",i+1,fEPRes[i]);
  }


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
		leg->AddEntry(fHists[j],fPlaneLabels[j].c_str(),"");
		fHists[j]->GetYaxis()->SetRangeUser(fCommonMin,fCommonMax);
    fHists[j]->UseCurrentStyle();
		fHists[j]->Draw();
		if (bIncludeFit) {
			fFits[j]->Draw("SAME");
		}	

		leg->SetTextSize(0.06);
		leg->SetBorderSize(0);
//		leg->SetFillColorAlpha(10,0);

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

  for (int k = 1; k < nRPFMethods; k++) {
    //Int_t nPar = 1 + fGlobalFit->GetNpar();
 //   Int_t nPar = fPyBkgParGraphs.size();
    // FIXME
    Int_t nPar = 6;
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
  lCmp->AddEntry(fChiSqGraph,"C++ Bkg-Only Fit","lp");
  lCmp->AddEntry(fPyBkgChiSqGraph,"Py Bkg-Only Fit","lp");
  if (nRPFMethods > 2) lCmp->AddEntry(fPyRPSChiSqGraph,"Py RP-Signal Fit","lp");

  SetGraphColorStyle(fChiSqGraph,kCMethodColor,kCMethodMarkerStyle);
  SetGraphColorStyle(fPyBkgChiSqGraph,kPyBkgColor,kPyBkgMarkerStyle);
  if (nRPFMethods > 2) SetGraphColorStyle(fPyRPSChiSqGraph,kPyRPDepColor,kPyRPDepMarkerStyle);

  mg1->Add(fChiSqGraph);
  mg1->Add(fPyBkgChiSqGraph);
  if (nRPFMethods > 2) mg1->Add(fPyRPSChiSqGraph);
  mg1->SetTitle("#chi/N_{dof}");
  mg1->GetXaxis()->SetTitle("z_{T}");
  mg1->GetYaxis()->SetTitle("#chi/N_{dof}");
  mg1->Draw("ALP");
  lCmp->Draw("SAME");

  cComparison->Print(Form("%s/Cmp_ChiSq.pdf",fOutputDir.Data()));
  cComparison->Print(Form("%s/CFiles/Cmp_ChiSq.C",fOutputDir.Data()));


  int nParams = 9; 
  int maxParamsFromPython = 6;

  printf("  C++ ParGraphs: %d\n",(int)fParGraphs.size());
  printf("  Py  ParGraphs: %d\n",(int)fPyBkgParGraphs.size());


  for (int i = 0; i < nParams; i++) {
 // for (int i = 0; i < (int) fPyBkgParGraphs.size(); i++) {
    printf("Drawing the comparison graphs for i = %d / %d\n",i,nParams);
//    printf("Drawing the comparison graphs for i = %d / %d\n",i,(int) fPyBkgParGraphs.size());
//    printf("    graph name: %s\n",fPyBkgParGraphs[i]->GetName());
//    mg1->Clear();
    cComparison->Clear();
    TMultiGraph * mg2 = new TMultiGraph();
    TLegend * lCmp2 = new TLegend(0.55,0.75,0.95,0.95);
    lCmp2->AddEntry(fParGraphs[i],"C++ Bkg-Only Fit","lp");
    if (i < maxParamsFromPython) {
      lCmp2->AddEntry(fPyBkgParGraphs[i],"Py Bkg-Only Fit","lp");
      if (nRPFMethods > 2) lCmp2->AddEntry(fPyRPSParGraphs[i],"Py RP-Signal Fit","lp");
    }

    SetGraphColorStyle(fParGraphs[i],kCMethodColor,kCMethodMarkerStyle);
    if (i < maxParamsFromPython) {
      SetGraphColorStyle(fPyBkgParGraphs[i],kPyBkgColor,kPyBkgMarkerStyle);
      if (nRPFMethods > 2) SetGraphColorStyle(fPyRPSParGraphs[i],kPyRPDepColor,kPyRPDepMarkerStyle);
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
      lCmp2->AddEntry(gTrack_V2,"Flow Fit Parameter","flp");
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
      gTrack_V6->Draw("SAME P");
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
			fParamGraph->SetPoint(j,fMeanObs,fRPFFits[0][j]->GetParameter(i+1));
			fParamGraph->SetPointError(j,fWidthObs,fRPFFits[0][j]->GetParError(i+1));
		}
    fParamGraph->UseCurrentStyle();
		fParamGraph->Draw("AP");
		fParGraphs.push_back(fParamGraph);
		cFitParams->Print(Form("%s/FitParam_Param%d.pdf",fOutputDir.Data(),i));
		cFitParams->Print(Form("%s/CFiles/FitParam_Param%d.C",fOutputDir.Data(),i));
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
      fV2T_Fixed = fRPFFits[0][i]->GetParameter(2);  
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
  // FIXME
  for (Int_t i = 0; i < RPF_Functor::kTotalNumberOfRn; i++) {
    fFitFunctor->SetEPRes(i,fEPRes[i]);
  }

  Int_t indexForTriggerFlow = iPtBin;
  if (fObservable == 0) {
    printf("ERROR, I didn't code for this yet\n"); // FIXME
  }

  double fV2T = gTrigger_V2->GetY()[indexForTriggerFlow];
  double fV2Te= gTrigger_V2->GetEY()[indexForTriggerFlow];

  double fV4T = gTrigger_V4->GetY()[indexForTriggerFlow];
  double fV4Te= gTrigger_V4->GetEY()[indexForTriggerFlow];

  double fV2A = 0;
  double fV2Ae = 0.2;

  double fV4A = 0;
  double fV4Ae = 0.2;

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
      fprintf(stderr,"MISSING Track ProjectionSE\n");
    }
    // FIXME simple using middle bin
    double fPtAValue = 0.5 * (fPtAMin + fPtAMax);
    printf("DEBUGFlow evaluating V2 at pT %f  (bin %d should be [%f , %f) )\n",fPtAValue,iObsBin,fPtAMin,fPtAMax);
    fV2A = gTrack_V2->Eval(fPtAValue);
    // get error from slope or something
    double fV2A_min = gTrack_V2->Eval(fPtAMin);
    double fV2A_max = gTrack_V2->Eval(fPtAMax);
    fV2Ae = 0.5 * TMath::Abs(fV2A_max - fV2A_min);

    fV4A = gTrack_V4->Eval(fPtAValue);
    double fV4A_min = gTrack_V4->Eval(fPtAMin);
    double fV4A_max = gTrack_V4->Eval(fPtAMax);
    fV4Ae = 0.5 * TMath::Abs(fV4A_max - fV4A_min);

  }
  printf("DEBUGFlow Found Measured V2T = %f +- %f\n",fV2T,fV2Te);
  printf("DEBUGFlow Found Measured V4T = %f +- %f\n",fV4T,fV4Te);
  printf("DEBUGFlow Found Measured V2A = %f +- %f\n",fV2A,fV2Ae);
  printf("DEBUGFlow Found Measured V4A = %f +- %f\n",fV4A,fV4Ae);

  // Pass the FitFunctor information from the Vn graphs
    // FIXME need to get index for trigger pt, track pt
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
      break; 
    default:
    case 0:
      break; 
  }

  switch (iFlowV5Mode) {

    default:
    case 0:
      fFitFunctor->SetFixedV5(0);
  }

  switch (iFlowV6TMode) {

    default:
    case 0:
      fFitFunctor->SetFixedV6T(0);
  }

  switch (iFlowV6AMode) {

    default:
    case 0:
      fFitFunctor->SetFixedV6A(0);
  }


  fFitFunctor->DebugPrint();   
	// Fit with RPF
	TF1 * fit = FitRPF(fMergedHist,fFitFunctor,fLabel,fV2T_Fixed);

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

  TCanvas * cOmniSandwich = new TCanvas("cOmniSandwich","cOmniSandwich");
  cOmniSandwich->SetGridx(kEnableGridX);
  cOmniSandwich->SetGridy(kEnableGridY);
  cOmniSandwich->Divide(kNEPBins + 1,3,0,0);
  

  TLegend * legTop = new TLegend(0.4,0.45,1.0,0.95);

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

  // Drawing the (data - fit) / fit
  for (Int_t j = 0; j <= kNEPBins; j++) {
    cOmniSandwich->cd(j+1+(kNEPBins+1));
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
    cOmniSandwich->cd(j+1+2*(kNEPBins+1));

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
  if (fIsMCGenMode) fExtraScale = 1./fMCRescaleFactor;

  for (Int_t iEPBin = 0; iEPBin < kNEPBins; iEPBin++) {
    TH1D * hLocalDPhiProj_Rescale = (TH1D *) fRegionDPhiProj_Sub[iObsBin][iEPBin]->Clone(Form("%s_Rescale",fRegionDPhiProj_Sub[iObsBin][iEPBin]->GetName()));
    // Calculate the ReScale
    if (!fAllTriggerPt || !fEPBinTriggerPt[iEPBin]) {
      fprintf(stderr,"Missing a trigger distribution!!!");
      return ;
    }   
    Double_t fNumTriggersInclusive  = fAllTriggerPt->Integral();
    Double_t fNumTriggersInEvtPlane = fEPBinTriggerPt[iEPBin]->Integral();
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

void TaskEventPlane::Run_Part1() {
	cout<<"Beginning Task Event Plane"<<endl;
	if (fDebugLevel) Debug(0);

	LoadHistograms();

	InitArrays();

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


