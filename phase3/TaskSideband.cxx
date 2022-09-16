#if defined(__CINT__)
#define _SYS_TYPES_H_
#endif

#include <TMultiGraph.h>

#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>

#include "TaskSideband.h"

//using std::vector;

using namespace std;

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::flush;
using std::ios;
/// \cond CLASSIMP
ClassImp(TaskSideband);

TaskSideband::TaskSideband():TObject(),
fObservable(-1),fTriggerName("#pi^{0}"),fObservableName(),nObsBins(0),fObsBins(),
fFullDPhiPi0(), fFullDPhiSB(), fTriggerPt(0), fTriggerPtWithinEPBin(0), fMassPtBinPi0(), fMassPtBinAll(), fMassPtBinSB(),
fFullDPhiFinal(),
fBackgroundSelection(4),fScalingFitFunction(0), iPtBin(3) {
	fDebugLevel = 1;
	
	SetStyle();
}


Int_t fSBColor[4] = {kAzure,kAzure-2,kAzure-4,kAzure-9};	
Int_t fSBStyle[4] = {kFullSquare,kFullCircle,kFullDiamond,kOpenSquare};

void TaskSideband::SetStyle() {
//	gStyle->SetCanvasColor(kBlack);
//	gStyle->SetAxisColor(0);
  TGaxis::SetMaxDigits(2);

	gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);

//  gStyle->SetPadTopMargin(0.07);//0.05
//  gStyle->SetPadBottomMargin(0.18);//0.15
////  gStyle->SetPadRightMargin(0.045);
//  gStyle->SetPadRightMargin(0.08);
//  gStyle->SetPadLeftMargin(0.21);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.07);
  //gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);


}

void TaskSideband::PrintCanvas(TCanvas * canvas, TString name) {
	if (!canvas) return;
	canvas->Print(Form("%s/%s.pdf",fOutputDir.Data(),name.Data()));
	canvas->Print(Form("%s/%s.png",fOutputDir.Data(),name.Data()));
	//canvas->Print(Form("%s/%s.eps",fOutputDir.Data(),name.Data()));
	canvas->Print(Form("%s/CFiles/%s.C",fOutputDir.Data(),name.Data()));
}

/**
 * iterate over obs bins
 * get dphi histograms for obs bin for each sideband
 * create scale fit as a function of central mass
 *  note: need to map from sideband index to mass
 *      this depeds on pt? extract mean mass from accepted histgrams
 *     produce average mass 
 * produce the candidate histogram for subtraction
 * This function also loads v_n values from the same pi0Cand file
 */

void TaskSideband::LoadPurity() {
	TCanvas * cPurity = new TCanvas("cPurity","cPurity");

	//Pi0YieldTotalRatio
	Pi0YieldTotalRatio = (TGraphErrors *) fPi0PurityFile->Get("Pi0YieldTotalRatio");

	MCPi0YieldTotalRatio = (TGraphErrors *) fPi0PurityFile->Get("MCPi0YieldTotalRatio");

  MCPhase2Pi0YieldTotalRatio = (TH1D *) fPi0CorrFile->Get("fPhase2Purity");

	if (!Pi0YieldTotalRatio) {
		fprintf(stderr,"Error: Could not find Pi0YieldTotalRatio in file %s\n",fPi0PurityFile->GetName());
		return;
	} 

  if (MCPi0YieldTotalRatio) {
    printf("Found an MC True Yield/Total Ratio from the Phase 1 file.\n");
  }

  if (MCPhase2Pi0YieldTotalRatio) {
    printf("Found an MC True Yield/Total Ratio from the Phase 2 file.\n");
  }

  //TH1D * fMCTriggerDistPi0 = 0;
  //vector<TH1D *> fMCTriggerDistSBs;
  fMCTriggerDistPi0 = (TH1D *) fPi0CorrFile->Get("fMCTriggerDist");
  if (fMCTriggerDistPi0) {
    fMCTriggerDistPi0->SetName("fMCTriggerDistPi0");
    for (Int_t j = 0; j < kNSB; j++) {
      // Is kNSB correctly set by this point?
      TH1D * fMCTriggerDistSB = (TH1D *) fSidebandFile[j]->Get("fMCTriggerDist");
      if (fMCTriggerDistSB) fMCTriggerDistSB->SetName(Form("fMCTriggerDistSB%d",j));
      else printf("Did not find an fMCTriggerDist\n");
      fMCTriggerDistSBs.push_back(fMCTriggerDistSB);
    }
    TLegend * legMCDist = new TLegend(0.5,0.6,0.85,0.85);
    // Draw a comparison of the MC True information
    fMCTriggerDistPi0->SetMarkerStyle(kFullSquare);
    fMCTriggerDistPi0->SetMarkerColor(kBlack);
    fMCTriggerDistPi0->SetLineColor(kBlack);
    fMCTriggerDistPi0->Draw();
    TAxis * fXaxisMC = fMCTriggerDistPi0->GetXaxis();
    legMCDist->AddEntry(fMCTriggerDistPi0,"#pi^{0} Candidates","lp");
    fXaxisMC->SetBinLabel(1,"Pure Background");
    fXaxisMC->SetBinLabel(2,"MC #eta#rightarrow2#gamma");
    fXaxisMC->SetBinLabel(3,"MC #pi^{0}");
    for (Int_t j = 0; j < kNSB; j++) {
      TH1D * fSBMC = fMCTriggerDistSBs[j];
      if (!fSBMC) continue;
      fSBMC->SetLineColor(fSBColor[j]);
      fSBMC->SetMarkerColor(fSBColor[j]);
      fSBMC->SetMarkerStyle(fSBStyle[j]);
      fSBMC->Draw("SAME");
      legMCDist->AddEntry(fSBMC,Form("Sideband %d",j+1),"LP");
    }
    legMCDist->Draw("SAME");
    gPad->SetLogy(1);
    cPurity->Print(Form("%s/MCDist.pdf",fOutputDir.Data()));
    cPurity->Print(Form("%s/MCDist.png",fOutputDir.Data()));
    cPurity->Print(Form("%s/CFiles/MCDist.C",fOutputDir.Data()));
    gPad->SetLogy(0);
    cPurity->Clear();
  }
	// Throw out the points we don't need (3-4,4-5), and above 17
  // Also build the Purity Array that will be used

  for (int i = 0; i < Pi0YieldTotalRatio->GetN(); i++) {
    Double_t localX = Pi0YieldTotalRatio->GetX()[i];
    if (localX < 5. || localX > 17.) { 
      printf("Removing point %d with p_T value %f\n",i,localX);
      Pi0YieldTotalRatio->RemovePoint(i);
      if (MCPi0YieldTotalRatio) MCPi0YieldTotalRatio->RemovePoint(i);
      i--;
    } else {
      double fPurity     = Pi0YieldTotalRatio->GetY()[i];
      double fPurity_Err = Pi0YieldTotalRatio->GetEY()[i];

      if (fUseMCPurity == 1) { // Load Purity from MC Phase 1
        printf("Loading Purity from MC Phase1 ...\n");
        fPurity = MCPi0YieldTotalRatio->GetY()[i];
        fPurity_Err = 0;
      } else if (fUseMCPurity == 2) { // Load Purity from MC Phase 2
        printf("Loading Purity from MC Phase2 for Pt = %f\n",MCPhase2Pi0YieldTotalRatio->GetBinCenter(i+1));
        fPurity = MCPhase2Pi0YieldTotalRatio->GetBinContent(1); // This histogram will always have just one bin with the purity
        //fPurity = MCPhase2Pi0YieldTotalRatio->GetBinContent(i+1); // This index will always be 1,2,3,4,5 -> 5-7,7-9,...
        fPurity_Err = 0;
      }

      switch (iPurityChoice) {
        case 6:
          fPurity = fPurity + 2*fPurity_Err;
          break;
        case 5:
          fPurity = fPurity - 2*fPurity_Err;
          break;
        case 4:
          fPurity = fPurity + fPurity_Err;
          break;
        case 3:
          fPurity = fPurity - fPurity_Err;
          break;
        case 2:
          fPurity = 1;
          fPurity_Err = 0;
          break;
        case 0:
          fPurity = 0;
          fPurity_Err = 0;
          break;
        default:
        case 1:
          break;
      }
      // FIXME temp
 //     fPurity = 0.550;

      // FIXME add an option for MCMode = 1 (true background)
      // Should default to purity choice 0 (0 purity, full subtraction)
      // The alternative is just taking the time to set purity = 0 when doing background only
      if (iMCMode == 2) {
        // For MC true pi0s, just let purity = 1
        // and let the rest of the code run
        fPurity = 1.0;
        fPurity_Err = 0;
      }


      printf("Storing Purity Value %d as %f \\pm %f.\n",i,fPurity,fPurity_Err);

      fPurityArray.push_back(fPurity);
      fPurityArray_Err.push_back(fPurity_Err);
    }
  }


/*	Pi0YieldTotalRatio->RemovePoint(0);
	Pi0YieldTotalRatio->RemovePoint(0);
	Pi0YieldTotalRatio->RemovePoint(7);
	Pi0YieldTotalRatio->RemovePoint(6);
*/
	Pi0YieldTotalRatio->Draw();
  TLegend * legPurity = new TLegend(0.5,0.3,0.85,0.5);
  if (MCPi0YieldTotalRatio) {
    MCPi0YieldTotalRatio->Draw("SAME LP");
    legPurity->AddEntry(Pi0YieldTotalRatio,"Reconstructed Purity","lp");
    legPurity->AddEntry(MCPi0YieldTotalRatio,"True (MC) Purity","lp");

    if (MCPhase2Pi0YieldTotalRatio) {
      MCPhase2Pi0YieldTotalRatio->SetLineColor(kAzure);
      MCPhase2Pi0YieldTotalRatio->Draw("SAME");
      legPurity->AddEntry(MCPhase2Pi0YieldTotalRatio,"True (MC) Purity (Phase 2)","lp");
    }
    legPurity->Draw("SAME");

  }
	Pi0YieldTotalRatio->GetYaxis()->SetRangeUser(0.,1.);
	
	cPurity->Print(Form("%s/Purity.pdf",fOutputDir.Data()));
	cPurity->Print(Form("%s/Purity.png",fOutputDir.Data()));
	cPurity->Print(Form("%s/CFiles/Purity.C",fOutputDir.Data()));

  // FIXME Also load the trigger/track vN from the Pi0Cand calculation
  
	gTrigger_Bv = (TGraphErrors *) fPi0PurityFile->Get("Trigger_Bv");
	gTrigger_V2 = (TGraphErrors *) fPi0PurityFile->Get("Trigger_V2");
	gTrigger_V4 = (TGraphErrors *) fPi0PurityFile->Get("Trigger_V4");
	gTrigger_V6 = (TGraphErrors *) fPi0PurityFile->Get("Trigger_V6");
	gTrack_Bv = (TGraphErrors *) fPi0PurityFile->Get("Track_Bv");
	gTrack_V2 = (TGraphErrors *) fPi0PurityFile->Get("Track_V2");
	gTrack_V4 = (TGraphErrors *) fPi0PurityFile->Get("Track_V4");
	gTrack_V6 = (TGraphErrors *) fPi0PurityFile->Get("Track_V6");

  if (!gTrigger_Bv) fprintf(stderr,"Missing trigger vn  graphs from phase1\n");
  if (!gTrack_Bv) fprintf(stderr,"Missing track vn graphs from phase1\n");


}

void TaskSideband::LoadHistograms() {

  TH1F * fRenorm = 0;
  fRenorm = (TH1F *) fPi0CorrFile->Get("Renorm");
  if (fRenorm != 0) {
    bNeedToRenormalize = true;
    printf("Pi0-hadron correlations will be renormalized\n");
    fRenorm->SetName("Renorm_Pi0");
  } else {
    printf("Pi0-hadron correlations will NOT be renormalized\n");
  }



  vector<TH1F *> fRenormSB = {};
  for (Int_t j = 0; j < kNSB; j++) {
    TH1F * fRenormSB_local = 0;
    fRenormSB_local = (TH1F *) fSidebandFile[j]->Get("Renorm");
    bNeedToRenormalizeSB.push_back((fRenormSB_local != 0));
    if ((fRenormSB_local != 0)) {
      printf("SB%d-hadron correlations will be renormalized\n",j+1);
      fRenormSB_local->SetName(Form("%s_SB%d",fRenormSB_local->GetName(),j+1));
    } else {
      printf("SB%d-hadron correlations will NOT be renormalized\n",j+1);
    }
    fRenormSB.push_back(fRenormSB_local);
  }

  fHistEventHash = (TH1F *) fPi0CorrFile->Get("HistEventHash");
 // Loading our observable settings:
  VariableInfo = (TH1D *) fPi0CorrFile->Get("VariableInfo");
  if (VariableInfo) {
    fObservable = VariableInfo->GetBinContent(1);
  } else {
    cout<<"No Variable Input TH1D found in Input File "<<fPi0CorrFile->GetName()<<endl;
    cout<<"Using default values for Observable Info and name."<<endl;
    fObservable=0; // the P_t
  }
  // This has to be set manually because of merging of the VariableInfo histogram
  fObservable=2;

  if (fObservable==0)      fObservableName = "#it{p}_{T}";
  else if (fObservable==1) fObservableName = "z_{T}";
  else if (fObservable==2) fObservableName = "#xi";

  fTriggerPt = (TH1D *) fPi0CorrFile->Get("fTriggerPt");
  fTriggerPtWithinEPBin = (TH1D *) fPi0CorrFile->Get("fTriggerPtWithinEPBin");

  if (!fTriggerPt) {
    fprintf(stderr,"Error: missing pi0 fTriggerPtBin histogram\n");
    return;
  }
  for (Int_t j = 0; j < kNSB; j++) {
    TH1D * fTriggerPtSB_local = (TH1D *) fSidebandFile[j]->Get("fTriggerPt");
    fTriggerPtSB_local->SetName(Form("%s_SB%d",fTriggerPtSB_local->GetName(),j+1));
    fTriggerPtSB.push_back(fTriggerPtSB_local);
  }

  double fRenormScalePi0 = 1;
  if (bNeedToRenormalize) {
    // Will this histogram still have the same range set in phase2?
    // Yes, as that was limited before producing this histogram as a projection
    fRenormScalePi0 = fTriggerPt->Integral();
    if (fRenormScalePi0 != 0) fRenormScalePi0 = 1./fRenormScalePi0;
  }
  vector<double> fRenormScaleSB = {};
  for (Int_t j = 0; j < kNSB; j++) {
    double fRenormScaleSB_local = 1;
    if (fRenormSB[j]) {
      fRenormScaleSB_local = fTriggerPtSB[j]->Integral();
      if (fRenormScaleSB_local != 0) fRenormScaleSB_local = 1./fRenormScaleSB_local;
    }
    fRenormScaleSB.push_back(fRenormScaleSB_local);
  }
  fTrackPtProjectionSE = (TH1D *) fPi0CorrFile->Get("TrackPtProjectionSE");
  fTrackPtProjectionME = (TH1D *) fPi0CorrFile->Get("TrackPtProjectionME");

  fTrackPtFromTrackPsi = (TH1D *) fPi0CorrFile->Get("TrackPtFromTrackPsi");



	// Loading Histograms from Pi0 Corr File (and getting nObs)

  hHistTrackPsiEPPtCent = (TH3F *) fPi0CorrFile->Get("fHistTrackPsiEPPtCent");
  if (!hHistTrackPsiEPPtCent) {
    fprintf(stderr,"Could not find hHistTrackPsiEPPtCent in Pi0-Corr file\n");
  }
  hHistTrackPsiEPPtCent->SetDirectory(0);
  // 3rd and 4th order event planes. May be missing in older files.
  hHistTrackPsiEP3PtCent = (TH3F *) fPi0CorrFile->Get("fHistTrackPsiEP3PtCent");
  if (!hHistTrackPsiEP3PtCent) {
    fprintf(stderr,"Could not find hHistTrackPsiEP3PtCent\n");
  } else {
    hHistTrackPsiEP3PtCent->SetDirectory(0);
  }
  hHistTrackPsiEP4PtCent = (TH3F *) fPi0CorrFile->Get("fHistTrackPsiEP4PtCent");
  if (!hHistTrackPsiEP4PtCent) {
    fprintf(stderr,"Could not find hHistTrackPsiEP4PtCent\n");
  } else {
    hHistTrackPsiEP4PtCent->SetDirectory(0);
  }



	for (Int_t i = 0; i < 13; i++) {
		TH1D * fLocal = 0;
		TString fLocalName = Form("dPhi_ObsBin%d_Full",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) break;
    if (bNeedToRenormalize) fLocal->Scale(fRenormScalePi0);
		fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
		fFullDPhiPi0.push_back(fLocal);		

		fLocalName = Form("dPhi_ObsBin%d_NearEta",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) break;
    if (bNeedToRenormalize) fLocal->Scale(fRenormScalePi0);
		fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
		fNearEtaDPhiPi0.push_back(fLocal);		
  
		fLocalName = Form("dPhi_ObsBin%d_FarEta",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) break;
    if (bNeedToRenormalize) fLocal->Scale(fRenormScalePi0);
		fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
		fFarEtaDPhiPi0.push_back(fLocal);		

		nObsBins++;
	}
  printf("Finished loading histograms from Pi0 file\n");
  printf("Found number of observable bins = %d\n",nObsBins);

	// Getting FullDPhi Projections from Sidebands
	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fLocalVector = {};
		for (Int_t j = 0; j < kNSB; j++) {
			TH1D * fLocal = 0;
			TString fLocalName = Form("dPhi_ObsBin%d_Full",i);
			fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);
			if (!fLocal) {
				fprintf(stderr,"Histo %s Not found for ObsBin = %d and Trigger Type = %d!!\n",fLocalName.Data(),i,j);
				return;
			}
      if (bNeedToRenormalizeSB[j]) fLocal->Scale(fRenormScaleSB[j]);
			fLocal->SetName(Form("%s_SB%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fFullDPhiSB.push_back(fLocalVector);
	}
  // Getting NearEtaDPhi Projections from Sidebands
	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fLocalVector = {};
		for (Int_t j = 0; j < kNSB; j++) {
			TH1D * fLocal = 0;
			TString fLocalName = Form("dPhi_ObsBin%d_NearEta",i);
			fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);
			if (!fLocal) {
				fprintf(stderr,"Histo %s Not found for ObsBin = %d and Trigger Type = %d!!\n",fLocalName.Data(),i,j);
				return;
			}
      if (bNeedToRenormalizeSB[j]) fLocal->Scale(fRenormScaleSB[j]);
			fLocal->SetName(Form("%s_SB%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fNearEtaDPhiSB.push_back(fLocalVector);
	}
  // Getting FarEtaDPhi Projections from Sidebands
	for (Int_t i = 0; i < nObsBins; i++) {
		vector<TH1D *> fLocalVector = {};
		for (Int_t j = 0; j < kNSB; j++) {
			TH1D * fLocal = 0;
			TString fLocalName = Form("dPhi_ObsBin%d_FarEta",i);
			fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);
			if (!fLocal) {
				fprintf(stderr,"Histo %s Not found for ObsBin = %d and Trigger Type = %d!!\n",fLocalName.Data(),i,j);
				return;
			}
      if (bNeedToRenormalizeSB[j]) fLocal->Scale(fRenormScaleSB[j]);
			fLocal->SetName(Form("%s_SB%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fFarEtaDPhiSB.push_back(fLocalVector);
	}

  // 2nd Order Event Plane
  for (int i = 0; i < kNPtBins; i++) {
    TH1D * hPtEPAnglePionAcc_Pion_Indiv = 0;
    // The old histogram that covered all Cents
    //hPtEPAnglePionAcc_Pion_Indiv = (TH1D*) fPi0CorrFile->Get(Form("PtEPAnglePionAcc_Proj_%d",i));
    hPtEPAnglePionAcc_Pion_Indiv = (TH1D*) fPi0CorrFile->Get(Form("PtEPAnglePionAccCent_Proj_%d",i));
    if (!hPtEPAnglePionAcc_Pion_Indiv) {
      printf("Could not find PtEPAnglePionAccCent_Proj_%d\n",i);
      //fprintf(stderr,"Could not find hPtEPAnglePionAcc_Pion_Indiv %d\n",i);
      // Try the centrality integrated one. Useful for T53
      hPtEPAnglePionAcc_Pion_Indiv = (TH1D *) fPi0CorrFile->Get(Form("PtEPAnglePionAcc_Proj_%d",i));
      if (!hPtEPAnglePionAcc_Pion_Indiv) {
        fprintf(stderr,"Could not find PtEPAnglePionAcc_Proj_%d either\n",i);
        break;
      } else {
        printf("Using centrality integrated flow histogram\n");
      }
    }
    hPtEPAnglePionAcc_Pion_Indiv->SetName(Form("PtEPAnglePionAcc_Pi0_Proj_%d",i));
    hPtEPAnglePionAcc_Proj_Pion.push_back(hPtEPAnglePionAcc_Pion_Indiv);

    TH1D * fLocal = 0;
    vector<TH1D*> fLocalVector = {}; 
    // Loop over sidebands
		for (Int_t j = 0; j < kNSB; j++) {
      fLocal = (TH1D *) fSidebandFile[j]->Get(Form("PtEPAnglePionAccCent_Proj_%d",i));
      if (!fLocal) {
        fprintf(stderr,"Missing an EP2 angle proj %d, trying centrality integrated hist.\n",j);
        fLocal = (TH1D *) fSidebandFile[j]->Get(Form("PtEPAnglePionAcc_Proj_%d",i));
        if (!fLocal) {
          fprintf(stderr,"Also missing trying centrality integrated hist.\n");
          break;
        }
      }
      fLocal->SetName(Form("PtEPAnglePionAcc_SB%d_Proj_%d",j,i));
      fLocal->SetTitle(Form("Sideband %d",j+1));
      fLocalVector.push_back(fLocal);
    }
    hPtEPAnglePionAcc_Proj_SB.push_back(fLocalVector);
  }


  // 3rd Order Event Plane
  for (int i = 0; i < kNPtBins; i++) {
    TH1D * hPtEP3AnglePionAcc_Pion_Indiv = 0;
    // The old histogram that covered all Cents
    hPtEP3AnglePionAcc_Pion_Indiv = (TH1D*) fPi0CorrFile->Get(Form("PtEP3AnglePionAccCent_Proj_%d",i));
    if (!hPtEP3AnglePionAcc_Pion_Indiv) {
      fprintf(stderr,"Could not find PtEP3AnglePionAccCent_Proj_%d\n",i);
      break;
    }
    hPtEP3AnglePionAcc_Pion_Indiv->SetName(Form("PtEP3AnglePionAcc_Pi0_Proj_%d",i));
    hPtEP3AnglePionAcc_Proj_Pion.push_back(hPtEP3AnglePionAcc_Pion_Indiv);

    TH1D * fLocal = 0;
    vector<TH1D*> fLocalVector = {}; 
    // Loop over sidebands
		for (Int_t j = 0; j < kNSB; j++) {
      fLocal = (TH1D *) fSidebandFile[j]->Get(Form("PtEP3AnglePionAccCent_Proj_%d",i));
      if (!fLocal) {
        fprintf(stderr,"Missing an EP3 angle proj %d\n",j);
        break;
      }
      fLocal->SetName(Form("PtEP3AnglePionAcc_SB%d_Proj_%d",j,i));
      fLocal->SetTitle(Form("Sideband %d",j+1));
      fLocalVector.push_back(fLocal);
    }
    hPtEP3AnglePionAcc_Proj_SB.push_back(fLocalVector);
  }


  // 4th Order Event Plane
  for (int i = 0; i < kNPtBins; i++) {
    TH1D * hPtEP4AnglePionAcc_Pion_Indiv = 0;
    // The old histogram that covered all Cents
    hPtEP4AnglePionAcc_Pion_Indiv = (TH1D*) fPi0CorrFile->Get(Form("PtEP4AnglePionAccCent_Proj_%d",i));
    if (!hPtEP4AnglePionAcc_Pion_Indiv) {
      fprintf(stderr,"Could not find PtEP4AnglePionAccCent_Proj_%d\n",i);
      break;
    }
    hPtEP4AnglePionAcc_Pion_Indiv->SetName(Form("PtEP4AnglePionAcc_Pi0_Proj_%d",i));
    hPtEP4AnglePionAcc_Proj_Pion.push_back(hPtEP4AnglePionAcc_Pion_Indiv);

    TH1D * fLocal = 0;
    vector<TH1D*> fLocalVector = {}; 
    // Loop over sidebands
		for (Int_t j = 0; j < kNSB; j++) {
      fLocal = (TH1D *) fSidebandFile[j]->Get(Form("PtEP4AnglePionAccCent_Proj_%d",i));
      if (!fLocal) {
        fprintf(stderr,"Missing an EP4 angle proj %d\n",j);
        break;
      }
      fLocal->SetName(Form("PtEP4AnglePionAcc_SB%d_Proj_%d",j,i));
      fLocal->SetTitle(Form("Sideband %d",j+1));
      fLocalVector.push_back(fLocal);
    }
    hPtEP4AnglePionAcc_Proj_SB.push_back(fLocalVector);
  }



  // fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);

  // Load the delta eta projections
  // Awayside-subtracted nearside

  // Pi0 Candidates
  for (Int_t i = 0; i < nObsBins; i++) {
    TH1D * fLocal = 0;
    TString fLocalName = Form("fDeta_NearSideProjSub_%d",i);
    fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
    if (!fLocal) {
      fprintf(stderr,"Pi0 Histo %s Not found for ObsBin = %d!!\n",fLocalName.Data(),i);
      break;
    }
    if (bNeedToRenormalize) fLocal->Scale(fRenormScalePi0);
    fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
    fLocal->GetYaxis()->SetTitle("1/N^{#pi^{0}} dN/d#Delta#eta");
    fNearSideSubDEtaPi0.push_back(fLocal);
  }

  // Sidebands
  for (Int_t i = 0; i < nObsBins; i++) {
    vector<TH1D *> fLocalVector = {};
    for (Int_t j = 0; j < kNSB; j++) {
      TH1D * fLocal = 0;
      TString fLocalName = Form("fDeta_NearSideProjSub_%d",i);
			fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);
      if (!fLocal) {
				fprintf(stderr,"Sideband Histo %s Not found for ObsBin = %d and Trigger Type = %d!!\n",fLocalName.Data(),i,j);
				break;
      }
      if (bNeedToRenormalizeSB[j]) fLocal->Scale(fRenormScaleSB[j]);
			fLocal->SetName(Form("%s_SB%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
    }
    fNearSideSubDEtaSB.push_back(fLocalVector);
  }


  // Raw awayside



	// Getting the Mass Distributions
	// NTF: merge in pt if fObs is not pT
	// FIXME need to do something about pt bin choice 
	//   in case doing Z_T, or Xi

  printf("Loading mass information\n");

	// Pi0 Candidate Masses
	for (Int_t i = 0; i < kGammaNBINS; i++) {
    printf("Searching for mass information for gamma bin %d\n",i);
		TH1D * fLocal = 0;
		TString fLocalName = Form("fMassPtPionAccProj_%d",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) {
      printf("Could not find an acceptance projection!!!\n");
      continue;
    }
		fMassPtBinPi0.push_back(fLocal);
		fLocalName = Form("fMassPtPionRejProj_%d",i);
    TH1D * fLocalAll = (TH1D *) fPi0CorrFile->Get(fLocalName);
    if (!fLocalAll) {
      printf("Could not find a reject\n");
      continue;
    }
    fLocalAll->Add(fLocal);
    fLocalAll->SetName(Form("fMassPt_All_%d",i));
    fMassPtBinAll.push_back(fLocalAll);
	}
	// SB Masses 
	for (Int_t i = 0; i < kGammaNBINS; i++) {
		vector<TH1D *> fLocalVector = {};
		TString fLocalName = Form("fMassPtPionAccProj_%d",i);
		for (Int_t j = 0; j < kNSB; j++) {
			TH1D * fLocal = 0;
			fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);
			if (!fLocal) continue;
			fLocalVector.push_back(fLocal);
		}	
		if (fLocalVector.size() > 0) fMassPtBinSB.push_back(fLocalVector);
	}

}

void TaskSideband::InitArrays() {
  cout<<"Initializing Arrays ..."<<endl;

  Double_t fZtStep = 1.0/(7 - 1.0);
  Double_t fXiStep = 2.5/(8 - 1.0);

  Double_t * ObsArray = 0;

  Double_t array_G_BinsValue[kGammaNBINS+1] ={5,7,9,11,14,17,20,23,30,60};
  Double_t array_ZT_BinsValue[kZtNBINS+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
  Double_t array_XI_BinsValue[kXiNBINS+1]   ={-100,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};

  const int kNTrackPtBins = 8;
//  std::vector <double> fTrackPtBins ={0.2,0.4,0.8,1.5,2.5,4,7,11,17};
  Double_t array_PTA_BinsValue[kNTrackPtBins+1] ={0.2,0.4,0.8,1.5,2.5,4,7,11,17};

  if (fObservable == 0) {
    ObsArray = array_G_BinsValue;
    fGlobalMinPt = 5;
    fGlobalMaxPt = 17;
  }
  else if (fObservable == 1) ObsArray = array_ZT_BinsValue;
  //else if (fObservable == 2) ObsArray = array_XI_BinsValue;
  else if (fObservable == 2) {
    //ObsArray = array_XI_BinsValue;
    ObsArray = array_PTA_BinsValue;
    nObsBins = kNTrackPtBins;
  }

  for (Int_t i = 0; i <= nObsBins; i++) {
    fObsBins.push_back(ObsArray[i]);
  }
  
  memcpy(fPtBins,array_G_BinsValue,kGammaNBINS+1);
  cout<<"Finished InitArrays()"<<endl;
}

void TaskSideband::DrawAlicePerf(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size)
{ 
  const char *kMonthList[12] = {"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec."};
  const char *kCentList[5] = {"0-90%","0-10%","10-30%","30-50%","50-80%"}; // index=fCent+1
  TLegend * leg  = new TLegend(x,y,x+x_size,y+y_size);
  TDatime * time = new TDatime();
  const char * month = kMonthList[time->GetMonth()-1];
  
  //leg->SetHeader(Form("ALICE Performance - %d %s %d",time->GetDay(),month,time->GetYear()));
  //leg->AddEntry(Histo,Form("ALICE Performance %d %s %d",time->GetDay()-1,month,time->GetYear()),"");
  leg->AddEntry(Histo,"ALICE Performance Not","");
  //else leg->AddEntry(Histo,Form("Work in Progress - %d %s %d",time->GetDay()-1,month,time->GetYear()),"");
 // leg->AddEntry(Histo,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 0-90%","");
  leg->AddEntry(Histo,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
 // leg->AddEntry(Histo,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, %s Centrality",kCentList[fCent+1]),"");
  if (sLabel.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel.Data()),"");
  if (sLabel2.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel2.Data()),"");
  leg->SetTextSize(0.035); // 0.04
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(10,0);

  leg->Draw("SAME");
}


void TaskSideband::DrawWIP(TH1 *Histo, Float_t x, Float_t y, Float_t x_size, Float_t y_size) {
  const char *kMonthList[12] = {"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec."};
  const char *kCentList[5] = {"0-90%","0-10%","10-30%","30-50%","50-80%"}; // index=fCent+1              
  TLegend * leg  = new TLegend(x,y,x+x_size,y+y_size);                                    
  TDatime * time = new TDatime();                                                                        
  const char * month = kMonthList[time->GetMonth()-1];

 // leg->SetHeader("ALICE");
  leg->AddEntry(Histo,Form("Work in Progress %d %s %d",time->GetDay()-1,month,time->GetYear()),"");  
  leg->AddEntry(Histo,"Pb-Pb #sqrt{s_{NN}} = 5.02 TeV","");
  if (sLabel.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel.Data()),"");
  if (sLabel2.Length() > 0) leg->AddEntry(Histo,Form("%s",sLabel2.Data()),"");
  leg->SetTextSize(0.03); //0.04
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(10,0);

  leg->Draw("SAME");
}


void TaskSideband::SetSidebandFitMask(Int_t input) {
// fSidebandFitMask[4]
  if (fSidebandMode == 0) {
    switch (input) {
      default:
      case 0:
        break;
      case 1: // Exclude 1 (3 in old label) // This should be default in use
        fNSBFit=3;
        fSidebandFitMask[0] = 0;
        break;
      case 2: // Exclude 2 (4 in old label)
        fNSBFit=3;
        fSidebandFitMask[1] = 0;
        break;
      case 3: // Exclude 3 (5 in old label)
        fNSBFit=3;
        fSidebandFitMask[2] = 0;
        break;
      case 4: // Exclude 4 (6 in old label)
        fNSBFit=3;
        fSidebandFitMask[3] = 0;
        break;
      
      case 5: // Exclude 14
        fNSBFit=2;
        fSidebandFitMask[0] = 0;
        fSidebandFitMask[3] = 0;
        break;
      case 6: // Exclude 12
        fNSBFit=2;
        fSidebandFitMask[0] = 0;
        fSidebandFitMask[1] = 0;
        break;
    }
  } else if (fSidebandMode == 1) { // 3 Sidebands in total
    switch (input) {
      default:
      case 0:
        fNSBFit=3;
        break;
    }
  }
}

// Extract Mean mass 
void TaskSideband::MassAnalysis() {
// fMassPtBinPi0, fMassPtBinSB
	cout<<"Analyzing Mass Information"<<endl;
	// Getting the Mean Mass values for the pi0 candidates
	fMeanMassPi0Val = {};
	fMeanMassPi0Val_Un = {};

	Int_t nMassValues = 6;  
	//Int_t nMassValues = fMassPtBinPi0.size();  
	if (fObservable == 0) { // PT
		nMassValues = nObsBins;
  } else {
		nMassValues = 1;
	}

	for (Int_t i = 0; i < nMassValues; i++) {
    printf("Test 1\n");
		if (fMassPtBinPi0[i]) {
      printf("%s\n",fMassPtBinPi0[i]->GetName());
			fMeanMassPi0Val.push_back(fMassPtBinPi0[i]->GetMean());
			fMeanMassPi0Val_Un.push_back(fMassPtBinPi0[i]->GetMeanError());
		} else {
			printf("Did not find fMassPtBinPi0[%d]\n",i);
		}
	}
	printf("Mean masses for Pi0 Cands extracted\n");
	printf("  nMassValues = %d\n",nMassValues);
	for (Int_t i = 0; i < nMassValues; i++) {
		vector<Double_t> fLocalVector    = {};
		vector<Double_t> fLocalVector_Un = {};
		for (Int_t j = 0; j < kNSB; j++) {
			fLocalVector.push_back(fMassPtBinSB[i][j]->GetMean());
			fLocalVector_Un.push_back(fMassPtBinSB[i][j]->GetMeanError());
		}
		fMeanMassSBVal.push_back(fLocalVector);
		fMeanMassSBVal_Un.push_back(fLocalVector_Un);
	}
	cout<<"Finished extracting mean masses"<<endl;
	
	// Could make pt vs mass
	// Make TGraphErrors, Plot them
	//TGraphErrors * 
	TCanvas * cMassIntegral = new TCanvas("MassIntegral","MassIntegral");	

//	if (fObservable == 0) { // PT
//		nMassValues = nObsBins;
//  } else {
//		nMassValues = 1;
//	}
	cout<<"Building TGraphs of integrals vs masses (in PT bins)"<<endl;;
//	for (Int_t i = 0; i < nMassValues; i++) {
	for (Int_t i = 0; i < nObsBins; i++) {
		Int_t iMassIndex = i;
		if (fObservable != 0) iMassIndex = 0;
		printf("Mass Index, Pt Bin  = %d, %d\n",i,iMassIndex);
		TGraphErrors * fMassVsIntegralPt_Full = new TGraphErrors(fNSBFit+1); // Currently doing one for Pi0 and one for each SB
		// Adding Pi0 Candidate Integral, mass		

		Double_t fLocalIntegral_Un = 0;
		Double_t fLocalIntegral = 0;
		fLocalIntegral = fFullDPhiPi0[i]->IntegralAndError(1,fFullDPhiPi0[i]->GetNbinsX(),fLocalIntegral_Un);
    // FIXME Testing just one bin
		//fLocalIntegral = fFullDPhiPi0[i]->IntegralAndError(fFullDPhiPi0[i]->FindFixBin(1.0),fFullDPhiPi0[i]->FindFixBin(1.0),fLocalIntegral_Un);
//		fLocalIntegral = fFullDPhiPi0[i]->IntegralAndError(fFullDPhiPi0[i]->FindFixBin(0.0),fFullDPhiPi0[i]->FindFixBin(0.0),fLocalIntegral_Un);



		fMassVsIntegralPt_Full->SetPoint(0,fMeanMassPi0Val[iMassIndex],fLocalIntegral);
		fMassVsIntegralPt_Full->SetPointError(0,fMeanMassPi0Val_Un[iMassIndex],fLocalIntegral_Un);

    TGraphErrors * fMassVsIntegralPt_Pi0 = new TGraphErrors(1);
		fMassVsIntegralPt_Pi0->SetPoint(0,fMeanMassPi0Val[iMassIndex],fLocalIntegral);
		fMassVsIntegralPt_Pi0->SetPointError(0,fMeanMassPi0Val_Un[iMassIndex],fLocalIntegral_Un);
    
    fMassVsIntegralPt_Pi0->SetName(Form("Pi0Point_%d",i));
    fMassVsIntegralPt_Pi0->SetTitle("#pi^{0}");
    fMassVsIntegralPt_Pi0->SetMarkerStyle(kFullStar);
    fMassVsIntegralPt_Pi0->SetMarkerColor(kRed+2);
    fMassVsIntegralPt_Pi0->SetLineColor(kRed+2);
    fMassVsIntegralPt_Pi0->SetMarkerSize(2);

  printf("Number of sidebands for fitting (fNSBFit) = %d\n",fNSBFit);

		// Adding SB Integrals, masses
    Int_t graph_counter = 0;
    double fAverageSBInt = 0; // Average value of the sideband integrals
		for (Int_t j = 0; j < fNSBFit; j++) {
      if (!fSidebandFitMask[j]) continue; 
      graph_counter++;

			Double_t fLocalIntegral_Un = 0;
			Double_t fLocalIntegral = fFullDPhiSB[i][j]->IntegralAndError(1,fFullDPhiSB[i][j]->GetNbinsX(),fLocalIntegral_Un);
      // FIXME Testing just one bin
			//Double_t fLocalIntegral = fFullDPhiSB[i][j]->IntegralAndError(fFullDPhiSB[i][j]->FindFixBin(1.0),fFullDPhiSB[i][j]->FindFixBin(1.0),fLocalIntegral_Un);
			//Double_t fLocalIntegral = fFullDPhiSB[i][j]->IntegralAndError(fFullDPhiSB[i][j]->FindFixBin(0.0),fFullDPhiSB[i][j]->FindFixBin(0.0),fLocalIntegral_Un);

			fMassVsIntegralPt_Full->SetPoint(graph_counter,fMeanMassSBVal[iMassIndex][j],fLocalIntegral);
			fMassVsIntegralPt_Full->SetPointError(graph_counter,fMeanMassSBVal_Un[iMassIndex][j],fLocalIntegral_Un);

			//fMassVsIntegralPt_Full->SetPoint(j+1,fMeanMassSBVal[iMassIndex][j],fLocalIntegral);
			//fMassVsIntegralPt_Full->SetPointError(j+1,fMeanMassSBVal_Un[iMassIndex][j],fLocalIntegral_Un);
      fAverageSBInt += fLocalIntegral;
		}
    fAverageSBInt = fAverageSBInt / graph_counter;

		fMassVsIntegralPt_Full->SetName(Form("MassVsIntegral_ObsBin_%d",i));
		fMassVsIntegralPt_Full->SetTitle(Form("Mass Vs Corr. Integral Obs. Bin %d",i));
		fMassVsIntegralPt_Full->GetXaxis()->SetTitle("<m> (GeV/c^{2})");
		fMassVsIntegralPt_Full->GetYaxis()->SetTitle("#frac{1}{2 #pi N_{trig}} #frac{dN^{assoc}}{d#Delta#eta}");
//		fMassVsIntegralPt_Full->GetYaxis()->SetTitle("#int C(#Delta#phi)");

		// Fitting.  Could move this elsewhere
		//   N2S: Could define fit function so par0 is precisely the value at the pi0 peak location for this bin
		//TF1 * fMassEffectFit_l = new TF1(Form("MassEffectFit_%d",i),"[0]+[1]*x",0.1,0.5);

		TF1 * fMassEffectFit_l = 0;
    
    switch(fScalingFitFunction) {
      case 0:
        printf("Fitting the mass scaling with a constant fit\n");
        fMassEffectFit_l = new TF1(Form("MassEffectFit_%d",i),"[0]",0.1,0.5);
        break;
      default:
      case 1:
        printf("Fitting the mass scaling with a linear fit\n");
        fMassEffectFit_l = new TF1(Form("MassEffectFit_%d",i),"[0]+[1]*(x - [2])",0.1,0.5);
        fMassEffectFit_l->SetParameter(1, 0.);
        fMassEffectFit_l->FixParameter(2,fMeanMassPi0Val[iMassIndex]);
        break;
      case 2:
        printf("Fitting the mass scaling with a quadratic fit\n");
        // General quadratic, written so that f([2]) is always = [0]
        fMassEffectFit_l = new TF1(Form("MassEffectFit_%d",i),"[0]+[1]*(x - [2])*(x - [3])",0.1,0.5);
        fMassEffectFit_l->SetParameter(1, 0.);
        fMassEffectFit_l->SetParameter(3, 0.);
        fMassEffectFit_l->FixParameter(2,fMeanMassPi0Val[iMassIndex]);
		  break;
    }
    // Guess value for parameter 0
    fMassEffectFit_l->SetParameter(0, fAverageSBInt);


		fMassEffectFit_l->SetLineColor(kViolet);
		fMassVsIntegralPt_Full->Fit(fMassEffectFit_l,"","",fMassVsIntegralPt_Full->GetX()[0]+0.03,0.5);
//		fMassVsIntegralPt_Full->Fit(fMassEffectFit_l);
		
		fMassVsIntegral_Full.push_back(fMassVsIntegralPt_Full);
		fMassEffectFit.push_back(fMassEffectFit_l);
		
		fMassVsIntegralPt_Full->SetMarkerStyle(kOpenSquare);
		fMassVsIntegralPt_Full->Draw("AP");
//    fMassVsIntegralPt_Full->GetYaxis()->SetLabelSize(0.034);
//    fMassVsIntegralPt_Full->GetYaxis()->SetTitleOffset(0.9);

    fMassVsIntegralPt_Pi0->Draw("SAME P");
//    fMassVsIntegralPt_Pi0->Draw("SAME AP");
    fMassVsIntegralPt_Full->RemovePoint(0);


    TLegend * leg = new TLegend(0.26,0.18,0.47,0.39); //  0.25,0.21,0.52,0.38
    leg->AddEntry(fMassVsIntegralPt_Pi0,"#pi^{0}_{Cand}","p");
    leg->AddEntry(fMassVsIntegralPt_Full,"Sidebands","p");

		// Moving location of fit box
		gPad->Update();
		TPaveStats *st = (TPaveStats*)fMassVsIntegralPt_Full->FindObject("stats");
		st->SetTextColor(kViolet);
		st->SetBorderSize(0);
		st->SetFillStyle(0);
		st->SetX1NDC(0.60);
		st->SetX2NDC(0.85);
		st->SetY1NDC(0.2);
		st->SetY2NDC(0.55);
		fMassEffectFit_l->Draw("SAME");
    leg->Draw("SAME");

		PrintCanvas(cMassIntegral,TString::Format("MassVsIntegral_Index_%d",i));
		//PrintCanvas(cMassIntegral,TString::Format("MassVsIntegral_Pt_%d",i));
		cMassIntegral->Clear();
	}
	delete cMassIntegral;
}

// Uses the integral of each histogram to make a plot vs 
// mean mass
void TaskSideband::SimpleNormFit() {
	cout<<"Starting Simple Normalization Fit: using correlation integrals."<<endl;


}

// Produces a nice looking plot of the sideband regions
void TaskSideband::ProduceSidebandFigure() {
//	Int_t fSBColor[4] = {kAzure,kAzure-2,kAzure-4,kAzure-9};	
//	Int_t fSBColor[4] = {kViolet,kAzure,kCyan,kOrange};	

	TCanvas * cSBFigure = new TCanvas("SBFigure","SBFigure");
	TLegend * lSBFigure = new TLegend(0.43,0.55,0.88,0.88);
  Int_t nPtBinsToPlot = fMassPtBinAll.size();
//	if (fObservable == 0) nPtBinsToPlot = nObsBins?
//{ // PtBins
//		for (Int_t i = 0; i < nObsBins; i++) {
		for (Int_t i = 0; i < nPtBinsToPlot; i++) {
			TString fName = Form("SidebandFigure%d",i);
//      fMassPtBinPi0[i]->Draw();
      fMassPtBinAll[i]->Draw("E");
      fMassPtBinAll[i]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/#it{c}^{2})");
      if (i == 0) fMassPtBinAll[i]->GetYaxis()->SetRangeUser(0,1.414*fMassPtBinAll[i]->GetBinContent(fMassPtBinAll[i]->GetMaximumBin()));
			fMassPtBinPi0[i]->Draw("SAME BAR");
      fMassPtBinPi0[i]->SetMarkerColor(kRed+2);
      fMassPtBinPi0[i]->SetFillColor(fMassPtBinPi0[i]->GetMarkerColor());
    //  fMassPtBinPi0[i]->SetFillStyle(10);
			lSBFigure->AddEntry(fMassPtBinPi0[i],"#pi^{0}-Cand","flp");
			for (Int_t j = 0; j < kNSB; j++) {
//			for (Int_t j = 0; j < fNSB; j++) {
   //     if (!fSidebandFitMask[j]) continue;
				fMassPtBinSB[i][j]->Draw("SAME BAR");
				fMassPtBinSB[i][j]->SetLineColor(fSBColor[j]);
				fMassPtBinSB[i][j]->SetFillColor(fSBColor[j]);
        fMassPtBinSB[i][j]->SetFillColorAlpha(fSBColor[j],0.7);
				fMassPtBinSB[i][j]->SetMarkerColor(fSBColor[j]);
				lSBFigure->AddEntry(fMassPtBinSB[i][j],Form("SB %d",j+1),"flp");
			}	
			lSBFigure->Draw("SAME");
      lSBFigure->SetLineColor(0);
      fMassPtBinAll[i]->Draw("SAME E");
			PrintCanvas(cSBFigure,fName);
			lSBFigure->Clear();
		}
//	}


  // Do the same thing, but just for the selected sidebands
		for (Int_t i = 0; i < nPtBinsToPlot; i++) {
			TString fName = Form("SidebandFigureSelection%d",i);
//      fMassPtBinPi0[i]->Draw();
      fMassPtBinAll[i]->Draw("E");
      fMassPtBinAll[i]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/#it{c}^{2})");
      if (i == 0) fMassPtBinAll[i]->GetYaxis()->SetRangeUser(0,1.414*fMassPtBinAll[i]->GetBinContent(fMassPtBinAll[i]->GetMaximumBin()));
			fMassPtBinPi0[i]->Draw("SAME BAR");
      fMassPtBinPi0[i]->SetMarkerColor(kRed+2);
      fMassPtBinPi0[i]->SetFillColor(fMassPtBinPi0[i]->GetMarkerColor());
    //  fMassPtBinPi0[i]->SetFillStyle(10);
			lSBFigure->AddEntry(fMassPtBinPi0[i],"#pi^{0}-Cand","flp");
			for (Int_t j = 0; j < kNSB; j++) {
//			for (Int_t j = 0; j < fNSB; j++) {
        if (!fSidebandMask[j]) continue;
        //if (!fSidebandFitMask[j]) continue;
				fMassPtBinSB[i][j]->Draw("SAME BAR");
				fMassPtBinSB[i][j]->SetLineColor(fSBColor[j]);
				fMassPtBinSB[i][j]->SetFillColor(fSBColor[j]);
        fMassPtBinSB[i][j]->SetFillColorAlpha(fSBColor[j],0.7);
				fMassPtBinSB[i][j]->SetMarkerColor(fSBColor[j]);

        bool bDoThingForPerformance=false;
        if (bDoThingForPerformance) {
          fMassPtBinAll[i]->GetYaxis()->SetTitle("Counts");
          lSBFigure->AddEntry(fMassPtBinSB[i][j],Form("SB %d",j+1),"flp");
          //lSBFigure->AddEntry(fMassPtBinSB[i][j],Form("SB %d",j),"flp");
        } else {
          lSBFigure->AddEntry(fMassPtBinSB[i][j],Form("SB %d",j+1),"flp");
        }

			}	
			lSBFigure->Draw("SAME");
      lSBFigure->SetLineColor(0);
      fMassPtBinAll[i]->Draw("SAME E");
			PrintCanvas(cSBFigure,fName);
			lSBFigure->Clear();
		}


	PrintCanvas(cSBFigure,"SidebandFigure");
}

// Make Graphs comparing the DPhi correlations in the different sidebands
// In
void TaskSideband::ProduceSidebandComparison() {
 // return; // FIXME Something crashes here in data and sometimes MC
  printf("Starting the sideband comparison\n");
//	Int_t fSBColor[4] = {kAzure,kAzure-2,kAzure-4,kAzure-9};
//  Int_t fSBStyle[4] = {kFullSquare,kFullCircle,kFullDiamond,kOpenSquare};
  TCanvas * cSBCmp = new TCanvas("SBCmp","SBCmp");
  TLegend * lSBCmp = new TLegend(0.55,0.55,0.88,0.88);
  cSBCmp->Divide(1,2,0.0,0.0);

  // Full DEta

  for (Int_t i = 0; i < nObsBins; i++) {
    printf("Doing Sideband Comparison for FullDEta ObsBin %d\n",i);
    cSBCmp->Clear("D");
    //cSBCmp->Clear();
   // cSBCmp->Divide(1,2,0.0,0.0);
    cSBCmp->cd(1);
    TString fName = Form("dPhi_SB_ObsBin%d_Full",i);
    //TString fName = Form("dPhi_ObsBin%d_Full",i);

    Float_t fMin = 1e9;
    Float_t fMax = 0;

    fFullDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}");
    //fFullDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{Trig}} #frac{d^{2}N^{assoc}}{d#Delta#eta d#Delta#phi}");
    fFullDPhiSB[i][0]->GetYaxis()->SetTitleOffset(0.7);

    TH1D * hSum = (TH1D *) fFullDPhiSB[i][0]->Clone(Form("Sum_%d",i));
    printf("Debug fNSB = %d\n",fNSB);
    for (Int_t j = 0; j < fNSB; j++) {
      if (j == 0) {
        fFullDPhiSB[i][j]->Draw();
      }
      else {
        hSum->Add(fFullDPhiSB[i][j]);
        fFullDPhiSB[i][j]->Draw("SAME");
      }
      fFullDPhiSB[i][j]->SetLineColor(fSBColor[j]);
      fFullDPhiSB[i][j]->SetMarkerColor(fSBColor[j]);
      fFullDPhiSB[i][j]->SetMarkerStyle(fSBStyle[j]);
      fFullDPhiSB[i][j]->SetMarkerSize(1.0);

      Float_t fLocalMin = fFullDPhiSB[i][j]->GetBinContent(fFullDPhiSB[i][j]->GetMinimumBin());
      Float_t fLocalMax = fFullDPhiSB[i][j]->GetBinContent(fFullDPhiSB[i][j]->GetMaximumBin());
      fMin = min(fMin,fLocalMin);
      fMax = max(fMax,fLocalMax);


      lSBCmp->AddEntry(fFullDPhiSB[i][j],Form("SB %d",j+1),"lp");
    }
    fFullDPhiSB[i][0]->GetYaxis()->SetRangeUser(fMin,fMax);

    lSBCmp->SetHeader(Form("%.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}, %.1f #leq p_{T}^{a} < %.1f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt,fObsBins[i],fObsBins[i+1]));

    lSBCmp->Draw("SAME");
    // Drawing the mean ratio
    hSum->Scale(1.0/fNSB);
    vector<TH1D *> hRatios = {};
    for (Int_t j = 0; j < fNSB; j++) {
      TH1D * hRatio = (TH1D *) fFullDPhiSB[i][j]->Clone(Form("%s_Ratio",fFullDPhiSB[i][j]->GetName()));
      printf("  Dividing %d\n",j);
      hRatio->Divide(hSum);
      hRatio->GetYaxis()->SetRangeUser(2.*fMin/(fMin+fMax),2*fMax/(fMin+fMax));
      hRatio->SetLineColor(fSBColor[j]);
      hRatio->SetMarkerColor(fSBColor[j]);
      hRatio->SetMarkerStyle(fSBStyle[j]);

      hRatios.push_back(hRatio);
    }
    TF1 * fUnity = new TF1("unity","1",hRatios[0]->GetXaxis()->GetXmin(),hRatios[0]->GetXaxis()->GetXmax());
    fUnity->SetLineStyle(2);
    fUnity->SetLineWidth(4);
    fUnity->SetLineColor(kGray);
    cSBCmp->cd(2);
    for (Int_t j = 0; j < fNSB; j++) {
      if (j == 0) hRatios[j]->Draw();
      else hRatios[j]->Draw("SAME");
    }
    fUnity->Draw("SAME");
    printf("Trying to save the plots\n");
    PrintCanvas(cSBCmp,fName);
    lSBCmp->Clear();
 //   cSBCmp->Clear("D");
  //  cSBCmp->Clear(); // This was causing an error
    delete hSum; // This seemed to fix a major error. But maybe not
    delete fUnity;
    for (auto thing : hRatios) delete thing;
    printf(" ... Finished FullDEta Sideband Comparison for ObsBin %d\n",i);
  }

  // Near DEta

  for (Int_t i = 0; i < nObsBins; i++) {
    printf("Doing Sideband Comparison for NearDEta ObsBin %d\n",i);
    cSBCmp->Clear("D");
    cSBCmp->cd(1);
    TString fName = Form("dPhi_SB_ObsBin%d_NearEta",i);

    Float_t fMin = 1e9;
    Float_t fMax = 0;

    fNearEtaDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}");
    //fNearEtaDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{Trig}} #frac{d^{2}N^{assoc}}{d#Delta#eta d#Delta#phi}");
    fNearEtaDPhiSB[i][0]->GetYaxis()->SetTitleOffset(0.7);

    TH1D * hSum = (TH1D *) fNearEtaDPhiSB[i][0]->Clone(Form("Sum_%d",i));
    printf("Debug fNSB = %d\n",fNSB);
    for (Int_t j = 0; j < fNSB; j++) {
      if (j == 0) {
        fNearEtaDPhiSB[i][j]->Draw();
      }
      else {
        hSum->Add(fNearEtaDPhiSB[i][j]);
        fNearEtaDPhiSB[i][j]->Draw("SAME");
      }
      fNearEtaDPhiSB[i][j]->SetLineColor(fSBColor[j]);
      fNearEtaDPhiSB[i][j]->SetMarkerColor(fSBColor[j]);
      fNearEtaDPhiSB[i][j]->SetMarkerStyle(fSBStyle[j]);

      Float_t fLocalMin = fNearEtaDPhiSB[i][j]->GetBinContent(fNearEtaDPhiSB[i][j]->GetMinimumBin());
      Float_t fLocalMax = fNearEtaDPhiSB[i][j]->GetBinContent(fNearEtaDPhiSB[i][j]->GetMaximumBin());
      fMin = min(fMin,fLocalMin);
      fMax = max(fMax,fLocalMax);


      lSBCmp->AddEntry(fNearEtaDPhiSB[i][j],Form("SB %d",j+1),"lp");
    }
    fNearEtaDPhiSB[i][0]->GetYaxis()->SetRangeUser(fMin,fMax);

    lSBCmp->SetHeader(Form("%.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}, %.1f #leq p_{T}^{a} < %.1f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt,fObsBins[i],fObsBins[i+1]));

    lSBCmp->Draw("SAME");
    // Drawing the mean ratio
    hSum->Scale(1.0/fNSB);
    vector<TH1D *> hRatios = {};
    for (Int_t j = 0; j < fNSB; j++) {
      TH1D * hRatio = (TH1D *) fNearEtaDPhiSB[i][j]->Clone(Form("%s_Ratio",fNearEtaDPhiSB[i][j]->GetName()));
      printf("  Dividing %d\n",j);
      hRatio->Divide(hSum);
      hRatio->GetYaxis()->SetRangeUser(2.*fMin/(fMin+fMax),2*fMax/(fMin+fMax));
      hRatio->SetLineColor(fSBColor[j]);
      hRatio->SetMarkerColor(fSBColor[j]);
      hRatio->SetMarkerStyle(fSBStyle[j]);

      hRatios.push_back(hRatio);
    }
    TF1 * fUnity = new TF1("unity","1",hRatios[0]->GetXaxis()->GetXmin(),hRatios[0]->GetXaxis()->GetXmax());
    fUnity->SetLineStyle(2);
    fUnity->SetLineWidth(4);
    fUnity->SetLineColor(kGray);
    cSBCmp->cd(2);
    for (Int_t j = 0; j < fNSB; j++) {
      if (j == 0) hRatios[j]->Draw();
      else hRatios[j]->Draw("SAME");
    }
    fUnity->Draw("SAME");
    printf("Trying to save the plots\n");
    PrintCanvas(cSBCmp,fName);
    lSBCmp->Clear();
 //   cSBCmp->Clear("D");
  //  cSBCmp->Clear(); // This was causing an error
    delete hSum; // This seemed to fix a major error. But maybe not
    delete fUnity;
    for (auto thing : hRatios) delete thing;
    printf(" ... Finished NearEta Sideband Comparison for ObsBin %d\n",i);
  }

  // Far Eta

  for (Int_t i = 0; i < nObsBins; i++) {
    printf("Doing Sideband Comparison for FarDEta ObsBin %d\n",i);
    cSBCmp->Clear("D");
    cSBCmp->cd(1);
    TString fName = Form("dPhi_SB_ObsBin%d_FarEta",i);

    Float_t fMin = 1e9;
    Float_t fMax = 0;

    fFarEtaDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}");
    //fFarEtaDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{Trig}} #frac{d^{2}N^{assoc}}{d#Delta#eta d#Delta#phi}");
    fFarEtaDPhiSB[i][0]->GetYaxis()->SetTitleOffset(0.7);

    TH1D * hSum = (TH1D *) fFarEtaDPhiSB[i][0]->Clone(Form("Sum_%d",i));
    printf("Debug fNSB = %d\n",fNSB);
    for (Int_t j = 0; j < fNSB; j++) {
      if (j == 0) {
        fFarEtaDPhiSB[i][j]->Draw();
      }
      else {
        hSum->Add(fFarEtaDPhiSB[i][j]);
        fFarEtaDPhiSB[i][j]->Draw("SAME");
      }
      fFarEtaDPhiSB[i][j]->SetLineColor(fSBColor[j]);
      fFarEtaDPhiSB[i][j]->SetMarkerColor(fSBColor[j]);
      fFarEtaDPhiSB[i][j]->SetMarkerStyle(fSBStyle[j]);

      Float_t fLocalMin = fFarEtaDPhiSB[i][j]->GetBinContent(fFarEtaDPhiSB[i][j]->GetMinimumBin());
      Float_t fLocalMax = fFarEtaDPhiSB[i][j]->GetBinContent(fFarEtaDPhiSB[i][j]->GetMaximumBin());
      fMin = min(fMin,fLocalMin);
      fMax = max(fMax,fLocalMax);


      lSBCmp->AddEntry(fFarEtaDPhiSB[i][j],Form("SB %d",j+1),"lp");
    }
    fFarEtaDPhiSB[i][0]->GetYaxis()->SetRangeUser(fMin,fMax);

    lSBCmp->SetHeader(Form("%.0f #leq p_{T}^{#pi^{0}} < %.0f GeV/#it{c}, %.1f #leq p_{T}^{a} < %.1f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt,fObsBins[i],fObsBins[i+1]));

    lSBCmp->Draw("SAME");
    // Drawing the mean ratio
    hSum->Scale(1.0/fNSB);
    vector<TH1D *> hRatios = {};
    for (Int_t j = 0; j < fNSB; j++) {
      TH1D * hRatio = (TH1D *) fFarEtaDPhiSB[i][j]->Clone(Form("%s_Ratio",fFarEtaDPhiSB[i][j]->GetName()));
      printf("  Dividing %d\n",j);
      hRatio->Divide(hSum);
      hRatio->GetYaxis()->SetRangeUser(2.*fMin/(fMin+fMax),2*fMax/(fMin+fMax));
      hRatio->SetLineColor(fSBColor[j]);
      hRatio->SetMarkerColor(fSBColor[j]);
      hRatio->SetMarkerStyle(fSBStyle[j]);

      hRatios.push_back(hRatio);
    }
    TF1 * fUnity = new TF1("unity","1",hRatios[0]->GetXaxis()->GetXmin(),hRatios[0]->GetXaxis()->GetXmax());
    fUnity->SetLineStyle(2);
    fUnity->SetLineWidth(4);
    fUnity->SetLineColor(kGray);
    cSBCmp->cd(2);
    for (Int_t j = 0; j < fNSB; j++) {
      if (j == 0) hRatios[j]->Draw();
      else hRatios[j]->Draw("SAME");
    }
    fUnity->Draw("SAME");
    printf("Trying to save the plots\n");
    PrintCanvas(cSBCmp,fName);
    lSBCmp->Clear();
    delete hSum; // This seemed to fix a major error. But maybe not
    delete fUnity;
    for (auto thing : hRatios) delete thing;
    printf(" ... Finished NearEta Sideband Comparison for ObsBin %d\n",i);
  }


}


// Use the chosen background histograms to produce the predicted background correlations
//   for dphi.
//   Uses the Sideband mass effect to produce the mass scale, and applies the 1-p scale
void TaskSideband::ProduceBackground() {
	cout<<"Producing Background Prediction"<<endl;
	fFullPredBkgDPhi = {};

  // Delta Phi Projections
	// Iterate over full, near, far
	for (Int_t i = 0; i < nObsBins; i++) {
		fFullPredBkgDPhi.push_back(MergeAndScaleBkg(i,0));
	}
	for (Int_t i = 0; i < nObsBins; i++) {
		fNearEtaPredBkgDPhi.push_back(MergeAndScaleBkg(i,1));
	}
	for (Int_t i = 0; i < nObsBins; i++) {
		fFarEtaPredBkgDPhi.push_back(MergeAndScaleBkg(i,2));
	}
  // Delta Eta Projections
	for (Int_t i = 0; i < nObsBins; i++) {
		fNearSideSubPredBkgDEta.push_back(MergeAndScaleDEtaBkg(i,0));
	}





	TCanvas * cPredBkg = new TCanvas("PredBkg","PredBkg");
	Int_t nResults = (Int_t) fFullPredBkgDPhi.size(); // being lazy
	cPredBkg->Divide(3,3);
	for (Int_t i = 0; i < nResults; i++) {
		cPredBkg->cd(i+1);
		if (fFullPredBkgDPhi[i]) fFullPredBkgDPhi[i]->Draw();
    // Adding the far and near eta variations
		if (fNearEtaPredBkgDPhi[i]) {
      fNearEtaPredBkgDPhi[i]->SetMarkerStyle(kOpenCircle);
      fNearEtaPredBkgDPhi[i]->SetMarkerColor(kBlue);
      fNearEtaPredBkgDPhi[i]->SetLineColor(kBlue);
      fNearEtaPredBkgDPhi[i]->Draw("SAME");
    }
		if (fFarEtaPredBkgDPhi[i]) {
      fFarEtaPredBkgDPhi[i]->SetMarkerStyle(kOpenSquare);
      fFarEtaPredBkgDPhi[i]->Draw("SAME");
    }
	}
	PrintCanvas(cPredBkg,"PredBkg");

  // Now draw this for the thing.
  cPredBkg->Clear();
  cPredBkg->Divide(3,3);
  for (Int_t i = 0; i < nResults; i++) {
    cPredBkg->cd(i+1);
    if (fNearSideSubPredBkgDEta[i]) fNearSideSubPredBkgDEta[i]->Draw();
  }
  PrintCanvas(cPredBkg,"PredBkgDEta");
}

void TaskSideband::SetSidebandMode(Int_t input) {
  // Reset Mask
  fSidebandMask[0]=1;
  fSidebandMask[1]=1;
  fSidebandMask[2]=1;
  fSidebandMask[3]=1;

  fSidebandMode  = input;
  if (fSidebandMode  == 0) {

    switch (fBackgroundSelection) {
      case 0: // Merge All 4
      default:
        fNSB = 4; // stays 4
        break;
      case 1: // Merge SB3+SB4
        fNSB = 2;
        fSidebandMask[2]=0;
        fSidebandMask[3]=0;
        break;
      case 2: // Merge SB5+SB6
        fNSB = 2;
        fSidebandMask[0]=0;
        fSidebandMask[1]=0;
        break;
      case 3: //Merge SB4+SB5
        fNSB = 2;
        fSidebandMask[0]=0;
        fSidebandMask[3]=0;
        break;
      case 4: //Merge SB4+SB5+SB6
        fNSB = 3;
        fSidebandMask[0]=0;
        break;
      case 5: // SB3 Only
        fNSB = 1;
        fSidebandMask[1]=0;
        fSidebandMask[2]=0;
        fSidebandMask[3]=0;
        break;
      case 6: // SB4 Only
        fNSB = 1;
        fSidebandMask[0]=0;
        fSidebandMask[2]=0;
        fSidebandMask[3]=0;
        break;
      case 7: // SB5 Only
        fNSB = 1;
        fSidebandMask[0]=0;
        fSidebandMask[1]=0;
        fSidebandMask[3]=0;
        break;
      case 8: // SB6 Only
        fNSB = 1;
        fSidebandMask[0]=0;
        fSidebandMask[1]=0;
        fSidebandMask[2]=0;
        break;
    }
  } else if (fSidebandMode==1) {
    fSidebandMask[3]=0; // unnecessary?
    kNSB = 3;
    fNSBFit = 3;
    switch (fBackgroundSelection) {
      case 0: // Merge All 3
      default:
        fNSB = 3; // stays 3
        break;
      case 1: // Merge SB1+SB2
        fNSB = 2;
        fSidebandMask[2]=0;
        break;
      case 2: // Merge SB2+SB3
        fNSB = 2;
        fSidebandMask[0]=0;
        break;
      case 3: //Merge SB1+SB3
        fNSB = 2;
        fSidebandMask[1]=0;
        break;
      case 4: // Use SB1 only
        fNSB = 1;
        fSidebandMask[1]=0;
        fSidebandMask[2]=0;
        break;
      case 5: // Use SB2 only
        fNSB = 1;
        fSidebandMask[0]=0;
        fSidebandMask[2]=0;
        break;
      case 6: // Use SB3 only
        fNSB = 1;
        fSidebandMask[0]=0;
        fSidebandMask[1]=0;
        break;
    }
  }


}

// Merges the selected DeltaPhi histograms and applies (1-P) scale
// Add index to change between full, near, far
TH1D * TaskSideband::MergeAndScaleBkg(Int_t index, Int_t iType = 0) {
	cout<<"Merging and Scaling background for Obs Bin "<<index<<endl;
	// fBackgroundSelection

//  Bool_t fSidebandMask[4] = {1,1,1,1};

  printf("Number of sidebands for summing (fNSB) = %d\n",fNSB);
 
	// FIXME for scaling with the mass effect, i need the final mean mass.  This should be weighted between the four sidebands based on number of triggers

	TH1D * fPredBkg = 0;
	TString fName = ""; //Form("PredBkgDPhi_%d",index);

  // Concept: change the next step to just create the fPredBkg object and array of used sideband histograms

  vector<TH1D *> fSBHists = {};
  vector<double> fSBMasses = {};

  printf("Debug: fMeanMassSBVal has size %d\n",(int) fMeanMassSBVal.size());
  printf("Debug: fMeanMassSBVal[0] has size %d\n",(int) fMeanMassSBVal[0].size());

  printf("Just made the temporary SB mass array\n");
  switch (iType) {
    case 2:
      fName = Form("PredBkgDPhi_FarEta_%d",index);
      fPredBkg = (TH1D *) fFarEtaDPhiSB[index][0]->Clone(fName);
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
        fPredBkg->Add(fFarEtaDPhiSB[index][j]);
        fSBHists.push_back(fFarEtaDPhiSB[index][j]);
      }
    break;
    case 1:
      fName = Form("PredBkgDPhi_NearEta_%d",index);
      fPredBkg = (TH1D *) fNearEtaDPhiSB[index][0]->Clone(fName);
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
        fPredBkg->Add(fNearEtaDPhiSB[index][j]);
        fSBHists.push_back(fNearEtaDPhiSB[index][j]);
      }
    break;
    case 0:
    default:
      fName = Form("PredBkgDPhi_Full_%d",index);
      fPredBkg = (TH1D *) fFullDPhiSB[index][0]->Clone(fName);
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
        printf("Adding to predicted background histogram %s\n",fFullDPhiSB[index][j]->GetName());
        fPredBkg->Add(fFullDPhiSB[index][j]);
        fSBHists.push_back(fFullDPhiSB[index][j]);
      }
  }
  fPredBkg->GetXaxis()->SetTitle("#Delta#phi");

  printf("Just did the sideband steps\n");

  if (iExtrapolatorMode == 0) {
    printf("Integral before scale = %f\n",fPredBkg->Integral());

    printf("Scaling background prediction by 1. / %d\n",fNSB);
    fPredBkg->Scale(1./(fNSB));

    printf("Integral after scale = %f\n",fPredBkg->Integral());
  }
	// Add the background, then scale with Sideband mass effect
//	switch (fBackgroundSelection) {

//		case 0: // Merge SB3,4,5,6
//		default:
//		fPredBkg = (TH1D *) fFullDPhiSB[index][0]->Clone(fName);
//		fPredBkg->GetXaxis()->SetTitle("#Delta#phi");
//		for (Int_t j = 1; j < fNSB; j++) fPredBkg->Add(fFullDPhiSB[index][j]);
		// scale by 1/4?
//		fPredBkg->Scale(1./(fNSB));
//	}


	Int_t iMassIndex = index;
	Int_t iPurityIndex = index;
	if (fObservable != 0) {
		iMassIndex = 0;
		iPurityIndex = iPtBin-1;
	}


  for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
    fSBMasses.push_back(fMeanMassSBVal[iMassIndex][j]);
  }
  if (iExtrapolatorMode == 1) {
    RunSidebandExtrapolator(fPredBkg,fSBHists,fSBMasses,fMeanMassPi0Val[iMassIndex]);
  }

	Double_t fMassScale = 1;
	// vector<Double_t> fMassScale_Err
	Double_t fEffectiveMass = 0.;  // Effective mass of sideband region
	// FIXME For now averaging the sideband mean masses.  Should be weighted mean based on statistics
	// reference: fMeanMassPi0Val, fMeanMassSBVal[iObs][jSB]
	//          fit: fMassEffectFit
	printf("Obs Index = %d\n",index);
	for (Int_t j = 0; j < kNSB; j++ ) {
    if (fSidebandMask[j]) fEffectiveMass += fMeanMassSBVal[iMassIndex][j];
	}
	fEffectiveMass = fEffectiveMass / fNSB;
	TF1 * fMassFit = 0;
	fMassFit = fMassEffectFit[index];
	if (!fMassFit) {
		printf("Problem missing fMassFit \n");
  }
	printf("Effective mass = %f\n",fEffectiveMass);
	printf("Want to scale with par0 = %f\n",fMassFit->GetParameter(0));
	printf("                 effVal = %f\n",fMassFit->Eval(fEffectiveMass));

	fMassScale = fMassFit->GetParameter(0) / fMassFit->Eval(fEffectiveMass);
	// FIXME this index is only correct when Obs=0 (Trigger Pt).

  // possible fixed with iPtBin - 1
  printf("Using Purity index %d\n",iPurityIndex);
  //Double_t fPurity = Pi0YieldTotalRatio->GetY()[iPurityIndex]
  Double_t fPurity     = fPurityArray[iPurityIndex];
  Double_t fPurity_Err = fPurityArray[iPurityIndex];

  // These options are already applied in the fPurity Array
 /* switch (iPurityChoice) {
    case 4:
      fPurity = fPurity + fPurity_Err;
      break;
    case 3:
      fPurity = fPurity - fPurity_Err;
      break;
    case 2:
      fPurity = 1;
      break;
    case 0:
      fPurity = 0;
      break;
    default:
    case 1:
      break;
  }*/

	Double_t fPurityScale = 1. - fPurity ; // (1 - purity)

//	Double_t fPurityScale = 1. - Pi0YieldTotalRatio->GetY()[3]; // (1 - purity) // FIXME temp
//	Double_t fPurityScale = 1. - Pi0YieldTotalRatio->GetY()[index]; // (1 - purity)


  if (iExtrapolatorMode != 0) fMassScale = 1; //clunky way of disabling the default SB scale correction

	printf("Using mass scaling   = %f\n",fMassScale);
  printf("Found pi0 purity = %f\n",fPurity);
	printf("Using Purity scaling (1-p) = %f\n",fPurityScale);

	fPredBkg->Scale(fMassScale * fPurityScale);	

  printf("Done producing background estimate %s\n\n",fPredBkg->GetName());

  fPredBkg->SetMarkerSize(1); // Not sure where this gets messed up

	return fPredBkg;
}


// Merges the selected DeltaEta histograms and applies (1-P) scale
// Add index to change between different DeltaPhi regions
// Code mostly copied from the DeltaPhi version. Could be merged
TH1D * TaskSideband::MergeAndScaleDEtaBkg(Int_t index, Int_t iType = 0) {
	cout<<"Merging and Scaling DeltaEta background for Obs Bin "<<index<<endl;

	TH1D * fPredBkg = 0;
	TString fName = ""; //Form("PredBkgTitle_%d",index);

  vector<TH1D *> fSBHists = {};
  vector<double> fSBMasses = {};

  printf("Debug: fMeanMassSBVal has size %d\n",(int) fMeanMassSBVal.size());
  printf("Debug: fMeanMassSBVal[0] has size %d\n",(int) fMeanMassSBVal[0].size());

  switch (iType) {

    case 0: // Nearside -pi/2 to pi/2 Minus Awayside
    default:
      cout<<"Creating background predictions for the nearside (with awayside subtraction)"<<endl;
      fName = Form("PredBkgDEta_NearsideSub_%d",index);

      cout<<"NearSideDEtaSub vector has size "<<fNearSideSubDEtaSB.size()<<endl;
      cout<<"NearSideDEtaSub[index] vector has size "<<fNearSideSubDEtaSB[index].size()<<endl;

      fPredBkg = (TH1D *)  fNearSideSubDEtaSB[index][0]->Clone(fName);
      if (!fPredBkg) {
        cout<<"Error: missing "<<fName.Data()<<endl;
      }
      cout<<"Built the template histogram"<<endl;
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
        if (!fNearSideSubDEtaSB[index][j]) cout<<"Missing NearSideSub DEta for SB "<<j<<endl;
        fPredBkg->Add(fNearSideSubDEtaSB[index][j]);
        fSBHists.push_back(fNearSideSubDEtaSB[index][j]);
      }
  }
  fPredBkg->GetXaxis()->SetTitle("#Delta#eta");

  
  if (iExtrapolatorMode == 0) {
    fPredBkg->Scale(1./fNSB);
    cout<<"Finished averaging sidebands"<<endl;
  }


	Int_t iMassIndex = index;
	Int_t iPurityIndex = index;
	if (fObservable != 0) {
		iMassIndex = 0;
		iPurityIndex = iPtBin-1;
	}

  for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
    fSBMasses.push_back(fMeanMassSBVal[iMassIndex][j]);
  }
  if (iExtrapolatorMode == 1) {
    RunSidebandExtrapolator(fPredBkg,fSBHists,fSBMasses,fMeanMassPi0Val[iMassIndex]);
  }

	Double_t fMassScale = 1;
	Double_t fEffectiveMass = 0.;  // Effective mass of sideband region

	printf("Obs Index = %d\n",index);
	for (Int_t j = 0; j < kNSB; j++ ) {
    if (fSidebandMask[j]) fEffectiveMass += fMeanMassSBVal[iMassIndex][j];
	}
	fEffectiveMass = fEffectiveMass / fNSB;
	TF1 * fMassFit = 0;
	fMassFit = fMassEffectFit[index];
	if (!fMassFit) {
		printf("Problem missing fMassFit \n");
  }

	fMassScale = fMassFit->GetParameter(0) / fMassFit->Eval(fEffectiveMass);


  if (iExtrapolatorMode != 0) fMassScale = 1; //clunky way of disabling the default SB scale correction


  Double_t fPurity     = fPurityArray[iPurityIndex];
  Double_t fPurity_Err = fPurityArray[iPurityIndex];

	Double_t fPurityScale = 1. - fPurity ; // (1 - purity)
	fPredBkg->Scale(fMassScale * fPurityScale);	
  printf("Done producing background estimate %s\n\n",fPredBkg->GetName());

  fPredBkg->SetMarkerSize(1); // Not sure where this gets messed up

	return fPredBkg;

}

void TaskSideband::PlotBkgAndSignal() {
	cout<<"Plotting Signal and Background"<<endl;

	TCanvas * cBkgSignal = new TCanvas("BkgSignal","BkgSignal");
	Int_t nX = 3;
	Int_t nY = 3;
	cBkgSignal->Divide(nX,nY); // FIXME what if nObs different?
	for (Int_t i = 0; i < nObsBins; i++) {
		Double_t fMin = 0;
		Double_t fMax = max(fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMaximumBin()),fFullPredBkgDPhi[i]->GetBinContent(fFullPredBkgDPhi[i]->GetMaximumBin()));
		fMax = 1.1 * fMax;

		cBkgSignal->cd(i+1);

    fFullDPhiPi0[i]->SetMarkerColor(kRawPi0Color);
    fFullDPhiPi0[i]->SetLineColor(kRawPi0Color);
		fFullDPhiPi0[i]->SetMarkerStyle(kFullSquare);
    fFullDPhiPi0[i]->SetMarkerSize(1.0);
    fFullDPhiPi0[i]->GetXaxis()->SetLabelSize(0.04);
    fFullDPhiPi0[i]->GetYaxis()->SetLabelSize(0.042);
    fFullDPhiPi0[i]->GetXaxis()->SetTitleSize(0.05);
    fFullDPhiPi0[i]->GetYaxis()->SetTitleSize(0.053);
    fFullDPhiPi0[i]->GetYaxis()->SetTitleOffset(0.88);
    fFullDPhiPi0[i]->GetXaxis()->SetTitleOffset(0.7); //0.8
    fFullDPhiPi0[i]->GetXaxis()->CenterTitle(kFALSE);
    fFullDPhiPi0[i]->GetXaxis()->SetTitle("#Delta#varphi"); // Removing the pi0-h part, to be inclusive to sideband
    fFullDPhiPi0[i]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}");

		fFullDPhiPi0[i]->Draw();
		fFullDPhiPi0[i]->GetYaxis()->SetRangeUser(fMin,fMax);
		fFullPredBkgDPhi[i]->SetLineColor(kScaledSBColor);
		fFullPredBkgDPhi[i]->SetMarkerColor(kScaledSBColor);
		fFullPredBkgDPhi[i]->SetMarkerStyle(kFullSquare);
		fFullPredBkgDPhi[i]->Draw("SAME");	



	}
	TLegend * leg = new TLegend(0.40,0.38,0.85,0.7); //(0.44,0.4,0.89,0.7)
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->AddEntry((TObject *)0,Form("%.0f #leq #it{p}_{T}^{trigger} < %.0f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt),"");
//  leg->AddEntry((TObject *) 0,Form("%s GeV/#it{c}",fFullDPhiPi0[0]->GetTitle()),"");
	leg->AddEntry(fFullDPhiPi0[0],"#pi^{0}_{cand} Corr.","lp");
  if (fFixedPurity < 0) leg->AddEntry(fFullPredBkgDPhi[0],"Scaled SB Corr.","lp");
  else leg->AddEntry(fFullPredBkgDPhi[0],Form("Scaled SB Corr. \n(Assuming %.0f%% purity)",fFixedPurity*100.));
	cBkgSignal->cd(nX*nY);
	leg->Draw("SAME");

	PrintCanvas(cBkgSignal,"BkgSigCmp");


  // Now Draw individual versions
  cBkgSignal->Clear();
	for (Int_t i = 0; i < nObsBins; i++) {
    leg->Clear();

    leg->SetY1NDC(0.49);
    leg->SetY2NDC(0.73);

    Float_t fArbScale = 4.5;
    Float_t kArbScaleColor = kScaledSBColor;


    Float_t fRatioSBtoPi0 = fFullDPhiPi0[i]->Integral()/fFullPredBkgDPhi[i]->Integral();
    if (fRatioSBtoPi0 > 5) fArbScale = 9.;


    /*
    // FIXME adding an arbitrary scaling just for the performance plot 
    Float_t fMaxValue  = fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMaximumBin());
    Float_t fMoreArbScale = 2. / fMaxValue;
    if (fMaxValue < 1.8) {
      fFullDPhiPi0[i]->Scale(fMoreArbScale);
      fFullPredBkgDPhi[i]->Scale(fMoreArbScale);
    }
    */

		fFullPredBkgDPhi_ArbScaled.push_back((TH1D *) fFullPredBkgDPhi[i]->Clone(Form("%s_ArbScaled",fFullPredBkgDPhi[i]->GetName())));	
    fFullPredBkgDPhi_ArbScaled[i]->Scale(fArbScale);
    fFullPredBkgDPhi_ArbScaled[i]->SetLineColor(kArbScaleColor);
    fFullPredBkgDPhi_ArbScaled[i]->SetMarkerColor(kArbScaleColor);
    fFullPredBkgDPhi_ArbScaled[i]->SetMarkerStyle(kOpenSquare);

    leg->AddEntry((TObject *)0,Form("%.0f #leq #it{p}_{T}^{trigger} < %.0f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt),"");
    leg->AddEntry((TObject *) 0,Form("%s GeV/#it{c}",fFullDPhiPi0[i]->GetTitle()),"");
    leg->AddEntry(fFullDPhiPi0[0],"#pi^{0}_{cand} Corr.","lp");
    if (fFixedPurity < 0) leg->AddEntry(fFullPredBkgDPhi[i],"Scaled SB Corr.","lp");
    else leg->AddEntry(fFullPredBkgDPhi[i],Form("Scaled SB Corr. \n(Assuming %.0f%% purity)",fFixedPurity*100.));
    leg->AddEntry(fFullPredBkgDPhi_ArbScaled[i],Form("Scaled SB Corr (#times%.1f for shape comparison)",fArbScale),"lp");

		Double_t fMin = 0;
		Double_t fMax = max(fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMaximumBin()),fFullPredBkgDPhi[i]->GetBinContent(fFullPredBkgDPhi[i]->GetMaximumBin()));

    double fPi0HistMin = fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMinimumBin());
    double fPi0HistMax = fMax;
    double fPi0AwaySideVal = fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->FindFixBin(TMath::Pi()));

    double fSBHistMin = fFullPredBkgDPhi[i]->GetBinContent(fFullPredBkgDPhi[i]->GetMinimumBin());
    double fSBHistMax = fFullPredBkgDPhi[i]->GetBinContent(fFullPredBkgDPhi[i]->GetMaximumBin());


		fMax = 1.4 * fMax;
		fFullDPhiPi0[i]->GetYaxis()->SetRangeUser(fMin,fMax);

    float yLegMin = 0.25; //0.33
    float yLegMax = 0.48; //0.55
    float yPerfMin = 0.74;
    float yPerfHeight = 0.15;

    if ((fMax - fPi0HistMax) < 0.8*(fPi0HistMin - fSBHistMin)) { // Cent 10-30% , where the legend goes under the pi0-can correlations
      leg->SetY1NDC(yLegMin);
      leg->SetY2NDC(yLegMax);
 //     yPerfMin = 0.55;
    } else { 
      yPerfMin = 0.74; //0.76
    }


    if (bNoYLabel) {
      fFullDPhiPi0[i]->GetYaxis()->SetLabelColor(kWhite);
      fFullDPhiPi0[i]->GetYaxis()->SetLabelSize(0.1);
      fFullDPhiPi0[i]->GetYaxis()->SetTitleOffset(0.6);
      //fFullDPhiPi0[i]->GetYaxis()->SetTitle("");
    }

		fFullDPhiPi0[i]->Draw();
//		fFullDPhiPi0[i]->GetYaxis()->SetRangeUser(fMin,fMax);
    fFullPredBkgDPhi_ArbScaled[i]->Draw("SAME");
		fFullPredBkgDPhi[i]->Draw("SAME");

    gPad->SetTickx();
    gPad->SetTicky();

    leg->Draw("SAME");
    if (bEnablePerformance) DrawAlicePerf(fFullDPhiPi0[i],0.41,yPerfMin,0.24,0.15);
    else DrawWIP(fFullDPhiPi0[i],0.28,0.44,0.33,0.25);
    PrintCanvas(cBkgSignal,Form("BkgSigCmp_Indiv_%d",i));
  }

  // Plotting these with an aribtrary scaling, and individually
	for (Int_t i = 0; i < nObsBins; i++) {

    TLegend * leg2 = new TLegend(0.45,0.48,0.85,0.70);
    //leg2->SetHeader(Form("%.0f #leq #it{p}_{T}^{trigger} < %.0f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt));
    leg2->AddEntry((TObject *)0,Form("%.0f #leq #it{p}_{T}^{trigger} < %.0f GeV/#it{c}",fGlobalMinPt,fGlobalMaxPt),"");
    leg2->AddEntry((TObject *) 0,Form("%s GeV/#it{c}",fFullDPhiPi0[i]->GetTitle()),"");
    leg2->AddEntry(fFullDPhiPi0[0],"#pi^{0}_{cand}-h corr.","lp");
    leg2->AddEntry(fFullPredBkgDPhi[0],"SB-h corr.","lp");
    leg2->SetTextSize(0.035);

    cBkgSignal->Clear(); // I hope that gets rid of the division 

		fFullPredBkgDPhi_Unscaled.push_back((TH1D *) fFullPredBkgDPhi[i]->Clone(Form("%s_Unscaled",fFullPredBkgDPhi[i]->GetName())));	

    double fAwaysideScale = fFullDPhiPi0[i]->Integral(fFullDPhiPi0[i]->FindFixBin(TMath::Pi()/2.),fFullDPhiPi0[i]->FindFixBin(3.*TMath::Pi()/2.));
    double fAwaysideScaleBkg = fFullPredBkgDPhi_Unscaled[i]->Integral(fFullDPhiPi0[i]->FindFixBin(TMath::Pi()/2.),fFullDPhiPi0[i]->FindFixBin(3.*TMath::Pi()/2.));



  //// disabling rescale
  //  fFullPredBkgDPhi_Unscaled[i]->Scale(0.5*fAwaysideScale/fAwaysideScaleBkg);

		Double_t fMin = fmin(fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMinimumBin()),fFullPredBkgDPhi_Unscaled[i]->GetBinContent(fFullPredBkgDPhi_Unscaled[i]->GetMinimumBin()));
		Double_t fMax = fmax(fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMaximumBin()),fFullPredBkgDPhi_Unscaled[i]->GetBinContent(fFullPredBkgDPhi[i]->GetMaximumBin()));
		fMax = 1.1 * fMax;
    fMin -= 0.1*(fMax - fMin);
		fFullDPhiPi0[i]->GetYaxis()->SetRangeUser(fMin,fMax);
		fFullDPhiPi0[i]->Draw();
    fFullDPhiPi0[i]->SetMarkerSize(1.0);

    if (bNoYLabel) {
      fFullDPhiPi0[i]->GetYaxis()->SetLabelColor(kWhite);
      fFullDPhiPi0[i]->GetYaxis()->SetLabelSize(0.1);
      fFullDPhiPi0[i]->GetYaxis()->SetTitleOffset(0.6);
      //fFullDPhiPi0[i]->GetYaxis()->SetTitle("");
    }

		fFullPredBkgDPhi_Unscaled[i]->Draw("SAME");	

    printf("Drawing a comparison of signal \"%s\" and background \"%s\"\n",fFullDPhiPi0[i]->GetTitle(),fFullPredBkgDPhi_Unscaled[i]->GetTitle());

    leg2->Draw("SAME");
    //DrawAlicePerf(fFullDPhiPi0[i],0.35,0.72,0.33,0.15);
    //DrawAlicePerf(fFullDPhiPi0[i],0.35,0.40,0.38,0.2);
    PrintCanvas(cBkgSignal,Form("BkgSigCmp_Unscaled_%d",i));
  }


  // FIXME Do similar for the Delta Eta Projections.


}


// May want a mode to use this to replace result of MergeAndScaleBkg(index,iType)
void TaskSideband::RunSidebandExtrapolator(TH1D * fPredBkg, vector<TH1D *> fSBHists, vector<double> fSBMasses, double fLocalPi0Mass) {
  printf("Running Sideband extrapolator\n");
  // could check for full vs near vs deta. Assume that some parts of the correction (alpha, x) will be the same for all. Purely point-by-point extrapolation would not worry.

  TCanvas * cExtrap = new TCanvas("Extrap","Extrap");

  TString sName = fPredBkg->GetName();

  // Use Vn information from gTriggerFlowSidebands_V2 ?

  // First step: linear extrapolation point-to-point.
  // Can try more fancy determination of g,h functions + alpha, x constants later

  TString fXAxisTitle = fPredBkg->GetXaxis()->GetTitle();

  TH1D * fOldSidebandSum = (TH1D *) fPredBkg->Clone("OldSidebandSum");

  fPredBkg->Reset();

  TF1 * fExtrapolatorFit = 0;

  int nPar = 2; // Default for the constant fit. Includes the Pi0 Mass parameter

  switch(fScalingFitFunction) {
    case 4:
      // Average the results of a linear fit and the closest sideband
      // Linear Fit
      fExtrapolatorFit = new TF1("ExtrapolatorFit","[0]+[2]*(x - [1])",0.1,0.5);
      nPar = 3;
      break;
    case 3:
      // Fractional Fit
      //TF1 * fExtrapolatorFit = new TF1("ExtrapolatorFit","(1+[1]*(x-[2])) / ((1/[0]) + ([1] + [3])*(x - [2]))",0.1,0.5);
      fExtrapolatorFit = new TF1("ExtrapolatorFit","([0]+[2]*(x-[1])) / (1 + ([3])*(x - [1]))",0.1,0.5);
      //  fExtrapolatorFit->SetParameter(0,0.5);
      // For fractional fit, need limits
      fExtrapolatorFit->SetParLimits(3,0,10);
      fExtrapolatorFit->SetParameter(3, 0.);
      nPar = 4;
      break;
    case 2:
      // Quadratic Fit
      fExtrapolatorFit = new TF1("ExtrapolatorFit","[0]+[2]*(x-[1])*(x - [3])",0.1,0.5);
      fExtrapolatorFit->SetParameter(3, 0.);
      nPar = 4;
      break;
    case 1:
      // Linear Fit
      fExtrapolatorFit = new TF1("ExtrapolatorFit","[0]+[2]*(x - [1])",0.1,0.5);
      nPar = 3;
      break;
    case 0:
    default:
      fExtrapolatorFit = new TF1("ExtrapolatorFit","[0]+0*(x-[1])",0.1,0.5);
      break;
  }

  // Need pi0 mass
  fExtrapolatorFit->FixParameter(1,fLocalPi0Mass);

//  fExtrapolatorFit->SetParameter(2, 0.);

  // remember: kNSB is 3 or 4, fNSB counts the selected sidebands
  TGraphErrors * gLocalSidebandPoints = new TGraphErrors(fNSB);

  int nBinsX = fPredBkg->GetNbinsX();
  double fMinX = fPredBkg->GetXaxis()->GetXmin();
  double fMaxX = fPredBkg->GetXaxis()->GetXmax();

  vector<TH1D *> ExtrapParameterHistograms = {};

  TH1D * hChiSquare = new TH1D(Form("%s_ExtrapChiSquare",sName.Data()),"",nBinsX,fMinX,fMaxX);

  for (int i = 0; i < nPar; i++) {
    TH1D * fParHist = new TH1D(Form("%s_ExtrapPar_%d",sName.Data(),i),"",nBinsX,fMinX,fMaxX);
    ExtrapParameterHistograms.push_back(fParHist);
  }


  for (int i = 1; i <= nBinsX; i++) {
    // FIXME what the hell, why was the pi0 mass being changed???
    //fExtrapolatorFit->SetParameter(1, 0.);
    //fExtrapolatorFit->SetParameter(3, 0.);

    fExtrapolatorFit->SetParameter(0, 0.); // Reseting the parameter for this next bin.
    if (fScalingFitFunction > 1) fExtrapolatorFit->SetParameter(2,0.);
    if (fScalingFitFunction > 2) fExtrapolatorFit->SetParameter(3,0.);

    std::vector<int> iEmptyBins = {};

    double fLocalAngleMin = fSBHists[0]->GetXaxis()->GetBinLowEdge(i);
    double fLocalAngleMax = fSBHists[0]->GetXaxis()->GetBinUpEdge(i);

    double fLocalMaxYErr = 0;
    for (int j = 0; j < (int) fSBHists.size(); j++) {
      double fLocalX = fSBMasses[j];
      double fLocalX_Err = 0;
      double fLocalY = fSBHists[j]->GetBinContent(i);
      double fLocalY_Err = fSBHists[j]->GetBinError(i);


      // Check if a point has value 0
      if (fLocalY == 0) iEmptyBins.push_back(j);
      if (fLocalY_Err > fLocalMaxYErr) fLocalMaxYErr = fLocalY_Err;

      gLocalSidebandPoints->SetPoint(j,fLocalX,fLocalY);
      gLocalSidebandPoints->SetPointError(j,fLocalX_Err,fLocalY_Err);
    }
    int iNumEmpty = iEmptyBins.size();
    bool isAllEmpty = false;
    switch (iNumEmpty) {
      case 3:
        isAllEmpty = true;
        printf("Extrapolation has a bin with all 3 Sidebands empty!\n");
        break;
      case 2:
      case 1:
        for (int l : iEmptyBins) {
          gLocalSidebandPoints->SetPointError(l,gLocalSidebandPoints->GetErrorX(l),fLocalMaxYErr);
        }
        break;
      case 0:
      default:
        // we are happy
        break;
    }

    double FinalY = 0;
    double FinalYErr = 0;
    
    if (!isAllEmpty)  {
      gLocalSidebandPoints->Fit(fExtrapolatorFit,"Q");

      FinalY = fExtrapolatorFit->GetParameter(0);
      FinalYErr = fExtrapolatorFit->GetParError(0);
    }

    if (fScalingFitFunction == 4) {
      double FirstSBVal = fSBHists[0]->GetBinContent(i);
      double FirstSBValErr = fSBHists[0]->GetBinError(i);

      double OldFinal = FinalY;
      FinalY = 0.5 * (FinalY + FirstSBVal);
      FinalYErr = 0.5 * TMath::Sqrt(TMath::Power(FinalYErr,2) + TMath::Power(FirstSBValErr,2));

    }

    fPredBkg->SetBinContent(i,FinalY);
    fPredBkg->SetBinError(i,FinalYErr);

    // Note: not using chisquare/ndf
    hChiSquare->SetBinContent(i,fExtrapolatorFit->GetChisquare());
    for (int k = 0; k < nPar; k++) {
      ExtrapParameterHistograms[k]->SetBinContent(i,fExtrapolatorFit->GetParameter(k));
      ExtrapParameterHistograms[k]->SetBinError(i,fExtrapolatorFit->GetParError(k));
    }
    if (fDebugLevel > 1) {
      TMultiGraph * mgDebug = new TMultiGraph();
      mgDebug->Add(gLocalSidebandPoints);


//      gLocalSidebandPoints->Draw("ALP");
      //gLocalSidebandPoints->GetXaxis()->SetRangeUser(0.1,0.45);
      //gLocalSidebandPoints->GetXaxis()->SetTitle("m_{#gamma#gamma}");
      //gLocalSidebandPoints->GetYaxis()->SetTitle("1/N_{trig} d^{2}/d#Delta#eta#Delta#phi");
      cExtrap->SetRightMargin(0.1);
      TGraphErrors * fPi0ExtrapolationPoint = new TGraphErrors(1);
      fPi0ExtrapolationPoint->SetMarkerStyle(kOpenCircle);
      fPi0ExtrapolationPoint->SetMarkerColor(kRed);
      fPi0ExtrapolationPoint->SetLineColor(kRed);
      fPi0ExtrapolationPoint->SetPoint(0,fLocalPi0Mass,FinalY);
      fPi0ExtrapolationPoint->SetPointError(0,0,FinalYErr);
      
      mgDebug->Add(fPi0ExtrapolationPoint);
      //fPi0ExtrapolationPoint->Draw("LP SAME");
      mgDebug->Draw("ALP");
      mgDebug->GetXaxis()->SetTitle("m_{#gamma#gamma}");
      mgDebug->GetYaxis()->SetTitle("1/N_{trig} d^{2}/d#Delta#eta#Delta#phi");

      TLegend * legDebug = new TLegend(0.7,0.3,0.85,0.4);

      legDebug->SetHeader(Form("%.2f#leq %s <%.2f",fLocalAngleMin,fXAxisTitle.Data(),fLocalAngleMax),"c");

      cExtrap->Print(Form("%s/QA/Indiv/Extrap_%s_Bin%d.png",fOutputDir.Data(),fPredBkg->GetName(),i));
    }
  }

  hChiSquare->Draw();
  cExtrap->Print(Form("%s/QA/Extrap_%s_ChiSquare.png",fOutputDir.Data(),sName.Data()));
  for (int i = 0; i < nPar; i++) {
    if (i == 1) continue; // Skip the fixed pi0 mass parameter  
    ExtrapParameterHistograms[i]->Draw();
    cExtrap->Print(Form("%s/QA/Extrap_%s_Param%d.png",fOutputDir.Data(),sName.Data(),i));
  }
  delete cExtrap;

}


void TaskSideband::Subtract() {
	cout<<"Subtracting Scaled background from Pi0Cand Correlations"<<endl;

	fFullDPhiFinal = {};
	fNearEtaDPhiFinal = {};
	fFarEtaDPhiFinal = {};

  fNearSideSubDEtaFinal = {};
  fAwaySideSubDEtaFinal = {};

	for (Int_t i = 0; i < nObsBins; i++) {


    // Delta Phi
		TString fName = Form("SBSub_FullDPhi_ObsBin%d",i);
		TH1D * fFullDPhiFinal_Local = (TH1D *) fFullDPhiPi0[i]->Clone(fName.Data());
		fFullDPhiFinal_Local->Add(fFullPredBkgDPhi[i],-1);

		fName = Form("SBSub_NearEtaDPhi_ObsBin%d",i);
		TH1D * fNearEtaDPhiFinal_Local = (TH1D *) fNearEtaDPhiPi0[i]->Clone(fName.Data());
		fNearEtaDPhiFinal_Local->Add(fNearEtaPredBkgDPhi[i],-1);

		fName = Form("SBSub_FarEtaDPhi_ObsBin%d",i);
		TH1D * fFarEtaDPhiFinal_Local = (TH1D *) fFarEtaDPhiPi0[i]->Clone(fName.Data());
		fFarEtaDPhiFinal_Local->Add(fFarEtaPredBkgDPhi[i],-1);

    // Delta Eta
    fName = Form("SBSub_NearSideDEta_ObsBin%d",i);
    TH1D * fNearSideSubDEtaFinal_Local = (TH1D *) fNearSideSubDEtaPi0[i]->Clone(fName.Data());
    fNearSideSubDEtaFinal_Local->Add(fNearSideSubPredBkgDEta[i],-1);


    // Scale by 1 / Purity
    Int_t iPurityIndex = i;
    if (fObservable != 0) {
      iPurityIndex = iPtBin - 1;
    }
    // Should this be using the fPurity Array?
  printf("Using Purity index %d\n",iPurityIndex);
    Double_t fPurity     = fPurityArray[iPurityIndex];
    Double_t fPurity_Err = fPurityArray[iPurityIndex]; // incorrect. Is this value ever used?
   // Double_t fPurity = Pi0YieldTotalRatio->GetY()[iPurityIndex];
    // Only the statistical uncertainty
  //  Double_t fPurity_Err = Pi0YieldTotalRatio->GetEY()[iPurityIndex];
    if (iMCMode == 2) {
      // For MC true pi0s, just let purity = 1
      // and let the rest of the code run
      fPurity = 1.0;
      fPurity_Err = 0;
    }
    if (fPurity > 0) {
      printf("Normalizing histograms by 1/Purity = 1 / %f\n",fPurity);
      fFullDPhiFinal_Local->Scale(1./fPurity);
      fNearEtaDPhiFinal_Local->Scale(1./fPurity);
      fFarEtaDPhiFinal_Local->Scale(1./fPurity);
      fNearSideSubDEtaFinal_Local->Scale(1./fPurity);
    } else {
      printf("Not Normalizing histograms, as purity <= 0\n");
    }

    //fFullDPhiFinal_Local->SetTitle(""); // FIXME I hope this doesn't affect anything later
    fFullDPhiFinal_Local->GetXaxis()->SetTitle("#Delta#phi^{#pi^{0}-h}");
    fFullDPhiFinal_Local->GetYaxis()->SetTitle("1/N^{#pi^{0}} dN/d#Delta#varphi");
    //fFullDPhiFinal_Local->GetYaxis()->SetTitle("1/N^{#pi^{0}} d^{2}N/d#Delta#etad#Delta#phi");

    //fNearSideSubDEtaFinal_Loca

		fFullDPhiFinal.push_back(fFullDPhiFinal_Local);
		fNearEtaDPhiFinal.push_back(fNearEtaDPhiFinal_Local);
		fFarEtaDPhiFinal.push_back(fFarEtaDPhiFinal_Local);

    fNearSideSubDEtaFinal.push_back(fNearSideSubDEtaFinal_Local);
	}

  Debug(3);

  // Draw a nice plot of the subtracted and rescaled correlations
  // Drawing Full range
  TCanvas * cSub = new TCanvas ("Sub","Sub",fCanvasWidth,fCanvasHeight);
  cSub->Divide(3,3,0.001,0.0012);
  Int_t nResults = fFullDPhiFinal.size();
  for (Int_t i = 0; i < nResults; i++) {
    cSub->cd(i+1);
    TH1D * hLocal = fFullDPhiFinal[i];
    //hLocal->SetMarkerStyle(4);
    hLocal->GetXaxis()->SetLabelSize(0.04);
    hLocal->GetYaxis()->SetLabelSize(0.042);
    hLocal->GetXaxis()->SetTitleSize(0.05);
    hLocal->GetYaxis()->SetTitleSize(0.053);
    hLocal->GetYaxis()->SetTitleOffset(0.88);
    hLocal->GetXaxis()->SetTitleOffset(0.7); //0.8
    hLocal->GetXaxis()->CenterTitle(kFALSE);

    TLegend* leg4=nullptr;
    if(fObservable==0)leg4 = new TLegend(0.42,0.78,0.7,0.8); //..Bkg subtracted
    if(fObservable!=0)leg4 = new TLegend(0.42,0.78,0.7,0.8); //..Bkg subtracted
    //leg4->AddEntry(PprojXFull,Form("Proj. over #Delta#eta [-1.2,1.2]"),"");

    //if (fObservable==0) leg4->AddEntry(hLocal,Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fObsBins[i],fTriggerName.Data(),fObsBins[i+1]),"pe");
    if(fObservable==1)leg4->AddEntry(hLocal,Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]),"pe");
 //   if(fObservable==2)leg4->AddEntry(hLocal,Form("%0.2f < #it{#xi} < %0.2f",ObsArray[i],ObsArray[i+1]),"pe");
  //  if(fObservable!=0)leg4->AddEntry(hLocal,"0.15 < #it{p}_{T}^{h} < 30 GeV/#it{c}","");
 //   if(fObservable!=0)leg4->AddEntry(hLocal,Form("%0.1f < #it{p}_{T}^{%s} < %0.1f GeV/#it{c}",fMinPt,fTriggerName.Data(),fMaxPt),"");
    leg4->SetTextColor(kBlack);
    leg4->SetTextSize(0.045);
    leg4->SetBorderSize(0);
    leg4->SetFillColorAlpha(10, 0);


    fFullDPhiFinal[i]->Draw();
    if (i==2) {
      if (bEnablePerformance) DrawAlicePerf(hLocal,0.28,0.5,0.33,0.25);
      else DrawWIP(hLocal,0.28,0.5,0.33,0.25);
    }
    leg4->Draw("same");

  }
  PrintCanvas(cSub,"SBSub");

  Int_t fOrigColor = kRed+1;// kOrange+1;
  Int_t fSubColor = kBlack;

  Int_t fOrigMarker = kOpenSquare;
  Int_t fSubMarker = kFullSquare;

  Float_t fOrigSize = 1.0;
  Float_t fSubSize = 1.0;


  // Draw a comparison to the input correlations
  //TCanvas * cSubCmp = new TCanvas ("SubCmp","SubCmp",fCanvasWidth,fCanvasHeight);
  cSub->Clear();
  cSub->Divide(3,3,0.001,0.0012);
  for (Int_t i = 0; i < nResults; i++) {
    cSub->cd(i+1);
    TH1D * hLocalSub = fFullDPhiFinal[i];
    TH1D * hLocalOrig = fFullDPhiPi0[i];
  
    hLocalOrig->SetLineColor(fOrigColor);
    hLocalOrig->SetMarkerColor(fOrigColor);
    hLocalSub->SetLineColor(fSubColor);
    hLocalSub->SetMarkerColor(fSubColor);


    Double_t fLocalYMin = min(hLocalSub->GetBinContent(hLocalSub->GetMinimumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMinimumBin()));
    Double_t fLocalYMax = max(hLocalSub->GetBinContent(hLocalSub->GetMaximumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMaximumBin()));
    fLocalYMax += 0.15*(fLocalYMax-fLocalYMin);
    fLocalYMin -= 0.05*(fLocalYMax-fLocalYMin);
    
    TLegend* leg4=nullptr;
    if (i==2) {
      if(fObservable==0)leg4 = new TLegend(0.38,0.57,0.87,0.8);
      if(fObservable!=0)leg4 = new TLegend(0.38,0.57,0.87,0.8);
    } else {
      if(fObservable==0)leg4 = new TLegend(0.38,0.75,0.87,0.8);
      if(fObservable!=0)leg4 = new TLegend(0.38,0.75,0.87,0.8);
    }
    //leg4->AddEntry(PprojXFull,Form("Proj. over #Delta#eta [-1.2,1.2]"),"");

    //if (fObservable==0) leg4->AddEntry(hLocal,Form("%0.0f < #it{p}_{T}^{%s} < %0.0f GeV/#it{c}",fObsBins[i],fTriggerName.Data(),fObsBins[i+1]),"pe");
    if(fObservable==1)leg4->SetHeader(Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]));
    if (i==2) {
      if(fObservable==1)leg4->AddEntry(hLocalOrig,"Raw #pi^{0}_{Cand.}-h Corr.","pe");
      if(fObservable==1)leg4->AddEntry(hLocalSub,"Sideband-Subtracted Corr.","pe");
    }
    hLocalSub->GetYaxis()->SetRangeUser(fLocalYMin,fLocalYMax);

    hLocalSub->Draw();
    hLocalOrig->Draw("SAME");
    if (i==3) {
      if (bEnablePerformance) DrawAlicePerf(hLocalSub,0.28,0.44,0.33,0.25);
      else DrawWIP(hLocalSub,0.28,0.44,0.33,0.25);
    }
    leg4->Draw("SAME");
    hLocalSub->Draw("SAME");
  
  }
  PrintCanvas(cSub,"SBSubCmp");
  // Do the same for Delta Eta plots
  cSub->Clear();
  cSub->Divide(3,3,0.001,0.0012);
  for (Int_t i = 0; i < nResults; i++) {
    cSub->cd(i+1);
    TH1D * hLocalSub = fNearSideSubDEtaFinal[i];
    TH1D * hLocalOrig = fNearSideSubDEtaPi0[i];


    hLocalOrig->SetLineColor(fOrigColor);
    hLocalOrig->SetMarkerColor(fOrigColor);
    hLocalSub->SetLineColor(fSubColor);
    hLocalSub->SetMarkerColor(fSubColor);
    Double_t fLocalYMin = min(hLocalSub->GetBinContent(hLocalSub->GetMinimumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMinimumBin()));
    Double_t fLocalYMax = max(hLocalSub->GetBinContent(hLocalSub->GetMaximumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMaximumBin()));
    fLocalYMax += 0.15*(fLocalYMax-fLocalYMin);
    fLocalYMin -= 0.05*(fLocalYMax-fLocalYMin);
    
    TLegend* leg4=nullptr;
    if (i==2) {
      if(fObservable==0)leg4 = new TLegend(0.38,0.57,0.87,0.8);
      if(fObservable!=0)leg4 = new TLegend(0.38,0.57,0.87,0.8);
    } else {
      if(fObservable==0)leg4 = new TLegend(0.38,0.75,0.87,0.8);
      if(fObservable!=0)leg4 = new TLegend(0.38,0.75,0.87,0.8);
    }
    if(fObservable==1)leg4->SetHeader(Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]));
    if (i==2) {
      if(fObservable==1)leg4->AddEntry(hLocalOrig,"Raw #pi^{0}_{Cand.}-h Corr.","pe");
      if(fObservable==1)leg4->AddEntry(hLocalSub,"Sideband-Subtracted Corr.","pe");
    }
    hLocalSub->GetYaxis()->SetRangeUser(fLocalYMin,fLocalYMax);


    hLocalSub->Draw("");
    hLocalOrig->Draw("SAME");
    if (i==3) {
      if (bEnablePerformance) DrawAlicePerf(hLocalSub,0.28,0.44,0.33,0.25);
      else DrawWIP(hLocalSub,0.28,0.44,0.33,0.25);
    }
    leg4->Draw("SAME");
    hLocalSub->Draw("SAME");

  }

  PrintCanvas(cSub,"SBSubCmpDEta");




  // Draw individual plots
  cSub->Clear();
  cSub->SetWindowSize(400,600);
  //cSub->SetWidth(400);
  //cSub->SetHeight(600);

  // Delta Phi (Full Delta Eta)
  for (Int_t i = 0; i < nResults; i++) {
    TH1D * hLocalSub = fFullDPhiFinal[i];
    TH1D * hLocalOrig = fFullDPhiPi0[i];

    hLocalOrig->SetLineColor(fOrigColor);
    hLocalOrig->SetMarkerColor(fOrigColor);
    hLocalOrig->SetMarkerStyle(fOrigMarker);
    hLocalOrig->SetMarkerSize(fOrigSize);
    hLocalSub->SetLineColor(fSubColor);
    hLocalSub->SetMarkerColor(fSubColor);
    hLocalSub->SetMarkerStyle(fSubMarker);
    hLocalSub->SetMarkerSize(fSubSize);

    Double_t fLocalYMin = min(hLocalSub->GetBinContent(hLocalSub->GetMinimumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMinimumBin()));
    Double_t fLocalYMax = max(hLocalSub->GetBinContent(hLocalSub->GetMaximumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMaximumBin()));
    fLocalYMax += 0.15*(fLocalYMax-fLocalYMin);
    fLocalYMin -= 0.05*(fLocalYMax-fLocalYMin);

    TLegend* leg4= new TLegend(0.38,0.70,0.87,0.87);
/*
    if (i==2) {
      if(fObservable==0)leg4 = new TLegend(0.38,0.57,0.87,0.8);
      if(fObservable!=0)leg4 = new TLegend(0.38,0.57,0.87,0.8);
    } else {
      if(fObservable==0)leg4 = new TLegend(0.38,0.75,0.87,0.8);
      if(fObservable!=0)leg4 = new TLegend(0.38,0.75,0.87,0.8);
    }
  */
    if(fObservable==1)leg4->SetHeader(Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]));
    if(fObservable==2)leg4->SetHeader(Form("%0.1f < #it{p}^{a}_{T} < %0.1f GeV/#it{c}",fObsBins[i],fObsBins[i+1]),"c");
    leg4->AddEntry(hLocalOrig,"Raw #pi^{0}_{Cand.}-h Corr.","pe");
    leg4->AddEntry(hLocalSub,"Sideband-Subtracted Corr.","pe");
    hLocalSub->GetYaxis()->SetRangeUser(fLocalYMin,fLocalYMax);

    hLocalSub->Draw();
    hLocalOrig->Draw("SAME");
    if (bEnablePerformance) DrawAlicePerf(hLocalSub,0.28,0.44,0.33,0.25);
    else DrawWIP(hLocalSub,0.28,0.44,0.33,0.25);
    leg4->Draw("SAME");
    hLocalSub->Draw("SAME");
    PrintCanvas(cSub,Form("SBSubCmp_Bin%d",i));
  }
  // Delta Phi (Near Delta Eta)
  for (Int_t i = 0; i < nResults; i++) {
    TH1D * hLocalSub = fNearEtaDPhiFinal[i];
    TH1D * hLocalOrig = fNearEtaDPhiPi0[i];

    hLocalOrig->SetLineColor(fOrigColor);
    hLocalOrig->SetMarkerColor(fOrigColor);
    hLocalOrig->SetMarkerStyle(fOrigMarker);
    hLocalOrig->SetMarkerSize(fOrigSize);
    hLocalSub->SetLineColor(fSubColor);
    hLocalSub->SetMarkerColor(fSubColor);
    hLocalSub->SetMarkerStyle(fSubMarker);
    hLocalSub->SetMarkerSize(fSubSize);

    Double_t fLocalYMin = min(hLocalSub->GetBinContent(hLocalSub->GetMinimumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMinimumBin()));
    Double_t fLocalYMax = max(hLocalSub->GetBinContent(hLocalSub->GetMaximumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMaximumBin()));
    fLocalYMax += 0.15*(fLocalYMax-fLocalYMin);
    fLocalYMin -= 0.05*(fLocalYMax-fLocalYMin);

    TLegend* leg4= new TLegend(0.38,0.70,0.87,0.87);
    if(fObservable==1)leg4->SetHeader(Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]));
    if(fObservable==2)leg4->SetHeader(Form("%0.1f < #it{p}^{a}_{T} < %0.1f GeV/#it{c}",fObsBins[i],fObsBins[i+1]),"c");
    leg4->AddEntry(hLocalOrig,"Raw #pi^{0}_{Cand.}-h Corr.","pe");
    leg4->AddEntry(hLocalSub,"Sideband-Subtracted Corr.","pe");
    hLocalSub->GetYaxis()->SetRangeUser(fLocalYMin,fLocalYMax);

    hLocalSub->Draw();
    hLocalOrig->Draw("SAME");
    if (bEnablePerformance) DrawAlicePerf(hLocalSub,0.28,0.44,0.33,0.25);
    else DrawWIP(hLocalSub,0.28,0.44,0.33,0.25);
    leg4->Draw("SAME");
    hLocalSub->Draw("SAME");
    PrintCanvas(cSub,Form("SBSubCmpNearDeltaEta_Bin%d",i));
  }


  // Delta Phi (Far Delta Eta)
  for (Int_t i = 0; i < nResults; i++) {
    TH1D * hLocalSub = fFarEtaDPhiFinal[i];
    TH1D * hLocalOrig = fFarEtaDPhiPi0[i];

    hLocalOrig->SetLineColor(fOrigColor);
    hLocalOrig->SetMarkerColor(fOrigColor);
    hLocalOrig->SetMarkerStyle(fOrigMarker);
    hLocalOrig->SetMarkerSize(fOrigSize);
    hLocalSub->SetLineColor(fSubColor);
    hLocalSub->SetMarkerColor(fSubColor);
    hLocalSub->SetMarkerStyle(fSubMarker);
    hLocalSub->SetMarkerSize(fSubSize);

    Double_t fLocalYMin = min(hLocalSub->GetBinContent(hLocalSub->GetMinimumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMinimumBin()));
    Double_t fLocalYMax = max(hLocalSub->GetBinContent(hLocalSub->GetMaximumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMaximumBin()));
    fLocalYMax += 0.15*(fLocalYMax-fLocalYMin);
    fLocalYMin -= 0.05*(fLocalYMax-fLocalYMin);

    TLegend* leg4= new TLegend(0.38,0.70,0.87,0.87);
    if(fObservable==1)leg4->SetHeader(Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]));
    if(fObservable==2)leg4->SetHeader(Form("%0.1f < #it{p}^{a}_{T} < %0.1f GeV/#it{c}",fObsBins[i],fObsBins[i+1]),"c");
    leg4->AddEntry(hLocalOrig,"Raw #pi^{0}_{Cand.}-h Corr.","pe");
    leg4->AddEntry(hLocalSub,"Sideband-Subtracted Corr.","pe");
    hLocalSub->GetYaxis()->SetRangeUser(fLocalYMin,fLocalYMax);

    hLocalSub->Draw();
    hLocalOrig->Draw("SAME");
    if (bEnablePerformance) DrawAlicePerf(hLocalSub,0.28,0.44,0.33,0.25);
    else DrawWIP(hLocalSub,0.28,0.44,0.33,0.25);
    leg4->Draw("SAME");
    hLocalSub->Draw("SAME");
    PrintCanvas(cSub,Form("SBSubCmpFarDeltaEta_Bin%d",i));
  }





  //Delta Eta
  for (Int_t i = 0; i < nResults; i++) {
    TH1D * hLocalSub = fNearSideSubDEtaFinal[i];
    TH1D * hLocalOrig = fNearSideSubDEtaPi0[i];

    Double_t fLocalYMin = min(hLocalSub->GetBinContent(hLocalSub->GetMinimumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMinimumBin()));
    Double_t fLocalYMax = max(hLocalSub->GetBinContent(hLocalSub->GetMaximumBin()),hLocalOrig->GetBinContent(hLocalSub->GetMaximumBin()));
    fLocalYMax += 0.15*(fLocalYMax-fLocalYMin);
    fLocalYMin -= 0.05*(fLocalYMax-fLocalYMin);

    TLegend* leg4= new TLegend(0.58,0.70,0.87,0.87);
    if(fObservable==1)leg4->SetHeader(Form("%0.2f < #it{z}_{T} < %0.2f",fObsBins[i],fObsBins[i+1]));
    if(fObservable==2)leg4->SetHeader(Form("%0.1f < #it{p}^{a}_{T} < %0.1f GeV/#it{c}",fObsBins[i],fObsBins[i+1]),"c");
    leg4->AddEntry(hLocalOrig,"Raw #pi^{0}_{Cand.}-h Corr.","pe");
    leg4->AddEntry(hLocalSub,"Sideband-Subtracted Corr.","pe");
    hLocalSub->GetYaxis()->SetRangeUser(fLocalYMin,fLocalYMax);

    hLocalSub->GetXaxis()->SetLabelSize(0.04);
    hLocalSub->GetYaxis()->SetLabelSize(0.036);
    hLocalSub->GetXaxis()->SetTitleSize(0.05);
    hLocalSub->GetYaxis()->SetTitleSize(0.053);
    hLocalSub->GetXaxis()->SetTitleOffset(0.7);
    hLocalSub->GetYaxis()->SetTitleOffset(0.99);

    hLocalSub->Draw();
    hLocalOrig->Draw("SAME");
    if (bEnablePerformance) DrawAlicePerf(hLocalSub,0.58,0.44,0.33,0.25);
    else DrawWIP(hLocalSub,0.58,0.44,0.33,0.25);
    leg4->Draw("SAME");
    hLocalSub->Draw("SAME");
    PrintCanvas(cSub,Form("SBSubCmpDEta_Bin%d",i));
  }










	cout<<"It is done."<<endl;
}

void TaskSideband::ProcessFlow() {
	Int_t fSBColor[4] = {kAzure,kAzure-2,kAzure-4,kAzure-9};
  Int_t fSBStyle[4] = {kFullSquare,kFullCircle,kFullDiamond,kOpenSquare};

  int nRebinFlow = 8;

  cout<<"Processing Flow"<<endl;




  if (hPtEPAnglePionAcc_Proj_Pion.size() == 0) {
    fprintf(stderr,"Missing Pion EP Angle Hist\n");
    return;
  }


  // Format the flow histograms
  for (int i = 0; i < kNPtBins; i++) {
    hPtEPAnglePionAcc_Proj_Pion[i]->Rebin(nRebinFlow);
    hPtEPAnglePionAcc_Proj_Pion[i]->SetLineColor(kBlack);
    hPtEPAnglePionAcc_Proj_Pion[i]->SetMarkerStyle(kFullSquare);
    //for (int j = 0; j < fNSB; j++) {
    for (int j = 0; j < kNSB; j++) {
      if (!hPtEPAnglePionAcc_Proj_SB[i][j]) {
        fprintf(stderr,"Missing sideband %d projection %d\n",j,i);
        return;
      }
      hPtEPAnglePionAcc_Proj_SB[i][j]->Rebin(nRebinFlow);
      hPtEPAnglePionAcc_Proj_SB[i][j]->SetLineColor(fSBColor[j]);
      hPtEPAnglePionAcc_Proj_SB[i][j]->SetMarkerColor(fSBColor[j]);
      hPtEPAnglePionAcc_Proj_SB[i][j]->SetMarkerStyle(fSBStyle[j]);
    }
  }
  if (hPtEP3AnglePionAcc_Proj_Pion.size() != 0) {
    for (int i = 0; i < kNPtBins; i++) {
      hPtEP3AnglePionAcc_Proj_Pion[i]->Rebin(nRebinFlow);
      hPtEP3AnglePionAcc_Proj_Pion[i]->SetLineColor(kBlack);
      hPtEP3AnglePionAcc_Proj_Pion[i]->SetMarkerStyle(kFullSquare);
      //for (int j = 0; j < fNSB; j++) {
      for (int j = 0; j < kNSB; j++) {
        if (!hPtEP3AnglePionAcc_Proj_SB[i][j]) {
          fprintf(stderr,"Missing sideband %d projection %d\n",j,i);
          return;
        }
        hPtEP3AnglePionAcc_Proj_SB[i][j]->Rebin(nRebinFlow);
        hPtEP3AnglePionAcc_Proj_SB[i][j]->SetLineColor(fSBColor[j]);
        hPtEP3AnglePionAcc_Proj_SB[i][j]->SetMarkerColor(fSBColor[j]);
        hPtEP3AnglePionAcc_Proj_SB[i][j]->SetMarkerStyle(fSBStyle[j]);
      }
    }
  }
  printf("Test 2\n");

  TCanvas * cFlowCanvas = new TCanvas("FlowCanvas","FlowCanvas");
  cFlowCanvas->cd();
  TLegend * legFlowCanvas = new TLegend(0.70,0.70,0.9,0.95); // TLegend(0.70,0.85,0.9,0.95);

  TLegend * leg = new TLegend(0.65,0.67,0.92,0.92);
  if (hPtEPAnglePionAcc_Proj_Pion.size() != 0) { 
    for (int i = 0; i < kNPtBins; i++) {
      leg->Clear();
      hPtEPAnglePionAcc_Proj_Pion[i]->Draw("MIN0");
      double yMax = hPtEPAnglePionAcc_Proj_Pion[i]->GetBinContent(hPtEPAnglePionAcc_Proj_Pion[i]->GetMaximumBin());
      leg->AddEntry(hPtEPAnglePionAcc_Proj_Pion[i],"#pi^{0}_{cand}","lp");
      for (int j = 0; j < kNSB; j++) {
        hPtEPAnglePionAcc_Proj_SB[i][j]->Draw("SAME");
        yMax = fmax(yMax,hPtEPAnglePionAcc_Proj_SB[i][j]->GetBinContent(hPtEPAnglePionAcc_Proj_SB[i][j]->GetMaximumBin()));
        leg->AddEntry(hPtEPAnglePionAcc_Proj_SB[i][j],Form("Sideband %d",j+1),"lp");
      }
      leg->Draw("SAME");
      hPtEPAnglePionAcc_Proj_Pion[i]->GetYaxis()->SetRangeUser(0.,1.2*yMax);
      cFlowCanvas->Print(Form("%s/FlowAll_Obs%d.pdf",fOutputDir.Data(),i)); 
      cFlowCanvas->Print(Form("%s/FlowAll_Obs%d.png",fOutputDir.Data(),i)); 
      cFlowCanvas->Print(Form("%s/CFiles/FlowAll_Obs%d.C",fOutputDir.Data(),i)); 
    }
  }
  if (hPtEP3AnglePionAcc_Proj_Pion.size() != 0) {
    for (int i = 0; i < kNPtBins; i++) {
      leg->Clear();
      hPtEP3AnglePionAcc_Proj_Pion[i]->Draw("MIN0");
      double yMax = hPtEP3AnglePionAcc_Proj_Pion[i]->GetBinContent(hPtEP3AnglePionAcc_Proj_Pion[i]->GetMaximumBin());
      leg->AddEntry(hPtEP3AnglePionAcc_Proj_Pion[i],"#pi^{0}_{cand}","lp");
      for (int j = 0; j < kNSB; j++) {
        hPtEP3AnglePionAcc_Proj_SB[i][j]->Draw("SAME");
        yMax = fmax(yMax,hPtEP3AnglePionAcc_Proj_SB[i][j]->GetBinContent(hPtEP3AnglePionAcc_Proj_SB[i][j]->GetMaximumBin()));
        leg->AddEntry(hPtEP3AnglePionAcc_Proj_SB[i][j],Form("Sideband %d",j+1),"lp");
      }
      leg->Draw("SAME");
      hPtEP3AnglePionAcc_Proj_Pion[i]->GetYaxis()->SetRangeUser(0.,1.2*yMax);
      cFlowCanvas->Print(Form("%s/FlowEP3All_Obs%d.pdf",fOutputDir.Data(),i)); 
      cFlowCanvas->Print(Form("%s/FlowEP3All_Obs%d.png",fOutputDir.Data(),i)); 
      cFlowCanvas->Print(Form("%s/CFiles/FlowEP3All_Obs%d.C",fOutputDir.Data(),i)); 
    }
  printf("Test 3\n");
  }

  printf("Preparing to calculate vn before subtraction\n");

  printf("Normalizing ... \n");
  printf("DEBUG: fNSB = %d\n",fNSB);
  // Normalize the 5 histograms
  // Or their clones?
  for (int i = 0; i < kNPtBins; i++) {
    double scale = hPtEPAnglePionAcc_Proj_Pion[i]->Integral("width");
    if (scale > 0) hPtEPAnglePionAcc_Proj_Pion[i]->Scale(1./scale);

    printf("  Pi0 hist normalization: %f\n",hPtEPAnglePionAcc_Proj_Pion[i]->Integral("width"));
    for (int j = 0; j < kNSB; j++) { // FIXME fNSB is not the right variable here
      scale = hPtEPAnglePionAcc_Proj_SB[i][j]->Integral("width");
      if(scale > 0) { 
        hPtEPAnglePionAcc_Proj_SB[i][j]->Scale(1./scale);
      } else printf("Debug: scale was negative or zero!\n");

      printf("  SB hist normalization: %f\n",hPtEPAnglePionAcc_Proj_SB[i][j]->Integral("width"));

    }
  }
  printf("Done normalizing the EP2 flow histograms\n");
  
  if (hPtEPAnglePionAcc_Proj_Pion.size() != 0) {
    for (int i = 0; i < kNPtBins; i++) {
      leg->Clear();
      hPtEPAnglePionAcc_Proj_Pion[i]->Draw("MIN0");
      double yMax = hPtEPAnglePionAcc_Proj_Pion[i]->GetBinContent(hPtEPAnglePionAcc_Proj_Pion[i]->GetMaximumBin());
      leg->AddEntry(hPtEPAnglePionAcc_Proj_Pion[i],"#pi^{0}_{cand}","lp");
      for (int j = 0; j < kNSB; j++) {
        hPtEPAnglePionAcc_Proj_SB[i][j]->Draw("SAME");
        yMax = fmax(yMax,hPtEPAnglePionAcc_Proj_SB[i][j]->GetBinContent(hPtEPAnglePionAcc_Proj_SB[i][j]->GetMaximumBin()));
        leg->AddEntry(hPtEPAnglePionAcc_Proj_SB[i][j],Form("Sideband %d",j+1),"lp");
      }
      leg->Draw("SAME");
      hPtEPAnglePionAcc_Proj_Pion[i]->GetYaxis()->SetRangeUser(0.,1.2*yMax);
      cFlowCanvas->Print(Form("%s/FlowAllNorm_Obs%d.pdf",fOutputDir.Data(),i)); 
      cFlowCanvas->Print(Form("%s/FlowAllNorm_Obs%d.png",fOutputDir.Data(),i)); 
      cFlowCanvas->Print(Form("%s/CFiles/FlowAllNorm_Obs%d.C",fOutputDir.Data(),i)); 
    }
  }

  // Preparing the EP3 histograms
  if (hPtEP3AnglePionAcc_Proj_Pion.size() != 0) {
    for (int i = 0; i < kNPtBins; i++) {
      double scale = hPtEP3AnglePionAcc_Proj_Pion[i]->Integral("width");
      if (scale > 0) hPtEP3AnglePionAcc_Proj_Pion[i]->Scale(1./scale);

      printf("  Pi0 hist normalization: %f\n",hPtEP3AnglePionAcc_Proj_Pion[i]->Integral("width"));
      for (int j = 0; j < kNSB; j++) { // FIXME fNSB is not the right variable here
        scale = hPtEP3AnglePionAcc_Proj_SB[i][j]->Integral("width");
        if(scale > 0) { 
          hPtEP3AnglePionAcc_Proj_SB[i][j]->Scale(1./scale);
        } else printf("Debug: scale was negative or zero!\n");

        printf("  SB hist normalization: %f\n",hPtEP3AnglePionAcc_Proj_SB[i][j]->Integral("width"));

      }
    }
    printf("Done normalizing the EP3 flow histograms\n");
    if (hPtEP3AnglePionAcc_Proj_Pion.size() != 0) {
      for (int i = 0; i < kNPtBins; i++) {
        leg->Clear();
        hPtEP3AnglePionAcc_Proj_Pion[i]->Draw("MIN0");
        double yMax = hPtEP3AnglePionAcc_Proj_Pion[i]->GetBinContent(hPtEP3AnglePionAcc_Proj_Pion[i]->GetMaximumBin());
        leg->AddEntry(hPtEP3AnglePionAcc_Proj_Pion[i],"#pi^{0}_{cand}","lp");
        for (int j = 0; j < kNSB; j++) {
          hPtEP3AnglePionAcc_Proj_SB[i][j]->Draw("SAME");
          yMax = fmax(yMax,hPtEP3AnglePionAcc_Proj_SB[i][j]->GetBinContent(hPtEP3AnglePionAcc_Proj_SB[i][j]->GetMaximumBin()));
          leg->AddEntry(hPtEP3AnglePionAcc_Proj_SB[i][j],Form("Sideband %d",j+1),"lp");
        }
        leg->Draw("SAME");
        hPtEP3AnglePionAcc_Proj_Pion[i]->GetYaxis()->SetRangeUser(0.,1.2*yMax);
        cFlowCanvas->Print(Form("%s/FlowEP3AllNorm_Obs%d.pdf",fOutputDir.Data(),i)); 
        cFlowCanvas->Print(Form("%s/FlowEP3AllNorm_Obs%d.png",fOutputDir.Data(),i)); 
        cFlowCanvas->Print(Form("%s/CFiles/FlowEP3AllNorm_Obs%d.C",fOutputDir.Data(),i)); 
      }
    }
  }

  // Compute Vn prior to subtraction
  CalculateVnSub(0);
  CalculateV3Sub(0);
  printf("Finished calculating Vn before subtraction\n");

  // FIXME FIXME disabling here in MC
  if (fMCTriggerDistPi0 != 0) return;


  // Build the background estimate using the same subtraction scheme as the sideband correction

  // EP2
  for (int i = 0; i < kNPtBins; i++) {
    printf("Calculating background for V2 Trigger flow.\n");
    // Start with an emptied out copy of the first sideband.
    printf("hPtEPAnglePionAcc_Proj_SB has size %d\n", (int) hPtEPAnglePionAcc_Proj_SB.size());
    printf("hPtEPAnglePionAcc_Proj_SB[i] has size %d\n", (int) hPtEPAnglePionAcc_Proj_SB[i].size());
    TH1D * hPtEPAngleBackgroundEstimate_PtBin = (TH1D *) hPtEPAnglePionAcc_Proj_SB[i][0]->Clone(Form("PtEPAngleBackgroundEstimate_PtBin%d",i+1));
    hPtEPAngleBackgroundEstimate_PtBin->Reset();
    printf("SB hist normalization at start = %f\n",hPtEPAngleBackgroundEstimate_PtBin->Integral("width"));

    for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
      printf("Adding a histogram to background estimate\n");
      printf("    The histogram being added has normalization %f\n",hPtEPAnglePionAcc_Proj_SB[i][j]->Integral("width"));
      hPtEPAngleBackgroundEstimate_PtBin->Add(hPtEPAnglePionAcc_Proj_SB[i][j]);
      printf("    After the addition, the total histogram has normalization %f\n",hPtEPAngleBackgroundEstimate_PtBin->Integral("width"));
    }

    printf("SB hist normalization = %f\n",hPtEPAngleBackgroundEstimate_PtBin->Integral("width"));

    printf("Scaling background prediction flow by 1. / %d\n",fNSB);
    hPtEPAngleBackgroundEstimate_PtBin->Scale(1./(fNSB));
    printf("SB hist normalization = %f\n",hPtEPAngleBackgroundEstimate_PtBin->Integral("width"));

    hPtEPAngleBackgroundEstimate.push_back(hPtEPAngleBackgroundEstimate_PtBin);
  }
  // EP3
  if (hPtEP3AnglePionAcc_Proj_SB[0].size() != 0) {
    printf("Calculating background for V3 Trigger flow.\n");
  //if (hPtEP3AnglePionAcc_Proj_SB[0][0]) {
    for (int i = 0; i < kNPtBins; i++) {
      // Start with an emptied out copy of the first sideband.
      TH1D * hPtEP3AngleBackgroundEstimate_PtBin = (TH1D *) hPtEP3AnglePionAcc_Proj_SB[i][0]->Clone(Form("PtEP3AngleBackgroundEstimate_PtBin%d",i+1));
      hPtEP3AngleBackgroundEstimate_PtBin->Reset();
      printf("SB hist normalization at start = %f\n",hPtEP3AngleBackgroundEstimate_PtBin->Integral("width"));

      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) {
        printf("Adding a histogram to background estimate\n");
        printf("    The histogram being added has normalization %f\n",hPtEP3AnglePionAcc_Proj_SB[i][j]->Integral("width"));
        hPtEP3AngleBackgroundEstimate_PtBin->Add(hPtEP3AnglePionAcc_Proj_SB[i][j]);
        printf("    After the addition, the total histogram has normalization %f\n",hPtEP3AngleBackgroundEstimate_PtBin->Integral("width"));
      }

      printf("SB hist normalization = %f\n",hPtEP3AngleBackgroundEstimate_PtBin->Integral("width"));

      printf("Scaling background prediction flow by 1. / %d\n",fNSB);
      hPtEP3AngleBackgroundEstimate_PtBin->Scale(1./(fNSB));
      printf("SB hist normalization = %f\n",hPtEP3AngleBackgroundEstimate_PtBin->Integral("width"));

      hPtEP3AngleBackgroundEstimate.push_back(hPtEP3AngleBackgroundEstimate_PtBin);
    }
  }

  // Subtract
  printf("About to subtract Pion for EP\n");
  // EP2 
  for (int i = 0; i < kNPtBins; i++) {
    TH1D * hPtEPAnglePionAcc_PionPostSub_PtBin = (TH1D *) hPtEPAnglePionAcc_Proj_Pion[i]->Clone(Form("PtEPAnglePionAcc_PionPostSub_PtBin%d",i+1));

    // FIXME make sure this is correct, but I think it is.
    Int_t iPurityIndex = i; //iPtBin-1;

    Double_t fPurity     = fPurityArray[iPurityIndex];
    printf("Using Purity index %d, found value %f\n",iPurityIndex,fPurity);

    Double_t fPurityScale = 1. - fPurity ; // (1 - purity)

    printf("Scaling background estimate by (1-p) = %f\n",fPurityScale);
    hPtEPAngleBackgroundEstimate[i]->Scale(fPurityScale);
    hPtEPAnglePionAcc_PionPostSub_PtBin->Add(hPtEPAngleBackgroundEstimate[i],-1);
    // scaling result by 1 / purity
    if (fPurity > 0) {
      printf("Normalizing v2 flow histograms by 1/Purity = 1 / %f\n",fPurity);
      hPtEPAnglePionAcc_PionPostSub_PtBin->Scale(1./fPurity);
    } else {
      printf("Not Normalizing flow histogram, as purity <= 0\n");
    }  

    legFlowCanvas->Clear();
    cFlowCanvas->Clear();
    
    hPtEPAnglePionAcc_Proj_Pion[i]->SetMarkerStyle(kOpenSquare);
    hPtEPAnglePionAcc_Proj_Pion[i]->Draw("MIN0");
    legFlowCanvas->AddEntry(hPtEPAnglePionAcc_Proj_Pion[i],"#pi^{0}_{Cand}","lp");
    hPtEPAngleBackgroundEstimate[i]->Draw("SAME");
    legFlowCanvas->AddEntry(hPtEPAngleBackgroundEstimate[i],"Background Estimate","lp");

    hPtEPAnglePionAcc_PionPostSub_PtBin->Draw("SAME");
    legFlowCanvas->AddEntry(hPtEPAnglePionAcc_PionPostSub_PtBin,"#pi^{0} Corrected","lp");

    legFlowCanvas->Draw("SAME");

    cFlowCanvas->Print(Form("%s/FlowPi0_Subtraction_PtBin%d.pdf",fOutputDir.Data(),i+1)); 
    cFlowCanvas->Print(Form("%s/FlowPi0_Subtraction_PtBin%d.png",fOutputDir.Data(),i+1)); 
    cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_Subtraction_PtBin%d.C",fOutputDir.Data(),i+1)); 
    hPtEPAnglePionAcc_PionPostSub.push_back(hPtEPAnglePionAcc_PionPostSub_PtBin);
  }

  // EP3 
  if (hPtEP3AnglePionAcc_Proj_Pion[0]) {
    for (int i = 0; i < kNPtBins; i++) {
      TH1D * hPtEP3AnglePionAcc_PionPostSub_PtBin = (TH1D *) hPtEP3AnglePionAcc_Proj_Pion[i]->Clone(Form("PtEP3AnglePionAcc_PionPostSub_PtBin%d",i+1));

      Int_t iPurityIndex = i; //iPtBin-1;

      Double_t fPurity     = fPurityArray[iPurityIndex];
      printf("Using Purity index %d, found value %f\n",iPurityIndex,fPurity);

      Double_t fPurityScale = 1. - fPurity ; // (1 - purity)

      printf("Scaling background estimate by (1-p) = %f\n",fPurityScale);
      hPtEP3AngleBackgroundEstimate[i]->Scale(fPurityScale);
      hPtEP3AnglePionAcc_PionPostSub_PtBin->Add(hPtEP3AngleBackgroundEstimate[i],-1);
      // scaling result by 1 / purity
      if (fPurity > 0) {
        printf("Normalizing v3 flow histograms by 1/Purity = 1 / %f\n",fPurity);
        hPtEP3AnglePionAcc_PionPostSub_PtBin->Scale(1./fPurity);
      } else {
        printf("Not Normalizing flow histogram, as purity <= 0\n");
      }  

      legFlowCanvas->Clear();
      cFlowCanvas->Clear();
      
      hPtEP3AnglePionAcc_Proj_Pion[i]->SetMarkerStyle(kOpenSquare);
      hPtEP3AnglePionAcc_Proj_Pion[i]->Draw("MIN0");
      legFlowCanvas->AddEntry(hPtEP3AnglePionAcc_Proj_Pion[i],"#pi^{0}_{Cand}","lp");
      hPtEP3AngleBackgroundEstimate[i]->Draw("SAME");
      legFlowCanvas->AddEntry(hPtEP3AngleBackgroundEstimate[i],"Background Estimate","lp");

      hPtEP3AnglePionAcc_PionPostSub_PtBin->Draw("SAME");
      legFlowCanvas->AddEntry(hPtEP3AnglePionAcc_PionPostSub_PtBin,"#pi^{0} Corrected","lp");

      legFlowCanvas->Draw("SAME");

      cFlowCanvas->Print(Form("%s/FlowPi0_EP3_Subtraction_PtBin%d.pdf",fOutputDir.Data(),i+1)); 
      cFlowCanvas->Print(Form("%s/FlowPi0_EP3_Subtraction_PtBin%d.png",fOutputDir.Data(),i+1)); 
      cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_EP3_Subtraction_PtBin%d.C",fOutputDir.Data(),i+1)); 
      hPtEP3AnglePionAcc_PionPostSub.push_back(hPtEP3AnglePionAcc_PionPostSub_PtBin);
    }
  }





  printf("Preparing to calculate v2 after subtraction\n");
  CalculateVnSub(1);
  printf("Preparing to calculate v3 after subtraction\n");
  CalculateV3Sub(1);

  printf("Finished calculating Vn after subtraction\n");

  // Compare Vn before and after subtraction

  TMultiGraph * mgFlow = new TMultiGraph();

  cFlowCanvas->cd();

  Int_t fRawColor = kBlack;
  Int_t fCorrColor = kOrange+1;

  // Settings for sidebands
  //for (int i = 0; i < kNSB; i++) {
  for (int i = 0; i < kNSB; i++) {
    gTriggerFlowSidebands_V2[i]->SetLineColor(fSBColor[i]);
    gTriggerFlowSidebands_V2[i]->SetMarkerColor(fSBColor[i]);
    gTriggerFlowSidebands_V2[i]->SetMarkerStyle(fSBStyle[i]);

    gTriggerFlowSidebands_V3[i]->SetLineColor(fSBColor[i]);
    gTriggerFlowSidebands_V3[i]->SetMarkerColor(fSBColor[i]);
    gTriggerFlowSidebands_V3[i]->SetMarkerStyle(fSBStyle[i]);

    gTriggerFlowSidebands_V4[i]->SetLineColor(fSBColor[i]);
    gTriggerFlowSidebands_V4[i]->SetMarkerColor(fSBColor[i]);
    gTriggerFlowSidebands_V4[i]->SetMarkerStyle(fSBStyle[i]);
  }

  printf("Done setting the Sideband visual settings\n");

  // V2
  cFlowCanvas->Clear();
  legFlowCanvas->Clear();
  gTriggerFlowPreSub_V2->SetLineColor(fRawColor);
  gTriggerFlowPreSub_V2->SetMarkerColor(fRawColor);
  gTriggerFlowPreSub_V2->SetMarkerStyle(kOpenSquare);
  gTriggerFlowPreSub_V2->Draw("ALP");
  legFlowCanvas->AddEntry(gTriggerFlowPreSub_V2,"Raw","lp");
  gTriggerFlowPostSub_V2->SetLineColor(fCorrColor);
  gTriggerFlowPostSub_V2->SetMarkerColor(fCorrColor);
  gTriggerFlowPostSub_V2->SetMarkerStyle(kFullSquare);
  gTriggerFlowPostSub_V2->Draw("LP SAME");
  legFlowCanvas->AddEntry(gTriggerFlowPostSub_V2,"Corrected","lp");
  legFlowCanvas->Draw("SAME");
  cFlowCanvas->Print(Form("%s/FlowPi0_V2_Cmp.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowPi0_V2_Cmp.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_V2_Cmp.C",fOutputDir.Data()));


  cFlowCanvas->Clear();
  mgFlow->Add(gTriggerFlowPreSub_V2,"lp");
  mgFlow->GetXaxis()->SetTitle(gTriggerFlowPreSub_V2->GetXaxis()->GetTitle());
  mgFlow->GetYaxis()->SetTitle(gTriggerFlowPreSub_V2->GetYaxis()->GetTitle());
  mgFlow->Add(gTriggerFlowPostSub_V2,"lp");

  double fMaxVn = fmax(gTriggerFlowPreSub_V2->GetMaximum(),gTriggerFlowPostSub_V2->GetMaximum());
  double fMinVn = fmin(gTriggerFlowPreSub_V2->GetMinimum(),gTriggerFlowPostSub_V2->GetMinimum());

  printf("About to make the Flow V2 plot (kNSB=%d)\n",kNSB);

  printf("gTriggerFlowSidebands_V2 has size %d\n",(int) gTriggerFlowSidebands_V2.size());

  for (int i = 0; i < kNSB; i++) {
    if (fSidebandMask[i]) 
    {
      fMaxVn = fmax(fMaxVn,gTriggerFlowSidebands_V2[i]->GetMaximum());
      fMinVn = fmin(fMinVn,gTriggerFlowSidebands_V2[i]->GetMinimum());
      mgFlow->Add(gTriggerFlowSidebands_V2[i],"lp");
      //gTriggerFlowSidebands_V2[i]->Draw("SAME LP");
      legFlowCanvas->AddEntry(gTriggerFlowSidebands_V2[i],Form("Sideband %d",i+1),"lp");
    }
  }

  printf("Done doing the thing\n");

  // The multigraph should be taking care of this
 // mgFlow->GetYaxis()->SetRangeUser(-0.2,0.2);

  mgFlow->Draw("ALP");
  //fMaxVn += 0.1 * (fMaxVn - fMinVn);
  mgFlow->GetYaxis()->SetRangeUser(fMinVn,fMaxVn);
  legFlowCanvas->Draw("SAME");

  cFlowCanvas->Print(Form("%s/FlowPi0_V2_SB_Cmp.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowPi0_V2_SB_Cmp.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_V2_SB_Cmp.C",fOutputDir.Data())); 
 



  // V4
  cFlowCanvas->Clear();
  legFlowCanvas->Clear();
  gTriggerFlowPreSub_V4->SetLineColor(fRawColor);
  gTriggerFlowPreSub_V4->SetMarkerColor(fRawColor);
  gTriggerFlowPreSub_V4->SetMarkerStyle(kOpenSquare);
  gTriggerFlowPreSub_V4->Draw("ALP");
  legFlowCanvas->AddEntry(gTriggerFlowPreSub_V4,"Raw","lp");
  gTriggerFlowPostSub_V4->SetLineColor(fCorrColor);
  gTriggerFlowPostSub_V4->SetMarkerColor(fCorrColor);
  gTriggerFlowPostSub_V4->SetMarkerStyle(kFullSquare);
  gTriggerFlowPostSub_V4->Draw("LP SAME");
  legFlowCanvas->AddEntry(gTriggerFlowPostSub_V4,"Corrected","lp");
  legFlowCanvas->Draw("SAME");
  cFlowCanvas->Print(Form("%s/FlowPi0_V4_Cmp.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowPi0_V4_Cmp.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_V4_Cmp.C",fOutputDir.Data())); 
  // V6
  cFlowCanvas->Clear();
  legFlowCanvas->Clear();
  gTriggerFlowPreSub_V6->SetLineColor(fRawColor);
  gTriggerFlowPreSub_V6->SetMarkerColor(fRawColor);
  gTriggerFlowPreSub_V6->SetMarkerStyle(kOpenSquare);
  gTriggerFlowPreSub_V6->Draw("ALP");
  legFlowCanvas->AddEntry(gTriggerFlowPreSub_V6,"Raw","lp");
  gTriggerFlowPostSub_V6->SetLineColor(fCorrColor);
  gTriggerFlowPostSub_V6->SetMarkerColor(fCorrColor);
  gTriggerFlowPostSub_V6->SetMarkerStyle(kFullSquare);
  gTriggerFlowPostSub_V6->Draw("LP SAME");
  legFlowCanvas->AddEntry(gTriggerFlowPostSub_V6,"Corrected","lp");
  legFlowCanvas->Draw("SAME");
  cFlowCanvas->Print(Form("%s/FlowPi0_V6_Cmp.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowPi0_V6_Cmp.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_V6_Cmp.C",fOutputDir.Data())); 

  // V3
  cFlowCanvas->Clear();
  legFlowCanvas->Clear();
  gTriggerFlowPreSub_V3->SetLineColor(fRawColor);
  gTriggerFlowPreSub_V3->SetMarkerColor(fRawColor);
  gTriggerFlowPreSub_V3->SetMarkerStyle(kOpenSquare);
  gTriggerFlowPreSub_V3->Draw("ALP");


  legFlowCanvas->AddEntry(gTriggerFlowPreSub_V3,"Raw","lp");
  gTriggerFlowPostSub_V3->SetLineColor(fCorrColor);
  gTriggerFlowPostSub_V3->SetMarkerColor(fCorrColor);
  gTriggerFlowPostSub_V3->SetMarkerStyle(kFullSquare);
  gTriggerFlowPostSub_V3->Draw("LP SAME");
  legFlowCanvas->AddEntry(gTriggerFlowPostSub_V3,"Corrected","lp");
  legFlowCanvas->Draw("SAME");
  cFlowCanvas->Print(Form("%s/FlowEP3Pi0_V3_Cmp.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowEP3Pi0_V3_Cmp.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowEP3Pi0_V3_Cmp.C",fOutputDir.Data())); 

  // Could also show the sideband vn, since I think i have those

//  cFlowCanvas->Clear();
//  legFlowCanvas->Clear();
//  gTriggerFlowPreSub_V3->Draw("ALP");
//  gTriggerFlowPreSub_V3->Draw("ALP");

//  for (int i = 0; i < kNSB; i++) {
//    gTriggerFlowSidebands_V3[i]->SetLineColor(fSBColor[i]);
//    gTriggerFlowSidebands_V3[i]->SetMarkerColor(fSBColor[i]);
//    gTriggerFlowSidebands_V3[i]->SetMarkerStyle(fSBStyle[i]);
//  }

  TMultiGraph * mgFlowV3 = new TMultiGraph();

  cFlowCanvas->Clear();
  mgFlowV3->Add(gTriggerFlowPreSub_V3,"lp");
  mgFlowV3->GetXaxis()->SetTitle(gTriggerFlowPreSub_V3->GetXaxis()->GetTitle());
  mgFlowV3->GetYaxis()->SetTitle(gTriggerFlowPreSub_V3->GetYaxis()->GetTitle());
  mgFlowV3->Add(gTriggerFlowPostSub_V3,"lp");


  printf("About to make the Flow V3 plot\n");

  fMaxVn = 0;
  fMinVn = 0;

  fMaxVn = fmax(gTriggerFlowPreSub_V3->GetMaximum(),gTriggerFlowPostSub_V3->GetMaximum());
  fMinVn = fmin(gTriggerFlowPreSub_V3->GetMinimum(),gTriggerFlowPostSub_V3->GetMinimum());

  for (int i = 0; i < kNSB; i++) {
    if (fSidebandMask[i]) 
    {
      //gTriggerFlowSidebands_V3[i]->Draw("SAME LP");
      mgFlowV3->Add(gTriggerFlowSidebands_V3[i]);
      legFlowCanvas->AddEntry(gTriggerFlowSidebands_V3[i],Form("Sideband %d",i+1),"lp");
    }
  }

  mgFlowV3->Draw("ALP");
  mgFlowV3->GetYaxis()->SetRangeUser(fMinVn,fMaxVn);
  legFlowCanvas->Draw("SAME");

  printf("Done making the Flow V3 plot\n");

  cFlowCanvas->Print(Form("%s/FlowEP3Pi0_V3_SB_Cmp.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowEP3Pi0_V3_SB_Cmp.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowEP3Pi0_V3_SB_Cmp.C",fOutputDir.Data())); 



  printf("Done with ProcessFlow\n");

}


// Note: these are not corrected for Event Plane Resolution
// May copy this to phase 4, where the EPR is easily applied
void TaskSideband::CalculateVnSub(int iPostSub = 0) {
  if (iPostSub == 0) printf("Calculating Vn Before Subtraction\n");
  if (iPostSub == 1) printf("Calculating Vn After Subtraction\n");

  // FIXME add a switch for turning off v6.

  TH1D * fHist = 0;
  TF1 * fFlowFunction = 0;

  TGraphErrors * gTriggerFlow_V2 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V2->SetName("TriggerFlow_V2");
  TGraphErrors * gTriggerFlow_V4 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V4->SetName("TriggerFlow_V4");
  TGraphErrors * gTriggerFlow_V6 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V6->SetName("TriggerFlow_V6");

  if (iPostSub == 0) {

    printf("Debug: fNSB = %d\n",fNSB);

    for (int j = 0; j < kNSB; j++) { // Do this for all Sidebands, even if not included in final result
      TGraphErrors * gTriggerFlowSB_V2 = new TGraphErrors(kNPtBins);
      gTriggerFlowSB_V2->SetName(Form("TriggerFlowSB_V2_SB%d",j));
      TGraphErrors * gTriggerFlowSB_V4 = new TGraphErrors(kNPtBins);
      gTriggerFlowSB_V4->SetName(Form("TriggerFlowSB_V4_SB%d",j));
      TGraphErrors * gTriggerFlowSB_V6 = new TGraphErrors(kNPtBins);
      gTriggerFlowSB_V6->SetName(Form("TriggerFlowSB_V6_SB%d",j));
      gTriggerFlowSidebands_V2.push_back(gTriggerFlowSB_V2);
      gTriggerFlowSidebands_V4.push_back(gTriggerFlowSB_V4);
      gTriggerFlowSidebands_V6.push_back(gTriggerFlowSB_V6);
    }
  }
  printf("gTriggerFlowSidebands_V2 has size %d\n", (int) gTriggerFlowSidebands_V2.size());

  for (int i = 0; i < kNPtBins; i++) {
    double fMinPt = fPtBins[i];
    double fMaxPt = fPtBins[i+1];

    if (iPostSub == 0) {
      fHist = hPtEPAnglePionAcc_Proj_Pion[i];
    } else {
      fHist = hPtEPAnglePionAcc_PionPostSub[i];
    }
    if (fHist == 0) {
      fprintf(stderr,"Missing Pion vs EP2 histogram\n");
    }
    double average_value = fHist->Integral("width") / TMath::PiOver2();


    fFlowFunction = new TF1(Form("TriggerFlowFit_ObsBin%d",i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x)+2*[3]*TMath::Cos(6*x))",fHist->GetXaxis()->GetXmin(),fHist->GetXaxis()->GetXmax());
    fFlowFunction->SetParName(0,"B");
    for (int j = 1; j <= 3; j++) fFlowFunction->SetParName(j,Form("#tilde{v}_{%d}",j*2));
//    fFlowFunction->SetParName(2,"#tilde{v}_{4}");
 //   fFlowFunction->SetParName(3,"#tilde{v}_{6}");
    fFlowFunction->SetParameter(0,average_value);
    fFlowFunction->SetParameter(1,0.1);
    fFlowFunction->SetParameter(2,0.05);
    fFlowFunction->SetParameter(3,0.01);


    if (iPostSub == 0) {
      fFlowFunction->SetName(Form("%s_PreSub",fFlowFunction->GetName()));
    } else {
      fFlowFunction->SetName(Form("%s_PostSub",fFlowFunction->GetName()));
    }

    for (int j = 1; j <= 3; j++) fFlowFunction->SetParLimits(j,-0.5,0.5);

    if (!bEnableV6) fFlowFunction->FixParameter(3,0.0);

    fHist->Fit(fFlowFunction,"Q");

    gTriggerFlow_V2->SetPoint(i,(fMinPt+fMaxPt)/2.,fFlowFunction->GetParameter(1));
    gTriggerFlow_V4->SetPoint(i,(fMinPt+fMaxPt)/2.,fFlowFunction->GetParameter(2));
    gTriggerFlow_V6->SetPoint(i,(fMinPt+fMaxPt)/2.,fFlowFunction->GetParameter(3));

    gTriggerFlow_V2->SetPointError(i,(fMaxPt-fMinPt)/2.,fFlowFunction->GetParError(1));
    gTriggerFlow_V4->SetPointError(i,(fMaxPt-fMinPt)/2.,fFlowFunction->GetParError(2));
    gTriggerFlow_V6->SetPointError(i,(fMaxPt-fMinPt)/2.,fFlowFunction->GetParError(3));

    if (iPostSub == 0) { // For pre sub only, do the sidebands
      
      for (int j = 0; j < fNSB; j++) {

 //       TF1 * fSBFlowFunction = (TF1 *) fFlowFunction->Clone(Form("%s_SB%d",fFlowFunction->GetName(),j));
        TF1 * fSBFlowFunction = new TF1(Form("TriggerFlowFit_SB%d_ObsBin%d",j,i),"[0]*(1+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(4*x)+2*[3]*TMath::Cos(6*x))",fHist->GetXaxis()->GetXmin(),fHist->GetXaxis()->GetXmax());

        for (int j = 1; j <= 3; j++) fSBFlowFunction->SetParName(j,Form("#tilde{v}_{%d}",j*2));
        fSBFlowFunction->SetParameter(0,0.8);
        fSBFlowFunction->SetParameter(1,0.1);
        fSBFlowFunction->SetParameter(2,0.05);
        fSBFlowFunction->SetParameter(3,0.01);
        for (int j = 1; j <= 3; j++) fSBFlowFunction->SetParLimits(j,-0.5,0.5);
        if (!bEnableV6) fSBFlowFunction->FixParameter(3,0.0);

        TH1D * fSBFlowHist = hPtEPAnglePionAcc_Proj_SB[i][j];

        fSBFlowHist->Fit(fSBFlowFunction,"Q");

        gTriggerFlowSidebands_V2[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(1));
        gTriggerFlowSidebands_V4[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(2));
        gTriggerFlowSidebands_V6[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(3));
        gTriggerFlowSidebands_V2[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(1));
        gTriggerFlowSidebands_V4[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(2));
        gTriggerFlowSidebands_V6[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(3));
      }
    }
  }

  if (iPostSub == 0) {
    gTriggerFlow_V2->SetName("TriggerFlowPreSub_V2");
    gTriggerFlow_V4->SetName("TriggerFlowPreSub_V4");
    gTriggerFlow_V6->SetName("TriggerFlowPreSub_V6");
    gTriggerFlow_V2->SetTitle("Raw Trigger #tilde{v}_{2};p_{T} (GeV/#it{c});#tilde{v}_{2}");
    gTriggerFlow_V4->SetTitle("Raw Trigger #tilde{v}_{4};p_{T} (GeV/#it{c});#tilde{v}_{4}");
    gTriggerFlow_V6->SetTitle("Raw Trigger #tilde{v}_{6};p_{T} (GeV/#it{c});#tilde{v}_{6}");


    gTriggerFlowPreSub_V2 = gTriggerFlow_V2;
    gTriggerFlowPreSub_V4 = gTriggerFlow_V4;
    gTriggerFlowPreSub_V6 = gTriggerFlow_V6;
  } else {
    gTriggerFlow_V2->SetName("TriggerFlowPostSub_V2");
    gTriggerFlow_V4->SetName("TriggerFlowPostSub_V4");
    gTriggerFlow_V6->SetName("TriggerFlowPostSub_V6");
    gTriggerFlow_V2->SetTitle("Corrected Trigger #tilde{v}_{2};p_{T} (GeV/#it{c});#tilde{v}_{2}");
    gTriggerFlow_V4->SetTitle("Corrected Trigger #tilde{v}_{4};p_{T} (GeV/#it{c});#tilde{v}_{4}");
    gTriggerFlow_V6->SetTitle("Corrected Trigger #tilde{v}_{6};p_{T} (GeV/#it{c});#tilde{v}_{6}");

    gTriggerFlowPostSub_V2 = gTriggerFlow_V2;
    gTriggerFlowPostSub_V4 = gTriggerFlow_V4;
    gTriggerFlowPostSub_V6 = gTriggerFlow_V6;
  }
}



// Note: these are not corrected for Event Plane Resolution
// Just For V3. This is a separate function so it can be easily run without
// many checks for the existing of v3 histograms
void TaskSideband::CalculateV3Sub(int iPostSub = 0) {
  if (iPostSub == 0) printf("Calculating V3 Before Subtraction\n");
  if (iPostSub == 1) printf("Calculating V3 After Subtraction\n");

  TH1D * fHist = 0;
  TF1 * fFlowFunction = 0;

  TGraphErrors * gTriggerFlow_V3 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V3->SetName("TriggerFlow_V3");

  for (int j = 0; j < kNSB; j++) {
    TGraphErrors * gTriggerFlowSB_V3 = new TGraphErrors(kNPtBins);
    gTriggerFlowSB_V3->SetName(Form("TriggerFlowSB_V3_SB%d",j));
    gTriggerFlowSidebands_V3.push_back(gTriggerFlowSB_V3);
  }

  printf("test1\n");

  for (int i = 0; i < kNPtBins; i++) {
    double fMinPt = fPtBins[i];
    double fMaxPt = fPtBins[i+1];

    if (iPostSub == 0) {
      printf("hPtEP3AnglePionAcc_Proj_Pion has size %d\n",(int) hPtEP3AnglePionAcc_Proj_Pion.size());
      if (hPtEP3AnglePionAcc_Proj_Pion.size() == 0) return;
      fHist = hPtEP3AnglePionAcc_Proj_Pion[i];
    } else {
      printf("hPtEP3AnglePionAcc_PionPostSub has size %d\n", (int) hPtEP3AnglePionAcc_PionPostSub.size());
      if (hPtEP3AnglePionAcc_PionPostSub.size() == 0) return;
      fHist = hPtEP3AnglePionAcc_PionPostSub[i];
    }
    if (!fHist) {
      printf("Missing the EP3 Trigger Histogram\n");
      return;
    }
    double average_value = fHist->Integral("width") / TMath::PiOver2();

    fFlowFunction = new TF1(Form("TriggerFlowV3Fit_ObsBin%d",i),"[0]*(1+2*[1]*TMath::Cos(3*x))",fHist->GetXaxis()->GetXmin(),fHist->GetXaxis()->GetXmax());
    fFlowFunction->SetParName(0,"B");
    fFlowFunction->SetParName(1,"#tilde{v}_{3}");
   // for (int j = 1; j <= 3; j++) fFlowFunction->SetParName(j,Form("#tilde{v}_{%d}",j*2));
//    fFlowFunction->SetParName(2,"#tilde{v}_{4}");
 //   fFlowFunction->SetParName(3,"#tilde{v}_{6}");
    fFlowFunction->SetParameter(0,average_value);
    fFlowFunction->SetParameter(1,0.1);

    if (iPostSub == 0) {
      fFlowFunction->SetName(Form("%s_PreSub",fFlowFunction->GetName()));
    } else {
      fFlowFunction->SetName(Form("%s_PostSub",fFlowFunction->GetName()));
    }

    fFlowFunction->SetParLimits(1,-0.5,0.5);
    //for (int j = 1; j <= 3; j++) fFlowFunction->SetParLimits(j,-0.5,0.5);

    printf("Trying the V3 fit on histogram %s (%s)\n",fHist->GetName(),fHist->GetTitle());

    fHist->Fit(fFlowFunction,"Q");

    gTriggerFlow_V3->SetPoint(i,(fMinPt+fMaxPt)/2.,fFlowFunction->GetParameter(1));

    gTriggerFlow_V3->SetPointError(i,(fMaxPt-fMinPt)/2.,fFlowFunction->GetParError(1));


    if (iPostSub == 0) { // For pre sub only, do the sidebands
      
      printf("There are %d sideband histograms to analyze here\n",(int) hPtEP3AnglePionAcc_Proj_SB[i].size());

      //for (int j = 0; j < fNSB; j++) {
      for (int j = 0; j < kNSB; j++) {

 //       TF1 * fSBFlowFunction = (TF1 *) fFlowFunction->Clone(Form("%s_SB%d",fFlowFunction->GetName(),j));
        TF1 * fSBFlowFunction = new TF1(Form("TriggerFlowV3Fit_SB%d_ObsBin%d",j,i),"[0]*(1+2*[1]*TMath::Cos(3*x))",fHist->GetXaxis()->GetXmin(),fHist->GetXaxis()->GetXmax());

        //for (int j = 1; j <= 3; j++) fSBFlowFunction->SetParName(j,Form("#tilde{v}_{%d}",j*2));
        fSBFlowFunction->SetParameter(0,0.8);
        fSBFlowFunction->SetParName(0,"B");
        fSBFlowFunction->SetParameter(1,0.1);
        fSBFlowFunction->SetParName(1,"#tilde{v}_{3}");
        fSBFlowFunction->SetParLimits(1,-0.5,0.5);

        //TH1D * fSBFlowHist = hPtEPAnglePionAcc_Proj_SB[i][j];
        printf(" hPtEPAnglePionAcc_Proj_SB[i] has size %d\n",(int) hPtEPAnglePionAcc_Proj_SB[i].size());

        TH1D * fSBFlowHist = hPtEP3AnglePionAcc_Proj_SB[i][j];

        printf("Trying the V3 fit on histogram %s (%s)\n",fSBFlowHist->GetName(),fSBFlowHist->GetTitle());

        fSBFlowHist->Fit(fSBFlowFunction,"Q");
        printf(" gTriggerFlowSidebands_V3 has size %d\n", (int) gTriggerFlowSidebands_V3.size());
        gTriggerFlowSidebands_V3[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(1));
        gTriggerFlowSidebands_V3[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(1));
      }

      printf("Finished the V3 fits for pt bin %d\n",i);
    }
  }

  printf("test2\n");

  if (iPostSub == 0) {
    gTriggerFlow_V3->SetName("TriggerFlowPreSub_V3");
    gTriggerFlow_V3->SetTitle("Raw Trigger #tilde{v}_{3};p_{T} (GeV/#it{c});#tilde{v}_{3}");

    gTriggerFlowPreSub_V3 = gTriggerFlow_V3;
  }
  else {
    gTriggerFlow_V3->SetName("TriggerFlowPostSub_V3");
    gTriggerFlow_V3->SetTitle("Corrected Trigger #tilde{v}_{3};p_{T} (GeV/#it{c});#tilde{v}_{3}");
    gTriggerFlowPostSub_V3 = gTriggerFlow_V3;
  }

}

void TaskSideband::SaveResults() {
	cout<<"Saving Results"<<endl;

  TFile * fOutputFile = TFile::Open(Form("output/%s",fOutputFileName.Data()),"RECREATE");
	if (!fOutputFile) return;

  if (fHistEventHash) fOutputFile->Add(fHistEventHash);
  if (VariableInfo) fOutputFile->Add(VariableInfo);

  if (fTriggerPt) fOutputFile->Add(fTriggerPt); 
  if (fTriggerPtWithinEPBin) fOutputFile->Add(fTriggerPtWithinEPBin); 

  if (hHistTrackPsiEPPtCent)  fOutputFile->Add(hHistTrackPsiEPPtCent);
  if (hHistTrackPsiEP3PtCent) fOutputFile->Add(hHistTrackPsiEP3PtCent);
  if (hHistTrackPsiEP4PtCent) fOutputFile->Add(hHistTrackPsiEP4PtCent);

  if (fTrackPtProjectionSE) fOutputFile->Add(fTrackPtProjectionSE);
  if (fTrackPtProjectionME) fOutputFile->Add(fTrackPtProjectionME);
  if (fTrackPtFromTrackPsi) fOutputFile->Add(fTrackPtFromTrackPsi);

	for (Int_t i = 0; i < (Int_t) fMassVsIntegral_Full.size(); i++) fOutputFile->Add(fMassVsIntegral_Full[i]);

	for (Int_t i = 0; i < (Int_t) fFullPredBkgDPhi.size(); i++) fOutputFile->Add(fFullPredBkgDPhi[i]);
	for (Int_t i = 0; i < (Int_t) fNearEtaPredBkgDPhi.size(); i++) fOutputFile->Add(fNearEtaPredBkgDPhi[i]);
	for (Int_t i = 0; i < (Int_t) fFarEtaPredBkgDPhi.size(); i++) fOutputFile->Add(fFarEtaPredBkgDPhi[i]);

//	for (Int_t i = 0; i < (Int_t) fFullDPhiPi0.size(); i++) fOutputFile->Add(fFullDPhiFinal[i]);;
	for (Int_t i = 0; i < (Int_t) fFullDPhiFinal.size(); i++) fOutputFile->Add(fFullDPhiFinal[i]);;
	for (Int_t i = 0; i < (Int_t) fNearEtaDPhiFinal.size(); i++) fOutputFile->Add(fNearEtaDPhiFinal[i]);;
	for (Int_t i = 0; i < (Int_t) fFarEtaDPhiFinal.size(); i++) fOutputFile->Add(fFarEtaDPhiFinal[i]);;

  // Delta Eta
  for (Int_t i = 0; i < (Int_t) fNearSideSubPredBkgDEta.size(); i++) fOutputFile->Add(fNearSideSubPredBkgDEta[i]);
  for (Int_t i = 0; i < (Int_t) fAwaySideSubPredBkgDEta.size(); i++) fOutputFile->Add(fAwaySideSubPredBkgDEta[i]);

  for (Int_t i = 0; i < (Int_t) fNearSideSubDEtaFinal.size(); i++) fOutputFile->Add(fNearSideSubDEtaFinal[i]);
  for (Int_t i = 0; i < (Int_t) fAwaySideSubDEtaFinal.size(); i++) fOutputFile->Add(fAwaySideSubDEtaFinal[i]);


  if (gTrigger_Bv) fOutputFile->Add(gTrigger_Bv);
  if (gTrigger_V2) fOutputFile->Add(gTrigger_V2);
  if (gTrigger_V4) fOutputFile->Add(gTrigger_V4);
  if (gTrigger_V6) fOutputFile->Add(gTrigger_V6);

  if (gTrack_Bv) fOutputFile->Add(gTrack_Bv);
  if (gTrack_V2) fOutputFile->Add(gTrack_V2);
  if (gTrack_V4) fOutputFile->Add(gTrack_V4);
  if (gTrack_V6) fOutputFile->Add(gTrack_V6);

  if (gTriggerFlowPreSub_V2) fOutputFile->Add(gTriggerFlowPreSub_V2);
  if (gTriggerFlowPreSub_V3) fOutputFile->Add(gTriggerFlowPreSub_V3);
  if (gTriggerFlowPreSub_V4) fOutputFile->Add(gTriggerFlowPreSub_V4);
  if (gTriggerFlowPreSub_V6) fOutputFile->Add(gTriggerFlowPreSub_V6);

  for (auto gTriggerFlowSidebands_V2_Indiv : gTriggerFlowSidebands_V2) fOutputFile->Add(gTriggerFlowSidebands_V2_Indiv);
  for (auto gTriggerFlowSidebands_V4_Indiv : gTriggerFlowSidebands_V4) fOutputFile->Add(gTriggerFlowSidebands_V4_Indiv);
  for (auto gTriggerFlowSidebands_V6_Indiv : gTriggerFlowSidebands_V6) fOutputFile->Add(gTriggerFlowSidebands_V6_Indiv);

  if (gTriggerFlowPostSub_V2) fOutputFile->Add(gTriggerFlowPostSub_V2);
  if (gTriggerFlowPostSub_V3) fOutputFile->Add(gTriggerFlowPostSub_V3);
  if (gTriggerFlowPostSub_V4) fOutputFile->Add(gTriggerFlowPostSub_V4);
  if (gTriggerFlowPostSub_V6) fOutputFile->Add(gTriggerFlowPostSub_V6);

  for (auto hPtEPAnglePionAcc_Proj_Pion_Indiv : hPtEPAnglePionAcc_Proj_Pion) fOutputFile->Add(hPtEPAnglePionAcc_Proj_Pion_Indiv);
  for (int i = 0; i < (int) hPtEPAnglePionAcc_Proj_SB.size(); i++) {
    for (auto hPtEPAnglePionAcc_Proj_SB_Indiv : hPtEPAnglePionAcc_Proj_SB[i]) {
      fOutputFile->Add(hPtEPAnglePionAcc_Proj_SB_Indiv);
    }
  }
  for (auto hPtEPAnglePionAcc_Proj_PionPostSub_Indiv : hPtEPAnglePionAcc_PionPostSub) fOutputFile->Add(hPtEPAnglePionAcc_Proj_PionPostSub_Indiv);

  for (auto hIndiv : hPtEP3AnglePionAcc_PionPostSub) fOutputFile->Add(hIndiv);
  for (auto hIndiv : hPtEP4AnglePionAcc_PionPostSub) fOutputFile->Add(hIndiv);

  cout<<"Finished adding things to outfile"<<endl;

  fOutputFile->Write();
  cout<<"Wrote file"<<endl;
  fOutputFile->Close();
  cout<<"Closed file"<<endl;
}


void TaskSideband::Debug(Int_t input) {
	cout<<"Debug stage "<<input<<endl;	
	cout<<"NObs bins = "<<nObsBins<<endl;
	cout<<"fObservable = "<<fObservable<<endl;

	cout<<"ObsBins: {";
	for (Int_t i = 0; i <= nObsBins; i++) {
		cout<<fObsBins[i]; if (i != nObsBins) cout<<",";
	}
	cout<<"}"<<endl;

  TCanvas * cDebugCanvas = new TCanvas("DebugCanvas","DebugCanvas");
  cDebugCanvas->Divide(3,3);
  TLegend * legend = new TLegend(0.1,0.1,0.9,0.9);
  switch(input) {
    default:
    case 1:
    case 2:
      for (int i = 0; i < nObsBins; i++) {
        cDebugCanvas->cd(i+1);
        fFullDPhiPi0[i]->Draw();
        fNearEtaDPhiPi0[i]->SetMarkerStyle(kOpenCircle);
        fFarEtaDPhiPi0[i]->SetMarkerStyle(kOpenSquare);


        fNearEtaDPhiPi0[i]->Draw("SAME");
        fFarEtaDPhiPi0[i]->Draw("SAME");
      }
      legend->AddEntry(fFullDPhiPi0[0],"Full Region","lp");
      legend->AddEntry(fNearEtaDPhiPi0[0],"Near #Delta#eta Region","lp");
      legend->AddEntry(fFarEtaDPhiPi0[0],"Far #Delta#eta Region","lp");


    break;
    case 3:
    case 4:

      for (int i = 0; i < nObsBins; i++) {
        cDebugCanvas->cd(i+1);
        fFullDPhiFinal[i]->Draw();
        fNearEtaDPhiFinal[i]->Draw("SAME");
        fFarEtaDPhiFinal[i]->Draw("SAME");
      }
      legend->AddEntry(fFullDPhiFinal[0],"Full Region","lp");
      legend->AddEntry(fNearEtaDPhiFinal[0],"Near #Delta#eta Region","lp");
      legend->AddEntry(fFarEtaDPhiFinal[0],"Far #Delta#eta Region","lp");

    break;
  }
  cDebugCanvas->cd(3*3);
  legend->Draw();
  cDebugCanvas->Print(Form("%s/Debug_%d.pdf",fOutputDir.Data(),input));

  delete cDebugCanvas;
}

void TaskSideband::Run() {
	cout<<"Beginning Sideband Subtraction Task"<<endl;

  gErrorIgnoreLevel = kWarning;	

	LoadPurity();
	LoadHistograms();
	InitArrays();
	Debug(1);

	MassAnalysis();
	ProduceSidebandFigure();
  printf("About to ProduceSidebandComparison\n");
  ProduceSidebandComparison();
  printf("Done with ProduceSidebandComparison\n");
//	SimpleNormFit();
	ProduceBackground();

	PlotBkgAndSignal();

  Debug(2);

	Subtract();

	Debug(4);

  ProcessFlow();

	SaveResults();

  // Debugging masks
  if (kNSB == 4) {
    printf("  SidebandFitMask = %d %d %d %d\n",fSidebandFitMask[0],fSidebandFitMask[1],fSidebandFitMask[2],fSidebandFitMask[3]);
    printf("  SidebandMask = %d %d %d %d\n",fSidebandMask[0],fSidebandMask[1],fSidebandMask[2],fSidebandMask[3]);
  } else if (kNSB == 3) {
    printf("  SidebandFitMask = %d %d %d\n",fSidebandFitMask[0],fSidebandFitMask[1],fSidebandFitMask[2]);
    printf("  SidebandMask = %d %d %d\n",fSidebandMask[0],fSidebandMask[1],fSidebandMask[2]);
  }
}



