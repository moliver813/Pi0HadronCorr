#if defined(__CINT__)
#define _SYS_TYPES_H_
#endif


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
	fDebugLevel = 0;
	
	SetStyle();
}

void TaskSideband::SetStyle() {
//	gStyle->SetCanvasColor(kBlack);
//	gStyle->SetAxisColor(0);
  TGaxis::SetMaxDigits(2);

	gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
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
	canvas->Print(Form("%s/%s.eps",fOutputDir.Data(),name.Data()));
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
	//Pi0YieldTotalRatio
	Pi0YieldTotalRatio = (TGraphErrors *) fPi0PurityFile->Get("Pi0YieldTotalRatio");
	if (!Pi0YieldTotalRatio) {
		fprintf(stderr,"Error: Could not find Pi0YieldTotalRatio in file %s\n",fPi0PurityFile->GetName());
		return;
	} 
	// Throw out the points we don't need (3-4,4-5), and above 17
  // Also build the Purity Array that will be used

  for (int i = 0; i < Pi0YieldTotalRatio->GetN(); i++) {
    Double_t localX = Pi0YieldTotalRatio->GetX()[i];
    if (localX < 5. || localX > 17.) { 
      printf("Removing point %d with p_T value %f\n",i,localX);
      Pi0YieldTotalRatio->RemovePoint(i);
      i--;
    } else {
      double fPurity     = Pi0YieldTotalRatio->GetY()[i];
      double fPurity_Err = Pi0YieldTotalRatio->GetEY()[i];
      switch (iPurityChoice) {
        case 2:
          fPurity = fPurity - fPurity_Err;
          break;
        case 0:
          fPurity = 0;
          break;
        default:
        case 1:
          break;
      }
      fPurityArray.push_back(fPurity);
      fPurityArray_Err.push_back(fPurity_Err);
    }
  }


/*	Pi0YieldTotalRatio->RemovePoint(0);
	Pi0YieldTotalRatio->RemovePoint(0);
	Pi0YieldTotalRatio->RemovePoint(7);
	Pi0YieldTotalRatio->RemovePoint(6);
*/
	TCanvas * cPurity = new TCanvas("cPurity","cPurity");
	Pi0YieldTotalRatio->Draw();
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

  if (!gTrigger_Bv || !gTrack_Bv) {
    fprintf(stderr,"Missing some flow graphs from phase1\n");
  } else {
    printf("Successfully loaded the flow graphs from phase1\n");
  }

}

void TaskSideband::LoadHistograms() {

 // Loading our observable settings:
  VariableInfo = (TH1D *) fPi0CorrFile->Get("VariableInfo");
  if (VariableInfo) {
    fObservable = VariableInfo->GetBinContent(1);
  } else {
    cout<<"No Variable Input TH1D found in Input File "<<fPi0CorrFile->GetName()<<endl;
    cout<<"Using default values for Observable Info and name."<<endl;
    fObservable=0; // the P_t
  }
  if (fObservable==0)      fObservableName = "#it{p}_{T}";
  else if (fObservable==1) fObservableName = "z_{T}";
  else if (fObservable==2) fObservableName = "#xi";

  fTriggerPt = (TH1D *) fPi0CorrFile->Get("fTriggerPt");
  fTriggerPtWithinEPBin = (TH1D *) fPi0CorrFile->Get("fTriggerPtWithinEPBin");

  fTrackPtProjectionSE = (TH1D *) fPi0CorrFile->Get("TrackPtProjectionSE");
  fTrackPtProjectionME = (TH1D *) fPi0CorrFile->Get("TrackPtProjectionME");

  fTrackPtFromTrackPsi = (TH1D *) fPi0CorrFile->Get("TrackPtFromTrackPsi");



	// Loading Histograms from Pi0 Corr File (and getting nObs)
	for (Int_t i = 0; i < 13; i++) {
		TH1D * fLocal = 0;
		TString fLocalName = Form("dPhi_ObsBin%d_Full",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) break;
		fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
		fFullDPhiPi0.push_back(fLocal);		

		fLocalName = Form("dPhi_ObsBin%d_NearEta",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) break;
		fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
		fNearEtaDPhiPi0.push_back(fLocal);		
  
		fLocalName = Form("dPhi_ObsBin%d_FarEta",i);
		fLocal = (TH1D *) fPi0CorrFile->Get(fLocalName);
		if (!fLocal) break;
		fLocal->SetName(Form("%s_Pi0",fLocalName.Data()));
		fFarEtaDPhiPi0.push_back(fLocal);		

		nObsBins++;
	}

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
			fLocal->SetName(Form("%s_SB%d",fLocalName.Data(),j));
			fLocalVector.push_back(fLocal);
		}
		fFarEtaDPhiSB.push_back(fLocalVector);
	}

  for (int i = 0; i < kNPtBins; i++) {
    TH1D * hPtEPAnglePionAcc_Pion_Indiv = 0;
    hPtEPAnglePionAcc_Pion_Indiv = (TH1D*) fPi0CorrFile->Get(Form("PtEPAnglePionAcc_Proj_%d",i));
    if (!hPtEPAnglePionAcc_Pion_Indiv) {
      fprintf(stderr,"Could not find hPtEPAnglePionAcc_Pion_Indiv %d\n",i);
      break;
    }
    hPtEPAnglePionAcc_Pion_Indiv->SetName(Form("PtEPAnglePionAcc_Pi0_Proj_%d",i));
    hPtEPAnglePionAcc_Proj_Pion.push_back(hPtEPAnglePionAcc_Pion_Indiv);

    TH1D * fLocal = 0;
    vector<TH1D*> fLocalVector = {}; 
    // Loop over sidebands
		for (Int_t j = 0; j < kNSB; j++) {
      fLocal = (TH1D *) fSidebandFile[j]->Get(Form("PtEPAnglePionAcc_Proj_%d",i));
      if (!fLocal) {
        fprintf(stderr,"Missing an ep angle proj %d\n",j);
        break;
      }
      fLocal->SetName(Form("PtEPAnglePionAcc_SB%d_Proj_%d",j,i));
      fLocal->SetTitle(Form("Sideband %d",j+1));
      fLocalVector.push_back(fLocal);
    }
    hPtEPAnglePionAcc_Proj_SB.push_back(fLocalVector);
  }

  // fLocal = (TH1D *) fSidebandFile[j]->Get(fLocalName);

	// Getting the Mass Distributions
	// NTF: merge in pt if fObs is not pT
	// FIXME need to do something about pt bin choice 
	//   in case doing Z_T, or Xi

  printf("Loading mass information\n");

	// Pi0 Candidate Masses
	for (Int_t i = 0; i < kGammaNBINS; i++) {
    printf("Seaching for mass information for gamma bin %d\n",i);
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

  if (fObservable == 0) {
    ObsArray = array_G_BinsValue;
    fGlobalMinPt = 5;
    fGlobalMaxPt = 17;
  }
  else if (fObservable == 1) ObsArray = array_ZT_BinsValue;
  else if (fObservable == 2) ObsArray = array_XI_BinsValue;

  for (Int_t i = 0; i <= nObsBins; i++) {
    fObsBins.push_back(ObsArray[i]);
  }
  
  memcpy(fPtBins,array_G_BinsValue,kGammaNBINS+1);
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
  leg->AddEntry(Histo,"ALICE Performance","");
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

		// Adding SB Integrals, masses
    Int_t graph_counter = 0;
		for (Int_t j = 0; j < fNSB; j++) {
      if (!fSidebandFitMask[j]) continue; 
      graph_counter++;
			Double_t fLocalIntegral_Un = 0;
			Double_t fLocalIntegral = fFullDPhiSB[i][j]->IntegralAndError(1,fFullDPhiSB[i][j]->GetNbinsX(),fLocalIntegral_Un);
			fMassVsIntegralPt_Full->SetPoint(graph_counter,fMeanMassSBVal[iMassIndex][j],fLocalIntegral);
			fMassVsIntegralPt_Full->SetPointError(graph_counter,fMeanMassSBVal_Un[iMassIndex][j],fLocalIntegral_Un);
			//fMassVsIntegralPt_Full->SetPoint(j+1,fMeanMassSBVal[iMassIndex][j],fLocalIntegral);
			//fMassVsIntegralPt_Full->SetPointError(j+1,fMeanMassSBVal_Un[iMassIndex][j],fLocalIntegral_Un);
		}

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
        fMassEffectFit_l->FixParameter(2,fMeanMassPi0Val[iMassIndex]);
        break;
      case 2:
        printf("Fitting the mass scaling with a quadratic fit\n");
        // General quadratic, written so that f([2]) is always = [0]
        fMassEffectFit_l = new TF1(Form("MassEffectFit_%d",i),"[0]+[1]*(x - [2])*(x - [3])",0.1,0.5);
        fMassEffectFit_l->FixParameter(2,fMeanMassPi0Val[iMassIndex]);
		  break;
    }


		fMassEffectFit_l->SetLineColor(kViolet);
		fMassVsIntegralPt_Full->Fit(fMassEffectFit_l,"Q","",fMassVsIntegralPt_Full->GetX()[0]+0.03,0.5);
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
	Int_t fSBColor[4] = {kAzure,kAzure-2,kAzure-4,kAzure-9};	
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
        if (!fSidebandFitMask[j]) continue;
				fMassPtBinSB[i][j]->Draw("SAME BAR");
				fMassPtBinSB[i][j]->SetLineColor(fSBColor[j]);
				fMassPtBinSB[i][j]->SetFillColor(fSBColor[j]);
        fMassPtBinSB[i][j]->SetFillColorAlpha(fSBColor[j],0.7);
				fMassPtBinSB[i][j]->SetMarkerColor(fSBColor[j]);

        bool bDoThingForPerformance=true;
        if (bDoThingForPerformance) {
          fMassPtBinAll[i]->GetYaxis()->SetTitle("Counts");
          lSBFigure->AddEntry(fMassPtBinSB[i][j],Form("SB %d",j),"flp");
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
	Int_t fSBColor[4] = {kAzure,kAzure-2,kAzure-4,kAzure-9};
  Int_t fSBStyle[4] = {kFullSquare,kFullCircle,kFullDiamond,kOpenSquare};
  TCanvas * cSBCmp = new TCanvas("SBCmp","SBCmp");
  TLegend * lSBCmp = new TLegend(0.43,0.55,0.88,0.88);
  for (Int_t i = 0; i < nObsBins; i++) {
    cSBCmp->Divide(1,2,0.0,0.0);
    cSBCmp->cd(1);
    TString fName = Form("dPhi_ObsBin%d_Full",i);

    Float_t fMin = 1e9;
    Float_t fMax = 0;

    fFullDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN^{assoc}}{d#Delta#varphi}");
    //fFullDPhiSB[i][0]->GetYaxis()->SetTitle("#frac{1}{N_{Trig}} #frac{d^{2}N^{assoc}}{d#Delta#eta d#Delta#phi}");
    fFullDPhiSB[i][0]->GetYaxis()->SetTitleOffset(0.7);

    TH1D * hSum = (TH1D *) fFullDPhiSB[i][0]->Clone(Form("Sum_%d",i));

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

      Float_t fLocalMin = fFullDPhiSB[i][j]->GetBinContent(fFullDPhiSB[i][j]->GetMinimumBin());
      Float_t fLocalMax = fFullDPhiSB[i][j]->GetBinContent(fFullDPhiSB[i][j]->GetMaximumBin());
      fMin = min(fMin,fLocalMin);
      fMax = max(fMax,fLocalMax);


      lSBCmp->AddEntry(fFullDPhiSB[i][j],Form("SB %d",j+1),"lp");
    }
    fFullDPhiSB[i][0]->GetYaxis()->SetRangeUser(fMin,fMax);
    lSBCmp->Draw("SAME");
    // Drawing the mean ratio
    hSum->Scale(1.0/fNSB);
    vector<TH1D *> hRatios = {};
    for (Int_t j = 0; j < fNSB; j++) {
      TH1D * hRatio = (TH1D *) fFullDPhiSB[i][j]->Clone(Form("%s_Ratio",fFullDPhiSB[i][j]->GetName()));
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

    PrintCanvas(cSBCmp,fName);
    lSBCmp->Clear();
    cSBCmp->Clear();
  }
}


// Use the chosen background histograms to produce the predicted background correlations
//   for dphi.
//   Uses the Sideband mass effect to produce the mass scale, and applies the 1-p scale
void TaskSideband::ProduceBackground() {
	cout<<"Producing Background Prediction"<<endl;
	fFullPredBkgDPhi = {};

	// nObsBins
	// FIXME iterate over full, near, far
	for (Int_t i = 0; i < nObsBins; i++) {
		fFullPredBkgDPhi.push_back(MergeAndScaleBkg(i,0));
	}
	for (Int_t i = 0; i < nObsBins; i++) {
		fNearEtaPredBkgDPhi.push_back(MergeAndScaleBkg(i,1));
	}
	for (Int_t i = 0; i < nObsBins; i++) {
		fFarEtaPredBkgDPhi.push_back(MergeAndScaleBkg(i,2));
	}

	TCanvas * cPredBkg = new TCanvas("PredBkg","PredBkg");
	Int_t nResults = (Int_t) fFullPredBkgDPhi.size(); // being lazy
	cPredBkg->Divide(3,2);
	for (Int_t i = 0; i < nResults; i++) {
		cPredBkg->cd(i+1);
		if (fFullPredBkgDPhi[i]) fFullPredBkgDPhi[i]->Draw();
	}
	PrintCanvas(cPredBkg,"PredBkg");
}

// Add index to change between full, near, far
TH1D * TaskSideband::MergeAndScaleBkg(Int_t index, Int_t iType = 0) {
	cout<<"Merging and Scaling background for Obs Bin "<<index<<endl;
	// fBackgroundSelection

//  Bool_t fSidebandMask[4] = {1,1,1,1};

  // FIXME move this part elsewhere, so it is only done once

	switch (fBackgroundSelection) {

		case 0: // Merge All 4
		default:
      fNSB = 4; // stays 4

//		fPredBkg = (TH1D *) fFullDPhiSB[index][0]->Clone(fName);
//		fPredBkg->GetXaxis()->SetTitle("#Delta#phi");
//		for (Int_t j = 1; j < fNSB; j++) fPredBkg->Add(fFullDPhiSB[index][j]);
		// scale by 1/4?
//		fPredBkg->Scale(1./(fNSB));
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

 

	// FIXME for scaling with the mass effect, i need the final mean mass.  This should be weighted between the four sidebands based on number of triggers

	TH1D * fPredBkg = 0;
	TString fName = ""; //Form("PredBkgDPhi_%d",index);

  switch (iType) {
    case 2:
      fName = Form("PredBkgDPhi_FarEta_%d",index);
      fPredBkg = (TH1D *) fFarEtaDPhiSB[index][0]->Clone(fName);
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) fPredBkg->Add(fFarEtaDPhiSB[index][j]);
    break;
    case 1:
      fName = Form("PredBkgDPhi_NearEta_%d",index);
      fPredBkg = (TH1D *) fNearEtaDPhiSB[index][0]->Clone(fName);
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) fPredBkg->Add(fNearEtaDPhiSB[index][j]);
    break;
    case 0:
    default:
      fName = Form("PredBkgDPhi_Full_%d",index);
      fPredBkg = (TH1D *) fFullDPhiSB[index][0]->Clone(fName);
      fPredBkg->Reset();
      for (Int_t j = 0; j < kNSB; j++) if (fSidebandMask[j]) fPredBkg->Add(fFullDPhiSB[index][j]);
  }
  fPredBkg->GetXaxis()->SetTitle("#Delta#phi");
  printf("Scaling background prediction by 1. / %d\n",fNSB);
  fPredBkg->Scale(1./(fNSB));
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
		iPurityIndex = iPtBin;
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


  //Double_t fPurity = Pi0YieldTotalRatio->GetY()[iPurityIndex]
  Double_t fPurity     = fPurityArray[iPurityIndex];
  Double_t fPurity_Err = fPurityArray[iPurityIndex];

//  if (iPurityChoice == 0) fPurity = 0;




	Double_t fPurityScale = 1. - fPurity ; // (1 - purity)

 
    

//	Double_t fPurityScale = 1. - Pi0YieldTotalRatio->GetY()[3]; // (1 - purity) // FIXME temp
//	Double_t fPurityScale = 1. - Pi0YieldTotalRatio->GetY()[index]; // (1 - purity)
	printf("Using mass scaling   = %f\n",fMassScale);
  printf("Found pi0 purity = %f\n",fPurity);
	printf("Using Purity scaling = %f\n",fPurityScale);
	fPredBkg->Scale(fMassScale * fPurityScale);	


	return fPredBkg;
}

void TaskSideband::PlotBkgAndSignal() {
	cout<<"Plotting Signal and Background"<<endl;

	TCanvas * cBkgSignal = new TCanvas("BkgSignal","BkgSignal");
	Int_t nX = 3;
	Int_t nY = 2;
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

    // FIXME adding an arbitrary scaling just for the performance plot 
    Float_t fMaxValue  = fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMaximumBin());
    Float_t fMoreArbScale = 2. / fMaxValue;
    if (fMaxValue < 1.8) {
      fFullDPhiPi0[i]->Scale(fMoreArbScale);
      fFullPredBkgDPhi[i]->Scale(fMoreArbScale);
    }


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
    DrawAlicePerf(fFullDPhiPi0[i],0.41,yPerfMin,0.24,0.15);
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



  // disabling rescale
    fFullPredBkgDPhi_Unscaled[i]->Scale(0.5*fAwaysideScale/fAwaysideScaleBkg);

		Double_t fMin = fmin(fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMinimumBin()),fFullPredBkgDPhi_Unscaled[i]->GetBinContent(fFullPredBkgDPhi_Unscaled[i]->GetMinimumBin()));
		Double_t fMax = fmax(fFullDPhiPi0[i]->GetBinContent(fFullDPhiPi0[i]->GetMaximumBin()),fFullPredBkgDPhi_Unscaled[i]->GetBinContent(fFullPredBkgDPhi[i]->GetMaximumBin()));
		fMax = 1.1 * fMax;
    fMin -= 0.1*(fMax - fMin);
		fFullDPhiPi0[i]->GetYaxis()->SetRangeUser(fMin,fMax);
		fFullDPhiPi0[i]->Draw();
//    fFullDPhiPi0[i]->SetMarkerSize(1.0);

    if (bNoYLabel) {
      fFullDPhiPi0[i]->GetYaxis()->SetLabelColor(kWhite);
      fFullDPhiPi0[i]->GetYaxis()->SetLabelSize(0.1);
      fFullDPhiPi0[i]->GetYaxis()->SetTitleOffset(0.6);
      //fFullDPhiPi0[i]->GetYaxis()->SetTitle("");
    }

		fFullPredBkgDPhi_Unscaled[i]->Draw("SAME");	

    leg2->Draw("SAME");
    //DrawAlicePerf(fFullDPhiPi0[i],0.35,0.72,0.33,0.15);
    //DrawAlicePerf(fFullDPhiPi0[i],0.35,0.40,0.38,0.2);
    PrintCanvas(cBkgSignal,Form("BkgSigCmp_Unscaled_%d",i));
  }

}

void TaskSideband::Subtract() {
	cout<<"Subtracting Scaled background from Pi0Cand Correlations"<<endl;

	fFullDPhiFinal = {};
	fNearEtaDPhiFinal = {};
	fFarEtaDPhiFinal = {};
	for (Int_t i = 0; i < nObsBins; i++) {
		TString fName = Form("SBSub_FullDPhi_ObsBin%d",i);
		TH1D * fFullDPhiFinal_Local = (TH1D *) fFullDPhiPi0[i]->Clone(fName.Data());
		fFullDPhiFinal_Local->Add(fFullPredBkgDPhi[i],-1);

		fName = Form("SBSub_NearEtaDPhi_ObsBin%d",i);
		TH1D * fNearEtaDPhiFinal_Local = (TH1D *) fNearEtaDPhiPi0[i]->Clone(fName.Data());
		fNearEtaDPhiFinal_Local->Add(fNearEtaPredBkgDPhi[i],-1);

		fName = Form("SBSub_FarEtaDPhi_ObsBin%d",i);
		TH1D * fFarEtaDPhiFinal_Local = (TH1D *) fFarEtaDPhiPi0[i]->Clone(fName.Data());
		fFarEtaDPhiFinal_Local->Add(fFarEtaPredBkgDPhi[i],-1);

    // Scale by 1 / Purity
    Int_t iPurityIndex = i;
    if (fObservable != 0) {
      iPurityIndex = iPtBin;
    }
    Double_t fPurity = Pi0YieldTotalRatio->GetY()[iPurityIndex];
    // Only the statistical uncertainty
    Double_t fPurity_Err = Pi0YieldTotalRatio->GetEY()[iPurityIndex];

    if (fPurity > 0) {
      fFullDPhiFinal_Local->Scale(1./fPurity);
      fNearEtaDPhiFinal_Local->Scale(1./fPurity);
      fFarEtaDPhiFinal_Local->Scale(1./fPurity);
    }

    fFullDPhiFinal_Local->SetTitle(""); // FIXME I hope this doesn't affect anything later
    fFullDPhiFinal_Local->GetXaxis()->SetTitle("#Delta#phi^{#pi^{0}-h}");
    fFullDPhiFinal_Local->GetYaxis()->SetTitle("1/N^{#pi^{0}} dN/d#Delta#varphi");
    //fFullDPhiFinal_Local->GetYaxis()->SetTitle("1/N^{#pi^{0}} d^{2}N/d#Delta#etad#Delta#phi");

		fFullDPhiFinal.push_back(fFullDPhiFinal_Local);
		fNearEtaDPhiFinal.push_back(fNearEtaDPhiFinal_Local);
		fFarEtaDPhiFinal.push_back(fFarEtaDPhiFinal_Local);

	}

  // Draw a nice plot of the subtracted and rescaled correlations
  // Drawing Full range
  TCanvas * cSub = new TCanvas ("Sub","Sub",fCanvasWidth,fCanvasHeight);
  cSub->Divide(3,2,0.001,0.0012);
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

  Int_t fOrigColor = kBlack;// kOrange+1;
  Int_t fSubColor = kOrange+1;

  // Draw a comparison to the input correlations
  //TCanvas * cSubCmp = new TCanvas ("SubCmp","SubCmp",fCanvasWidth,fCanvasHeight);
  cSub->Clear();
  cSub->Divide(3,2,0.001,0.0012);
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

  // Draw individual plots
  cSub->Clear();
  cSub->SetWindowSize(400,600);
  //cSub->SetWidth(400);
  //cSub->SetHeight(600);
  for (Int_t i = 0; i < nResults; i++) {
    TH1D * hLocalSub = fFullDPhiFinal[i];
    TH1D * hLocalOrig = fFullDPhiPi0[i];

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
    hPtEPAnglePionAcc_Proj_Pion[i]->SetMarkerStyle(kFullSquare);
    for (int j = 0; j < fNSB; j++) {
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

  TCanvas * cFlowCanvas = new TCanvas("FlowCanvas","FlowCanvas");
  cFlowCanvas->cd();
  


  TLegend * leg = new TLegend(0.65,0.67,0.92,0.92);
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
  
  printf("Preparing to calculate vn before subtraction\n");

  printf("Normalizing ... \n");
  // Normalize the 5 histograms
  // Or their clones?
  for (int i = 0; i < kNPtBins; i++) {
    double scale = hPtEPAnglePionAcc_Proj_Pion[i]->Integral("width");
    if (scale > 0) hPtEPAnglePionAcc_Proj_Pion[i]->Scale(1./scale);
    for (int j = 0; j < fNSB; j++) {
      scale = hPtEPAnglePionAcc_Proj_SB[i][j]->Integral("width");
      if(scale > 0) hPtEPAnglePionAcc_Proj_SB[i][j]->Scale(1./scale);
    
    }
  }

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



  // Compute Vn prior to subtraction
  CalculateVnSub(0);

  gTriggerPreSub_V2 = new TGraphErrors(kNPtBins);
  gTriggerPreSub_V2->SetName("TriggerPreSub_V2");

  // Compare Vn Calculations
//  TLegend * leg1 = new TLegend(0.7,0.6,0.87,0.87);
  
//  leg1->Draw("SAME");
//  cFlowCanvas->Print(Form("%s/FlowAll_V2T.pdf",fOutputDir.Data())); 
//  cFlowCanvas->Print(Form("%s/CFiles/FlowAll_V2T.C",fOutputDir.Data())); 





  // Subtract


  cFlowCanvas->Clear();
  cFlowCanvas->Print(Form("%s/FlowPi0_Subtraction.pdf",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/FlowPi0_Subtraction.png",fOutputDir.Data())); 
  cFlowCanvas->Print(Form("%s/CFiles/FlowPi0_Subtraction.C",fOutputDir.Data())); 




//  printf("Preparing to calculate vn after subtraction\n");
//  CalculateVnSub(1);
  

}


// Note: these are not corrected for Event Plane Resolution
// May copy this to phase 4, where the EPR is easily applied
void TaskSideband::CalculateVnSub(int iPostSub = 0) {

  TH1D * fHist = 0;
  TF1 * fFlowFunction = 0;

  TGraphErrors * gTriggerFlow_V2 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V2->SetName("TriggerFlow_V2");
  TGraphErrors * gTriggerFlow_V4 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V4->SetName("TriggerFlow_V4");
  TGraphErrors * gTriggerFlow_V6 = new TGraphErrors(kNPtBins);
  gTriggerFlow_V6->SetName("TriggerFlow_V6");

  for (int j = 0; j < fNSB; j++) {
    TGraphErrors * gTriggerFlowSB_V2 = new TGraphErrors(kNPtBins);
    gTriggerFlowSB_V2->SetName(Form("TriggerFlowSB_V2_SB%d",j));
    TGraphErrors * gTriggerFlowSB_V4 = new TGraphErrors(kNPtBins);
    gTriggerFlowSB_V4->SetName(Form("TriggerFlowSB_V4_SB%d",j));
    TGraphErrors * gTriggerFlowSB_V6 = new TGraphErrors(kNPtBins);
    gTriggerFlowSB_V6->SetName(Form("TriggerFlowSB_V6_SB%d",j));
    gTriggerSidebands_V2.push_back(gTriggerFlowSB_V2);
    gTriggerSidebands_V4.push_back(gTriggerFlowSB_V4);
    gTriggerSidebands_V6.push_back(gTriggerFlowSB_V6);
  }

  for (int i = 0; i < kNPtBins; i++) {
    double fMinPt = fPtBins[i];
    double fMaxPt = fPtBins[i+1];

    if (iPostSub == 0) {
      fHist = hPtEPAnglePionAcc_Proj_Pion[i];
    } else {
      fHist = hPtEPAnglePionAcc_Proj_PionPostSub[i];
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

    fHist->Fit(fFlowFunction);

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

        fSBFlowHist->Fit(fSBFlowFunction);
        gTriggerSidebands_V2[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(1));
        gTriggerSidebands_V4[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(2));
        gTriggerSidebands_V6[j]->SetPoint(i,(fMinPt+fMaxPt)/2.,fSBFlowFunction->GetParameter(3));
        gTriggerSidebands_V2[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(1));
        gTriggerSidebands_V4[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(2));
        gTriggerSidebands_V6[j]->SetPointError(i,(fMaxPt-fMinPt)/2.,fSBFlowFunction->GetParError(3));
      }
    }
  }

  if (iPostSub == 0) {
    gTriggerFlow_V2->SetName("TriggerPreSub_V2");
    gTriggerFlow_V4->SetName("TriggerPreSub_V4");
    gTriggerFlow_V6->SetName("TriggerPreSub_V6");

    gTriggerPreSub_V2 = gTriggerFlow_V2;
    gTriggerPreSub_V4 = gTriggerFlow_V4;
    gTriggerPreSub_V6 = gTriggerFlow_V6;
  } else {
    gTriggerFlow_V2->SetName("TriggerPostSub_V2");
    gTriggerFlow_V4->SetName("TriggerPostSub_V4");
    gTriggerFlow_V6->SetName("TriggerPostSub_V6");

    gTriggerPostSub_V2 = gTriggerFlow_V2;
    gTriggerPostSub_V4 = gTriggerFlow_V4;
    gTriggerPostSub_V6 = gTriggerFlow_V6;
  }


}





void TaskSideband::SaveResults() {
	cout<<"Saving Results"<<endl;

  TFile * fOutputFile = TFile::Open(Form("output/%s",fOutputFileName.Data()),"RECREATE");
	if (!fOutputFile) return;

  if (VariableInfo) fOutputFile->Add(VariableInfo);

  if (fTriggerPt) fOutputFile->Add(fTriggerPt); 
  if (fTriggerPtWithinEPBin) fOutputFile->Add(fTriggerPtWithinEPBin); 

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


  if (gTrigger_Bv) fOutputFile->Add(gTrigger_Bv);
  if (gTrigger_V2) fOutputFile->Add(gTrigger_V2);
  if (gTrigger_V4) fOutputFile->Add(gTrigger_V4);
  if (gTrigger_V6) fOutputFile->Add(gTrigger_V6);

  if (gTrack_Bv) fOutputFile->Add(gTrack_Bv);
  if (gTrack_V2) fOutputFile->Add(gTrack_V2);
  if (gTrack_V4) fOutputFile->Add(gTrack_V4);
  if (gTrack_V6) fOutputFile->Add(gTrack_V6);

  if (gTriggerPreSub_V2) fOutputFile->Add(gTriggerPreSub_V2);
  if (gTriggerPreSub_V4) fOutputFile->Add(gTriggerPreSub_V4);
  if (gTriggerPreSub_V6) fOutputFile->Add(gTriggerPreSub_V6);

  for (auto gTriggerSidebands_V2_Indiv : gTriggerSidebands_V2) fOutputFile->Add(gTriggerSidebands_V2_Indiv);
  for (auto gTriggerSidebands_V4_Indiv : gTriggerSidebands_V4) fOutputFile->Add(gTriggerSidebands_V4_Indiv);
  for (auto gTriggerSidebands_V6_Indiv : gTriggerSidebands_V6) fOutputFile->Add(gTriggerSidebands_V6_Indiv);

  if (gTriggerPostSub_V2) fOutputFile->Add(gTriggerPostSub_V2);
  if (gTriggerPostSub_V4) fOutputFile->Add(gTriggerPostSub_V4);
  if (gTriggerPostSub_V6) fOutputFile->Add(gTriggerPostSub_V6);

  for (auto hPtEPAnglePionAcc_Proj_Pion_Indiv : hPtEPAnglePionAcc_Proj_Pion) fOutputFile->Add(hPtEPAnglePionAcc_Proj_Pion_Indiv);
  for (int i = 0; i < (int) hPtEPAnglePionAcc_Proj_SB.size(); i++) {
    for (auto hPtEPAnglePionAcc_Proj_SB_Indiv : hPtEPAnglePionAcc_Proj_SB[i]) {
      fOutputFile->Add(hPtEPAnglePionAcc_Proj_SB_Indiv);
    }
  }
  for (auto hPtEPAnglePionAcc_Proj_PionPostSub_Indiv : hPtEPAnglePionAcc_Proj_PionPostSub) fOutputFile->Add(hPtEPAnglePionAcc_Proj_PionPostSub_Indiv);

  fOutputFile->Write();
  fOutputFile->Close();
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
}

void TaskSideband::Run() {
	cout<<"Beginning Sideband Subtraction Task"<<endl;
	
	LoadPurity();
	LoadHistograms();
	InitArrays();
	Debug(1);

	MassAnalysis();
	ProduceSidebandFigure();
  ProduceSidebandComparison();
//	SimpleNormFit();
	ProduceBackground();


	PlotBkgAndSignal();

	Subtract();

  ProcessFlow();

	SaveResults();

  // Debugging masks
  printf("  SidebandFitMask = %d %d %d %d\n",fSidebandFitMask[0],fSidebandFitMask[1],fSidebandFitMask[2],fSidebandFitMask[3]);
  printf("  SidebandMask = %d %d %d %d\n",fSidebandMask[0],fSidebandMask[1],fSidebandMask[2],fSidebandMask[3]);

}



